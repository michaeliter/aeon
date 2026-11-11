//! GTDB-Tk-compatible genome classification: ANI screen, AAI, then a
//! marker-gene search.

pub(crate) mod aai;
pub(crate) mod ani;
pub(crate) mod external;
pub(crate) mod genes;
pub(crate) mod log;
pub(crate) mod markers;
pub(crate) mod refdata;
pub(crate) mod refdb;
pub(crate) mod run_record;
pub(crate) mod summary;

use rayon::prelude::*;
use std::collections::{BTreeMap, HashSet};
use std::error::Error;
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

use crate::classify::aai::{AaiDecision, AaiHit, compute_aai};
use crate::classify::ani::{
    AniDecision, AniHit, Candidate, MIN_SPECIES_RADIUS, Screen, decide, refine_candidates, screen,
    sketch_ani,
};
use crate::classify::external::{find_program, find_tool};
use crate::classify::genes::{GeneCaller, TRANSLATION_TABLE};
use crate::classify::log::ClassifyLog;
use crate::classify::markers::{
    CopyNumber, MARKER_SETS, MarkerTools, PfamInfo, align_markers, apply_mask, filtered_reason,
    genome_domain, identify_markers, marker_hmms, marker_sequence, multihit_warning, read_mask,
    write_filtered, write_marker_info, write_markers_summary, write_msa_gz,
};
use crate::classify::refdata::{marker_set_for_taxonomy, standardise_taxonomy};
use crate::classify::refdb::{ClassifyDb, Index};
use crate::classify::run_record::{ALIGN, ANI_SCREEN, CLASSIFY, IDENTIFY, RunInfo, RunRecord};
use crate::classify::summary::{
    AniScreenRow, NONE_VALUE, RelatedReference, SummaryRow, format_related_references, py_round,
    write_ani_summary, write_summary, write_unassigned_hits,
};
use crate::database::search_vectors;
use crate::sketch::sketch_from_params;
use crate::utils::read_list_file;
use crate::warmup::{Preload, Warmup};

pub(crate) enum QueryInput {
    /// Every file in `dir` ending in `.{extension}`; IDs drop the extension.
    GenomeDir { dir: String, extension: String },
    /// GTDB-Tk batchfile: `path[\tgenome_id]` per line.
    Batchfile(String),
}

pub(crate) struct ClassifyOptions {
    pub(crate) db: String,
    pub(crate) input: QueryInput,
    pub(crate) out_dir: String,
    pub(crate) prefix: String,
    pub(crate) min_af: f64,
    pub(crate) candidates: usize,
    pub(crate) beam_width: usize,
    pub(crate) min_sketch_ani: f64,
    pub(crate) min_marker_identity: f64,
    pub(crate) genus_aai: f64,
    pub(crate) family_aai: f64,
    pub(crate) aai_margin: f64,
    pub(crate) fgs: String,
    pub(crate) fgs_train_dir: Option<String>,
    pub(crate) diamond: String,
    pub(crate) min_perc_aa: f64,
    pub(crate) bin_dir: Option<String>,
    /// Read each index into the page cache before searching it.
    pub(crate) preload: bool,
    pub(crate) threads: usize,
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct QueryGenome {
    pub(crate) id: String,
    pub(crate) path: String,
}

pub(crate) fn resolve_queries(input: &QueryInput) -> Result<Vec<QueryGenome>, Box<dyn Error>> {
    let mut genomes = match input {
        QueryInput::GenomeDir { dir, extension } => {
            let suffix = format!(".{}", extension.trim_start_matches('.'));
            let mut genomes = Vec::new();
            for entry in fs::read_dir(dir)? {
                let path = entry?.path();
                let Some(name) = path.file_name().and_then(|name| name.to_str()) else {
                    continue;
                };
                if let Some(id) = name.strip_suffix(&suffix) {
                    if path.is_file() && !id.is_empty() {
                        genomes.push(QueryGenome {
                            id: id.to_owned(),
                            path: path.to_string_lossy().into_owned(),
                        });
                    }
                }
            }
            genomes.sort_by(|left, right| left.id.cmp(&right.id));
            genomes
        }
        QueryInput::Batchfile(path) => read_list_file(path)?
            .into_iter()
            .map(|line| {
                let mut fields = line.split('\t');
                let path = fields.next().unwrap_or_default().to_owned();
                let id = fields
                    .next()
                    .map(str::to_owned)
                    .unwrap_or_else(|| file_stem(&path));
                QueryGenome { id, path }
            })
            .collect(),
    };
    if genomes.is_empty() {
        return Err("no query genomes were found".into());
    }
    let mut ids = HashSet::with_capacity(genomes.len());
    for genome in &genomes {
        if !ids.insert(genome.id.as_str()) {
            return Err(format!("duplicate query genome ID: {}", genome.id).into());
        }
        if !Path::new(&genome.path).is_file() {
            return Err(format!("query genome does not exist: {}", genome.path).into());
        }
    }
    genomes.shrink_to_fit();
    Ok(genomes)
}

/// File name without directories or FASTA/compression extensions.
fn file_stem(path: &str) -> String {
    refdata::file_stem(path, &[".fna", ".fa", ".fasta", ".fas"]).to_owned()
}

/// Everything known about one query after each tier it reached.
struct QueryResult<'a> {
    genome: &'a QueryGenome,
    /// The ANI screen's hits, reported for every genome.
    hits: Vec<AniHit>,
    stage: Stage,
}

/// How far down the cascade a genome went.
enum Stage {
    /// The ANI screen assigned a species.
    Species(AniHit),
    /// The ANI screen left it unresolved and genes are not yet called.
    Unresolved,
    /// FragGeneScanRs called no genes.
    NoGenes,
    /// Genes were called; the AAI and marker tiers fill in their results
    /// as each runs.
    Called {
        protein_faa: PathBuf,
        aai: Option<AaiResult>,
        markers: Option<MarkerResult>,
    },
}

struct AaiResult {
    hits: Vec<AaiHit>,
    decision: AaiDecision,
}

/// Where the marker pipeline left a genome.
#[derive(Debug)]
enum MarkerOutcome {
    /// In the masked user MSA and searched against the marker index.
    Aligned {
        msa_percent: f64,
        closest: Option<MarkerHit>,
    },
    /// No single-copy marker of either set.
    NoMarkers,
    /// Too few amino acids in the masked MSA; carries GTDB-Tk's reason.
    Filtered(String),
}

struct MarkerResult {
    marker_set: &'static str,
    outcome: MarkerOutcome,
    /// Domain-questionable and multiple-hit warnings.
    warnings: Vec<String>,
}

/// The closest reference in the marker index.
#[derive(Debug, Clone, PartialEq)]
struct MarkerHit {
    accession: String,
    /// Sketch-estimated identity of the concatenated markers, 0-100.
    identity: f64,
}

impl Stage {
    /// The proteins of a genome whose genes were called.
    fn protein_faa(&self) -> Option<&Path> {
        match self {
            Self::Called { protein_faa, .. } => Some(protein_faa),
            _ => None,
        }
    }

    /// Genes called, and unresolved by AAI where it ran.
    fn needs_markers(&self) -> bool {
        matches!(self, Self::Called { aai, markers: None, .. }
            if aai.as_ref().is_none_or(|aai| aai.decision == AaiDecision::Unresolved))
    }
}

pub(crate) fn classify_genomes(options: &ClassifyOptions) -> Result<(), Box<dyn Error>> {
    let started = Instant::now();
    let out_dir = PathBuf::from(&options.out_dir);
    fs::create_dir_all(&out_dir)?;
    let log = ClassifyLog::create(&out_dir, &options.prefix)?;
    log.info(format!("aeon v{} classify", env!("CARGO_PKG_VERSION")));

    let queries = resolve_queries(&options.input)?;
    let db = ClassifyDb::open(&options.db)?;
    log.info(format!(
        "Using {} reference database {} ({} reference genomes)",
        db.manifest.gtdb_release,
        options.db,
        db.nt.paths.len()
    ));
    log.info(format!("Identified {} genomes to classify", queries.len()));
    // Warm all three indexes in search order while the queries are processed.
    let warmup = Warmup::start(
        &[&db.nt.prefix, &db.aa.prefix, &db.markers.prefix],
        options.preload,
    )?;
    if !options.preload {
        log.info(Preload::Disabled.to_string());
    }

    let mut record = RunRecord::start();
    let mut results = run_ani_tier(&queries, &db, &warmup, options, &out_dir, &log, &mut record)?;
    record.mark(ANI_SCREEN);
    call_query_genes(&mut results, options, &out_dir, &log)?;
    if results
        .iter()
        .any(|result| !matches!(result.stage, Stage::Species(_)))
    {
        record.mark(IDENTIFY);
    }
    run_aai_tier(
        &mut results,
        &db,
        &warmup,
        options,
        &out_dir,
        &log,
        &mut record,
    )?;
    record.mark(CLASSIFY);
    run_marker_tier(
        &mut results,
        &db,
        &warmup,
        options,
        &out_dir,
        &log,
        &mut record,
    )?;

    write_identify_files(&results, &out_dir, &options.prefix, &mut record)?;
    write_summaries(&results, &db, options, &out_dir, &log, &mut record)?;
    record.mark(CLASSIFY);
    let (genome_dir, extension, batchfile) = match &options.input {
        QueryInput::GenomeDir { dir, extension } => {
            (Some(dir.clone()), Some(extension.clone()), None)
        }
        QueryInput::Batchfile(path) => (None, None, Some(path.clone())),
    };
    record.write(&RunInfo {
        command_line: std::iter::once("aeon".to_owned())
            .chain(std::env::args().skip(1))
            .collect::<Vec<_>>()
            .join(" "),
        database_version: db.manifest.gtdb_release.clone(),
        database_path: PathBuf::from(&options.db),
        out_dir: out_dir.clone(),
        prefix: options.prefix.clone(),
        genome_dir,
        batchfile,
        extension,
        min_af: options.min_af,
        min_perc_aa: options.min_perc_aa,
    })?;
    log.info(format!(
        "Done. Classified {} genomes in {:.3}s",
        queries.len(),
        started.elapsed().as_secs_f64()
    ));
    Ok(())
}

const NO_GENES: &str = "No genes were called by FragGeneScanRs";
const NO_MARKERS: &str = "No bacterial or archaeal marker";

fn marker_genes_dir(out_dir: &Path) -> PathBuf {
    out_dir
        .join("identify")
        .join("intermediate_results")
        .join("marker_genes")
}

/// Screen every query against the nucleotide index: a species from the
/// sketch alone when the candidates name one species, else TurboANI against
/// the candidates and GTDB-Tk's radius rule.
fn run_ani_tier<'a>(
    queries: &'a [QueryGenome],
    db: &ClassifyDb,
    warmup: &Warmup,
    options: &ClassifyOptions,
    out_dir: &Path,
    log: &ClassifyLog,
    record: &mut RunRecord,
) -> Result<Vec<QueryResult<'a>>, Box<dyn Error>> {
    log.info("Running ANI screen with DiskANN candidate search and TurboANI");
    let query_paths = queries
        .iter()
        .map(|genome| genome.path.clone())
        .collect::<Vec<_>>();
    let vectors = sketch_from_params(&query_paths, &db.nt.params);
    // Per query, the sketch-only species hit or the candidates for TurboANI.
    let (estimated, exact): (Vec<_>, Vec<_>) =
        search_index(&db.nt, warmup, &vectors, options, log)?
            .into_iter()
            .map(|neighbors| {
                let candidates = neighbors
                    .into_iter()
                    .map(|(reference, distance)| Candidate {
                        reference,
                        sketch_ani: sketch_ani(distance, db.nt.params.kmer_size),
                    })
                    .collect::<Vec<_>>();
                match screen(
                    &candidates,
                    &db.nt.accessions,
                    &db.data,
                    options.min_sketch_ani,
                ) {
                    Screen::Estimated(hit) => (Some(hit), Vec::new()),
                    Screen::Exact(references) => (None, references),
                }
            })
            .unzip();
    let turboani = exact
        .iter()
        .filter(|references| !references.is_empty())
        .count();
    log.info(format!(
        "TurboANI checks candidates for {turboani} genome(s); {} matched one species by sketch alone.",
        estimated.iter().flatten().count()
    ));
    let exact_hits = refine_candidates(&query_paths, &exact, &db.nt.paths)?;
    let results = queries
        .iter()
        .zip(estimated)
        .zip(exact_hits)
        .map(|((genome, estimated), exact_hits)| {
            let hits = estimated.map_or(exact_hits, |hit| vec![hit]);
            let decision = decide(&hits, &db.nt.accessions, &db.data, options.min_af);
            ::log::debug!(
                "ani genome={} decision={decision:?} hits={}",
                genome.id,
                hits.iter()
                    .map(|hit| format!(
                        "{}:{:.2}/{:?}",
                        db.nt.accessions[hit.reference as usize], hit.ani, hit.af
                    ))
                    .collect::<Vec<_>>()
                    .join(",")
            );
            let stage = match decision {
                AniDecision::Species { hit } => Stage::Species(hit),
                AniDecision::Unresolved => Stage::Unresolved,
            };
            QueryResult {
                genome,
                hits,
                stage,
            }
        })
        .collect::<Vec<_>>();
    let species = results
        .iter()
        .filter(|result| matches!(result.stage, Stage::Species(_)))
        .count();
    log.info(format!(
        "{species} genome(s) have been classified using the ANI pre-screening step."
    ));

    // The hits of unassigned genomes, kept by GTDB-Tk for its near-miss
    // warnings; sketch estimates carry no AF and only ever assign species.
    let mut unassigned = results
        .iter()
        .filter(|result| !matches!(result.stage, Stage::Species(_)))
        .collect::<Vec<_>>();
    unassigned.sort_by(|left, right| left.genome.id.cmp(&right.genome.id));
    let path = out_dir.join("classify").join("ani_screen").join(format!(
        "{}.ani_screen_unassigned_hits.tsv.gz",
        options.prefix
    ));
    write_unassigned_hits(
        &path,
        unassigned.iter().flat_map(|result| {
            result.hits.iter().filter_map(|hit| {
                Some((
                    result.genome.id.as_str(),
                    db.nt.accessions[hit.reference as usize].as_str(),
                    hit.ani,
                    hit.af?,
                ))
            })
        }),
    )?;
    log.info(format!(
        "ANI hits of {} genome(s) not assigned by the ANI screen saved to: {}",
        unassigned.len(),
        path.display()
    ));
    record.output(ANI_SCREEN, "all", &path);
    Ok(results)
}

/// Write the `{prefix}.{marker_set}.summary.tsv` files, linked from the
/// output directory, and the ANI screen's summaries.
fn write_summaries(
    results: &[QueryResult],
    db: &ClassifyDb,
    options: &ClassifyOptions,
    out_dir: &Path,
    log: &ClassifyLog,
    record: &mut RunRecord,
) -> Result<(), Box<dyn Error>> {
    let mut rows: BTreeMap<&'static str, Vec<SummaryRow>> = BTreeMap::new();
    let mut ani_rows: BTreeMap<&'static str, Vec<AniScreenRow>> = BTreeMap::new();
    for result in results {
        let (marker_set, row) = summary_row(result, db, options);
        if row.warnings.iter().any(|warning| !warning.is_empty()) {
            log.warn(format!("{}: {}", row.user_genome, row.warnings.join(";")));
        }
        if let Stage::Species(hit) = &result.stage {
            ani_rows.entry(marker_set).or_default().push(AniScreenRow {
                user_genome: result.genome.id.clone(),
                reference_genome: db.nt.accessions[hit.reference as usize].clone(),
                ani: hit.ani,
                af: hit.af,
                reference_taxonomy: row.closest_genome_taxonomy.clone().unwrap_or_default(),
                other_related_references: row.other_related_references.clone(),
            });
        }
        rows.entry(marker_set).or_default().push(row);
    }

    for (marker_set, rows) in &ani_rows {
        let path = out_dir
            .join("classify")
            .join("ani_screen")
            .join(format!("{}.{marker_set}.ani_summary.tsv", options.prefix));
        write_ani_summary(&path, rows)?;
        record.output(ANI_SCREEN, marker_set, &path);
    }
    for (marker_set, rows) in &rows {
        let name = format!("{}.{marker_set}.summary.tsv", options.prefix);
        let path = out_dir.join("classify").join(&name);
        write_summary(&path, rows)?;
        record.output(CLASSIFY, marker_set, &path);
        link_summary(out_dir, &name)?;
        log.info(format!(
            "{} {marker_set} summary written to: {}",
            rows.len(),
            out_dir.join(&name).display()
        ));
    }
    Ok(())
}

/// Call genes with FragGeneScanRs for genomes the ANI tier left unresolved.
fn call_query_genes(
    results: &mut [QueryResult],
    options: &ClassifyOptions,
    out_dir: &Path,
    log: &ClassifyLog,
) -> Result<(), Box<dyn Error>> {
    let pending = results
        .iter()
        .enumerate()
        .filter(|(_, result)| matches!(result.stage, Stage::Unresolved))
        .map(|(index, _)| index)
        .collect::<Vec<_>>();
    if pending.is_empty() {
        return Ok(());
    }
    let caller = GeneCaller::locate(&options.fgs, options.fgs_train_dir.as_deref())?;
    log.info(format!(
        "Calling genes with FragGeneScanRs for {} genome(s) not resolved by ANI",
        pending.len()
    ));
    let marker_dir = marker_genes_dir(out_dir);
    // FragGeneScanRs splits a genome's contigs across its threads.
    let job_threads = (options.threads / pending.len()).max(1);
    let genes = pending
        .par_iter()
        .map(|&index| {
            let genome = results[index].genome;
            caller
                .call(&genome.path, &genome.id, &marker_dir, job_threads, true)
                .map_err(|error| format!("gene calling failed for {}: {error}", genome.id))
        })
        .collect::<Result<Vec<_>, _>>()?;
    for (index, protein_faa) in pending.into_iter().zip(genes) {
        results[index].stage = match protein_faa {
            Some(protein_faa) => Stage::Called {
                protein_faa,
                aai: None,
                markers: None,
            },
            None => Stage::NoGenes,
        };
    }
    Ok(())
}

/// Search `index` for each sketch's `--candidates` nearest references, once
/// the index's warm-up has finished.
fn search_index(
    index: &Index,
    warmup: &Warmup,
    vectors: &[Vec<u16>],
    options: &ClassifyOptions,
    log: &ClassifyLog,
) -> Result<Vec<Vec<(u32, f32)>>, Box<dyn Error>> {
    match warmup.wait(&index.prefix)? {
        Preload::Disabled => {}
        preload => log.info(preload.to_string()),
    }
    search_vectors(
        &index.prefix,
        index.len(),
        vectors,
        options.candidates,
        options.beam_width,
    )
}

/// Classify genomes the ANI tier left unresolved by AAI against candidate
/// reference proteomes.
fn run_aai_tier(
    results: &mut [QueryResult],
    db: &ClassifyDb,
    warmup: &Warmup,
    options: &ClassifyOptions,
    out_dir: &Path,
    log: &ClassifyLog,
    record: &mut RunRecord,
) -> Result<(), Box<dyn Error>> {
    let called = results
        .iter()
        .enumerate()
        .filter_map(|(index, result)| {
            result
                .stage
                .protein_faa()
                .map(|faa| (index, faa.to_path_buf()))
        })
        .collect::<Vec<_>>();
    if called.is_empty() {
        return Ok(());
    }
    let diamond = find_program(&options.diamond)?;
    log.info(format!(
        "Running AAI search with DiskANN candidate search and DIAMOND for {} genome(s)",
        called.len()
    ));
    let proteome_paths = called
        .iter()
        .map(|(_, faa)| faa.to_string_lossy().into_owned())
        .collect::<Vec<_>>();
    let vectors = sketch_from_params(&proteome_paths, &db.aa.params);
    let candidates = search_index(&db.aa, warmup, &vectors, options, log)?
        .into_iter()
        .map(|neighbors| neighbors.into_iter().map(|(id, _)| id).collect::<Vec<_>>())
        .collect::<Vec<_>>();
    let proteomes = called
        .iter()
        .map(|(_, faa)| faa.as_path())
        .collect::<Vec<_>>();
    let aai_hits = compute_aai(
        &diamond,
        &out_dir
            .join("classify")
            .join("aai")
            .join("intermediate_results"),
        &proteomes,
        &candidates,
        &db.aa.paths,
        options.threads,
    )?;

    let mut resolved = 0;
    for ((index, _), hits) in called.iter().zip(aai_hits) {
        let index = *index;
        let decision = aai::decide(
            &hits,
            &db.nt.accessions,
            &db.data,
            options.genus_aai,
            options.family_aai,
            options.aai_margin,
        );
        if decision != AaiDecision::Unresolved {
            resolved += 1;
        }
        ::log::debug!(
            "aai genome={} decision={:?} hits={}",
            results[index].genome.id,
            decision,
            hits.iter()
                .map(|hit| format!(
                    "{}:{:.2}/{}",
                    db.nt.accessions[hit.reference as usize], hit.aai, hit.reciprocal_hits
                ))
                .collect::<Vec<_>>()
                .join(",")
        );
        if let Stage::Called { aai, .. } = &mut results[index].stage {
            *aai = Some(AaiResult { hits, decision });
        }
    }
    log.info(format!(
        "{resolved} genome(s) have been classified using AAI."
    ));
    let path = write_aai_summary(results, db, out_dir, &options.prefix)?;
    record.output(CLASSIFY, "all", &path);
    Ok(())
}

/// GTDB-Tk's `identify` and `align` steps for genomes neither ANI nor AAI
/// resolved (marker search, copy number, domain call, and the masked
/// concatenated alignment), then a search of the marker index for the
/// closest reference.
fn run_marker_tier(
    results: &mut [QueryResult],
    db: &ClassifyDb,
    warmup: &Warmup,
    options: &ClassifyOptions,
    out_dir: &Path,
    log: &ClassifyLog,
    record: &mut RunRecord,
) -> Result<(), Box<dyn Error>> {
    let pending = results
        .iter()
        .enumerate()
        .filter(|(_, result)| result.stage.needs_markers())
        .map(|(index, result)| {
            let faa = result.stage.protein_faa().expect("genes").to_path_buf();
            (index, result.genome.id.clone(), faa)
        })
        .collect::<Vec<_>>();
    if pending.is_empty() {
        return Ok(());
    }
    let root = db.root().to_path_buf();
    let (hmmsearch, hmmalign) = (
        find_tool(options.bin_dir.as_deref(), "hmmsearch")?,
        find_tool(options.bin_dir.as_deref(), "hmmalign")?,
    );

    log.info(format!(
        "Identifying markers in {} genomes with {} threads.",
        thousands(pending.len()),
        options.threads
    ));
    let pfam_info = PfamInfo::load(&root)?;
    let tools = MarkerTools {
        hmmsearch: &hmmsearch,
        hmmalign: &hmmalign,
        root: &root,
        pfam_info: &pfam_info,
    };
    let proteomes = pending
        .iter()
        .map(|(_, id, faa)| (id.clone(), faa.clone()))
        .collect::<Vec<_>>();
    let searched = identify_markers(
        &tools,
        &proteomes,
        &marker_genes_dir(out_dir),
        options.threads,
    )?
    .into_iter()
    .zip(&pending)
    .map(|((id, copies), (index, _, _))| (*index, id, copies))
    .collect::<Vec<_>>();
    let mut summaries: BTreeMap<&str, BTreeMap<String, CopyNumber>> = BTreeMap::new();
    for (_, id, copies) in &searched {
        for (set, copies) in MARKER_SETS.iter().zip(copies) {
            summaries
                .entry(set)
                .or_default()
                .insert(id.clone(), copies.clone());
        }
    }
    let identify_dir = out_dir.join("identify");
    for (set, genomes) in &summaries {
        let path = identify_dir.join(format!("{}.{set}.markers_summary.tsv", options.prefix));
        write_markers_summary(&path, genomes)?;
        record.output(IDENTIFY, set, &path);
    }
    log.info("Done.");
    record.mark(IDENTIFY);

    let align_dir = out_dir.join("align");
    let align_tmp = align_dir.join("intermediate_results");
    fs::create_dir_all(&align_tmp)?;
    let mut hmms = BTreeMap::new();
    for set in MARKER_SETS {
        let set_hmms = marker_hmms(&root, set)?;
        write_marker_info(
            &align_tmp.join(format!("{}.{set}.marker_info.tsv", options.prefix)),
            &set_hmms,
        )?;
        hmms.insert(set, set_hmms);
    }
    let mut domains = BTreeMap::new();
    for (index, id, copies) in &searched {
        let [ar53, bac120] = copies;
        let (set, questionable) =
            genome_domain(bac120.single_copy_count(), ar53.single_copy_count());
        let mut warnings = Vec::new();
        if let Some((bacterial, archaeal)) = questionable {
            warnings.push(format!(
                "Genome domain questionable ( {}% Bacterial, {}% Archaeal)",
                summary::py_float(bacterial),
                summary::py_float(archaeal)
            ));
        }
        let domain_copies = if set == "bac120" { bac120 } else { ar53 };
        warnings.extend(multihit_warning(domain_copies, set));
        domains.insert(id.clone(), (*index, set, domain_copies, warnings));
    }

    log.info(format!(
        "Aligning markers in {} genomes with {} CPUs.",
        thousands(searched.len()),
        options.threads
    ));
    let mut outcomes: BTreeMap<String, MarkerOutcome> = BTreeMap::new();
    let mut aligned: BTreeMap<String, (String, f64, &'static str)> = BTreeMap::new();
    for set in MARKER_SETS {
        let genomes = domains
            .iter()
            .filter(|(_, (_, domain, _, _))| *domain == set)
            .map(|(id, (_, _, copies, _))| (id.clone(), (*copies).clone()))
            .collect::<BTreeMap<_, _>>();
        if genomes.is_empty() {
            continue;
        }
        let domain_name = domain_adjective(set);
        log.info(format!(
            "Processing {} genomes identified as {domain_name}.",
            thousands(genomes.len())
        ));
        let user = align_markers(&hmmalign, &hmms[set], &genomes, &align_tmp.join("hmmalign"))?;
        fs::remove_dir_all(align_tmp.join("hmmalign"))?;
        for id in genomes.keys().filter(|id| !user.contains_key(*id)) {
            outcomes.insert(id.clone(), MarkerOutcome::NoMarkers);
        }
        if user.is_empty() {
            log.warn(format!("Identified 0 single copy {domain_name} hits."));
            continue;
        }
        let mask = read_mask(&root, set)?;
        let unmasked_length = user.values().next().map_or(0, String::len);
        let (mut kept, mut pruned) = (BTreeMap::new(), BTreeMap::new());
        for (id, alignment) in user {
            let (masked, valid, msa_percent) = apply_mask(&alignment, mask.as_deref())?;
            if valid < options.min_perc_aa / 100.0 {
                pruned.insert(id, msa_percent);
            } else {
                kept.insert(id, (masked, msa_percent));
            }
        }
        log.info(match &mask {
            Some(mask) => format!(
                "Masked {domain_name} alignment from {} to {} AAs.",
                thousands(unmasked_length),
                thousands(mask.iter().filter(|keep| **keep).count())
            ),
            None => format!(
                "Kept all {} columns of the {domain_name} alignment (the database has no mask).",
                thousands(unmasked_length)
            ),
        });
        if options.min_perc_aa > 0.0 {
            log.info(format!(
                "{} {domain_name} user genomes have amino acids in <{:.1}% of columns in filtered MSA.",
                thousands(pruned.len()),
                options.min_perc_aa
            ));
        }
        let filtered_path = align_dir.join(format!("{}.{set}.filtered.tsv", options.prefix));
        write_filtered(&filtered_path, &pruned)?;
        record.output(ALIGN, set, &filtered_path);

        if !kept.is_empty() {
            log.info(format!(
                "Creating concatenated alignment for {} {domain_name} user genomes.",
                thousands(kept.len())
            ));
            let msa_path = align_dir.join(format!("{}.{set}.user_msa.fasta.gz", options.prefix));
            write_msa_gz(&msa_path, kept.iter().map(|(id, (masked, _))| (id, masked)))?;
            record.output(ALIGN, set, &msa_path);
        } else {
            log.info(format!(
                "All {domain_name} user genomes have been filtered out."
            ));
        }
        for (id, msa_percent) in pruned {
            outcomes.insert(id, MarkerOutcome::Filtered(filtered_reason(msa_percent)));
        }
        for (id, (masked, msa_percent)) in kept {
            aligned.insert(id, (masked, msa_percent, set));
        }
    }
    record.mark(ALIGN);
    let closest = search_marker_index(&aligned, db, warmup, options, out_dir, log)?;
    if !aligned.is_empty() {
        record.output(
            CLASSIFY,
            "all",
            &marker_summary_path(out_dir, &options.prefix),
        );
    }
    for (id, (_, msa_percent, _)) in aligned {
        outcomes.insert(
            id.clone(),
            MarkerOutcome::Aligned {
                msa_percent,
                closest: closest.get(&id).cloned(),
            },
        );
    }
    let no_markers = outcomes
        .iter()
        .filter(|(_, outcome)| matches!(outcome, MarkerOutcome::NoMarkers))
        .map(|(id, _)| id.as_str())
        .collect::<Vec<_>>();
    if !no_markers.is_empty() {
        let mut output = BufWriter::new(File::create(
            align_tmp.join(format!("{}.align.failed.tsv", options.prefix)),
        )?);
        for id in no_markers {
            writeln!(output, "{id}\t{NO_MARKERS}")?;
        }
        output.flush()?;
    }
    for (id, (index, set, _, warnings)) in domains {
        let outcome = outcomes.remove(&id).unwrap_or(MarkerOutcome::NoMarkers);
        if let Stage::Called { markers, .. } = &mut results[index].stage {
            *markers = Some(MarkerResult {
                marker_set: set,
                outcome,
                warnings,
            });
        }
    }
    Ok(())
}

/// Search the marker index with each aligned genome's marker residues, and
/// write the hits of the genome's own domain at or above
/// `--min-marker-identity` to `classify/markers/{prefix}.marker_summary.tsv`.
/// Returns the closest such reference per genome.
fn search_marker_index(
    aligned: &BTreeMap<String, (String, f64, &'static str)>,
    db: &ClassifyDb,
    warmup: &Warmup,
    options: &ClassifyOptions,
    out_dir: &Path,
    log: &ClassifyLog,
) -> Result<BTreeMap<String, MarkerHit>, Box<dyn Error>> {
    if aligned.is_empty() {
        return Ok(BTreeMap::new());
    }
    log.info(format!(
        "Searching the marker index for {} genome(s).",
        thousands(aligned.len())
    ));
    let summary_path = marker_summary_path(out_dir, &options.prefix);
    let sequence_dir = summary_path.with_file_name("intermediate_results");
    fs::create_dir_all(&sequence_dir)?;
    let mut paths = Vec::with_capacity(aligned.len());
    for (id, (masked, _, _)) in aligned {
        let path = sequence_dir.join(format!("{id}.faa"));
        fs::write(&path, format!(">{id}\n{}\n", marker_sequence(masked)))?;
        paths.push(path.to_string_lossy().into_owned());
    }
    let vectors = sketch_from_params(&paths, &db.markers.params);
    let neighbors = search_index(&db.markers, warmup, &vectors, options, log)?;
    let mut output = BufWriter::new(File::create(&summary_path)?);
    writeln!(
        output,
        "user_genome\treference_genome\testimated_marker_identity\treference_taxonomy"
    )?;
    let mut closest = BTreeMap::new();
    for ((id, (_, _, marker_set)), hits) in aligned.iter().zip(neighbors) {
        for (reference, distance) in hits {
            let accession = &db.markers.accessions[reference as usize];
            let identity = sketch_ani(distance, db.markers.params.kmer_size);
            // Unrelated marker sets share k-mers only by chance, and the
            // other domain's markers are different genes.
            let Some(taxonomy) = db.data.taxonomy_of(accession).filter(|taxonomy| {
                marker_set_for_taxonomy(taxonomy) == *marker_set
                    && identity >= options.min_marker_identity
            }) else {
                continue;
            };
            writeln!(
                output,
                "{id}\t{accession}\t{}\t{}",
                py_round(identity, 2),
                taxonomy.join(";")
            )?;
            closest.entry(id.clone()).or_insert_with(|| MarkerHit {
                accession: accession.clone(),
                identity,
            });
        }
    }
    output.flush()?;
    Ok(closest)
}

fn marker_summary_path(out_dir: &Path, prefix: &str) -> PathBuf {
    out_dir
        .join("classify")
        .join("markers")
        .join(format!("{prefix}.marker_summary.tsv"))
}

fn domain_adjective(marker_set: &str) -> &'static str {
    if marker_set == "ar53" {
        "archaeal"
    } else {
        "bacterial"
    }
}

/// Python's `{:,}` formatting.
fn thousands(value: usize) -> String {
    let digits = value.to_string();
    let mut out = String::new();
    for (index, digit) in digits.chars().enumerate() {
        if index > 0 && (digits.len() - index) % 3 == 0 {
            out.push(',');
        }
        out.push(digit);
    }
    out
}

/// `classify/aai/{prefix}.aai_summary.tsv`: every AAI hit per genome.
fn write_aai_summary(
    results: &[QueryResult],
    db: &ClassifyDb,
    out_dir: &Path,
    prefix: &str,
) -> Result<PathBuf, Box<dyn Error>> {
    let dir = out_dir.join("classify").join("aai");
    fs::create_dir_all(&dir)?;
    let path = dir.join(format!("{prefix}.aai_summary.tsv"));
    let mut output = BufWriter::new(File::create(&path)?);
    writeln!(
        output,
        "user_genome\treference_genome\taai\taai_sd\treciprocal_best_hits\tsmaller_proteome\treference_taxonomy"
    )?;
    let mut sorted = results.iter().collect::<Vec<_>>();
    sorted.sort_by(|left, right| left.genome.id.cmp(&right.genome.id));
    for result in sorted {
        let Stage::Called { aai: Some(aai), .. } = &result.stage else {
            continue;
        };
        for hit in &aai.hits {
            let accession = &db.nt.accessions[hit.reference as usize];
            writeln!(
                output,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                result.genome.id,
                accession,
                py_round(hit.aai, 2),
                py_round(hit.sd, 2),
                hit.reciprocal_hits,
                hit.min_proteins,
                db.data
                    .taxonomy_of(accession)
                    .map_or_else(|| NONE_VALUE.to_owned(), |taxonomy| taxonomy.join(";"))
            )?;
        }
    }
    output.flush()?;
    Ok(path)
}

/// GTDB-Tk's `identify/{prefix}.translation_table_summary.tsv` and
/// `identify/{prefix}.failed_genomes.tsv` for genomes that went through gene calling.
fn write_identify_files(
    results: &[QueryResult],
    out_dir: &Path,
    prefix: &str,
    record: &mut RunRecord,
) -> Result<(), Box<dyn Error>> {
    let mut tables = BTreeMap::new();
    let mut failed = BTreeMap::new();
    for result in results {
        let id = result.genome.id.as_str();
        match result.stage {
            Stage::Called { .. } => {
                tables.insert(id, TRANSLATION_TABLE);
            }
            Stage::NoGenes => {
                failed.insert(id, NO_GENES);
            }
            Stage::Species(_) | Stage::Unresolved => {}
        }
    }
    if tables.is_empty() && failed.is_empty() {
        return Ok(());
    }
    let dir = out_dir.join("identify");
    fs::create_dir_all(&dir)?;
    let failed_path = dir.join(format!("{prefix}.failed_genomes.tsv"));
    let tables_path = dir.join(format!("{prefix}.translation_table_summary.tsv"));
    let mut output = BufWriter::new(File::create(&tables_path)?);
    for (id, table) in tables {
        writeln!(output, "{id}\t{table}")?;
    }
    output.flush()?;
    let mut output = BufWriter::new(File::create(&failed_path)?);
    for (id, reason) in failed {
        writeln!(output, "{id}\t{reason}")?;
    }
    output.flush()?;
    record.output(IDENTIFY, "all", &failed_path);
    record.output(IDENTIFY, "all", &tables_path);
    Ok(())
}

/// GTDB-Tk symlinks `{out_dir}/{prefix}.{ms}.summary.tsv` to the copy in `classify/`.
fn link_summary(out_dir: &Path, name: &str) -> std::io::Result<()> {
    let link = out_dir.join(name);
    match fs::remove_file(&link) {
        Ok(()) => {}
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => {}
        Err(error) => return Err(error),
    }
    #[cfg(unix)]
    {
        std::os::unix::fs::symlink(Path::new("classify").join(name), &link)
    }
    #[cfg(not(unix))]
    {
        fs::copy(out_dir.join("classify").join(name), &link).map(|_| ())
    }
}

fn related(hits: &[AniHit], db: &ClassifyDb, min_af: f64) -> Vec<RelatedReference> {
    ani::reportable(hits)
        .filter(|hit| hit.passes_af(min_af))
        .map(|hit| RelatedReference {
            genome_id: db.nt.accessions[hit.reference as usize].clone(),
            ani: hit.ani,
            af: hit.af,
        })
        .collect()
}

/// Build the summary row and its marker set for one ANI-tier result,
/// following GTDB-Tk 2.8's field conventions.
fn summary_row(
    result: &QueryResult,
    db: &ClassifyDb,
    options: &ClassifyOptions,
) -> (&'static str, SummaryRow) {
    let mut row = SummaryRow::new(&result.genome.id);
    let related = related(&result.hits, db, options.min_af);
    let taxonomy_of = |hit: &AniHit| {
        db.data
            .taxonomy_of(&db.nt.accessions[hit.reference as usize])
            .map(<[String]>::to_vec)
    };
    let best_taxonomy = result.hits.first().and_then(taxonomy_of);
    let marker_set = best_taxonomy
        .as_deref()
        .map_or("bac120", marker_set_for_taxonomy);

    match &result.stage {
        Stage::Species(hit) => {
            let accession = &db.nt.accessions[hit.reference as usize];
            let taxonomy = taxonomy_of(hit).unwrap_or_default().join(";");
            row.classification = Some(standardise_taxonomy(&taxonomy));
            row.closest_genome_reference = Some(accession.clone());
            row.closest_genome_reference_radius = db.data.radius_of(accession);
            row.closest_genome_taxonomy = Some(taxonomy.clone());
            row.closest_genome_ani = Some(hit.ani);
            row.closest_genome_af = hit.af;
            row.classification_method = Some("ani_screen".into());
            row.note = Some("classification based on ANI only".into());
            row.other_related_references =
                format_related_references(&related, &[accession.as_str()], &db.data);
            let marker_set = taxonomy_of(hit)
                .as_deref()
                .map_or(marker_set, marker_set_for_taxonomy);
            (marker_set, row)
        }
        stage => {
            row.other_related_references = format_related_references(&related, &[], &db.data);
            row.warnings = radius_warnings(&result.hits, db, options.min_af);
            match stage {
                // GTDB-Tk reports genomes without genes in the bac120 summary.
                Stage::NoGenes => {
                    row.classification = Some("Unclassified".into());
                    row.warnings.push(NO_GENES.into());
                    ("bac120", row)
                }
                Stage::Called {
                    markers: Some(markers),
                    ..
                } => marker_summary_row(row, markers, db),
                Stage::Called { aai: Some(aai), .. } => aai_summary_row(row, aai, marker_set, db),
                _ => {
                    row.classification = Some(unclassified(marker_set, best_taxonomy.is_some()));
                    (marker_set, row)
                }
            }
        }
    }
}

fn unclassified(marker_set: &str, domain_known: bool) -> String {
    match marker_set {
        "ar53" if domain_known => "Unclassified Archaea".into(),
        _ if domain_known => "Unclassified Bacteria".into(),
        _ => "Unclassified".into(),
    }
}

/// Fill in the AAI-tier fields of a row the ANI tier could not resolve.
fn aai_summary_row(
    mut row: SummaryRow,
    aai: &AaiResult,
    ani_marker_set: &'static str,
    db: &ClassifyDb,
) -> (&'static str, SummaryRow) {
    row.translation_table = Some(TRANSLATION_TABLE);
    let best_taxonomy = aai.hits.first().and_then(|hit| {
        db.data
            .taxonomy_of(&db.nt.accessions[hit.reference as usize])
    });
    let marker_set = best_taxonomy.map_or(ani_marker_set, marker_set_for_taxonomy);
    let (call, rank) = match &aai.decision {
        AaiDecision::Genus(call) => (call, "genus"),
        AaiDecision::Family(call) => (call, "family"),
        AaiDecision::Unresolved => {
            row.classification = Some(unclassified(marker_set, best_taxonomy.is_some()));
            return (marker_set, row);
        }
    };
    row.classification = Some(standardise_taxonomy(&call.taxonomy));
    row.classification_method = Some("taxonomic classification defined by AAI".into());
    let lead = call.lead.map_or_else(
        || format!("no other {rank} among the AAI hits"),
        |lead| format!("{} points above the next {rank}", py_round(lead, 2)),
    );
    row.note = Some(format!(
        "{rank} of the closest reference by AAI ({}%, {lead})",
        py_round(call.aai, 2)
    ));
    (marker_set, row)
}

/// Fill in the marker-tier fields: GTDB-Tk's rows for genomes without
/// markers or filtered from the MSA, else the closest marker reference's
/// lineage down to order.
fn marker_summary_row(
    mut row: SummaryRow,
    markers: &MarkerResult,
    db: &ClassifyDb,
) -> (&'static str, SummaryRow) {
    match &markers.outcome {
        MarkerOutcome::NoMarkers => {
            row.classification = Some("Unclassified".into());
            row.warnings.push(NO_MARKERS.into());
            ("bac120", row)
        }
        MarkerOutcome::Filtered(reason) => {
            row.classification = Some(unclassified(markers.marker_set, true));
            row.warnings.push(reason.clone());
            (markers.marker_set, row)
        }
        MarkerOutcome::Aligned {
            msa_percent,
            closest,
        } => {
            row.translation_table = Some(TRANSLATION_TABLE);
            row.msa_percent = Some(*msa_percent);
            row.warnings.extend(markers.warnings.iter().cloned());
            let Some((hit, taxonomy)) = closest.as_ref().and_then(|hit| {
                db.data
                    .taxonomy_of(&hit.accession)
                    .map(|taxonomy| (hit, taxonomy))
            }) else {
                row.classification = Some(unclassified(markers.marker_set, true));
                return (markers.marker_set, row);
            };
            row.classification = Some(standardise_taxonomy(&taxonomy[..4].join(";")));
            row.closest_placement_reference = Some(hit.accession.clone());
            row.closest_placement_radius = db.data.radius_of(&hit.accession);
            row.closest_placement_taxonomy = Some(taxonomy.join(";"));
            row.classification_method =
                Some("taxonomic classification defined by marker gene search".into());
            row.note = Some(format!(
                "order of the closest marker gene match (estimated marker identity {}%)",
                py_round(hit.identity, 2)
            ));
            (markers.marker_set, row)
        }
    }
}

/// GTDB-Tk's warnings for near misses that did not receive a species call.
fn radius_warnings(hits: &[AniHit], db: &ClassifyDb, min_af: f64) -> Vec<String> {
    let mut warnings = Vec::new();
    if let Some(best) = hits.iter().find(|hit| hit.passes_af(min_af)) {
        let accession = &db.nt.accessions[best.reference as usize];
        if best.ani >= MIN_SPECIES_RADIUS {
            let radius = db
                .data
                .radius_of(accession)
                .map_or_else(|| "None".to_owned(), summary::py_float);
            warnings.push(format!(
                "Genome not assigned to closest species as it falls outside its pre-defined ANI radius (closest {accession} ANI={}%, AF={}, radius={radius}%)",
                py_round(best.ani, 2),
                best.af.map_or_else(|| NONE_VALUE.to_owned(), |af| py_round(af, 3))
            ));
        }
    } else if let Some(top) = hits.first() {
        if top.ani >= MIN_SPECIES_RADIUS {
            warnings.push(format!(
                "Closest reference {} has ANI={}% but alignment fraction AF={} below the {} threshold; no species assignment",
                db.nt.accessions[top.reference as usize],
                py_round(top.ani, 2),
                top.af.map_or_else(|| NONE_VALUE.to_owned(), |af| py_round(af, 3)),
                summary::py_float(min_af)
            ));
        }
    }
    warnings
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn genome_dir_ids_drop_the_extension_and_batchfiles_take_ids() {
        let dir = std::env::temp_dir().join(format!("aeon-queries-{}", std::process::id()));
        fs::create_dir_all(&dir).unwrap();
        for name in ["b.fna.gz", "a.fna.gz", "ignored.faa"] {
            fs::write(dir.join(name), ">x\nACGT\n").unwrap();
        }
        let genomes = resolve_queries(&QueryInput::GenomeDir {
            dir: dir.to_string_lossy().into_owned(),
            extension: "fna.gz".into(),
        })
        .unwrap();
        assert_eq!(
            genomes.iter().map(|g| g.id.as_str()).collect::<Vec<_>>(),
            ["a", "b"]
        );

        let batch = dir.join("batch.tsv");
        fs::write(
            &batch,
            format!(
                "{}\tsample_1\n{}\n",
                dir.join("a.fna.gz").display(),
                dir.join("b.fna.gz").display()
            ),
        )
        .unwrap();
        let genomes =
            resolve_queries(&QueryInput::Batchfile(batch.to_string_lossy().into_owned())).unwrap();
        fs::remove_dir_all(&dir).unwrap();
        assert_eq!(
            genomes.iter().map(|g| g.id.as_str()).collect::<Vec<_>>(),
            ["sample_1", "b"]
        );
    }
}
