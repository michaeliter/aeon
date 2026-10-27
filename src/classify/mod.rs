//! GTDB-Tk-compatible genome classification: ANI screen, AAI, then marker placement.

pub(crate) mod aai;
pub(crate) mod ani;
pub(crate) mod external;
pub(crate) mod genes;
pub(crate) mod log;
pub(crate) mod markers;
pub(crate) mod placement;
pub(crate) mod refdata;
pub(crate) mod refdb;
pub(crate) mod split;
pub(crate) mod summary;

use rayon::prelude::*;
use std::cmp::Reverse;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::error::Error;
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

use crate::classify::aai::{AaiDecision, AaiHit, compute_aai};
use crate::classify::ani::{AniDecision, AniHit, MIN_SPECIES_RADIUS, decide, refine_candidates};
use crate::classify::external::find_program;
use crate::classify::genes::{GeneCalls, call_genes};
use crate::classify::log::ClassifyLog;
use crate::classify::markers::{
    CopyNumber, MARKER_SETS, PfamInfo, align_markers, apply_mask, copy_number, filtered_reason,
    genome_domain, marker_hmms, multihit_warning, read_mask, read_msa, read_proteins,
    read_reference_msa, search_markers, write_filtered, write_marker_info, write_markers_summary,
    write_msa_gz,
};
use crate::classify::placement::{
    Placement, classify_tree, program_version, run_pplacer, write_pplacer_classification,
    write_red_dictionary,
};
use crate::classify::refdata::{marker_set_for_taxonomy, standardise_taxonomy, strip_db_prefix};
use crate::classify::refdb::ClassifyDb;
use crate::classify::split::{
    backbone_placement, class_tree, classify_backbone, phyla_of_classes, read_tree_mapping,
    reconcile, write_backbone_classification, write_tree_mapping,
};
use crate::classify::summary::{
    AniScreenRow, RelatedReference, SummaryRow, format_related_references, py_round,
    write_ani_summary, write_summary,
};
use crate::database::search_vectors;
use crate::sketch::sketch_from_params;
use crate::utils::read_list_file;

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
    pub(crate) genus_ani: f64,
    pub(crate) genus_aai: f64,
    pub(crate) family_aai: f64,
    pub(crate) aai_top: usize,
    pub(crate) pyrodigal: String,
    pub(crate) diamond: String,
    pub(crate) min_perc_aa: f64,
    pub(crate) bin_dir: Option<String>,
    pub(crate) full_tree: bool,
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
    let mut name = Path::new(path)
        .file_name()
        .and_then(|name| name.to_str())
        .unwrap_or(path);
    for suffix in [".gz", ".fna", ".fa", ".fasta", ".fas"] {
        name = name.strip_suffix(suffix).unwrap_or(name);
    }
    name.to_owned()
}

/// Everything known about one query after each tier it reached.
struct QueryResult<'a> {
    genome: &'a QueryGenome,
    hits: Vec<AniHit>,
    decision: AniDecision,
    /// Gene calls for genomes the ANI tier left unresolved; the inner
    /// `None` means pyrodigal called no genes.
    genes: Option<Option<GeneCalls>>,
    aai: Option<AaiResult>,
    markers: Option<MarkerResult>,
}

struct AaiResult {
    hits: Vec<AaiHit>,
    decision: AaiDecision,
}

/// Where the marker pipeline left a genome.
#[derive(Debug)]
enum MarkerOutcome {
    /// In the masked user MSA, ready for placement.
    Aligned { msa_percent: f64 },
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
    /// Set once pplacer has placed an aligned genome.
    placement: Option<Placement>,
}

impl QueryResult<'_> {
    fn called_genes(&self) -> Option<&GeneCalls> {
        self.genes.as_ref().and_then(Option::as_ref)
    }

    /// Unresolved by ANI and, where it ran, by AAI.
    fn needs_markers(&self) -> bool {
        self.decision == AniDecision::Unresolved
            && self.called_genes().is_some()
            && self
                .aai
                .as_ref()
                .is_none_or(|aai| aai.decision == AaiDecision::Unresolved)
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
        "Using GTDB {} reference database {} ({} representative genomes)",
        db.manifest.gtdb_release,
        options.db,
        db.nt_paths.len()
    ));
    log.info(format!("Identified {} genomes to classify", queries.len()));

    log.info("Running ANI screen with DiskANN candidate search and TurboANI");
    let query_paths = queries
        .iter()
        .map(|genome| genome.path.clone())
        .collect::<Vec<_>>();
    let vectors = sketch_from_params(&query_paths, &db.nt_params);
    let candidates = search_vectors(
        &db.nt_prefix(),
        db.nt_paths.len(),
        &vectors,
        options.candidates,
        options.beam_width,
    )?
    .into_iter()
    .map(|neighbors| neighbors.into_iter().map(|(id, _)| id).collect::<Vec<_>>())
    .collect::<Vec<_>>();
    let hits = refine_candidates(&query_paths, &candidates, &db.nt_paths)?;
    let results = queries
        .iter()
        .zip(hits)
        .map(|(genome, hits)| QueryResult {
            decision: decide(
                &hits,
                &db.nt_accessions,
                &db.data,
                options.min_af,
                options.genus_ani,
            ),
            genome,
            hits,
            genes: None,
            aai: None,
            markers: None,
        })
        .collect::<Vec<_>>();
    for result in &results {
        ::log::debug!(
            "ani genome={} decision={:?} hits={}",
            result.genome.id,
            result.decision,
            result
                .hits
                .iter()
                .map(|hit| format!(
                    "{}:{:.2}/{:.3}",
                    db.nt_accessions[hit.reference as usize], hit.ani, hit.af
                ))
                .collect::<Vec<_>>()
                .join(",")
        );
    }
    let species = results
        .iter()
        .filter(|result| matches!(result.decision, AniDecision::Species { .. }))
        .count();
    log.info(format!(
        "{species} genome(s) have been classified using the ANI pre-screening step."
    ));
    let mut results = results;
    call_query_genes(&mut results, options, &out_dir, &log)?;
    run_aai_tier(&mut results, &db, options, &out_dir, &log)?;
    run_marker_tier(&mut results, &db, options, &out_dir, &log)?;
    run_placement_tier(&mut results, &db, options, &out_dir, &log)?;

    write_identify_files(&results, &out_dir, &options.prefix)?;
    let mut rows: BTreeMap<&'static str, Vec<SummaryRow>> = BTreeMap::new();
    let mut ani_rows: BTreeMap<&'static str, Vec<AniScreenRow>> = BTreeMap::new();
    for result in &results {
        let (marker_set, row) = summary_row(result, &db, options);
        if row.warnings.iter().any(|warning| !warning.is_empty()) {
            log.warn(format!("{}: {}", row.user_genome, row.warnings.join(";")));
        }
        if let AniDecision::Species { hit } = &result.decision {
            ani_rows.entry(marker_set).or_default().push(AniScreenRow {
                user_genome: result.genome.id.clone(),
                reference_genome: db.nt_accessions[hit.reference as usize].clone(),
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
    }
    for (marker_set, rows) in &rows {
        let name = format!("{}.{marker_set}.summary.tsv", options.prefix);
        write_summary(&out_dir.join("classify").join(&name), rows)?;
        link_summary(&out_dir, &name)?;
        log.info(format!(
            "{} {marker_set} summary written to: {}",
            rows.len(),
            out_dir.join(&name).display()
        ));
    }
    log.info(format!(
        "Done. Classified {} genomes in {:.3}s",
        queries.len(),
        started.elapsed().as_secs_f64()
    ));
    Ok(())
}

fn marker_genes_dir(out_dir: &Path) -> PathBuf {
    out_dir
        .join("identify")
        .join("intermediate_results")
        .join("marker_genes")
}

/// Call genes with pyrodigal for genomes the ANI tier left unresolved.
fn call_query_genes(
    results: &mut [QueryResult],
    options: &ClassifyOptions,
    out_dir: &Path,
    log: &ClassifyLog,
) -> Result<(), Box<dyn Error>> {
    let pending = results
        .iter()
        .enumerate()
        .filter(|(_, result)| result.decision == AniDecision::Unresolved)
        .map(|(index, _)| index)
        .collect::<Vec<_>>();
    if pending.is_empty() {
        return Ok(());
    }
    let pyrodigal = find_program(&options.pyrodigal)?;
    log.info(format!(
        "Calling genes with pyrodigal for {} genome(s) not resolved by ANI",
        pending.len()
    ));
    let marker_dir = marker_genes_dir(out_dir);
    let genes = pending
        .par_iter()
        .map(|&index| {
            let genome = results[index].genome;
            call_genes(&pyrodigal, &genome.path, &genome.id, &marker_dir)
                .map_err(|error| format!("gene calling failed for {}: {error}", genome.id))
        })
        .collect::<Result<Vec<_>, _>>()?;
    for (index, calls) in pending.into_iter().zip(genes) {
        results[index].genes = Some(calls);
    }
    Ok(())
}

/// Classify genomes the ANI tier left unresolved by AAI against candidate
/// reference proteomes.
fn run_aai_tier(
    results: &mut [QueryResult],
    db: &ClassifyDb,
    options: &ClassifyOptions,
    out_dir: &Path,
    log: &ClassifyLog,
) -> Result<(), Box<dyn Error>> {
    let called = results
        .iter()
        .enumerate()
        .filter(|(_, result)| result.decision == AniDecision::Unresolved)
        .filter_map(|(index, result)| {
            result
                .called_genes()
                .map(|calls| (index, calls.protein_faa.clone()))
        })
        .collect::<Vec<_>>();
    if called.is_empty() {
        return Ok(());
    }
    let Some((aa_params, aa_paths)) = &db.aa else {
        log.warn(format!(
            "{} genome(s) need the AAI tier, but the database has no protein index; rebuild it with --call-genes or --protein-list",
            called.len()
        ));
        return Ok(());
    };
    let diamond = find_program(&options.diamond)?;
    log.info(format!(
        "Running AAI search with DiskANN candidate search and DIAMOND for {} genome(s)",
        called.len()
    ));
    let proteome_paths = called
        .iter()
        .map(|(_, faa)| faa.to_string_lossy().into_owned())
        .collect::<Vec<_>>();
    let vectors = sketch_from_params(&proteome_paths, aa_params);
    let candidates = search_vectors(
        &db.aa_prefix(),
        aa_paths.len(),
        &vectors,
        options.candidates,
        options.beam_width,
    )?
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
        aa_paths,
        options.threads,
    )?;

    let mut hits_by_index = called
        .iter()
        .map(|(index, _)| *index)
        .zip(aai_hits)
        .collect::<BTreeMap<_, _>>();
    let mut resolved = 0;
    let indices = called.iter().map(|(index, _)| *index).collect::<Vec<_>>();
    for index in indices {
        let hits = hits_by_index.remove(&index).unwrap_or_default();
        let decision = aai::decide(
            &hits,
            &db.nt_accessions,
            &db.data,
            options.genus_aai,
            options.family_aai,
            options.aai_top,
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
                    db.nt_accessions[hit.reference as usize], hit.aai, hit.reciprocal_hits
                ))
                .collect::<Vec<_>>()
                .join(",")
        );
        results[index].aai = Some(AaiResult { hits, decision });
    }
    log.info(format!(
        "{resolved} genome(s) have been classified using AAI."
    ));
    write_aai_summary(results, db, out_dir, &options.prefix)
}

/// An external program from `--bin-dir`, or from PATH.
fn tool(options: &ClassifyOptions, program: &str) -> Result<PathBuf, Box<dyn Error>> {
    match &options.bin_dir {
        Some(dir) => find_program(&Path::new(dir).join(program).to_string_lossy()),
        None => find_program(program),
    }
}

/// GTDB-Tk's `classify` step for aligned genomes: pplacer into the reference
/// tree, then RED and topology. Bacteria go through the backbone and
/// class-level trees unless `--full_tree` is set.
fn run_placement_tier(
    results: &mut [QueryResult],
    db: &ClassifyDb,
    options: &ClassifyOptions,
    out_dir: &Path,
    log: &ClassifyLog,
) -> Result<(), Box<dyn Error>> {
    let Some(root) = db.data.root.clone() else {
        return Ok(());
    };
    for set in MARKER_SETS {
        let pending = results
            .iter()
            .enumerate()
            .filter(|(_, result)| {
                result.markers.as_ref().is_some_and(|markers| {
                    markers.marker_set == set
                        && matches!(markers.outcome, MarkerOutcome::Aligned { .. })
                })
            })
            .map(|(index, result)| (result.genome.id.clone(), index))
            .collect::<BTreeMap<_, _>>();
        if pending.is_empty() {
            continue;
        }
        let placer = Placer {
            pplacer: tool(options, "pplacer")?,
            guppy: tool(options, "guppy")?,
            root: &root,
            db,
            options,
            out_dir,
            log,
            red_dict: if set == "ar53" {
                &db.data.red_dist_arc
            } else {
                &db.data.red_dist_bac
            },
        };
        fs::create_dir_all(placer.pplacer_dir())?;
        write_red_dictionary(
            &placer
                .intermediate()
                .join(format!("{}.{set}.red_dictionary.tsv", options.prefix)),
            placer.red_dict,
        )?;
        let users = pending.keys().cloned().collect::<HashSet<_>>();
        let placements = if set == "bac120" && !options.full_tree {
            placer.place_split(&users)?
        } else {
            placer.place_full_tree(set, &users)?
        };
        for (genome, placement) in placements {
            let index = pending[&genome];
            add_placement_ani(&mut results[index], &placement, db)?;
            if let Some(markers) = results[index].markers.as_mut() {
                markers.placement = Some(placement);
            }
        }
    }
    Ok(())
}

/// Runs pplacer for one marker set and classifies its placements.
struct Placer<'a> {
    pplacer: PathBuf,
    guppy: PathBuf,
    root: &'a Path,
    db: &'a ClassifyDb,
    options: &'a ClassifyOptions,
    out_dir: &'a Path,
    log: &'a ClassifyLog,
    red_dict: &'a HashMap<String, f64>,
}

impl Placer<'_> {
    fn classify_dir(&self) -> PathBuf {
        self.out_dir.join("classify")
    }

    fn intermediate(&self) -> PathBuf {
        self.classify_dir().join("intermediate_results")
    }

    fn pplacer_dir(&self) -> PathBuf {
        self.intermediate().join("pplacer")
    }

    fn user_msa(&self, set: &str) -> PathBuf {
        self.out_dir
            .join("align")
            .join(format!("{}.{set}.user_msa.fasta.gz", self.options.prefix))
    }

    fn threads(&self) -> usize {
        self.options.threads.clamp(1, 64)
    }

    fn run(
        &self,
        refpkg: &Path,
        msa: &Path,
        json_out: &Path,
        pplacer_out: &Path,
        tree_out: &Path,
    ) -> Result<(), Box<dyn Error>> {
        run_pplacer(
            &self.pplacer,
            &self.guppy,
            refpkg,
            msa,
            json_out,
            pplacer_out,
            tree_out,
            self.options.threads,
        )?;
        self.log
            .info("Calculating RED values based on reference tree.");
        Ok(())
    }

    fn warn_unplaced(&self, expected: usize, placed: usize) {
        if expected > placed {
            self.log.warn(format!(
                "{} out of {expected} have not been properly placed by pplacer. This is a known issue with pplacer but we do not have a solution currently.",
                expected - placed
            ));
        }
    }

    /// The whole reference tree of `set` (GTDB-Tk's full-tree mode).
    fn place_full_tree(
        &self,
        set: &str,
        users: &HashSet<String>,
    ) -> Result<BTreeMap<String, Placement>, Box<dyn Error>> {
        let prefix = &self.options.prefix;
        self.log.info(format!(
            "Placing {} {} genomes into reference tree with pplacer using {} CPUs (be patient).",
            thousands(users.len()),
            domain_adjective(set),
            self.threads()
        ));
        self.log.info(format!(
            "pplacer version: {}",
            program_version(&self.pplacer)
        ));
        let tree_path = self
            .classify_dir()
            .join(format!("{prefix}.{set}.classify.tree"));
        self.run(
            &self
                .root
                .join("pplacer")
                .join(format!("gtdb_{}_{set}.refpkg", self.db.data.version)),
            &self.user_msa(set),
            &self.pplacer_dir().join(format!("pplacer.{set}.json")),
            &self.pplacer_dir().join(format!("pplacer.{set}.out")),
            &tree_path,
        )?;
        let placements = classify_tree(
            &tree_path,
            &self
                .root
                .join("mrca_red")
                .join(format!("gtdbtk_{}_{set}.tsv", self.db.data.version)),
            set,
            users,
            &self.db.data,
            self.red_dict,
        )?;
        write_pplacer_classification(
            &self
                .intermediate()
                .join(format!("{prefix}.{set}.classification_pplacer.tsv")),
            &placements,
        )?;
        self.log
            .info("Traversing tree to determine classification method.");
        self.warn_unplaced(users.len(), placements.len());
        Ok(placements)
    }

    /// GTDB-Tk's split-tree mode for bacteria: the backbone tree assigns a
    /// class, the class picks a class-level tree, and rules 1-6 reconcile
    /// the two placements.
    fn place_split(
        &self,
        users: &HashSet<String>,
    ) -> Result<BTreeMap<String, Placement>, Box<dyn Error>> {
        let prefix = &self.options.prefix;
        let split = self.root.join("split");
        let class_level = split.join("class_level");
        self.log.info(format!(
            "Placing {} bacterial genomes into backbone reference tree with pplacer using {} CPUs (be patient).",
            thousands(users.len()),
            self.threads()
        ));
        self.log.info(format!(
            "pplacer version: {}",
            program_version(&self.pplacer)
        ));
        let backbone_tree = self
            .classify_dir()
            .join(format!("{prefix}.backbone.bac120.classify.tree"));
        self.run(
            &split
                .join("backbone")
                .join("pplacer")
                .join("gtdbtk_package_backbone.refpkg"),
            &self.user_msa("bac120"),
            &self.pplacer_dir().join("pplacer.backbone.bac120.json"),
            &self.pplacer_dir().join("pplacer.backbone.bac120.out"),
            &backbone_tree,
        )?;
        let backbone = classify_backbone(
            &backbone_tree,
            &split
                .join("backbone")
                .join("red")
                .join("backbone_red_value.tsv"),
            users,
            &self.db.data,
            self.red_dict,
        )?;
        write_backbone_classification(
            &self.intermediate().join(format!(
                "{prefix}.bac120.backbone.classification_pplacer.tsv"
            )),
            &backbone,
        )?;
        self.warn_unplaced(users.len(), backbone.len());

        let mapping = read_tree_mapping(&class_level.join("tree_mapping.tsv"))?;
        let mut placements = BTreeMap::new();
        let mut rules = BTreeMap::new();
        let mut trees: Vec<(&str, Vec<&str>)> = Vec::new();
        for (genome, call) in &backbone {
            match class_tree(call, &mapping)? {
                Some(tree) => match trees.iter_mut().find(|(id, _)| *id == tree) {
                    Some((_, genomes)) => genomes.push(genome),
                    None => trees.push((tree, vec![genome])),
                },
                None => {
                    placements.insert(genome.clone(), backbone_placement(call));
                    rules.insert(genome.clone(), ("backbone".to_owned(), "Rule 1"));
                }
            }
        }
        self.log.info(format!(
            "{} out of {} have an class assignments. Those genomes will be reclassified.",
            trees
                .iter()
                .map(|(_, genomes)| genomes.len())
                .sum::<usize>(),
            users.len()
        ));
        trees.sort_by_key(|(_, genomes)| Reverse(genomes.len()));

        let msa = read_msa(&self.user_msa("bac120"))?;
        let calls = backbone
            .iter()
            .map(|(genome, call)| (genome.as_str(), call))
            .collect::<HashMap<_, _>>();
        for (index, (tree_id, genomes)) in trees.iter().enumerate() {
            let tree_dir = self.pplacer_dir().join(format!("tree_{tree_id}"));
            fs::create_dir_all(&tree_dir)?;
            let sub_msa = tree_dir.join("user_msa_file.fasta");
            let mut output = BufWriter::new(File::create(&sub_msa)?);
            for genome in genomes {
                writeln!(
                    output,
                    ">{genome}\n{}",
                    msa.get(*genome).map_or("", String::as_str)
                )?;
            }
            output.flush()?;
            self.log.info(format!(
                "Placing {} bacterial genomes into class-level reference tree {tree_id} ({}/{}) with pplacer using {} CPUs (be patient).",
                thousands(genomes.len()),
                index + 1,
                trees.len(),
                self.threads()
            ));
            let tree_path = self
                .classify_dir()
                .join(format!("{prefix}.bac120.classify.tree.{tree_id}.tree"));
            self.run(
                &class_level
                    .join("pplacer")
                    .join(format!("gtdbtk.package.{tree_id}.refpkg")),
                &sub_msa,
                &tree_dir.join("pplacer.class_level.bac120.json"),
                &tree_dir.join("pplacer.class_level.bac120.out"),
                &tree_path,
            )?;
            let tree_users = genomes.iter().map(|g| (*g).to_owned()).collect();
            let class_placements = classify_tree(
                &tree_path,
                &class_level
                    .join("red")
                    .join(format!("red_value_{tree_id}.tsv")),
                "bac120",
                &tree_users,
                &self.db.data,
                self.red_dict,
            )?;
            write_pplacer_classification(
                &self.intermediate().join(format!(
                    "{prefix}.bac120.class_level.classification_pplacer_tree_{tree_id}.tsv"
                )),
                &class_placements,
            )?;
            self.log
                .info("Traversing tree to determine classification method.");
            self.warn_unplaced(genomes.len(), class_placements.len());
            let valid_classes = mapping
                .iter()
                .filter(|(_, tree)| tree == tree_id)
                .map(|(class, _)| class.clone())
                .collect::<HashSet<_>>();
            let valid_phyla = phyla_of_classes(&valid_classes, &self.db.data);
            for (genome, placement) in class_placements {
                let reconciled = reconcile(
                    calls[genome.as_str()],
                    placement,
                    &valid_classes,
                    &valid_phyla,
                );
                let tree = if reconciled.on_backbone {
                    "backbone".to_owned()
                } else {
                    (*tree_id).to_owned()
                };
                rules.insert(genome.clone(), (tree, reconciled.rule));
                placements.insert(genome, reconciled.placement);
            }
        }
        write_tree_mapping(
            &self
                .classify_dir()
                .join(format!("{prefix}.bac120.tree.mapping.tsv")),
            &rules,
        )?;
        Ok(placements)
    }
}

/// GTDB-Tk compares a placed genome that has ANI hits with every reference,
/// so its placement reference always has an ANI if one exists. aeon's hits
/// come from candidate search, so measure the placement reference directly
/// when it is missing.
fn add_placement_ani(
    result: &mut QueryResult,
    placement: &Placement,
    db: &ClassifyDb,
) -> Result<(), Box<dyn Error>> {
    let Some(reference) = &placement.reference else {
        return Ok(());
    };
    if ani::reportable(&result.hits).next().is_none() {
        return Ok(());
    }
    let accession = strip_db_prefix(reference);
    let Some(index) = db.nt_accessions.iter().position(|id| id == accession) else {
        return Ok(());
    };
    if result
        .hits
        .iter()
        .any(|hit| hit.reference as usize == index)
    {
        return Ok(());
    }
    let hits = refine_candidates(
        std::slice::from_ref(&result.genome.path),
        &[vec![index as u32]],
        &db.nt_paths,
    )?;
    result.hits.extend(hits.into_iter().flatten());
    ani::sort_hits(&mut result.hits);
    Ok(())
}

/// GTDB-Tk's `identify` and `align` steps for genomes neither ANI nor AAI
/// resolved: marker search, copy number, domain call, and the masked
/// concatenated alignment.
fn run_marker_tier(
    results: &mut [QueryResult],
    db: &ClassifyDb,
    options: &ClassifyOptions,
    out_dir: &Path,
    log: &ClassifyLog,
) -> Result<(), Box<dyn Error>> {
    let pending = results
        .iter()
        .enumerate()
        .filter(|(_, result)| result.needs_markers())
        .map(|(index, result)| {
            let faa = result.called_genes().expect("genes").protein_faa.clone();
            (index, result.genome.id.clone(), faa)
        })
        .collect::<Vec<_>>();
    if pending.is_empty() {
        return Ok(());
    }
    let Some(root) = db.data.root.clone() else {
        log.warn(format!(
            "{} genome(s) need the marker tier, but the database was not built from a GTDB-Tk data package (--gtdbtk-data)",
            pending.len()
        ));
        return Ok(());
    };
    let (hmmsearch, hmmalign) = (tool(options, "hmmsearch")?, tool(options, "hmmalign")?);

    log.info(format!(
        "Identifying markers in {} genomes with {} threads.",
        thousands(pending.len()),
        options.threads
    ));
    let pfam_info = PfamInfo::load(&root)?;
    let marker_dir = marker_genes_dir(out_dir);
    let searched = pending
        .par_iter()
        .map(|(index, id, faa)| {
            let (pfam, tigrfam) =
                search_markers(&hmmsearch, &root, &pfam_info, id, faa, &marker_dir.join(id))
                    .map_err(|error| format!("marker search failed for {id}: {error}"))?;
            let proteins = read_proteins(faa).map_err(|error| format!("{id}: {error}"))?;
            let copies = MARKER_SETS.map(|set| copy_number(set, &proteins, &pfam, &tigrfam));
            Ok::<_, String>((*index, id.clone(), copies))
        })
        .collect::<Result<Vec<_>, _>>()?;
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
        write_markers_summary(
            &identify_dir.join(format!("{}.{set}.markers_summary.tsv", options.prefix)),
            genomes,
        )?;
    }
    log.info("Done.");

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
        domains.insert(id.clone(), (*index, set, domain_copies.clone(), warnings));
    }

    log.info(format!(
        "Aligning markers in {} genomes with {} CPUs.",
        thousands(searched.len()),
        options.threads
    ));
    let mut outcomes: BTreeMap<String, MarkerOutcome> = BTreeMap::new();
    for set in MARKER_SETS {
        let genomes = domains
            .iter()
            .filter(|(_, (_, domain, _, _))| *domain == set)
            .map(|(id, (_, _, copies, _))| (id.clone(), copies.clone()))
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
        let mask = read_mask(&db.data, &root, set)?;
        let unmasked_length = user.values().next().map_or(0, String::len);
        let (mut kept, mut pruned) = (BTreeMap::new(), BTreeMap::new());
        for (id, alignment) in user {
            let (masked, valid, msa_percent) = apply_mask(&alignment, &mask)?;
            if valid < options.min_perc_aa / 100.0 {
                pruned.insert(id, masked);
            } else {
                kept.insert(id, (masked, msa_percent));
            }
        }
        log.info(format!(
            "Masked {domain_name} alignment from {} to {} AAs.",
            thousands(unmasked_length),
            thousands(mask.iter().filter(|keep| **keep).count())
        ));
        if options.min_perc_aa > 0.0 {
            log.info(format!(
                "{} {domain_name} user genomes have amino acids in <{:.1}% of columns in filtered MSA.",
                thousands(pruned.len()),
                options.min_perc_aa
            ));
        }
        write_filtered(
            &align_dir.join(format!("{}.{set}.filtered.tsv", options.prefix)),
            &pruned,
        )?;

        let user_msa = kept
            .iter()
            .map(|(id, (masked, _))| (id.clone(), masked.clone()))
            .collect::<BTreeMap<_, _>>();
        let mut msa = user_msa.clone();
        for (id, alignment) in read_reference_msa(&root, &db.data.version, set)? {
            msa.insert(id, apply_mask(&alignment, &mask)?.0);
        }
        log.info(format!(
            "Creating concatenated alignment for {} {domain_name} GTDB and user genomes.",
            thousands(msa.len())
        ));
        write_msa_gz(
            &align_dir.join(format!("{}.{set}.msa.fasta.gz", options.prefix)),
            &msa,
            &db.data,
        )?;
        drop(msa);
        if !user_msa.is_empty() {
            log.info(format!(
                "Creating concatenated alignment for {} {domain_name} user genomes.",
                thousands(user_msa.len())
            ));
            write_msa_gz(
                &align_dir.join(format!("{}.{set}.user_msa.fasta.gz", options.prefix)),
                &user_msa,
                &db.data,
            )?;
        } else {
            log.info(format!(
                "All {domain_name} user genomes have been filtered out."
            ));
        }
        for (id, masked) in &pruned {
            outcomes.insert(id.clone(), MarkerOutcome::Filtered(filtered_reason(masked)));
        }
        for (id, (_, msa_percent)) in kept {
            outcomes.insert(id, MarkerOutcome::Aligned { msa_percent });
        }
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
            writeln!(output, "{id}\tNo bacterial or archaeal marker")?;
        }
        output.flush()?;
    }
    for (id, (index, set, _, warnings)) in domains {
        let outcome = outcomes.remove(&id).unwrap_or(MarkerOutcome::NoMarkers);
        results[index].markers = Some(MarkerResult {
            marker_set: set,
            outcome,
            warnings,
            placement: None,
        });
    }
    Ok(())
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
) -> Result<(), Box<dyn Error>> {
    let dir = out_dir.join("classify").join("aai");
    fs::create_dir_all(&dir)?;
    let mut output = BufWriter::new(File::create(dir.join(format!("{prefix}.aai_summary.tsv")))?);
    writeln!(
        output,
        "user_genome\treference_genome\taai\treciprocal_best_hits\treference_taxonomy"
    )?;
    let mut sorted = results.iter().collect::<Vec<_>>();
    sorted.sort_by(|left, right| left.genome.id.cmp(&right.genome.id));
    for result in sorted {
        let Some(aai) = &result.aai else { continue };
        for hit in &aai.hits {
            let accession = &db.nt_accessions[hit.reference as usize];
            writeln!(
                output,
                "{}\t{}\t{}\t{}\t{}",
                result.genome.id,
                accession,
                py_round(hit.aai, 2),
                hit.reciprocal_hits,
                db.data
                    .taxonomy_of(accession)
                    .map_or_else(|| "N/A".to_owned(), |taxonomy| taxonomy.join(";"))
            )?;
        }
    }
    output.flush()?;
    Ok(())
}

/// GTDB-Tk's `identify/{prefix}.translation_table_summary.tsv` and
/// `identify/{prefix}.failed_genomes.tsv` for genomes that went through gene calling.
fn write_identify_files(
    results: &[QueryResult],
    out_dir: &Path,
    prefix: &str,
) -> Result<(), Box<dyn Error>> {
    let called = results
        .iter()
        .filter_map(|result| result.genes.as_ref().map(|genes| (result.genome, genes)))
        .collect::<Vec<_>>();
    if called.is_empty() {
        return Ok(());
    }
    let dir = out_dir.join("identify");
    fs::create_dir_all(&dir)?;
    let mut tables = BTreeMap::new();
    let mut failed = BTreeMap::new();
    for (genome, genes) in called {
        match genes {
            Some(genes) => {
                tables.insert(genome.id.as_str(), genes.translation_table);
            }
            None => {
                failed.insert(genome.id.as_str(), "No genes were called by Prodigal");
            }
        }
    }
    let mut output = BufWriter::new(File::create(
        dir.join(format!("{prefix}.translation_table_summary.tsv")),
    )?);
    for (id, table) in tables {
        writeln!(output, "{id}\t{table}")?;
    }
    output.flush()?;
    let mut output = BufWriter::new(File::create(
        dir.join(format!("{prefix}.failed_genomes.tsv")),
    )?);
    for (id, reason) in failed {
        writeln!(output, "{id}\t{reason}")?;
    }
    output.flush()?;
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
        .filter(|hit| hit.af >= min_af)
        .map(|hit| RelatedReference {
            genome_id: db.nt_accessions[hit.reference as usize].clone(),
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
            .taxonomy_of(&db.nt_accessions[hit.reference as usize])
            .map(<[String]>::to_vec)
    };
    let best_taxonomy = result.hits.first().and_then(taxonomy_of);
    let marker_set = best_taxonomy
        .as_deref()
        .map_or("bac120", marker_set_for_taxonomy);

    match &result.decision {
        AniDecision::Species { hit } => {
            let accession = &db.nt_accessions[hit.reference as usize];
            let taxonomy = taxonomy_of(hit).unwrap_or_default().join(";");
            row.classification = Some(standardise_taxonomy(&taxonomy, None));
            row.closest_genome_reference = Some(accession.clone());
            row.closest_genome_reference_radius = db.data.radius_of(accession);
            row.closest_genome_taxonomy = Some(taxonomy.clone());
            row.closest_genome_ani = Some(hit.ani);
            row.closest_genome_af = Some(hit.af);
            row.classification_method = Some("ani_screen".into());
            row.note = Some("classification based on ANI only".into());
            row.other_related_references =
                format_related_references(&related, &[accession.as_str()], &db.data);
            let marker_set = taxonomy_of(hit)
                .as_deref()
                .map_or(marker_set, marker_set_for_taxonomy);
            return (marker_set, row);
        }
        AniDecision::Genus {
            genus_taxonomy,
            supporting,
            considered,
        } => {
            row.classification = Some(standardise_taxonomy(genus_taxonomy, None));
            row.classification_method = Some("taxonomic classification defined by ANI".into());
            row.note = Some(format!(
                "genus assigned by ANI consensus ({supporting}/{considered} references >= {}%)",
                options.genus_ani
            ));
        }
        AniDecision::Unresolved => {
            row.other_related_references = format_related_references(&related, &[], &db.data);
            row.warnings = radius_warnings(&result.hits, db, options.min_af);
            match (&result.genes, &result.markers, &result.aai) {
                // GTDB-Tk reports genomes without genes in the bac120 summary.
                (Some(None), _, _) => {
                    row.classification = Some("Unclassified".into());
                    row.warnings.push("No genes were called by Prodigal".into());
                    return ("bac120", row);
                }
                (Some(Some(genes)), Some(markers), _) => {
                    return marker_summary_row(row, genes, markers, &result.hits, db);
                }
                (Some(Some(genes)), None, Some(aai)) => {
                    return aai_summary_row(row, genes, aai, marker_set, db, options);
                }
                _ => {
                    row.classification = Some(unclassified(marker_set, best_taxonomy.is_some()));
                    return (marker_set, row);
                }
            }
        }
    }
    row.other_related_references = format_related_references(&related, &[], &db.data);
    row.warnings = radius_warnings(&result.hits, db, options.min_af);
    (marker_set, row)
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
    genes: &GeneCalls,
    aai: &AaiResult,
    ani_marker_set: &'static str,
    db: &ClassifyDb,
    options: &ClassifyOptions,
) -> (&'static str, SummaryRow) {
    row.translation_table = Some(genes.translation_table);
    let best_taxonomy = aai.hits.first().and_then(|hit| {
        db.data
            .taxonomy_of(&db.nt_accessions[hit.reference as usize])
    });
    let marker_set = best_taxonomy.map_or(ani_marker_set, marker_set_for_taxonomy);
    let (taxonomy, rank, threshold, supporting, considered) = match &aai.decision {
        AaiDecision::Genus {
            taxonomy,
            supporting,
            considered,
        } => (taxonomy, "genus", options.genus_aai, supporting, considered),
        AaiDecision::Family {
            taxonomy,
            supporting,
            considered,
        } => (
            taxonomy,
            "family",
            options.family_aai,
            supporting,
            considered,
        ),
        AaiDecision::Unresolved => {
            row.classification = Some(unclassified(marker_set, best_taxonomy.is_some()));
            return (marker_set, row);
        }
    };
    row.classification = Some(standardise_taxonomy(taxonomy, None));
    row.classification_method = Some("taxonomic classification defined by AAI".into());
    row.note = Some(format!(
        "{rank} assigned by AAI consensus ({supporting}/{considered} references >= {threshold}%)"
    ));
    (marker_set, row)
}

/// Fill in the marker-tier fields, following GTDB-Tk's rows for genomes
/// without markers, filtered from the MSA, or (until placement) aligned.
fn marker_summary_row(
    mut row: SummaryRow,
    genes: &GeneCalls,
    markers: &MarkerResult,
    hits: &[AniHit],
    db: &ClassifyDb,
) -> (&'static str, SummaryRow) {
    match &markers.outcome {
        MarkerOutcome::NoMarkers => {
            row.classification = Some("Unclassified".into());
            row.warnings.push("No bacterial or archaeal marker".into());
            ("bac120", row)
        }
        MarkerOutcome::Filtered(reason) => {
            row.classification = Some(unclassified(markers.marker_set, true));
            row.warnings.push(reason.clone());
            (markers.marker_set, row)
        }
        MarkerOutcome::Aligned { msa_percent } => {
            row.translation_table = Some(genes.translation_table);
            row.msa_percent = Some(*msa_percent);
            row.warnings.extend(markers.warnings.iter().cloned());
            let Some(placement) = &markers.placement else {
                row.classification = Some(unclassified(markers.marker_set, true));
                return (markers.marker_set, row);
            };
            row.classification = Some(placement.classification.clone());
            row.pplacer_taxonomy = Some(placement.pplacer_taxonomy.clone());
            row.red_value = Some(placement.red_value);
            row.note = placement.note.map(Into::into);
            row.warnings.extend(placement.warning.map(Into::into));
            let Some(method) = placement.method else {
                // GTDB-Tk reports genomes the backbone tree alone classified
                // on a fresh row, without ANI details.
                row.other_related_references = None;
                row.warnings = markers.warnings.clone();
                return (markers.marker_set, row);
            };
            row.classification_method = Some(method.into());
            // GTDB-Tk reports the placement reference for genomes with ANI hits.
            if ani::reportable(hits).next().is_some() {
                row.classification_method =
                    Some("taxonomic classification defined by topology and ANI".into());
                if let Some(reference) = &placement.reference {
                    let accession = strip_db_prefix(reference).to_owned();
                    row.closest_placement_radius = db.data.radius_of(&accession);
                    row.closest_placement_taxonomy = db
                        .data
                        .taxonomy_of(&accession)
                        .map(|taxonomy| taxonomy.join(";"));
                    if let Some(hit) = hits
                        .iter()
                        .find(|hit| db.nt_accessions[hit.reference as usize] == accession)
                    {
                        row.closest_placement_ani = Some(hit.ani);
                        row.closest_placement_af = Some(hit.af);
                    }
                    row.closest_placement_reference = Some(accession);
                }
            }
            (markers.marker_set, row)
        }
    }
}

/// GTDB-Tk's warnings for near misses that did not receive a species call.
fn radius_warnings(hits: &[AniHit], db: &ClassifyDb, min_af: f64) -> Vec<String> {
    let mut warnings = Vec::new();
    if let Some(best) = hits.iter().find(|hit| hit.af >= min_af) {
        let accession = &db.nt_accessions[best.reference as usize];
        if best.ani >= MIN_SPECIES_RADIUS {
            let radius = db
                .data
                .radius_of(accession)
                .map_or_else(|| "None".to_owned(), summary::py_float);
            warnings.push(format!(
                "Genome not assigned to closest species as it falls outside its pre-defined ANI radius (closest {accession} ANI={}%, AF={}, radius={radius}%)",
                py_round(best.ani, 2),
                py_round(best.af, 3)
            ));
        }
    } else if let Some(top) = hits.first() {
        if top.ani >= MIN_SPECIES_RADIUS {
            warnings.push(format!(
                "Closest reference {} has ANI={}% but alignment fraction AF={} below the {} threshold; no species assignment",
                db.nt_accessions[top.reference as usize],
                py_round(top.ani, 2),
                py_round(top.af, 3),
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
