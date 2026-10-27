//! GTDB-Tk-compatible genome classification: ANI screen, AAI, then marker placement.

pub(crate) mod ani;
pub(crate) mod log;
pub(crate) mod refdata;
pub(crate) mod refdb;
pub(crate) mod summary;

use std::collections::{BTreeMap, HashSet};
use std::error::Error;
use std::fs;
use std::path::{Path, PathBuf};
use std::time::Instant;

use crate::classify::ani::{AniDecision, AniHit, MIN_SPECIES_RADIUS, decide, refine_candidates};
use crate::classify::log::ClassifyLog;
use crate::classify::refdata::{marker_set_for_taxonomy, standardise_taxonomy};
use crate::classify::refdb::ClassifyDb;
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

/// Everything known about one query after the ANI tier.
struct AniResult<'a> {
    genome: &'a QueryGenome,
    hits: Vec<AniHit>,
    decision: AniDecision,
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
        .map(|(genome, hits)| AniResult {
            decision: decide(
                &hits,
                &db.nt_accessions,
                &db.data,
                options.min_af,
                options.genus_ani,
            ),
            genome,
            hits,
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
    hits.iter()
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
    result: &AniResult,
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
            row.classification = Some(match marker_set {
                "ar53" if best_taxonomy.is_some() => "Unclassified Archaea".into(),
                _ if best_taxonomy.is_some() => "Unclassified Bacteria".into(),
                _ => "Unclassified".into(),
            });
        }
    }
    row.other_related_references = format_related_references(&related, &[], &db.data);
    row.warnings = radius_warnings(&result.hits, db, options.min_af);
    (marker_set, row)
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
