use anndists::dist::DistHamming;
use log::debug;
use rayon::prelude::*;
use rust_diskann::{DiskANN, DiskAnnParams, GuardedDeleteResult, RoutabilityAdmissionReport};
use std::collections::{HashMap, HashSet};
use std::error::Error;
use std::fs::{self, File};
use std::io::{self, BufWriter, Write};
use std::path::Path;
use std::time::Instant;

use crate::metadata::{
    PrefixParams, genomes_path, idmap_path, index_path, load_params, save_params,
};
use crate::sketch::{
    GenomeContentHash, SketchedGenome, genome_content_hashes_from_params,
    sketch_and_hash_from_params, sketch_from_params,
};
use crate::utils::{read_idmap_content_hashes, read_list_file, write_genome_list, write_idmap_tsv};
use crate::warmup::Warmup;

pub(crate) fn build_database(
    prefix: &str,
    reference_list: &str,
    params: PrefixParams,
) -> Result<(), Box<dyn Error>> {
    let started = Instant::now();
    let genomes = read_list_file(reference_list)?;
    if genomes.is_empty() {
        return Err("reference list is empty".into());
    }

    eprintln!(
        "Building database: prefix={}, genomes={}, seq_type={}, k={}, sketch_size={}, R={}, build_beam={}",
        prefix,
        genomes.len(),
        params.seq_type,
        params.kmer_size,
        params.sketch_size,
        params.max_degree,
        params.build_beam_width
    );
    for genome in &genomes {
        if !Path::new(genome).is_file() {
            return Err(
                format!("reference genome does not exist or is not a file: {genome}").into(),
            );
        }
    }

    let records = sketch_and_hash_from_params(&genomes, &params);
    let no_existing_hashes = HashSet::new();
    let CanonicalGenomes {
        paths: genomes,
        vectors,
        hashes,
        report,
    } = canonicalize_genomes(genomes, records, &no_existing_hashes);
    if report.skipped_total() > 0 {
        eprintln!(
            "Deduplicated reference genomes by content hash: requested={}, unique={}, skipped_within_input={}",
            report.requested, report.unique, report.skipped_within_input
        );
    }
    let diskann_params = DiskAnnParams {
        max_degree: params.max_degree,
        build_beam_width: params.build_beam_width,
        alpha: params.alpha,
        extra_seeds: params.extra_seeds,
    };
    let index = DiskANN::<u16, DistHamming>::build_index_with_params(
        &vectors,
        DistHamming,
        &index_path(prefix),
        diskann_params,
    )?;
    sanity_check_index_mapping(&index, &vectors, 3, params.build_beam_width);
    write_genome_list(&genomes_path(prefix), &genomes)?;
    write_idmap_tsv(&idmap_path(prefix), &genomes, Some(&hashes))?;
    save_params(prefix, &params)?;
    eprintln!(
        "Built database: prefix={}, genomes={}, elapsed={:.3}s",
        prefix,
        genomes.len(),
        started.elapsed().as_secs_f64()
    );
    Ok(())
}

pub(crate) fn update_database(
    prefix: &str,
    insert_paths: Vec<String>,
    delete_names: Vec<String>,
    beam_width: Option<usize>,
) -> Result<(), Box<dyn Error>> {
    let started = Instant::now();
    if insert_paths.is_empty() && delete_names.is_empty() {
        return Err("update contains neither insertions nor deletions".into());
    }
    let params = load_params(prefix)?;
    let old_genomes = read_list_file(&genomes_path(prefix))?;
    let delete_ids = resolve_delete_ids(&old_genomes, &delete_names)?;
    let deleted_ids = delete_ids.iter().copied().collect::<HashSet<_>>();
    let retained_names = old_genomes
        .iter()
        .enumerate()
        .filter(|(id, _)| !deleted_ids.contains(&(*id as u32)))
        .map(|(_, name)| name.as_str())
        .collect::<HashSet<_>>();
    for path in &insert_paths {
        if !Path::new(path).is_file() {
            return Err(format!("insert genome does not exist or is not a file: {path}").into());
        }
    }

    let requested_insertion_count = insert_paths.len();
    let insert_records = if insert_paths.is_empty() {
        Vec::new()
    } else {
        sketch_and_hash_from_params(&insert_paths, &params)
    };
    let source = DiskANN::<u16, DistHamming>::open_index_with(&index_path(prefix), DistHamming)?;
    if source.num_vectors != old_genomes.len() {
        return Err(format!(
            "index has {} vectors but genome map has {} entries",
            source.num_vectors,
            old_genomes.len()
        )
        .into());
    }
    let old_hashes = if insert_paths.is_empty() {
        read_idmap_content_hashes(&idmap_path(prefix), &old_genomes)?
    } else {
        Some(load_or_recompute_genome_hashes(
            prefix,
            &old_genomes,
            &params,
        )?)
    };
    let retained_hashes = old_hashes
        .as_ref()
        .map(|hashes| {
            hashes
                .iter()
                .enumerate()
                .filter(|(id, _)| !deleted_ids.contains(&(*id as u32)))
                .map(|(_, hash)| *hash)
                .collect::<HashSet<_>>()
        })
        .unwrap_or_default();
    let CanonicalGenomes {
        paths: insert_paths,
        vectors: insert_vectors,
        hashes: insert_hashes,
        report: insertion_dedup,
    } = canonicalize_genomes(insert_paths, insert_records, &retained_hashes);
    let insertion_count = insert_paths.len();
    if let Some(path) = insert_paths
        .iter()
        .find(|path| retained_names.contains(path.as_str()))
    {
        return Err(format!(
            "genome path is already present with different compressed genome content: {path}; delete it before inserting a replacement"
        )
        .into());
    }
    if insertion_dedup.skipped_total() > 0 {
        eprintln!(
            "Deduplicated insertion genomes by content hash: requested={}, unique={}, skipped_existing={}, skipped_within_input={}",
            insertion_dedup.requested,
            insertion_dedup.unique,
            insertion_dedup.skipped_existing,
            insertion_dedup.skipped_within_input
        );
    }
    if delete_ids.is_empty() && insert_vectors.is_empty() {
        eprintln!(
            "No unique genome-content hashes to insert; database remains unchanged (requested={requested_insertion_count})"
        );
        return Ok(());
    }

    let final_count = old_genomes.len() - delete_ids.len() + insertion_count;
    if final_count == 0 {
        return Err("an update cannot produce an empty database".into());
    }
    let capacity = old_genomes.len().max(final_count);
    let work_path = format!("{prefix}.aeon-work-{}", std::process::id());
    let temporary_index = format!("{}.aeon-next-{}", index_path(prefix), std::process::id());
    eprintln!(
        "Updating database: prefix={}, current={}, delete={}, insert_requested={}, insert_unique={}, final={}",
        prefix,
        old_genomes.len(),
        delete_ids.len(),
        requested_insertion_count,
        insertion_count,
        final_count
    );
    let mut update = source.begin_updates(capacity, params.alpha, &work_path)?;
    let mut slot_names = old_genomes.into_iter().map(Some).collect::<Vec<_>>();
    slot_names.resize(capacity, None);
    let mut slot_hashes = old_hashes.map(|hashes| {
        let mut slots = hashes.into_iter().map(Some).collect::<Vec<_>>();
        slots.resize(capacity, None);
        slots
    });

    if !delete_ids.is_empty() {
        let guard = update.routability_guard_status()?;
        if !guard.enabled {
            drop(update);
            discard_update_workspace(&work_path);
            return Err(
                "rust-diskann did not enable routability admission control for the update session"
                    .into(),
            );
        }
        eprintln!(
            "Checking deletion routability: landmarks={}, pool={}, beam={}",
            guard.landmark_count, guard.landmark_pool_size, guard.beam_width
        );
        let stats = match update.delete_batch_with_admission_control(&delete_ids)? {
            GuardedDeleteResult::Applied { stats, report } => {
                eprintln!(
                    "Deletion admission accepted: checked_landmarks={}, baseline_routes={}, residual_routes={}",
                    report.checked_landmarks.len(),
                    report.baseline_reachable_landmarks,
                    report.residual_reachable_landmarks
                );
                stats
            }
            GuardedDeleteResult::Deferred(report) => {
                drop(update);
                discard_update_workspace(&work_path);
                return Err(deferred_delete_error(&report));
            }
        };
        for id in &delete_ids {
            slot_names[*id as usize] = None;
            if let Some(hashes) = &mut slot_hashes {
                hashes[*id as usize] = None;
            }
        }
        let recovered: usize = stats.iter().map(|stat| stat.recovered_in_neighbors).sum();
        let candidates: usize = stats.iter().map(|stat| stat.repair_candidates).sum();
        debug!(
            "Deleted {} genomes; MERIT averages: {:.2} recovered in-neighbors, {:.2} repair candidates",
            stats.len(),
            recovered as f64 / stats.len() as f64,
            candidates as f64 / stats.len() as f64
        );
        eprintln!("Deleted {} genomes", stats.len());
    }
    if !insert_vectors.is_empty() {
        let beam = beam_width.unwrap_or(params.build_beam_width);
        let inserted_ids = update.insert_batch(insert_vectors, beam)?;
        for ((id, name), hash) in inserted_ids
            .into_iter()
            .zip(insert_paths)
            .zip(insert_hashes)
        {
            slot_names[id as usize] = Some(name);
            if let Some(hashes) = &mut slot_hashes {
                hashes[id as usize] = Some(hash);
            }
        }
        eprintln!("Inserted {insertion_count} unique genomes with beam {beam}");
    }

    let (committed, old_to_new) = update.commit_updates_to_static(&temporary_index)?;
    if committed.num_vectors != final_count {
        return Err(format!(
            "committed index has {} vectors; expected {final_count}",
            committed.num_vectors
        )
        .into());
    }
    let new_genomes = rebuild_genome_order(&slot_names, &old_to_new, final_count)?;
    let new_hashes = match slot_hashes {
        Some(slot_hashes) => Some(rebuild_hash_order(&slot_hashes, &old_to_new, final_count)?),
        None => None,
    };
    drop(committed);
    install_updated_files(
        prefix,
        &temporary_index,
        &new_genomes,
        new_hashes.as_deref(),
    )?;
    eprintln!(
        "Updated database: prefix={}, genomes={}, elapsed={:.3}s",
        prefix,
        new_genomes.len(),
        started.elapsed().as_secs_f64()
    );
    Ok(())
}

#[derive(Debug, Default, PartialEq, Eq)]
struct GenomeDedupReport {
    requested: usize,
    unique: usize,
    skipped_existing: usize,
    skipped_within_input: usize,
}

impl GenomeDedupReport {
    fn skipped_total(&self) -> usize {
        self.skipped_existing + self.skipped_within_input
    }
}

#[derive(Debug)]
struct CanonicalGenomes {
    paths: Vec<String>,
    vectors: Vec<Vec<u16>>,
    hashes: Vec<GenomeContentHash>,
    report: GenomeDedupReport,
}

/// Retains only the first appearance of a SHA-256 genome-content hash. Unlike
/// b-bit MinHash equality, this does not merge merely similar genomes that
/// happen to produce the same approximate sketch.
fn canonicalize_genomes(
    paths: Vec<String>,
    records: Vec<SketchedGenome>,
    existing_hashes: &HashSet<GenomeContentHash>,
) -> CanonicalGenomes {
    assert_eq!(
        paths.len(),
        records.len(),
        "each path needs one sketch record"
    );
    let requested = paths.len();
    let mut canonical_paths = Vec::with_capacity(requested);
    let mut canonical_vectors = Vec::with_capacity(requested);
    let mut canonical_hashes = Vec::with_capacity(requested);
    let mut seen_in_input = HashSet::with_capacity(requested);
    let mut report = GenomeDedupReport {
        requested,
        ..GenomeDedupReport::default()
    };

    for (path, record) in paths.into_iter().zip(records) {
        if existing_hashes.contains(&record.content_hash) {
            report.skipped_existing += 1;
            continue;
        }
        if !seen_in_input.insert(record.content_hash) {
            report.skipped_within_input += 1;
            continue;
        }
        canonical_paths.push(path);
        canonical_vectors.push(record.vector);
        canonical_hashes.push(record.content_hash);
    }
    report.unique = canonical_paths.len();
    CanonicalGenomes {
        paths: canonical_paths,
        vectors: canonical_vectors,
        hashes: canonical_hashes,
        report,
    }
}

fn load_or_recompute_genome_hashes(
    prefix: &str,
    genomes: &[String],
    params: &PrefixParams,
) -> Result<Vec<GenomeContentHash>, Box<dyn Error>> {
    if let Some(hashes) = read_idmap_content_hashes(&idmap_path(prefix), genomes)? {
        return Ok(hashes);
    }
    eprintln!(
        "Genome content hashes are absent from {}. Recomputing hashes from {} stored genomes",
        idmap_path(prefix),
        genomes.len()
    );
    Ok(genome_content_hashes_from_params(genomes, params))
}

fn discard_update_workspace(work_path: &str) {
    if let Err(error) = fs::remove_file(work_path) {
        if error.kind() != io::ErrorKind::NotFound {
            debug!(
                "failed to remove deferred update workspace path={} error={}",
                work_path, error
            );
        }
    }
}

fn deferred_delete_error(report: &RoutabilityAdmissionReport) -> Box<dyn Error> {
    format!(
        "deletion deferred by routability-aware admission control: requested={}, checked_landmarks={}, baseline_routes={}, residual_routes={}, lost_landmark_ids={:?}; the database was not changed. Split the delete list or rebuild the index.",
        report.requested,
        report.checked_landmarks.len(),
        report.baseline_reachable_landmarks,
        report.residual_reachable_landmarks,
        report.lost_landmarks
    )
    .into()
}

pub(crate) fn search_database(
    prefix: &str,
    query_list: &str,
    k: usize,
    beam: usize,
    output_path: Option<&str>,
    preload: bool,
) -> Result<(), Box<dyn Error>> {
    let started = Instant::now();
    let params = load_params(prefix)?;
    let references = read_list_file(&genomes_path(prefix))?;
    let queries = read_list_file(query_list)?;
    if queries.is_empty() {
        return Err("query list is empty".into());
    }
    eprintln!(
        "Searching database: prefix={}, queries={}, references={}, k={}, beam={}",
        prefix,
        queries.len(),
        references.len(),
        k,
        beam
    );
    // Warm the index while the queries are sketched.
    let warmup = Warmup::start(&[prefix], preload)?;
    let vectors = sketch_from_params(&queries, &params);
    eprintln!("{}", warmup.wait(prefix)?);
    let hits = search_vectors(prefix, references.len(), &vectors, k, beam)?
        .into_iter()
        .enumerate()
        .collect::<Vec<_>>();
    let output: Box<dyn Write> = match output_path {
        Some(path) => Box::new(BufWriter::new(File::create(path)?)),
        None => Box::new(BufWriter::new(io::stdout())),
    };
    write_search_results(
        output,
        &queries,
        &hits,
        &references,
        &params.seq_type,
        params.kmer_size,
    )?;
    eprintln!(
        "Searched {} queries in {:.3}s",
        queries.len(),
        started.elapsed().as_secs_f64()
    );
    Ok(())
}

/// Search sketch vectors against the static index at `prefix`, returning the
/// ranked `(vector_id, distance)` neighbors of each vector in input order.
pub(crate) fn search_vectors(
    prefix: &str,
    reference_count: usize,
    vectors: &[Vec<u16>],
    k: usize,
    beam: usize,
) -> Result<Vec<Vec<(u32, f32)>>, Box<dyn Error>> {
    let index = DiskANN::<u16, DistHamming>::open_index_with(&index_path(prefix), DistHamming)?;
    if index.num_vectors != reference_count {
        return Err("index and genome mapping have different lengths".into());
    }
    Ok(vectors
        .par_iter()
        .map(|vector| index.search_with_dists(vector, k, beam))
        .collect())
}

fn resolve_delete_ids(genomes: &[String], names: &[String]) -> Result<Vec<u32>, Box<dyn Error>> {
    let mut by_name = HashMap::with_capacity(genomes.len());
    let mut by_basename: HashMap<&str, Vec<u32>> = HashMap::with_capacity(genomes.len());
    for (id, name) in genomes.iter().enumerate() {
        if by_name.insert(name.as_str(), id as u32).is_some() {
            return Err(format!("database contains duplicate genome name: {name}").into());
        }
        let basename = Path::new(name)
            .file_name()
            .and_then(|filename| filename.to_str())
            .unwrap_or(name);
        by_basename.entry(basename).or_default().push(id as u32);
    }
    let mut seen_ids = HashSet::with_capacity(names.len());
    let mut ids = Vec::with_capacity(names.len());
    for name in names {
        let id = match by_name.get(name.as_str()).copied() {
            Some(id) => id,
            None => match by_basename.get(name.as_str()) {
                Some(candidates) if candidates.len() == 1 => candidates[0],
                Some(_) => {
                    return Err(format!(
                        "genome basename is ambiguous in the database: {name}; use its complete stored path"
                    )
                    .into());
                }
                None => {
                    return Err(
                        format!("genome name is not present in the database: {name}").into(),
                    );
                }
            },
        };
        if seen_ids.insert(id) {
            ids.push(id);
        }
    }
    Ok(ids)
}

fn rebuild_genome_order(
    slot_names: &[Option<String>],
    old_to_new: &[Option<u32>],
    expected: usize,
) -> Result<Vec<String>, Box<dyn Error>> {
    let mut genomes = vec![None; expected];
    for (old_id, new_id) in old_to_new.iter().enumerate() {
        if let Some(new_id) = new_id {
            let name = slot_names
                .get(old_id)
                .and_then(Option::as_ref)
                .ok_or_else(|| format!("missing genome name for live slot {old_id}"))?;
            genomes[*new_id as usize] = Some(name.clone());
        }
    }
    genomes
        .into_iter()
        .enumerate()
        .map(|(id, name)| name.ok_or_else(|| format!("missing genome name for new ID {id}").into()))
        .collect()
}

fn rebuild_hash_order(
    slot_hashes: &[Option<GenomeContentHash>],
    old_to_new: &[Option<u32>],
    expected: usize,
) -> Result<Vec<GenomeContentHash>, Box<dyn Error>> {
    let mut hashes = vec![None; expected];
    for (old_id, new_id) in old_to_new.iter().enumerate() {
        if let Some(new_id) = new_id {
            let hash = slot_hashes
                .get(old_id)
                .and_then(Option::as_ref)
                .ok_or_else(|| format!("missing genome content hash for live slot {old_id}"))?;
            hashes[*new_id as usize] = Some(*hash);
        }
    }
    hashes
        .into_iter()
        .enumerate()
        .map(|(id, hash)| {
            hash.ok_or_else(|| format!("missing genome content hash for new ID {id}").into())
        })
        .collect()
}

fn install_updated_files(
    prefix: &str,
    temporary_index: &str,
    genomes: &[String],
    hashes: Option<&[GenomeContentHash]>,
) -> Result<(), Box<dyn Error>> {
    let suffix = format!("aeon-next-{}", std::process::id());
    let temporary_genomes = format!("{}.{}", genomes_path(prefix), suffix);
    let temporary_idmap = format!("{}.{}", idmap_path(prefix), suffix);
    write_genome_list(&temporary_genomes, genomes)?;
    write_idmap_tsv(&temporary_idmap, genomes, hashes)?;
    fs::rename(&temporary_genomes, genomes_path(prefix))?;
    fs::rename(&temporary_idmap, idmap_path(prefix))?;
    fs::rename(temporary_index, index_path(prefix))?;
    Ok(())
}

fn write_search_results(
    mut output: Box<dyn Write>,
    query_paths: &[String],
    hits: &[(usize, Vec<(u32, f32)>)],
    references: &[String],
    seq_type: &str,
    kmer_size: usize,
) -> io::Result<()> {
    writeln!(output, "Query\tHit\tVectorID\tRawJaccard\tHammingDist\tANI")?;
    for (query, neighbors) in hits {
        for (id, distance) in neighbors {
            let raw_jaccard = (1.0 - *distance).clamp(0.0, 1.0);
            let ani = if seq_type.eq_ignore_ascii_case("dna") {
                format!("{:.6}", ani_from_jaccard(raw_jaccard, kmer_size))
            } else {
                "NA".to_owned()
            };
            writeln!(
                output,
                "{}\t{}\t{}\t{:.6}\t{:.6}\t{}",
                query_paths[*query], references[*id as usize], id, raw_jaccard, distance, ani
            )?;
        }
    }
    Ok(())
}

/// Estimate nucleotide ANI from a MinHash-style Jaccard estimate using Mash's
/// k-mer survival model.
pub(crate) fn ani_from_jaccard(jaccard: f32, kmer_size: usize) -> f32 {
    if kmer_size == 0 {
        return 0.0;
    }
    let jaccard = jaccard.clamp(0.0, 1.0);
    if jaccard == 0.0 {
        return 0.0;
    }
    let shared_kmer_probability = (2.0 * jaccard / (1.0 + jaccard)).clamp(f32::MIN_POSITIVE, 1.0);
    (1.0 + shared_kmer_probability.ln() / kmer_size as f32).clamp(0.0, 1.0)
}

fn sanity_check_index_mapping(
    index: &DiskANN<u16, DistHamming>,
    vectors: &[Vec<u16>],
    num_checks: usize,
    build_beam_width: usize,
) {
    let checks = vectors.len().min(num_checks.max(1));
    let beam = build_beam_width.max(64);
    for (id, vector) in vectors.iter().take(checks).enumerate() {
        match index.search_with_dists(vector, 1, beam).first().copied() {
            Some((found, distance)) => debug!(
                "index mapping check reference_id={id} nearest_id={found} distance={distance:.6}"
            ),
            None => debug!("index mapping check reference_id={id} returned no neighbors"),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn delete_names_accept_unique_basenames_deduplicate_and_reject_ambiguous_ones() {
        let genomes = vec![
            "/db/a.fna.gz".into(),
            "/db/b.fna.gz".into(),
            "/other/a.fna.gz".into(),
        ];
        assert_eq!(
            resolve_delete_ids(&genomes, &["/db/b.fna.gz".into()]).unwrap(),
            vec![1]
        );
        assert_eq!(
            resolve_delete_ids(&genomes, &["b.fna.gz".into()]).unwrap(),
            vec![1]
        );
        assert!(resolve_delete_ids(&genomes, &["a.fna.gz".into()]).is_err());
        assert_eq!(
            resolve_delete_ids(
                &genomes,
                &["/db/b.fna.gz".into(), "b.fna.gz".into(), "b.fna.gz".into(),]
            )
            .unwrap(),
            vec![1]
        );
    }

    #[test]
    fn rebuilds_names_from_compacted_ids() {
        let slots = vec![Some("a".into()), None, Some("c".into()), Some("d".into())];
        let mapping = vec![Some(0), None, Some(1), Some(2)];
        assert_eq!(
            rebuild_genome_order(&slots, &mapping, 3).unwrap(),
            vec!["a", "c", "d"]
        );
    }

    #[test]
    fn ani_is_one_minus_mash_distance() {
        let identity = 0.95_f32;
        let shared_kmer_probability = identity.powi(16);
        let jaccard = shared_kmer_probability / (2.0 - shared_kmer_probability);
        let expected = 1.0 + identity.ln();
        assert!((ani_from_jaccard(jaccard, 16) - expected).abs() < 1e-6);
    }

    #[test]
    fn ani_is_bounded_at_extreme_jaccards() {
        assert_eq!(ani_from_jaccard(0.0, 16), 0.0);
        assert_eq!(ani_from_jaccard(1.0, 16), 1.0);
    }

    fn content_hash(value: u8) -> GenomeContentHash {
        [value; 32]
    }

    fn record(vector: Vec<u16>, hash: u8) -> SketchedGenome {
        SketchedGenome {
            vector,
            content_hash: content_hash(hash),
        }
    }

    #[test]
    fn content_hash_dedup_preserves_first_without_merging_equal_sketches() {
        let result = canonicalize_genomes(
            vec![
                "first".into(),
                "different-content".into(),
                "duplicate".into(),
            ],
            vec![
                record(vec![1, 2, 3], 7),
                record(vec![1, 2, 3], 8),
                record(vec![9, 9, 9], 7),
            ],
            &HashSet::new(),
        );

        assert_eq!(result.paths, vec!["first", "different-content"]);
        assert_eq!(result.vectors, vec![vec![1, 2, 3], vec![1, 2, 3]]);
        assert_eq!(result.hashes, vec![content_hash(7), content_hash(8)]);
        assert_eq!(
            result.report,
            GenomeDedupReport {
                requested: 3,
                unique: 2,
                skipped_existing: 0,
                skipped_within_input: 1,
            }
        );
    }

    #[test]
    fn existing_content_hash_is_not_reinserted() {
        let existing = HashSet::from([content_hash(3)]);
        let result = canonicalize_genomes(
            vec!["existing-copy".into(), "new-genome".into()],
            vec![record(vec![4, 5], 3), record(vec![4, 5], 4)],
            &existing,
        );

        assert_eq!(result.paths, vec!["new-genome"]);
        assert_eq!(result.hashes, vec![content_hash(4)]);
        assert_eq!(
            result.report,
            GenomeDedupReport {
                requested: 2,
                unique: 1,
                skipped_existing: 1,
                skipped_within_input: 0,
            }
        );
    }
}
