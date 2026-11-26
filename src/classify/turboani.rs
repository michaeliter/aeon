//! Optional TurboANI check for the ANI tier (`--turboani`, cargo feature
//! `turboani`): when the sketch alone cannot name one species, TurboANI
//! re-scores the plausible candidates with exact ANI and AF.

use std::collections::{HashMap, HashSet};
use std::error::Error;
use std::path::PathBuf;

use turboani::{AniConfig, compare_paths};

use crate::classify::ClassifyOptions;
use crate::classify::ani::{AniHit, Candidate, sort_hits};
use crate::classify::log::ClassifyLog;
use crate::classify::refdata::GtdbtkData;
use crate::classify::refdb::ClassifyDb;

/// Largest number of distinct candidate references TurboANI indexes at once.
const MAX_BATCH_REFERENCES: usize = 256;

/// With TurboANI, a sketch estimate must clear the species radius by this
/// many points to decide on its own: the estimate carries sketching error,
/// and it cannot check AF, so borderline genomes (or chimeras) go to TurboANI.
const SKETCH_RADIUS_MARGIN: f64 = 1.0;

/// How the ANI tier settles one query from its sketch candidates.
#[derive(Debug, Clone, PartialEq)]
enum Screen {
    /// The only plausible candidate is within its species radius by the
    /// sketch estimate alone.
    Estimated(AniHit),
    /// These candidates need TurboANI; none means no candidate is plausible.
    Exact(Vec<u32>),
}

/// Candidates with a sketch ANI of at least `min_sketch_ani` could be the
/// query's species. When exactly one is, and its estimate clears its radius
/// by `SKETCH_RADIUS_MARGIN`, the estimate decides; otherwise every plausible
/// candidate goes to TurboANI. Incomplete genomes lower the estimate, so
/// `min_sketch_ani` sits a few points below the radii.
fn screen(
    candidates: &[Candidate],
    accessions: &[String],
    data: &GtdbtkData,
    min_sketch_ani: f64,
) -> Screen {
    let plausible = candidates
        .iter()
        .filter(|candidate| candidate.sketch_ani >= min_sketch_ani)
        .collect::<Vec<_>>();
    if let [only] = plausible[..] {
        if data
            .radius_of(&accessions[only.reference as usize])
            .is_some_and(|radius| only.sketch_ani >= radius + SKETCH_RADIUS_MARGIN)
        {
            return Screen::Estimated(AniHit {
                reference: only.reference,
                ani: only.sketch_ani,
                af: None,
            });
        }
    }
    Screen::Exact(
        plausible
            .iter()
            .map(|candidate| candidate.reference)
            .collect(),
    )
}

/// Per query, the sketch-only species hit when the candidates name one
/// species, else TurboANI's hits against the plausible candidates.
pub(crate) fn hits(
    query_paths: &[String],
    candidates: &[Vec<Candidate>],
    db: &ClassifyDb,
    options: &ClassifyOptions,
    log: &ClassifyLog,
) -> Result<Vec<Vec<AniHit>>, Box<dyn Error>> {
    let (estimated, exact): (Vec<_>, Vec<_>) = candidates
        .iter()
        .map(|candidates| {
            match screen(
                candidates,
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
    let exact_hits = refine_candidates(query_paths, &exact, &db.nt.paths)?;
    Ok(estimated
        .into_iter()
        .zip(exact_hits)
        .map(|(estimated, exact_hits)| estimated.map_or(exact_hits, |hit| vec![hit]))
        .collect())
}

fn turboani_config() -> AniConfig {
    AniConfig {
        // Report every mapped pair; aeon applies GTDB-Tk's AF threshold itself.
        min_fraction: 0.0,
        ..AniConfig::default()
    }
}

/// Exact ANI between each query and its candidate references, sorted by
/// descending ANI then AF. Pairs below TurboANI's 80% identity floor are absent.
fn refine_candidates(
    query_paths: &[String],
    candidates: &[Vec<u32>],
    reference_paths: &[String],
) -> Result<Vec<Vec<AniHit>>, Box<dyn Error>> {
    let config = turboani_config();
    let mut hits = vec![Vec::new(); query_paths.len()];
    for batch in reference_batches(candidates, MAX_BATCH_REFERENCES) {
        let queries = batch
            .iter()
            .map(|&query| PathBuf::from(&query_paths[query]))
            .collect::<Vec<_>>();
        let mut reference_ids = batch
            .iter()
            .flat_map(|&query| candidates[query].iter().copied())
            .collect::<Vec<_>>();
        reference_ids.sort_unstable();
        reference_ids.dedup();
        let references = reference_ids
            .iter()
            .map(|&id| PathBuf::from(&reference_paths[id as usize]))
            .collect::<Vec<_>>();
        let id_of = references
            .iter()
            .cloned()
            .zip(reference_ids.iter().copied())
            .collect::<HashMap<_, _>>();
        let query_of = queries
            .iter()
            .cloned()
            .zip(batch.iter().copied())
            .collect::<HashMap<_, _>>();

        // (query index, reference id) -> (ANI, AF) in each mapping direction.
        let mut forward = HashMap::new();
        for result in compare_paths(&queries, &references, &config)? {
            let fraction = result.mapped_fragments as f64 / result.total_query_fragments as f64;
            forward.insert(
                (query_of[&result.query], id_of[&result.reference]),
                (result.ani, fraction),
            );
        }
        let mut reverse = HashMap::new();
        for result in compare_paths(&references, &queries, &config)? {
            let fraction = result.mapped_fragments as f64 / result.total_query_fragments as f64;
            reverse.insert(
                (query_of[&result.reference], id_of[&result.query]),
                (result.ani, fraction),
            );
        }

        for &query in &batch {
            let wanted = candidates[query].iter().copied().collect::<HashSet<_>>();
            let mut query_hits = wanted
                .into_iter()
                .filter_map(|reference| {
                    combine_directions(
                        forward.get(&(query, reference)).copied(),
                        reverse.get(&(query, reference)).copied(),
                    )
                    .map(|(ani, af)| AniHit {
                        reference,
                        ani,
                        af: Some(af),
                    })
                })
                .collect::<Vec<_>>();
            sort_hits(&mut query_hits);
            hits[query] = query_hits;
        }
    }
    Ok(hits)
}

/// ANI and AF are each the larger of the two mapping directions, as with
/// GTDB-Tk's skani AF (`max(Align_fraction_ref, Align_fraction_query)`).
fn combine_directions(
    forward: Option<(f64, f64)>,
    reverse: Option<(f64, f64)>,
) -> Option<(f64, f64)> {
    match (forward, reverse) {
        (Some((ani_f, af_f)), Some((ani_r, af_r))) => Some((ani_f.max(ani_r), af_f.max(af_r))),
        (Some(hit), None) | (None, Some(hit)) => Some(hit),
        (None, None) => None,
    }
}

/// Group queries so that each group needs at most `limit` distinct references
/// (a single query with more candidates forms its own group).
fn reference_batches(candidates: &[Vec<u32>], limit: usize) -> Vec<Vec<usize>> {
    let mut batches = Vec::new();
    let mut batch = Vec::new();
    let mut references = HashSet::new();
    for (query, query_candidates) in candidates.iter().enumerate() {
        if query_candidates.is_empty() {
            continue;
        }
        let added = query_candidates
            .iter()
            .filter(|reference| !references.contains(*reference))
            .count();
        if !batch.is_empty() && references.len() + added > limit {
            batches.push(std::mem::take(&mut batch));
            references.clear();
        }
        references.extend(query_candidates.iter().copied());
        batch.push(query);
    }
    if !batch.is_empty() {
        batches.push(batch);
    }
    batches
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::classify::ani::tests::data;

    #[test]
    fn turboani_runs_only_when_the_sketch_is_ambiguous() {
        let (data, accessions) = data();
        let candidate = |reference, sketch_ani| Candidate {
            reference,
            sketch_ani,
        };
        // One plausible candidate within its radius: the estimate decides.
        assert_eq!(
            screen(
                &[candidate(0, 97.0), candidate(2, 85.0)],
                &accessions,
                &data,
                90.0
            ),
            Screen::Estimated(AniHit {
                reference: 0,
                ani: 97.0,
                af: None,
            })
        );
        // Two plausible species, or one without a clear margin over its
        // radius, need TurboANI.
        assert_eq!(
            screen(
                &[candidate(0, 97.0), candidate(1, 93.0)],
                &accessions,
                &data,
                90.0
            ),
            Screen::Exact(vec![0, 1])
        );
        assert_eq!(
            screen(&[candidate(1, 96.0)], &accessions, &data, 90.0),
            Screen::Exact(vec![1])
        );
        assert_eq!(
            screen(&[candidate(0, 95.5)], &accessions, &data, 90.0),
            Screen::Exact(vec![0])
        );
        assert_eq!(
            screen(&[candidate(0, 80.0)], &accessions, &data, 90.0),
            Screen::Exact(vec![])
        );
    }

    #[test]
    fn directions_combine_by_maximum() {
        assert_eq!(
            combine_directions(Some((97.0, 0.6)), Some((96.0, 0.8))),
            Some((97.0, 0.8))
        );
        assert_eq!(
            combine_directions(None, Some((90.0, 0.5))),
            Some((90.0, 0.5))
        );
        assert_eq!(combine_directions(None, None), None);
    }

    #[test]
    fn batches_respect_reference_limit() {
        let candidates = vec![vec![1, 2], vec![2, 3], vec![], vec![4, 5, 6]];
        assert_eq!(reference_batches(&candidates, 3), vec![vec![0, 1], vec![3]]);
        assert_eq!(reference_batches(&candidates, 10), vec![vec![0, 1, 3]]);
    }
}
