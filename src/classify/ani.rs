//! Tier 1: nucleotide ANI. DiskANN sketch-search candidates are re-scored
//! with TurboANI, then GTDB-Tk's ANI screen rule decides species calls.

use std::collections::{HashMap, HashSet};
use std::error::Error;
use std::path::PathBuf;

use turboani::{AniConfig, compare_paths};

use crate::classify::refdata::GtdbtkData;

/// Largest number of distinct candidate references TurboANI indexes at once.
const MAX_BATCH_REFERENCES: usize = 256;

/// Lowest ANI radius used by GTDB; references can have a raised radius.
pub(crate) const MIN_SPECIES_RADIUS: f64 = 95.0;

/// skani, which GTDB-Tk uses, screens out pairs below about 80% ANI, so
/// GTDB-Tk never reports such hits; TurboANI can return genome ANI just
/// below 80 from fragments above it.
pub(crate) const MIN_REPORTED_ANI: f64 = 80.0;

pub(crate) fn reportable(hits: &[AniHit]) -> impl Iterator<Item = &AniHit> {
    hits.iter().filter(|hit| hit.ani >= MIN_REPORTED_ANI)
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct AniHit {
    /// Reference vector ID in the nucleotide index.
    pub(crate) reference: u32,
    /// Percent identity, 0-100.
    pub(crate) ani: f64,
    /// Larger of the query and reference aligned fractions, 0-1.
    pub(crate) af: f64,
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
pub(crate) fn refine_candidates(
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
                    .map(|(ani, af)| AniHit { reference, ani, af })
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

pub(crate) fn sort_hits(hits: &mut [AniHit]) {
    hits.sort_by(|left, right| {
        right
            .ani
            .total_cmp(&left.ani)
            .then(right.af.total_cmp(&left.af))
            .then(left.reference.cmp(&right.reference))
    });
}

/// Group queries so that each group needs at most `limit` distinct references
/// (a single query with more candidates forms its own group).
pub(crate) fn reference_batches(candidates: &[Vec<u32>], limit: usize) -> Vec<Vec<usize>> {
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

#[derive(Debug, Clone, PartialEq)]
pub(crate) enum AniDecision {
    /// Within the closest AF-passing reference's species radius.
    Species { hit: AniHit },
    /// Outside every radius but close enough for a consensus genus.
    Genus {
        genus_taxonomy: String,
        supporting: usize,
        considered: usize,
    },
    /// Too distant for a nucleotide-level call.
    Unresolved,
}

/// Apply GTDB-Tk 2.8's ANI screen to sorted hits: only the closest hit with
/// `af >= min_af` is checked against its own radius. Failing that, assign the
/// best hit's genus when it is also the strict-majority genus among hits with
/// ANI at or above `genus_ani` and `af >= min_af`.
pub(crate) fn decide(
    hits: &[AniHit],
    accessions: &[String],
    data: &GtdbtkData,
    min_af: f64,
    genus_ani: f64,
) -> AniDecision {
    if let Some(top) = hits.iter().find(|hit| hit.af >= min_af) {
        let accession = &accessions[top.reference as usize];
        if data
            .radius_of(accession)
            .is_some_and(|radius| top.ani >= radius)
        {
            return AniDecision::Species { hit: *top };
        }
    }

    let genus_of = |hit: &AniHit| {
        data.taxonomy_of(&accessions[hit.reference as usize])
            .map(|ranks| ranks[..6].join(";"))
    };
    // A high ANI over a small aligned fraction (shared operons, mobile
    // elements) is a poor genus signal, so only hits passing the AF
    // threshold vote.
    let close = hits
        .iter()
        .filter(|hit| hit.ani >= genus_ani && hit.af >= min_af)
        .filter_map(genus_of)
        .collect::<Vec<_>>();
    let Some(best_genus) = close.first() else {
        return AniDecision::Unresolved;
    };
    let supporting = close.iter().filter(|genus| *genus == best_genus).count();
    if supporting * 2 > close.len() {
        AniDecision::Genus {
            genus_taxonomy: best_genus.clone(),
            supporting,
            considered: close.len(),
        }
    } else {
        AniDecision::Unresolved
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::classify::refdata::Radius;

    fn data() -> (GtdbtkData, Vec<String>) {
        let mut data = GtdbtkData::default();
        let references = [
            ("GCF_000000001.1", "g__A;s__A one", 95.0),
            ("GCF_000000002.1", "g__A;s__A two", 97.0),
            ("GCF_000000003.1", "g__B;s__B one", 95.0),
        ];
        for (accession, genus_species, radius) in references {
            let canonical = crate::classify::refdata::canonical_gid(accession);
            data.taxonomy.insert(
                canonical.clone(),
                format!("d__Bacteria;p__P;c__C;o__O;f__F;{genus_species}")
                    .split(';')
                    .map(str::to_owned)
                    .collect(),
            );
            data.radii.insert(
                canonical,
                Radius {
                    species: genus_species.split(';').nth(1).unwrap().into(),
                    ani: radius,
                },
            );
        }
        let accessions = references.iter().map(|(a, _, _)| a.to_string()).collect();
        (data, accessions)
    }

    fn hit(reference: u32, ani: f64, af: f64) -> AniHit {
        AniHit { reference, ani, af }
    }

    #[test]
    fn species_needs_closest_af_passing_hit_within_its_radius() {
        let (data, accessions) = data();
        let hits = [hit(0, 96.0, 0.8), hit(1, 95.5, 0.9)];
        assert_eq!(
            decide(&hits, &accessions, &data, 0.5, 80.0),
            AniDecision::Species { hit: hits[0] }
        );
        // Low AF on the closest hit moves the check to the next hit.
        let hits = [hit(0, 99.0, 0.3), hit(1, 97.5, 0.6)];
        assert_eq!(
            decide(&hits, &accessions, &data, 0.5, 80.0),
            AniDecision::Species { hit: hits[1] }
        );
    }

    #[test]
    fn lower_hit_within_its_radius_does_not_rescue() {
        let (data, accessions) = data();
        // Closest AF-passing hit (radius 97) fails; GCF_1 (radius 95) would pass.
        let hits = [hit(1, 96.5, 0.8), hit(0, 95.2, 0.8)];
        assert!(matches!(
            decide(&hits, &accessions, &data, 0.5, 80.0),
            AniDecision::Genus {
                supporting: 2,
                considered: 2,
                ..
            }
        ));
    }

    #[test]
    fn genus_requires_majority_agreement() {
        let (data, accessions) = data();
        let hits = [hit(0, 88.0, 0.6), hit(1, 85.0, 0.5), hit(2, 84.0, 0.5)];
        assert_eq!(
            decide(&hits, &accessions, &data, 0.5, 80.0),
            AniDecision::Genus {
                genus_taxonomy: "d__Bacteria;p__P;c__C;o__O;f__F;g__A".into(),
                supporting: 2,
                considered: 3,
            }
        );
        let split = [hit(0, 88.0, 0.6), hit(2, 85.0, 0.5)];
        assert_eq!(
            decide(&split, &accessions, &data, 0.5, 80.0),
            AniDecision::Unresolved
        );
        assert_eq!(
            decide(&[hit(0, 79.0, 0.6)], &accessions, &data, 0.5, 80.0),
            AniDecision::Unresolved
        );
        // Hits below the AF threshold do not vote.
        assert_eq!(
            decide(&[hit(0, 90.0, 0.1)], &accessions, &data, 0.5, 80.0),
            AniDecision::Unresolved
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
