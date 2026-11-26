//! Tier 1: nucleotide ANI. DiskANN sketch search proposes candidate species
//! representatives, and GTDB-Tk's ANI screen rule decides the species call
//! over their sketch estimates, or over TurboANI's (see `turboani`).

use crate::classify::refdata::GtdbtkData;
use crate::database::ani_from_jaccard;

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
    /// Larger of the query and reference aligned fractions, 0-1; `None` for
    /// a sketch estimate, which TurboANI did not check.
    pub(crate) af: Option<f64>,
}

impl AniHit {
    /// Hits without an AF are sketch estimates, whose k-mer Jaccard already
    /// falls with the shared fraction, so they pass the AF threshold.
    pub(crate) fn passes_af(&self, min_af: f64) -> bool {
        self.af.is_none_or(|af| af >= min_af)
    }
}

/// A DiskANN neighbour with its sketch-estimated ANI.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Candidate {
    pub(crate) reference: u32,
    pub(crate) sketch_ani: f64,
}

/// Sketch-estimated ANI (0-100) from a DiskANN Hamming distance.
pub(crate) fn sketch_ani(distance: f32, kmer_size: usize) -> f64 {
    f64::from(ani_from_jaccard(1.0 - distance, kmer_size)) * 100.0
}

/// The candidates at or above `MIN_REPORTED_ANI` as sketch-estimated hits,
/// sorted by descending ANI, for an ANI screen without TurboANI.
pub(crate) fn estimated_hits(candidates: &[Candidate]) -> Vec<AniHit> {
    let mut hits = candidates
        .iter()
        .filter(|candidate| candidate.sketch_ani >= MIN_REPORTED_ANI)
        .map(|candidate| AniHit {
            reference: candidate.reference,
            ani: candidate.sketch_ani,
            af: None,
        })
        .collect::<Vec<_>>();
    sort_hits(&mut hits);
    hits
}

pub(crate) fn sort_hits(hits: &mut [AniHit]) {
    hits.sort_by(|left, right| {
        right
            .ani
            .total_cmp(&left.ani)
            .then(right.af.unwrap_or(1.0).total_cmp(&left.af.unwrap_or(1.0)))
            .then(left.reference.cmp(&right.reference))
    });
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) enum AniDecision {
    /// Within the closest AF-passing reference's species radius.
    Species { hit: AniHit },
    /// No species call; left for the AAI tier.
    Unresolved,
}

/// Apply GTDB-Tk 2.8's ANI screen to sorted hits: only the closest hit with
/// `af >= min_af` is checked against its own radius.
pub(crate) fn decide(
    hits: &[AniHit],
    accessions: &[String],
    data: &GtdbtkData,
    min_af: f64,
) -> AniDecision {
    match hits.iter().find(|hit| hit.passes_af(min_af)) {
        Some(top)
            if data
                .radius_of(&accessions[top.reference as usize])
                .is_some_and(|radius| top.ani >= radius) =>
        {
            AniDecision::Species { hit: *top }
        }
        _ => AniDecision::Unresolved,
    }
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    use crate::classify::refdata::Radius;

    pub(crate) fn data() -> (GtdbtkData, Vec<String>) {
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
        AniHit {
            reference,
            ani,
            af: Some(af),
        }
    }

    #[test]
    fn species_needs_closest_af_passing_hit_within_its_radius() {
        let (data, accessions) = data();
        let hits = [hit(0, 96.0, 0.8), hit(1, 95.5, 0.9)];
        assert_eq!(
            decide(&hits, &accessions, &data, 0.5),
            AniDecision::Species { hit: hits[0] }
        );
        // Low AF on the closest hit moves the check to the next hit.
        let hits = [hit(0, 99.0, 0.3), hit(1, 97.5, 0.6)];
        assert_eq!(
            decide(&hits, &accessions, &data, 0.5),
            AniDecision::Species { hit: hits[1] }
        );
    }

    #[test]
    fn lower_hit_within_its_radius_does_not_rescue() {
        let (data, accessions) = data();
        // Closest AF-passing hit (radius 97) fails; GCF_1 (radius 95) would pass.
        let hits = [hit(1, 96.5, 0.8), hit(0, 95.2, 0.8)];
        assert_eq!(
            decide(&hits, &accessions, &data, 0.5),
            AniDecision::Unresolved
        );
    }

    #[test]
    fn sketch_only_screen_uses_the_closest_estimate() {
        let (data, accessions) = data();
        let candidate = |reference, sketch_ani| Candidate {
            reference,
            sketch_ani,
        };
        let estimated = |reference, ani| AniHit {
            reference,
            ani,
            af: None,
        };
        // Sorted by estimate, without candidates below 80%.
        let hits = estimated_hits(&[candidate(2, 79.0), candidate(1, 96.5), candidate(0, 95.5)]);
        assert_eq!(hits, vec![estimated(1, 96.5), estimated(0, 95.5)]);
        // The closest (radius 97) falls short; GCF_1 (radius 95) does not rescue.
        assert_eq!(
            decide(&hits, &accessions, &data, 0.5),
            AniDecision::Unresolved
        );
        let hits = estimated_hits(&[candidate(0, 95.2), candidate(2, 90.0)]);
        assert_eq!(
            decide(&hits, &accessions, &data, 0.5),
            AniDecision::Species { hit: hits[0] }
        );
    }
}
