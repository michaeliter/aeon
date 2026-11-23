//! Marker-gene rank calls: identity over the masked bac120/ar53 marker
//! alignment against every reference of the query's marker set, and a sibling
//! radius for each taxon from phylum to genus.
//!
//! A reference's sibling identity in a taxon is its highest identity to a
//! member of a different child taxon. The taxon's sibling radius is the lowest
//! sibling identity among its members, so a genome of a child taxon missing
//! from the database still reaches the radius at its closest reference
//! whenever its closest relatives in the taxon are as close as those of the
//! references.

use rayon::prelude::*;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::error::Error;
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

use crate::classify::refdata::{is_named, read_lines};

/// The ranks with radii, phylum to genus, as taxonomy indexes.
const PHYLUM: usize = 1;
const GENUS: usize = 5;
const SPECIES: usize = 6;
/// Fewer shared residue columns than this make a pair incomparable.
const MIN_OVERLAP: u32 = 500;
/// A reference ties with the best match when its identity is within this
/// many standard errors of the difference between the two.
const TIE_STANDARD_ERRORS: f64 = 1.96;
const MARKER_RADII_FILE: &str = "marker_radii.tsv";

/// A masked marker alignment row as bytes: the 20 standard amino acids in
/// upper case, 0 for gaps and anything else.
pub(crate) fn encode_row(masked: &[u8]) -> Vec<u8> {
    masked.iter().map(|&residue| encode_residue(residue)).collect()
}

fn encode_residue(residue: u8) -> u8 {
    match residue.to_ascii_uppercase() {
        residue @ (b'A' | b'C' | b'D' | b'E' | b'F' | b'G' | b'H' | b'I' | b'K' | b'L' | b'M'
        | b'N' | b'P' | b'Q' | b'R' | b'S' | b'T' | b'V' | b'W' | b'Y') => residue,
        _ => 0,
    }
}

/// Identical residues and the columns where both rows have a residue.
fn identity_counts(left: &[u8], right: &[u8]) -> (u32, u32) {
    let (mut matches, mut overlap) = (0, 0);
    // Byte counters over chunks too short to overflow them vectorise well;
    // a multiple of 64 leaves no scalar tail.
    for (left, right) in left.chunks(192).zip(right.chunks(192)) {
        let (mut chunk_matches, mut chunk_overlap) = (0u8, 0u8);
        for (&a, &b) in left.iter().zip(right) {
            let both = (a != 0) & (b != 0);
            chunk_overlap += both as u8;
            chunk_matches += (both & (a == b)) as u8;
        }
        matches += u32::from(chunk_matches);
        overlap += u32::from(chunk_overlap);
    }
    (matches, overlap)
}

/// Percent identity and overlap of a comparable pair.
pub(crate) fn percent_identity(left: &[u8], right: &[u8]) -> Option<(f64, u32)> {
    let (matches, overlap) = identity_counts(left, right);
    (overlap >= MIN_OVERLAP).then(|| (100.0 * f64::from(matches) / f64::from(overlap), overlap))
}

fn alignment_paths(dir: &Path, marker_set: &str) -> (PathBuf, PathBuf) {
    (
        dir.join(format!("{marker_set}.marker_msa.txt")),
        dir.join(format!("{marker_set}.marker_msa.ids")),
    )
}

/// Write one marker set's masked reference rows, one line each, and their
/// accessions in the same order. A set without references gets empty files,
/// replacing any from an earlier build.
pub(crate) fn write_marker_alignment(
    dir: &Path,
    marker_set: &str,
    rows: &[(String, Vec<u8>)],
) -> std::io::Result<()> {
    let (msa_path, ids_path) = alignment_paths(dir, marker_set);
    let mut msa = BufWriter::new(File::create(msa_path)?);
    let mut ids = BufWriter::new(File::create(ids_path)?);
    let mut line = Vec::new();
    for (accession, row) in rows {
        line.clear();
        line.extend(row.iter().map(|&residue| if residue == 0 { b'-' } else { residue }));
        line.push(b'\n');
        msa.write_all(&line)?;
        writeln!(ids, "{accession}")?;
    }
    msa.flush()?;
    ids.flush()
}

/// One marker set's reference rows, encoded in place in a single buffer.
pub(crate) struct MarkerAlignment {
    /// Reference accessions in row order.
    pub(crate) accessions: Vec<String>,
    /// Rows of `width` columns, each followed by its (encoded) line break.
    data: Vec<u8>,
    width: usize,
}

impl MarkerAlignment {
    /// Read one marker set's rows (none when the database has no reference
    /// of the set).
    pub(crate) fn read(dir: &Path, marker_set: &str) -> Result<Self, Box<dyn Error>> {
        let (msa_path, ids_path) = alignment_paths(dir, marker_set);
        let missing = |path: &Path, error: std::io::Error| {
            format!(
                "cannot read {}: {error}; rebuild the database with this aeon classify-db",
                path.display()
            )
        };
        let accessions = fs::read_to_string(&ids_path)
            .map_err(|error| missing(&ids_path, error))?
            .lines()
            .map(str::to_owned)
            .collect::<Vec<_>>();
        let mut data = fs::read(&msa_path).map_err(|error| missing(&msa_path, error))?;
        let width = data.iter().position(|&byte| byte == b'\n').unwrap_or(0);
        if data.len() != accessions.len() * (width + 1) {
            return Err(format!(
                "{} does not hold one row of equal width for each accession in {}",
                msa_path.display(),
                ids_path.display()
            )
            .into());
        }
        data.par_iter_mut()
            .for_each(|byte| *byte = encode_residue(*byte));
        Ok(Self {
            accessions,
            data,
            width,
        })
    }

    pub(crate) fn len(&self) -> usize {
        self.accessions.len()
    }

    pub(crate) fn row(&self, index: usize) -> &[u8] {
        let start = index * (self.width + 1);
        &self.data[start..start + self.width]
    }
}

/// Sibling and nearest radii of one taxon. Nearest identity is a member's
/// highest identity to any other member.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct TaxonRadii {
    /// Taxonomy index of the taxon's rank.
    pub(crate) rank: usize,
    pub(crate) sibling: Option<f64>,
    pub(crate) nearest: Option<f64>,
}

/// The sibling and nearest radii of every named taxon from phylum to genus
/// among `references` (taxonomy, encoded row). Only pairs within a phylum
/// share a taxon at these ranks, so only they are compared.
pub(crate) fn taxon_radii(references: &[(&[String], &[u8])]) -> BTreeMap<String, TaxonRadii> {
    // Taxon IDs from class to species, so pairs compare integers.
    let mut taxa: HashMap<&str, u32> = HashMap::new();
    let ids = references
        .iter()
        .map(|(taxonomy, _)| {
            let mut ids = [0; SPECIES - PHYLUM];
            for (id, taxon) in ids.iter_mut().zip(&taxonomy[PHYLUM + 1..=SPECIES]) {
                let next = taxa.len() as u32;
                *id = *taxa.entry(taxon).or_insert(next);
            }
            ids
        })
        .collect::<Vec<_>>();
    let mut phyla: BTreeMap<&str, Vec<usize>> = BTreeMap::new();
    for (index, (taxonomy, _)) in references.iter().enumerate() {
        phyla.entry(&taxonomy[PHYLUM]).or_default().push(index);
    }
    // Per reference, the highest identity to a reference whose deepest
    // shared rank is each of phylum to genus, then to one of the same
    // species: the sibling identity at that rank. Identities are stored as
    // f64 bits plus one (0 for none), as non-negative floats order like
    // their bits.
    const SLOTS: usize = SPECIES - PHYLUM + 1;
    let maxima = (0..references.len() * SLOTS)
        .map(|_| AtomicU64::new(0))
        .collect::<Vec<_>>();
    let raise = |slot: &AtomicU64, bits: u64| {
        if slot.load(Ordering::Relaxed) < bits {
            slot.fetch_max(bits, Ordering::Relaxed);
        }
    };
    for members in phyla.values() {
        members.par_iter().enumerate().for_each(|(position, &left)| {
            let mut own = [0; SLOTS];
            for &right in &members[position + 1..] {
                let Some((identity, _)) = percent_identity(references[left].1, references[right].1)
                else {
                    continue;
                };
                let slot = (0..SLOTS - 1)
                    .find(|&rank| ids[left][rank] != ids[right][rank])
                    .unwrap_or(SLOTS - 1);
                let bits = identity.to_bits() + 1;
                own[slot] = own[slot].max(bits);
                raise(&maxima[right * SLOTS + slot], bits);
            }
            for (slot, bits) in own.into_iter().enumerate() {
                raise(&maxima[left * SLOTS + slot], bits);
            }
        });
    }
    let value = |bits: u64| (bits > 0).then(|| f64::from_bits(bits - 1));
    let lowest = |radius: &mut Option<f64>, identity: Option<f64>| {
        if let Some(identity) = identity {
            *radius = Some(radius.map_or(identity, |radius: f64| radius.min(identity)));
        }
    };
    let mut radii: BTreeMap<&str, TaxonRadii> = BTreeMap::new();
    for (reference, (taxonomy, _)) in references.iter().enumerate() {
        let maxima = &maxima[reference * SLOTS..(reference + 1) * SLOTS];
        // Nearest identity at a rank covers every deeper shared rank.
        let mut nearest = maxima[SLOTS - 1].load(Ordering::Relaxed);
        for rank in (PHYLUM..=GENUS).rev() {
            let sibling = maxima[rank - PHYLUM].load(Ordering::Relaxed);
            nearest = nearest.max(sibling);
            if !is_named(&taxonomy[rank]) {
                continue;
            }
            let radius = radii.entry(&taxonomy[rank]).or_insert(TaxonRadii {
                rank,
                sibling: None,
                nearest: None,
            });
            lowest(&mut radius.sibling, value(sibling));
            lowest(&mut radius.nearest, value(nearest));
        }
    }
    radii
        .into_iter()
        .map(|(taxon, radius)| (taxon.to_owned(), radius))
        .collect()
}

/// The radius the rank call uses for each taxon: its own sibling radius, or
/// its nearest radius when it has a single child taxon or no comparable
/// sibling pairs, or else the median sibling radius of its rank over the
/// taxa that have their own.
pub(crate) fn resolve_radii(radii: &BTreeMap<String, TaxonRadii>) -> BTreeMap<String, f64> {
    let mut by_rank: HashMap<usize, Vec<f64>> = HashMap::new();
    for radius in radii.values() {
        if let Some(sibling) = radius.sibling {
            by_rank.entry(radius.rank).or_default().push(sibling);
        }
    }
    let medians = by_rank
        .into_iter()
        .map(|(rank, mut values)| {
            values.sort_by(f64::total_cmp);
            let middle = values.len() / 2;
            let median = if values.len() % 2 == 1 {
                values[middle]
            } else {
                (values[middle - 1] + values[middle]) / 2.0
            };
            (rank, median)
        })
        .collect::<HashMap<_, _>>();
    radii
        .iter()
        .filter_map(|(taxon, radius)| {
            let resolved = radius
                .sibling
                .or(radius.nearest)
                .or_else(|| medians.get(&radius.rank).copied())?;
            Some((taxon.clone(), resolved))
        })
        .collect()
}

pub(crate) fn write_marker_radii(
    dir: &Path,
    radii: &BTreeMap<String, TaxonRadii>,
    resolved: &BTreeMap<String, f64>,
) -> std::io::Result<()> {
    let value = |radius: Option<f64>| radius.map_or_else(String::new, |radius| radius.to_string());
    let mut output = BufWriter::new(File::create(dir.join(MARKER_RADII_FILE))?);
    writeln!(output, "taxon\tradius\tsibling_radius\tnearest_radius")?;
    for (taxon, radius) in radii {
        writeln!(
            output,
            "{taxon}\t{}\t{}\t{}",
            value(resolved.get(taxon).copied()),
            value(radius.sibling),
            value(radius.nearest)
        )?;
    }
    output.flush()
}

/// The resolved radius of each taxon, as `resolve_radii` gives it.
pub(crate) fn read_marker_radii(dir: &Path) -> Result<HashMap<String, f64>, Box<dyn Error>> {
    let path = dir.join(MARKER_RADII_FILE);
    let lines = read_lines(&path).map_err(|error| {
        format!("{error}; rebuild the database with this aeon classify-db")
    })?;
    let mut radii = HashMap::new();
    for line in lines.skip(1) {
        let line = line?;
        let mut fields = line.split('\t');
        let (Some(taxon), Some(radius)) = (fields.next(), fields.next()) else {
            return Err(format!("malformed line in {}: {line}", path.display()).into());
        };
        if !radius.is_empty() {
            radii.insert(taxon.to_owned(), radius.parse()?);
        }
    }
    Ok(radii)
}

/// A reference comparable with the query over the marker alignment.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct MarkerMatch {
    /// Position in the marker set's alignment.
    pub(crate) row: usize,
    /// Percent identity over the overlap.
    pub(crate) identity: f64,
    pub(crate) overlap: u32,
}

/// Every reference comparable with `query`, in alignment order.
pub(crate) fn match_references(query: &[u8], alignment: &MarkerAlignment) -> Vec<MarkerMatch> {
    (0..alignment.len())
        .into_par_iter()
        .filter_map(|row| {
            let (identity, overlap) = percent_identity(query, alignment.row(row))?;
            Some(MarkerMatch {
                row,
                identity,
                overlap,
            })
        })
        .collect()
}

/// The marker call for one query.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct MarkerCall {
    /// The best match: highest identity, first in alignment order.
    pub(crate) best: MarkerMatch,
    /// Matches tied with the best, itself included, in alignment order.
    pub(crate) tied: Vec<MarkerMatch>,
    /// The called lineage from the domain down.
    pub(crate) lineage: Vec<String>,
}

/// Call the rank from a query's matches (in alignment order) and the
/// references' taxonomy: take the lineage shared by the references tied with
/// the best match, from the domain down to the genus, and keep each rank from
/// the phylum while the best identity reaches the taxon's radius, stopping at
/// the first it misses.
pub(crate) fn call_rank<'a>(
    matches: &[MarkerMatch],
    taxonomy: impl Fn(usize) -> &'a [String],
    radii: &HashMap<String, f64>,
) -> Option<MarkerCall> {
    let best = matches
        .iter()
        .reduce(|best, hit| if hit.identity > best.identity { hit } else { best })?
        .clone();
    let variance = |hit: &MarkerMatch| {
        let p = hit.identity / 100.0;
        p * (1.0 - p) / f64::from(hit.overlap)
    };
    let tied = matches
        .iter()
        .filter(|hit| {
            let bound = 100.0 * TIE_STANDARD_ERRORS * (variance(&best) + variance(hit)).sqrt();
            hit.identity >= best.identity - bound
        })
        .cloned()
        .collect::<Vec<_>>();
    let best_taxonomy = taxonomy(best.row);
    let mut shared = GENUS + 1;
    for hit in &tied {
        let other = taxonomy(hit.row);
        shared = (0..shared)
            .take_while(|&rank| other[rank] == best_taxonomy[rank])
            .count();
    }
    let depth = (PHYLUM..shared)
        .take_while(|&rank| {
            is_named(&best_taxonomy[rank])
                && radii
                    .get(&best_taxonomy[rank])
                    .is_some_and(|radius| best.identity >= *radius)
        })
        .count();
    Some(MarkerCall {
        lineage: best_taxonomy[..=depth].to_vec(),
        best,
        tied,
    })
}

/// The `count` closest matches and every match tied with the best, closest
/// first, each with whether it is tied.
pub(crate) fn listed_matches<'m>(
    matches: &'m [MarkerMatch],
    call: &MarkerCall,
    count: usize,
) -> Vec<(&'m MarkerMatch, bool)> {
    let by_identity =
        |left: &&MarkerMatch, right: &&MarkerMatch| right.identity.total_cmp(&left.identity);
    let tied = call.tied.iter().map(|hit| hit.row).collect::<HashSet<_>>();
    let mut ranked = matches.iter().collect::<Vec<_>>();
    if count > 0 && ranked.len() > count {
        ranked.select_nth_unstable_by(count - 1, by_identity);
    }
    let (closest, rest) = ranked.split_at(count.min(ranked.len()));
    let mut listed = closest.to_vec();
    listed.extend(rest.iter().filter(|hit| tied.contains(&hit.row)));
    listed.sort_by(by_identity);
    listed
        .into_iter()
        .map(|hit| (hit, tied.contains(&hit.row)))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn lineage(taxa: &str) -> Vec<String> {
        taxa.split(';').map(str::to_owned).collect()
    }

    /// A row of `length` columns: `same` copies of A, then C.
    fn row(same: usize, length: usize) -> Vec<u8> {
        (0..length).map(|column| if column < same { b'A' } else { b'C' }).collect()
    }

    #[test]
    fn identity_counts_standard_residues_shared_by_both() {
        let left = encode_row(b"AC-DxFB*");
        let right = encode_row(b"ACG-XYBA");
        assert_eq!(identity_counts(&left, &right), (2, 3));
        let long = vec![b'A'; 1000];
        assert_eq!(identity_counts(&long, &long), (1000, 1000));
        assert_eq!(percent_identity(&long[..499], &long[..499]), None);
    }

    #[test]
    fn alignment_rows_round_trip_encoded() {
        let dir = std::env::temp_dir().join(format!("aeon-marker-msa-{}", std::process::id()));
        fs::create_dir_all(&dir).unwrap();
        let rows = [
            ("GCF_1.1".to_owned(), encode_row(b"ACx-D")),
            ("GCA_2.1".to_owned(), encode_row(b"-CGYW")),
        ];
        write_marker_alignment(&dir, "bac120", &rows).unwrap();
        write_marker_alignment(&dir, "ar53", &[]).unwrap();
        let alignment = MarkerAlignment::read(&dir, "bac120").unwrap();
        let empty = MarkerAlignment::read(&dir, "ar53").unwrap();
        fs::remove_dir_all(&dir).unwrap();
        assert_eq!(alignment.accessions, ["GCF_1.1", "GCA_2.1"]);
        assert_eq!((alignment.row(0), alignment.row(1)), (&rows[0].1[..], &rows[1].1[..]));
        assert_eq!(empty.len(), 0);
    }

    #[test]
    fn radii_take_the_lowest_per_member_sibling_identity() {
        let g1a = lineage("d__B;p__P;c__C;o__O;f__F;g__G1;s__a");
        let g1b = lineage("d__B;p__P;c__C;o__O;f__F;g__G1;s__b");
        let g2 = lineage("d__B;p__P;c__C;o__O;f__F;g__G2;s__c");
        let rows = [row(1000, 1000), row(800, 1000), row(700, 1000)];
        let refs = [(&g1a[..], &rows[0][..]), (&g1b[..], &rows[1][..]), (&g2[..], &rows[2][..])];
        let radii = taxon_radii(&refs);
        let pair = |taxon: &str| (radii[taxon].sibling, radii[taxon].nearest);
        // g1a-g1b 80%, g1a-g2 70%, g1b-g2 90%.
        assert_eq!(pair("g__G1"), (Some(80.0), Some(80.0)));
        assert_eq!(pair("f__F"), (Some(70.0), Some(80.0)));
        assert_eq!(pair("o__O"), (None, Some(80.0)));
        assert_eq!(pair("g__G2"), (None, None));
        assert_eq!(radii["g__G2"].rank, 5);
        let resolved = resolve_radii(&radii);
        assert_eq!(resolved["o__O"], 80.0);
        // The only genus with its own sibling radius gives the genus median.
        assert_eq!(resolved["g__G2"], 80.0);
    }

    #[test]
    fn rank_walk_stops_at_the_first_missed_radius() {
        let radii = HashMap::from([
            ("p__P".to_owned(), 40.0),
            ("c__C".to_owned(), 50.0),
            ("o__O".to_owned(), 60.0),
            ("f__F".to_owned(), 70.0),
            ("g__G".to_owned(), 85.0),
            ("g__H".to_owned(), 85.0),
        ]);
        let lineages = ["g__G", "g__H", "g__G"]
            .map(|genus| lineage(&format!("d__B;p__P;c__C;o__O;f__F;{genus};s__x")));
        let taxonomy = |row: usize| &lineages[row][..];
        let hit = |row: usize, identity: f64| MarkerMatch {
            row,
            identity,
            overlap: 5000,
        };
        let call = call_rank(&[hit(0, 75.0)], taxonomy, &radii).unwrap();
        assert_eq!(call.lineage, lineage("d__B;p__P;c__C;o__O;f__F"));
        let call = call_rank(&[hit(0, 65.0)], taxonomy, &radii).unwrap();
        assert_eq!(call.lineage, lineage("d__B;p__P;c__C;o__O"));
        // A tie across genera stops at their family.
        let call = call_rank(&[hit(1, 89.5), hit(0, 90.0)], taxonomy, &radii).unwrap();
        assert_eq!((call.best.row, call.lineage.len(), call.tied.len()), (0, 5, 2));
        let matches = [hit(0, 90.0), hit(1, 85.0)];
        let call = call_rank(&matches, taxonomy, &radii).unwrap();
        assert_eq!(call.lineage, lineage("d__B;p__P;c__C;o__O;f__F;g__G"));
        assert!(call_rank(&[], taxonomy, &radii).is_none());
        // Listing: the closest match, plus a tied one beyond the count.
        let matches = [hit(2, 80.0), hit(1, 90.0), hit(0, 89.0)];
        let call = call_rank(&matches, taxonomy, &radii).unwrap();
        let listed = listed_matches(&matches, &call, 1)
            .into_iter()
            .map(|(hit, tied)| (hit.identity, tied))
            .collect::<Vec<_>>();
        assert_eq!(listed, [(90.0, true), (89.0, true)]);
    }
}
