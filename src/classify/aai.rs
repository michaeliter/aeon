//! Tier 2: amino-acid identity. Amino-acid sketch search proposes candidate
//! reference proteomes, DIAMOND reciprocal best hits give AAI, and the
//! closest reference's genus or family is assigned when it leads every other
//! genus or family by a margin.
//!
//! AAI comes from DIAMOND `blastp` in both directions with DIAMOND's default
//! e-value: each protein's first hit of at least 20% identity is its best
//! hit, and AAI is the unweighted mean identity of reciprocal best hits,
//! reported only when there are at least 50 of them.

use std::collections::{HashMap, HashSet};
use std::error::Error;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;
use std::process::Command;

use needletail::parse_fastx_file;
use rayon::prelude::*;

use crate::classify::external::run;
use crate::classify::genes::write_fasta_record;
use crate::classify::refdata::GtdbtkData;

/// DIAMOND threads per query job; small searches scale poorly past a few.
const THREADS_PER_JOB: usize = 2;
/// Genome tag of the query's proteins, distinct from every reference ID.
const QUERY: u32 = u32::MAX;
/// DIAMOND's default e-value cutoff, applied as if each reference proteome
/// were its own database, as in a search of one genome pair.
const MAX_EVALUE: f64 = 0.001;
/// Best hits need this identity (%) ...
const MIN_IDENTITY: f64 = 20.0;
/// ... and an AAI needs this many reciprocal best hits.
const MIN_RECIPROCAL_HITS: usize = 50;

#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct AaiHit {
    /// Reference ID shared with the nucleotide index.
    pub(crate) reference: u32,
    /// Mean identity of reciprocal best hits, 0-100.
    pub(crate) aai: f64,
    /// Population standard deviation of those identities.
    pub(crate) sd: f64,
    pub(crate) reciprocal_hits: usize,
    /// Protein count of the smaller of the two proteomes.
    pub(crate) min_proteins: usize,
}

/// One DIAMOND tabular hit: (query protein, subject protein, identity, bitscore).
#[derive(Debug, Clone, PartialEq)]
struct ProteinHit {
    query: String,
    subject: String,
    identity: f64,
    bitscore: f64,
}

/// Tagged protein IDs ("{genome}|{index}") of one or more proteomes.
struct Proteomes {
    /// Tagged ID -> residues.
    proteins: HashMap<String, Vec<u8>>,
    /// Genome -> (protein count, residue count).
    sizes: HashMap<u32, (usize, usize)>,
}

/// AAI between each query proteome and its candidate reference proteomes,
/// sorted by descending AAI. Each query is searched against its own
/// candidates only, so small DIAMOND jobs run side by side with a few
/// threads each.
pub(crate) fn compute_aai(
    diamond: &Path,
    work_dir: &Path,
    query_proteomes: &[&Path],
    candidates: &[Vec<u32>],
    reference_proteomes: &[String],
    threads: usize,
) -> Result<Vec<Vec<AaiHit>>, Box<dyn Error>> {
    assert_eq!(query_proteomes.len(), candidates.len());
    fs::create_dir_all(work_dir)?;
    let threads = threads.max(1);
    let jobs = query_proteomes
        .len()
        .clamp(1, (threads / THREADS_PER_JOB).max(1));
    let pool = rayon::ThreadPoolBuilder::new().num_threads(jobs).build()?;
    pool.install(|| {
        query_proteomes
            .par_iter()
            .zip(candidates)
            .enumerate()
            .map(|(query, (proteome, references))| {
                // Workers split the threads evenly, the first few taking the
                // remainder.
                let worker = rayon::current_thread_index().unwrap_or(0);
                let job_threads = threads / jobs + usize::from(worker < threads % jobs);
                let dir = work_dir.join(format!("query_{query}"));
                let result = query_aai(
                    diamond,
                    &dir,
                    proteome,
                    references,
                    reference_proteomes,
                    job_threads,
                );
                let _ = fs::remove_dir_all(&dir);
                result.map_err(|error| {
                    format!("AAI search for {} failed: {error}", proteome.display())
                })
            })
            .collect::<Result<Vec<_>, _>>()
    })
    .map_err(Into::into)
}

/// Reciprocal-best-hit AAI of one query proteome against its candidates.
fn query_aai(
    diamond: &Path,
    dir: &Path,
    proteome: &Path,
    candidates: &[u32],
    reference_proteomes: &[String],
    threads: usize,
) -> Result<Vec<AaiHit>, Box<dyn Error>> {
    let mut references = candidates.to_vec();
    references.sort_unstable();
    references.dedup();
    if references.is_empty() {
        return Ok(Vec::new());
    }
    fs::create_dir_all(dir)?;
    // Protein IDs become "{genome}|{protein}" so hits map back to genomes;
    // references keep their database IDs.
    let query_fasta = dir.join("query.faa");
    let reference_fasta = dir.join("references.faa");
    let query = write_tagged_proteins(&[(QUERY, proteome)], &query_fasta)?;
    let reference_files = references
        .iter()
        .map(|&id| (id, Path::new(&reference_proteomes[id as usize])))
        .collect::<Vec<_>>();
    let reference = write_tagged_proteins(&reference_files, &reference_fasta)?;
    let Some(&(query_count, _)) = query.sizes.get(&QUERY) else {
        return Ok(Vec::new());
    };
    // All candidates share one database. Its size is pinned to their total
    // length and each hit's e-value rescaled to its own proteome's length,
    // reproducing the e-values of one database per reference.
    let total_letters = reference
        .sizes
        .values()
        .map(|&(_, letters)| letters)
        .sum::<usize>();
    let smallest = reference
        .sizes
        .values()
        .map(|&(_, letters)| letters)
        .min()
        .unwrap_or(1);
    let evalue_scale = |hit: &ProteinHit| {
        reference.sizes[&genome_of(&hit.subject)].1 as f64 / total_letters as f64
    };
    let forward = blastp(
        diamond,
        &[
            "--dbsize",
            &total_letters.to_string(),
            "--evalue",
            &(MAX_EVALUE * total_letters as f64 / smallest.max(1) as f64).to_string(),
            // Every reference needs its own best hit, not the 25 best overall.
            "--max-target-seqs",
            "0",
        ],
        &dir.join("forward.tsv"),
        &query_fasta,
        &reference_fasta,
        threads,
        |hit, evalue| evalue * evalue_scale(hit) <= MAX_EVALUE,
    )?;
    let forward_best = best_hits(&forward);
    // Reciprocity only needs the reverse best hits of reference proteins that
    // are some query protein's best hit, so only those are searched back
    // against the whole query proteome.
    let best_subjects = forward_best
        .values()
        .map(|hit| hit.subject.as_str())
        .collect::<HashSet<_>>();
    if best_subjects.is_empty() {
        return Ok(Vec::new());
    }
    let reverse_fasta = dir.join("best_hits.faa");
    write_selected_proteins(&reference.proteins, &best_subjects, &reverse_fasta)?;
    let reverse = blastp(
        diamond,
        &[],
        &dir.join("reverse.tsv"),
        &reverse_fasta,
        &query_fasta,
        threads,
        |_, _| true,
    )?;
    let pairs = reciprocal_aai(&forward_best, &reverse);
    let mut hits = references
        .into_iter()
        .filter_map(|id| {
            let identity = pairs.get(&id)?;
            (identity.count >= MIN_RECIPROCAL_HITS).then(|| AaiHit {
                reference: id,
                aai: identity.mean(),
                sd: identity.sd(),
                reciprocal_hits: identity.count,
                min_proteins: query_count.min(reference.sizes[&id].0),
            })
        })
        .collect::<Vec<_>>();
    hits.sort_by(|left, right| {
        right
            .aai
            .total_cmp(&left.aai)
            .then(right.reciprocal_hits.cmp(&left.reciprocal_hits))
            .then(left.reference.cmp(&right.reference))
    });
    Ok(hits)
}

/// Concatenate proteomes with `{genome}|{index}` IDs, keeping stop codons.
fn write_tagged_proteins(
    proteomes: &[(u32, &Path)],
    output: &Path,
) -> Result<Proteomes, Box<dyn Error>> {
    let mut writer = BufWriter::new(File::create(output)?);
    let mut proteins = HashMap::new();
    let mut sizes = HashMap::new();
    for &(genome, path) in proteomes {
        let mut reader = parse_fastx_file(path)
            .map_err(|error| format!("invalid protein FASTA {}: {error}", path.display()))?;
        let (mut count, mut letters) = (0usize, 0usize);
        while let Some(record) = reader.next() {
            let record = record?;
            let sequence = record.seq();
            if sequence.iter().all(|&residue| residue == b'*') {
                continue;
            }
            let id = format!("{genome}|{count}");
            count += 1;
            letters += sequence.len();
            write_fasta_record(&mut writer, id.as_bytes(), &sequence)?;
            proteins.insert(id, sequence.into_owned());
        }
        if count > 0 {
            sizes.insert(genome, (count, letters));
        }
    }
    writer.flush()?;
    Ok(Proteomes { proteins, sizes })
}

/// Write the proteins named in `ids`, in ID order.
fn write_selected_proteins(
    proteins: &HashMap<String, Vec<u8>>,
    ids: &HashSet<&str>,
    output: &Path,
) -> Result<(), Box<dyn Error>> {
    let mut writer = BufWriter::new(File::create(output)?);
    let mut ids = ids.iter().copied().collect::<Vec<_>>();
    ids.sort_unstable();
    for id in ids {
        write_fasta_record(&mut writer, id.as_bytes(), &proteins[id])?;
    }
    writer.flush()?;
    Ok(())
}

/// DIAMOND `blastp` of `queries` against the FASTA `subjects`
/// (DIAMOND builds the database in memory), keeping hits of at least
/// `MIN_IDENTITY` that `keep(hit, evalue)` accepts.
fn blastp(
    diamond: &Path,
    options: &[&str],
    output: &Path,
    queries: &Path,
    subjects: &Path,
    threads: usize,
    keep: impl Fn(&ProteinHit, f64) -> bool,
) -> Result<Vec<ProteinHit>, Box<dyn Error>> {
    run(Command::new(diamond)
        .args(["blastp", "--quiet", "--threads", &threads.to_string()])
        .args(options)
        .args([
            "--outfmt", "6", "qseqid", "sseqid", "pident", "bitscore", "evalue",
        ])
        .arg("-q")
        .arg(queries)
        .arg("-d")
        .arg(subjects)
        .arg("-o")
        .arg(output))?;
    let mut hits = Vec::new();
    for line in BufReader::new(File::open(output)?).lines() {
        let line = line?;
        let fields = line.split('\t').collect::<Vec<_>>();
        if fields.len() < 5 {
            continue;
        }
        let hit = ProteinHit {
            query: fields[0].to_owned(),
            subject: fields[1].to_owned(),
            identity: fields[2].parse()?,
            bitscore: fields[3].parse()?,
        };
        if hit.identity >= MIN_IDENTITY && keep(&hit, fields[4].parse()?) {
            hits.push(hit);
        }
    }
    Ok(hits)
}

fn genome_of(protein: &str) -> u32 {
    protein
        .split_once('|')
        .and_then(|(genome, _)| genome.parse().ok())
        .expect("tagged protein ID")
}

/// Best hit of each protein in each subject genome: (protein, genome) -> hit.
/// Ties keep DIAMOND's first hit.
fn best_hits(hits: &[ProteinHit]) -> HashMap<(&str, u32), &ProteinHit> {
    let mut best: HashMap<(&str, u32), &ProteinHit> = HashMap::new();
    for hit in hits {
        let key = (hit.query.as_str(), genome_of(&hit.subject));
        match best.get(&key) {
            Some(current) if current.bitscore >= hit.bitscore => {}
            _ => {
                best.insert(key, hit);
            }
        }
    }
    best
}

/// Running sums of reciprocal-best-hit identities.
#[derive(Debug, Default, Clone, Copy, PartialEq)]
struct Identities {
    sum: f64,
    sum_squares: f64,
    count: usize,
}

impl Identities {
    fn mean(&self) -> f64 {
        self.sum / self.count as f64
    }

    /// Population standard deviation.
    fn sd(&self) -> f64 {
        (self.sum_squares / self.count as f64 - self.mean().powi(2))
            .max(0.0)
            .sqrt()
    }
}

/// Identities of reciprocal best hits per reference, from the query's
/// forward best hits and the reverse search.
fn reciprocal_aai(
    forward_best: &HashMap<(&str, u32), &ProteinHit>,
    reverse: &[ProteinHit],
) -> HashMap<u32, Identities> {
    let reverse_best = best_hits(reverse);
    // The reverse search's identity is averaged; the forward one is the
    // same alignment up to DIAMOND's heuristics.
    let mut reciprocal = forward_best
        .iter()
        .filter_map(|(&(query_protein, reference), hit)| {
            let back = reverse_best.get(&(hit.subject.as_str(), genome_of(query_protein)))?;
            (back.subject == query_protein).then_some((reference, query_protein, back.identity))
        })
        .collect::<Vec<_>>();
    // Sum in a fixed order, not hash order, so rounded AAIs are reproducible.
    reciprocal.sort_unstable_by(|left, right| (left.0, left.1).cmp(&(right.0, right.1)));
    let mut sums: HashMap<u32, Identities> = HashMap::new();
    for (reference, _, identity) in reciprocal {
        let entry = sums.entry(reference).or_default();
        entry.sum += identity;
        entry.sum_squares += identity.powi(2);
        entry.count += 1;
    }
    sums
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) enum AaiDecision {
    Genus(AaiCall),
    Family(AaiCall),
    Unresolved,
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct AaiCall {
    pub(crate) taxonomy: String,
    /// AAI to the closest reference, which gives the taxon.
    pub(crate) aai: f64,
    /// Points by which it beats the closest reference of another taxon at
    /// the same rank; `None` when every hit shares the taxon.
    pub(crate) lead: Option<f64>,
}

/// Assign the closest reference's genus (or else family) when its AAI
/// reaches the rank's threshold and beats every other genus (family) among
/// the hits by at least `margin` points. Unlike a vote, the margin does not
/// depend on how many references a taxon has. `hits` are sorted by
/// descending AAI.
pub(crate) fn decide(
    hits: &[AaiHit],
    accessions: &[String],
    data: &GtdbtkData,
    genus_aai: f64,
    family_aai: f64,
    margin: f64,
) -> AaiDecision {
    let ranked = hits
        .iter()
        .filter_map(|hit| {
            data.taxonomy_of(&accessions[hit.reference as usize])
                .map(|taxonomy| (hit.aai, taxonomy))
        })
        .collect::<Vec<_>>();
    let call = |threshold: f64, ranks: usize| {
        let ((aai, taxonomy), rest) = ranked.split_first()?;
        if *aai < threshold {
            return None;
        }
        let taxon = &taxonomy[..ranks];
        // A reference without a name at this rank cannot give one.
        if taxon[ranks - 1].len() <= 3 {
            return None;
        }
        let lead = rest
            .iter()
            .find(|(_, other)| other[..ranks] != *taxon)
            .map(|(other, _)| aai - other);
        lead.is_none_or(|lead| lead >= margin).then(|| AaiCall {
            taxonomy: taxon.join(";"),
            aai: *aai,
            lead,
        })
    };
    if let Some(call) = call(genus_aai, 6) {
        return AaiDecision::Genus(call);
    }
    call(family_aai, 5).map_or(AaiDecision::Unresolved, AaiDecision::Family)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn hit(query: &str, subject: &str, identity: f64, bitscore: f64) -> ProteinHit {
        ProteinHit {
            query: query.into(),
            subject: subject.into(),
            identity,
            bitscore,
        }
    }

    #[test]
    fn aai_averages_reciprocal_best_hits_only() {
        // Query genome 0 vs reference genome 7.
        let forward = [
            hit("0|0", "7|0", 80.0, 300.0),
            hit("0|0", "7|1", 60.0, 200.0),
            hit("0|1", "7|1", 70.0, 250.0),
            hit("0|2", "7|2", 50.0, 100.0), // not reciprocal
        ];
        let reverse = [
            hit("7|0", "0|0", 80.0, 300.0),
            hit("7|1", "0|1", 70.0, 250.0),
            hit("7|2", "0|1", 40.0, 150.0),
        ];
        let pairs = reciprocal_aai(&best_hits(&forward), &reverse);
        assert_eq!(pairs[&7].count, 2);
        assert_eq!(pairs[&7].mean(), 75.0);
        assert_eq!(pairs[&7].sd(), 5.0);
    }

    #[test]
    fn selected_proteins_keep_tagged_ids() {
        let dir = std::env::temp_dir().join(format!("aeon-aai-{}", std::process::id()));
        fs::create_dir_all(&dir).unwrap();
        let proteome = dir.join("ref.faa");
        fs::write(&proteome, ">a\nMKV*\n>b\n*\n>c\nMLL\n").unwrap();
        let tagged = dir.join("tagged.faa");
        let proteomes = write_tagged_proteins(&[(7, &proteome)], &tagged).unwrap();
        assert_eq!(
            fs::read_to_string(&tagged).unwrap(),
            ">7|0\nMKV*\n>7|1\nMLL\n"
        );
        assert_eq!(proteomes.sizes[&7], (2, 7));
        let selected = dir.join("selected.faa");
        write_selected_proteins(&proteomes.proteins, &HashSet::from(["7|1"]), &selected).unwrap();
        assert_eq!(fs::read_to_string(&selected).unwrap(), ">7|1\nMLL\n");
        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn margin_decides_genus_then_family() {
        let mut data = GtdbtkData::default();
        let taxa = [
            ("GCF_000000001.1", "f__F;g__A"),
            ("GCF_000000002.1", "f__F;g__B"),
            ("GCF_000000003.1", "f__G;g__C"),
        ];
        for (accession, family_genus) in taxa {
            data.taxonomy.insert(
                crate::classify::refdata::canonical_gid(accession),
                format!("d__Bacteria;p__P;c__C;o__O;{family_genus};s__")
                    .split(';')
                    .map(str::to_owned)
                    .collect(),
            );
        }
        let accessions = taxa.iter().map(|(a, _)| a.to_string()).collect::<Vec<_>>();
        let aai = |reference, aai| AaiHit {
            reference,
            aai,
            sd: 10.0,
            reciprocal_hits: 500,
            min_proteins: 2000,
        };
        let decide = |hits: &[AaiHit]| decide(hits, &accessions, &data, 65.0, 45.0, 5.0);
        let call = |taxonomy: &str, aai, lead| AaiCall {
            taxonomy: format!("d__Bacteria;p__P;c__C;o__O;{taxonomy}"),
            aai,
            lead,
        };
        // A genus with one reference still wins on its lead.
        assert_eq!(
            decide(&[aai(0, 80.0), aai(1, 70.0), aai(1, 69.0), aai(2, 60.0)]),
            AaiDecision::Genus(call("f__F;g__A", 80.0, Some(10.0)))
        );
        // Too close to genus B for a genus, but well ahead of family G.
        assert_eq!(
            decide(&[aai(0, 72.0), aai(1, 70.0), aai(2, 60.0)]),
            AaiDecision::Family(call("f__F", 72.0, Some(12.0)))
        );
        assert_eq!(
            decide(&[aai(0, 72.0), aai(2, 70.0)]),
            AaiDecision::Unresolved
        );
        // An unopposed taxon needs only the threshold.
        assert_eq!(
            decide(&[aai(0, 60.0)]),
            AaiDecision::Family(call("f__F", 60.0, None))
        );
        assert_eq!(decide(&[aai(0, 40.0)]), AaiDecision::Unresolved);
    }
}
