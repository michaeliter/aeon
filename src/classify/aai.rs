//! Tier 2: amino-acid identity. Amino-acid sketch search proposes candidate
//! reference proteomes, DIAMOND reciprocal best hits give AAI, and a
//! majority vote among close references assigns a genus or family.

use std::collections::{HashMap, HashSet};
use std::error::Error;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;
use std::process::Command;

use needletail::parse_fastx_file;
use rayon::prelude::*;

use crate::classify::external::run;
use crate::classify::refdata::GtdbtkData;

/// DIAMOND threads per query job; small searches scale poorly past a few.
const THREADS_PER_JOB: usize = 2;
/// Fixed DIAMOND database size (letters), roughly ten candidate proteomes, so
/// e-values do not depend on how many references a query was searched against.
const DBSIZE: &str = "10000000";
/// Genome tag of the query's proteins, distinct from every reference ID.
const QUERY: u32 = u32::MAX;
/// Reciprocal best hits must reach this identity (%) ...
const MIN_IDENTITY: f64 = 30.0;
/// ... and align this fraction of the shorter protein (MiGA/EzAAI-style).
const MIN_ALIGNED_FRACTION: f64 = 0.7;

#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct AaiHit {
    /// Reference ID shared with the nucleotide index.
    pub(crate) reference: u32,
    /// Mean identity of reciprocal best hits, 0-100.
    pub(crate) aai: f64,
    pub(crate) reciprocal_hits: usize,
}

/// One DIAMOND tabular hit: (query protein, subject protein, identity, bitscore).
#[derive(Debug, Clone, PartialEq)]
struct ProteinHit {
    query: String,
    subject: String,
    identity: f64,
    bitscore: f64,
}

/// Tagged protein ID ("{genome}|{index}") -> residues.
type Proteins = HashMap<String, Vec<u8>>;

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
    let query_proteins = write_tagged_proteins(&[(QUERY, proteome)], &query_fasta)?;
    let reference_files = references
        .iter()
        .map(|&id| (id, Path::new(&reference_proteomes[id as usize])))
        .collect::<Vec<_>>();
    let reference_proteins = write_tagged_proteins(&reference_files, &reference_fasta)?;
    let forward = blastp(
        diamond,
        dir,
        "forward",
        &query_fasta,
        &reference_fasta,
        threads,
        &query_proteins,
        &reference_proteins,
    )?;
    let forward_best = best_hits(&forward);
    // Reciprocity only needs the reverse best hits of reference proteins that
    // are some query protein's best hit, so only those are searched back.
    let best_subjects = forward_best
        .values()
        .map(|hit| hit.subject.as_str())
        .collect::<HashSet<_>>();
    if best_subjects.is_empty() {
        return Ok(Vec::new());
    }
    let reverse_fasta = dir.join("best_hits.faa");
    write_selected_proteins(&reference_proteins, &best_subjects, &reverse_fasta)?;
    let reverse = blastp(
        diamond,
        dir,
        "reverse",
        &reverse_fasta,
        &query_fasta,
        threads,
        &reference_proteins,
        &query_proteins,
    )?;
    let pairs = reciprocal_aai(&forward_best, &reverse);
    let mut hits = references
        .into_iter()
        .filter_map(|reference| {
            pairs.get(&reference).map(|&(aai, reciprocal_hits)| AaiHit {
                reference,
                aai,
                reciprocal_hits,
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

/// Concatenate proteomes with `{genome}|{index}` IDs, dropping stop codons.
/// Returns the written proteins.
fn write_tagged_proteins(
    proteomes: &[(u32, &Path)],
    output: &Path,
) -> Result<Proteins, Box<dyn Error>> {
    let mut writer = BufWriter::new(File::create(output)?);
    let mut proteins = HashMap::new();
    for &(genome, path) in proteomes {
        let mut reader = parse_fastx_file(path)
            .map_err(|error| format!("invalid protein FASTA {}: {error}", path.display()))?;
        let mut index = 0usize;
        while let Some(record) = reader.next() {
            let record = record?;
            let sequence = record
                .seq()
                .iter()
                .copied()
                .filter(|&residue| residue != b'*')
                .collect::<Vec<_>>();
            if sequence.is_empty() {
                continue;
            }
            let id = format!("{genome}|{index}");
            index += 1;
            write_record(&mut writer, &id, &sequence)?;
            proteins.insert(id, sequence);
        }
    }
    writer.flush()?;
    Ok(proteins)
}

/// Write the proteins named in `ids`, in ID order.
fn write_selected_proteins(
    proteins: &Proteins,
    ids: &HashSet<&str>,
    output: &Path,
) -> Result<(), Box<dyn Error>> {
    let mut writer = BufWriter::new(File::create(output)?);
    let mut ids = ids.iter().copied().collect::<Vec<_>>();
    ids.sort_unstable();
    for id in ids {
        write_record(&mut writer, id, &proteins[id])?;
    }
    writer.flush()?;
    Ok(())
}

fn write_record(writer: &mut impl Write, id: &str, sequence: &[u8]) -> std::io::Result<()> {
    writeln!(writer, ">{id}")?;
    writer.write_all(sequence)?;
    writer.write_all(b"\n")
}

/// DIAMOND blastp of `queries` against the FASTA `subjects` (DIAMOND builds
/// the database in memory), keeping hits that cover `MIN_ALIGNED_FRACTION`
/// of the shorter protein.
#[allow(clippy::too_many_arguments)]
fn blastp(
    diamond: &Path,
    dir: &Path,
    name: &str,
    queries: &Path,
    subjects: &Path,
    threads: usize,
    query_proteins: &Proteins,
    subject_proteins: &Proteins,
) -> Result<Vec<ProteinHit>, Box<dyn Error>> {
    let output = dir.join(format!("{name}.tsv"));
    run(Command::new(diamond)
        .args(["blastp", "--quiet", "--threads", &threads.to_string()])
        .args(["--id", &MIN_IDENTITY.to_string()])
        .args(["--dbsize", DBSIZE])
        // Report every subject: best hits are needed per subject genome.
        .args(["--max-target-seqs", "0"])
        .args([
            "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "bitscore",
        ])
        .arg("-q")
        .arg(queries)
        .arg("-d")
        .arg(subjects)
        .arg("-o")
        .arg(&output))?;
    let length = |proteins: &Proteins, id: &str| proteins.get(id).map_or(0, Vec::len);
    let mut hits = Vec::new();
    for line in BufReader::new(File::open(&output)?).lines() {
        let line = line?;
        let fields = line.split('\t').collect::<Vec<_>>();
        if fields.len() < 5 {
            continue;
        }
        let shorter = length(query_proteins, fields[0]).min(length(subject_proteins, fields[1]));
        let aligned: usize = fields[3].parse()?;
        if shorter == 0 || (aligned as f64) < MIN_ALIGNED_FRACTION * shorter as f64 {
            continue;
        }
        hits.push(ProteinHit {
            query: fields[0].to_owned(),
            subject: fields[1].to_owned(),
            identity: fields[2].parse()?,
            bitscore: fields[4].parse()?,
        });
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

/// Mean identity and count of reciprocal best hits per reference, from the
/// query's forward best hits and the reverse search.
fn reciprocal_aai(
    forward_best: &HashMap<(&str, u32), &ProteinHit>,
    reverse: &[ProteinHit],
) -> HashMap<u32, (f64, usize)> {
    let reverse_best = best_hits(reverse);
    let mut sums: HashMap<u32, (f64, usize)> = HashMap::new();
    for (&(query_protein, reference), hit) in forward_best {
        let reciprocal = reverse_best
            .get(&(hit.subject.as_str(), genome_of(query_protein)))
            .is_some_and(|back| back.subject == query_protein);
        if reciprocal {
            let entry = sums.entry(reference).or_default();
            entry.0 += hit.identity;
            entry.1 += 1;
        }
    }
    sums.into_iter()
        .map(|(reference, (sum, count))| (reference, (sum / count as f64, count)))
        .collect()
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) enum AaiDecision {
    Genus {
        taxonomy: String,
        supporting: usize,
        considered: usize,
    },
    Family {
        taxonomy: String,
        supporting: usize,
        considered: usize,
    },
    Unresolved,
}

/// Among the top `top_n` hits, assign the best hit's genus (or else family)
/// when it is the strict-majority taxon among hits at or above the threshold.
pub(crate) fn decide(
    hits: &[AaiHit],
    accessions: &[String],
    data: &GtdbtkData,
    genus_aai: f64,
    family_aai: f64,
    top_n: usize,
) -> AaiDecision {
    let top = &hits[..hits.len().min(top_n)];
    let consensus = |threshold: f64, ranks: usize| {
        let taxa = top
            .iter()
            .filter(|hit| hit.aai >= threshold)
            .filter_map(|hit| {
                data.taxonomy_of(&accessions[hit.reference as usize])
                    .map(|taxonomy| taxonomy[..ranks].join(";"))
            })
            .collect::<Vec<_>>();
        let best = taxa.first()?.clone();
        let supporting = taxa.iter().filter(|taxon| **taxon == best).count();
        (supporting * 2 > taxa.len()).then_some((best, supporting, taxa.len()))
    };
    if let Some((taxonomy, supporting, considered)) = consensus(genus_aai, 6) {
        return AaiDecision::Genus {
            taxonomy,
            supporting,
            considered,
        };
    }
    if let Some((taxonomy, supporting, considered)) = consensus(family_aai, 5) {
        return AaiDecision::Family {
            taxonomy,
            supporting,
            considered,
        };
    }
    AaiDecision::Unresolved
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
        assert_eq!(pairs[&7], (75.0, 2));
    }

    #[test]
    fn selected_proteins_keep_tagged_ids() {
        let dir = std::env::temp_dir().join(format!("aeon-aai-{}", std::process::id()));
        fs::create_dir_all(&dir).unwrap();
        let proteome = dir.join("ref.faa");
        fs::write(&proteome, ">a\nMKV*\n>b\n*\n>c\nMLL\n").unwrap();
        let tagged = dir.join("tagged.faa");
        let proteins = write_tagged_proteins(&[(7, &proteome)], &tagged).unwrap();
        assert_eq!(
            fs::read_to_string(&tagged).unwrap(),
            ">7|0\nMKV\n>7|1\nMLL\n"
        );
        let selected = dir.join("selected.faa");
        write_selected_proteins(&proteins, &HashSet::from(["7|1"]), &selected).unwrap();
        assert_eq!(fs::read_to_string(&selected).unwrap(), ">7|1\nMLL\n");
        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn rank_consensus_falls_back_to_family() {
        let mut data = GtdbtkData::default();
        let taxa = [
            ("GCF_000000001.1", "f__F;g__A"),
            ("GCF_000000002.1", "f__F;g__B"),
            ("GCF_000000003.1", "f__F;g__C"),
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
            reciprocal_hits: 500,
        };
        assert_eq!(
            decide(
                &[aai(0, 80.0), aai(0, 76.0), aai(1, 70.0)],
                &accessions,
                &data,
                75.0,
                65.0,
                10
            ),
            AaiDecision::Genus {
                taxonomy: "d__Bacteria;p__P;c__C;o__O;f__F;g__A".into(),
                supporting: 2,
                considered: 2,
            }
        );
        assert_eq!(
            decide(
                &[aai(0, 74.0), aai(1, 70.0), aai(2, 66.0)],
                &accessions,
                &data,
                75.0,
                65.0,
                10
            ),
            AaiDecision::Family {
                taxonomy: "d__Bacteria;p__P;c__C;o__O;f__F".into(),
                supporting: 3,
                considered: 3,
            }
        );
        assert_eq!(
            decide(&[aai(0, 60.0)], &accessions, &data, 75.0, 65.0, 10),
            AaiDecision::Unresolved
        );
    }
}
