//! Tier 3 preparation: find the GTDB bac120/ar53 marker genes with HMMER,
//! classify their copy number, pick each genome's domain, and build the
//! masked concatenated marker alignment, following GTDB-Tk 2.8's `identify`
//! and `align` steps.

use flate2::Compression;
use flate2::write::GzEncoder;
use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::error::Error;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::Command;

use needletail::parse_fastx_file;
use rayon::prelude::*;

use crate::classify::external::run;
use crate::classify::refdata::GtdbtkData;
use crate::classify::summary::{py_float, py_round};

pub(crate) const BAC120_PFAM: [&str; 6] = [
    "PF00380.20",
    "PF00410.20",
    "PF00466.21",
    "PF01025.20",
    "PF02576.18",
    "PF03726.15",
];
pub(crate) const BAC120_TIGRFAM: [&str; 114] = [
    "TIGR00006",
    "TIGR00019",
    "TIGR00020",
    "TIGR00029",
    "TIGR00043",
    "TIGR00054",
    "TIGR00059",
    "TIGR00061",
    "TIGR00064",
    "TIGR00065",
    "TIGR00082",
    "TIGR00083",
    "TIGR00084",
    "TIGR00086",
    "TIGR00088",
    "TIGR00090",
    "TIGR00092",
    "TIGR00095",
    "TIGR00115",
    "TIGR00116",
    "TIGR00138",
    "TIGR00158",
    "TIGR00166",
    "TIGR00168",
    "TIGR00186",
    "TIGR00194",
    "TIGR00250",
    "TIGR00337",
    "TIGR00344",
    "TIGR00362",
    "TIGR00382",
    "TIGR00392",
    "TIGR00396",
    "TIGR00398",
    "TIGR00414",
    "TIGR00416",
    "TIGR00420",
    "TIGR00431",
    "TIGR00435",
    "TIGR00436",
    "TIGR00442",
    "TIGR00445",
    "TIGR00456",
    "TIGR00459",
    "TIGR00460",
    "TIGR00468",
    "TIGR00472",
    "TIGR00487",
    "TIGR00496",
    "TIGR00539",
    "TIGR00580",
    "TIGR00593",
    "TIGR00615",
    "TIGR00631",
    "TIGR00634",
    "TIGR00635",
    "TIGR00643",
    "TIGR00663",
    "TIGR00717",
    "TIGR00755",
    "TIGR00810",
    "TIGR00922",
    "TIGR00928",
    "TIGR00959",
    "TIGR00963",
    "TIGR00964",
    "TIGR00967",
    "TIGR01009",
    "TIGR01011",
    "TIGR01017",
    "TIGR01021",
    "TIGR01029",
    "TIGR01032",
    "TIGR01039",
    "TIGR01044",
    "TIGR01059",
    "TIGR01063",
    "TIGR01066",
    "TIGR01071",
    "TIGR01079",
    "TIGR01082",
    "TIGR01087",
    "TIGR01128",
    "TIGR01146",
    "TIGR01164",
    "TIGR01169",
    "TIGR01171",
    "TIGR01302",
    "TIGR01391",
    "TIGR01393",
    "TIGR01394",
    "TIGR01510",
    "TIGR01632",
    "TIGR01951",
    "TIGR01953",
    "TIGR02012",
    "TIGR02013",
    "TIGR02027",
    "TIGR02075",
    "TIGR02191",
    "TIGR02273",
    "TIGR02350",
    "TIGR02386",
    "TIGR02397",
    "TIGR02432",
    "TIGR02729",
    "TIGR03263",
    "TIGR03594",
    "TIGR03625",
    "TIGR03632",
    "TIGR03654",
    "TIGR03723",
    "TIGR03725",
    "TIGR03953",
];
pub(crate) const AR53_PFAM: [&str; 12] = [
    "PF04919.13",
    "PF07541.13",
    "PF01000.27",
    "PF00687.22",
    "PF00466.21",
    "PF00827.18",
    "PF01280.21",
    "PF01090.20",
    "PF01200.19",
    "PF01015.19",
    "PF00900.21",
    "PF00410.20",
];
pub(crate) const AR53_TIGRFAM: [&str; 41] = [
    "TIGR00037",
    "TIGR00064",
    "TIGR00111",
    "TIGR00134",
    "TIGR00279",
    "TIGR00291",
    "TIGR00323",
    "TIGR00335",
    "TIGR00373",
    "TIGR00405",
    "TIGR00448",
    "TIGR00483",
    "TIGR00491",
    "TIGR00522",
    "TIGR00967",
    "TIGR00982",
    "TIGR01008",
    "TIGR01012",
    "TIGR01018",
    "TIGR01020",
    "TIGR01028",
    "TIGR01046",
    "TIGR01052",
    "TIGR01171",
    "TIGR01213",
    "TIGR01952",
    "TIGR02236",
    "TIGR02338",
    "TIGR02389",
    "TIGR02390",
    "TIGR03626",
    "TIGR03627",
    "TIGR03628",
    "TIGR03629",
    "TIGR03670",
    "TIGR03671",
    "TIGR03672",
    "TIGR03673",
    "TIGR03674",
    "TIGR03676",
    "TIGR03680",
];

pub(crate) const MARKER_SETS: [&str; 2] = ["ar53", "bac120"];
const DEFAULT_MULTIHIT_THRESHOLD: f64 = 10.0;

/// Marker IDs of a marker set in GTDB-Tk's concatenation order (sorted).
pub(crate) fn marker_ids(marker_set: &str) -> Vec<&'static str> {
    let mut ids = match marker_set {
        "bac120" => BAC120_PFAM
            .iter()
            .chain(BAC120_TIGRFAM.iter())
            .copied()
            .collect::<Vec<_>>(),
        "ar53" => AR53_PFAM
            .iter()
            .chain(AR53_TIGRFAM.iter())
            .copied()
            .collect::<Vec<_>>(),
        other => panic!("unknown marker set {other}"),
    };
    ids.sort_unstable();
    ids
}

/// One marker's individual HMM and its header fields.
#[derive(Debug, Clone)]
pub(crate) struct MarkerHmm {
    pub(crate) id: &'static str,
    pub(crate) name: String,
    pub(crate) description: String,
    pub(crate) length: usize,
    pub(crate) path: PathBuf,
}

/// The marker set's HMMs in concatenation order.
pub(crate) fn marker_hmms(
    data_root: &Path,
    marker_set: &str,
) -> Result<Vec<MarkerHmm>, Box<dyn Error>> {
    marker_ids(marker_set)
        .into_iter()
        .map(|id| {
            let path = if id.starts_with("PF") {
                data_root
                    .join("markers/pfam/individual_hmms")
                    .join(format!("{id}.hmm"))
            } else {
                data_root
                    .join("markers/tigrfam/individual_hmms")
                    .join(format!("{id}.HMM"))
            };
            let (mut name, mut description, mut length) = (None, None, None);
            let file = File::open(&path)
                .map_err(|error| format!("cannot read {}: {error}", path.display()))?;
            for line in BufReader::new(file).lines() {
                let line = line?;
                let value = || line[4..].trim().to_owned();
                match line.get(..4) {
                    Some("NAME") => name = Some(value()),
                    Some("DESC") => description = Some(value()),
                    Some("LENG") => length = Some(value().parse()?),
                    Some("HMM ") => break,
                    _ => {}
                }
            }
            Ok(MarkerHmm {
                id,
                name: name.unwrap_or_default(),
                description: description.unwrap_or_default(),
                length: length.ok_or_else(|| format!("{} has no LENG line", path.display()))?,
                path,
            })
        })
        .collect()
}

/// `align/intermediate_results/{prefix}.{marker_set}.marker_info.tsv`.
pub(crate) fn write_marker_info(path: &Path, hmms: &[MarkerHmm]) -> std::io::Result<()> {
    let mut output = BufWriter::new(File::create(path)?);
    writeln!(output, "Marker ID\tName\tDescription\tLength (bp)")?;
    for hmm in hmms {
        writeln!(
            output,
            "{}\t{}\t{}\t{}",
            hmm.id, hmm.name, hmm.description, hmm.length
        )?;
    }
    output.flush()
}

/// An HMM hit; ordering matches GTDB-Tk's `Hit.__lt__` (larger is more significant).
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Hit {
    pub(crate) gene_id: String,
    pub(crate) hmm_id: String,
    pub(crate) evalue: f64,
    pub(crate) bitscore: f64,
}

impl Hit {
    fn significance(&self, other: &Self) -> Ordering {
        self.bitscore
            .total_cmp(&other.bitscore)
            .then(other.evalue.total_cmp(&self.evalue))
            .then(other.hmm_id.cmp(&self.hmm_id))
            .then(other.gene_id.cmp(&self.gene_id))
    }
}

/// Most significant hit per (gene, HMM); TIGRFAM keeps only one HMM per gene.
#[derive(Debug, Default)]
pub(crate) struct TopHits {
    hits: BTreeMap<String, BTreeMap<String, Hit>>,
}

impl TopHits {
    fn add(&mut self, hit: Hit) {
        let gene = self.hits.entry(hit.gene_id.clone()).or_default();
        match gene.get(&hit.hmm_id) {
            Some(current) if current.significance(&hit) != Ordering::Less => {}
            _ => {
                gene.insert(hit.hmm_id.clone(), hit);
            }
        }
    }

    fn add_single_family(&mut self, hit: Hit) {
        let replace = match self.top_hit(&hit.gene_id) {
            Some(top) => top.significance(&hit) == Ordering::Less,
            None => true,
        };
        if replace {
            self.hits.insert(hit.gene_id.clone(), BTreeMap::new());
            self.add(hit);
        }
    }

    fn top_hit(&self, gene_id: &str) -> Option<&Hit> {
        self.hits
            .get(gene_id)?
            .values()
            .max_by(|left, right| left.significance(right))
    }

    fn iter(&self) -> impl Iterator<Item = &Hit> {
        self.hits.values().flat_map(BTreeMap::values)
    }

    /// GTDB-Tk's `*_tophit.tsv` format.
    fn write(&self, path: &Path) -> std::io::Result<()> {
        let mut output = BufWriter::new(File::create(path)?);
        writeln!(output, "Gene Id\tTop hits (Family id,e-value,bitscore)")?;
        for (gene, hits) in &self.hits {
            let mut hits = hits.values().collect::<Vec<_>>();
            hits.sort_by(|left, right| right.significance(left));
            let joined = hits
                .iter()
                .map(|hit| {
                    format!(
                        "{},{},{}",
                        hit.hmm_id,
                        py_float(hit.evalue),
                        py_float(hit.bitscore)
                    )
                })
                .collect::<Vec<_>>()
                .join(";");
            writeln!(output, "{gene}\t{joined}")?;
        }
        output.flush()
    }
}

/// Pfam family annotations from `Pfam-A.hmm.dat`, keyed by family name as
/// PfamScan does.
#[derive(Debug, Default)]
pub(crate) struct PfamInfo {
    clan: HashMap<String, String>,
    family_type: HashMap<String, String>,
    /// PfamScan keeps only the last `NE` line seen for each family, in both directions.
    nested: HashMap<String, String>,
    /// Gathering thresholds (sequence, domain).
    gathering: HashMap<String, (f64, f64)>,
}

impl PfamInfo {
    pub(crate) fn load(data_root: &Path) -> Result<Self, Box<dyn Error>> {
        let path = data_root.join("markers/pfam/Pfam-A.hmm.dat");
        let mut info = Self::default();
        let mut name = String::new();
        for line in BufReader::new(File::open(&path)?).lines() {
            let line = line?;
            let Some(rest) = line.strip_prefix("#=GF ") else {
                continue;
            };
            let mut fields = rest.split_whitespace();
            let (Some(key), Some(value)) = (fields.next(), fields.next()) else {
                continue;
            };
            match key {
                "ID" => name = value.to_owned(),
                "TP" => {
                    info.family_type.insert(name.clone(), value.to_owned());
                }
                "CL" => {
                    info.clan.insert(name.clone(), value.to_owned());
                }
                "NE" => {
                    info.nested.insert(name.clone(), value.to_owned());
                    info.nested.insert(value.to_owned(), name.clone());
                }
                "GA" => {
                    let domain = fields.next().unwrap_or(value);
                    let parse = |text: &str| text.trim_end_matches(';').parse::<f64>();
                    info.gathering
                        .insert(name.clone(), (parse(value)?, parse(domain)?));
                }
                _ => {}
            }
        }
        Ok(info)
    }
}

#[derive(Debug, Clone)]
struct PfamDomain {
    gene_id: String,
    accession: String,
    name: String,
    ali: (usize, usize),
    env: (usize, usize),
    hmm: (usize, usize),
    hmm_length: usize,
    sequence_bitscore: f64,
    bitscore: f64,
    evalue: f64,
}

/// Parse `hmmsearch --domtblout` domains (score = domain score, E-value = i-Evalue).
fn read_pfam_domains(path: &Path) -> Result<Vec<PfamDomain>, Box<dyn Error>> {
    let mut domains = Vec::new();
    for line in BufReader::new(File::open(path)?).lines() {
        let line = line?;
        if line.starts_with('#') {
            continue;
        }
        let fields = line.split_whitespace().collect::<Vec<_>>();
        if fields.len() < 22 {
            continue;
        }
        let number = |index: usize| fields[index].parse::<usize>();
        domains.push(PfamDomain {
            gene_id: fields[0].to_owned(),
            name: fields[3].to_owned(),
            accession: fields[4].to_owned(),
            hmm_length: number(5)?,
            sequence_bitscore: fields[7].parse()?,
            evalue: fields[12].parse()?,
            bitscore: fields[13].parse()?,
            hmm: (number(15)?, number(16)?),
            ali: (number(17)?, number(18)?),
            env: (number(19)?, number(20)?),
        });
    }
    Ok(domains)
}

/// PfamScan's clan-overlap resolution and significance filter. Per gene,
/// domains are visited in E-value order (stable, so ties keep hmmsearch's
/// order) and dropped when they overlap a kept domain of the same clan that
/// is not nested in it. A kept domain is reported only when both its
/// sequence and domain scores reach the family's gathering thresholds.
fn resolve_clan_overlaps(domains: Vec<PfamDomain>, info: &PfamInfo) -> Vec<PfamDomain> {
    let mut by_gene: BTreeMap<String, Vec<PfamDomain>> = BTreeMap::new();
    for domain in domains {
        by_gene
            .entry(domain.gene_id.clone())
            .or_default()
            .push(domain);
    }
    let mut kept_all = Vec::new();
    for (_, mut domains) in by_gene {
        domains.sort_by(|left, right| left.evalue.total_cmp(&right.evalue));
        let mut kept: Vec<PfamDomain> = Vec::new();
        for domain in domains {
            let overlaps = kept.iter().any(|other| {
                let same_clan = match (info.clan.get(&domain.name), info.clan.get(&other.name)) {
                    (Some(left), Some(right)) => left == right,
                    _ => false,
                };
                let (first, second) = if domain.ali.0 <= other.ali.0 {
                    (&domain, other)
                } else {
                    (other, &domain)
                };
                same_clan
                    && second.ali.0 <= first.ali.1
                    && info.nested.get(&domain.name) != Some(&other.name)
            });
            if !overlaps {
                kept.push(domain);
            }
        }
        // Report in alignment-start order, as PfamScan writes them.
        kept.sort_by_key(|domain| domain.ali.0);
        kept_all.extend(kept.into_iter().filter(|domain| {
            info.gathering
                .get(&domain.name)
                .is_some_and(|&(sequence, single)| {
                    domain.sequence_bitscore >= sequence && domain.bitscore >= single
                })
        }));
    }
    kept_all
}

/// PfamScan's `write_results` layout. The banner names aeon rather than
/// PyPfam; header lines and the fixed-width body match it.
fn write_pfam_table(
    path: &Path,
    domains: &[PfamDomain],
    info: &PfamInfo,
    protein_faa: &Path,
    pfam_hmm: &Path,
    max_seqname: usize,
) -> std::io::Result<()> {
    let mut output = BufWriter::new(File::create(path)?);
    writeln!(output, "{}", "#".repeat(84))?;
    let title = format!(
        "aeon v{} Pfam marker search (PfamScan output format)",
        env!("CARGO_PKG_VERSION")
    );
    writeln!(output, "# {title:^80} #")?;
    writeln!(output, "# {} #", "=".repeat(80))?;
    for (key, value) in [
        (
            "Execution time:",
            chrono::Local::now()
                .format("%Y-%m-%d %H:%M:%S%.6f")
                .to_string(),
        ),
        ("Query sequence file:", protein_faa.display().to_string()),
        ("CPU threads:", "1".to_owned()),
        ("Searching against:", pfam_hmm.display().to_string()),
        ("Resolve clan overlaps:", "on".to_owned()),
    ] {
        writeln!(output, "# {key:<25} {value}")?;
    }
    writeln!(output, "# {} #", "=".repeat(80))?;
    writeln!(
        output,
        "# <seq id> <alignment start> <alignment end> <envelope start> <envelope end> <hmm acc> <hmm name> <type> <hmm start> <hmm end> <hmm length> <bit score> <E-value> <significance> <clan> "
    )?;
    let width = if max_seqname < 1 { 20 } else { max_seqname };
    for domain in domains {
        writeln!(
            output,
            "{:<width$} {:>6} {:>6} {:>6} {:>6} {:<11} {:<16} {:>7} {:>5} {:>5} {:>5} {:>8} {:>9} {:>3} {:<8} ",
            domain.gene_id,
            domain.ali.0,
            domain.ali.1,
            domain.env.0,
            domain.env.1,
            domain.accession,
            domain.name,
            info.family_type
                .get(&domain.name)
                .map_or("Family", String::as_str),
            domain.hmm.0,
            domain.hmm.1,
            domain.hmm_length,
            py_float(domain.bitscore),
            py_float(domain.evalue),
            1,
            info.clan
                .get(&domain.name)
                .map_or("No_clan", String::as_str)
        )?;
    }
    output.flush()
}

/// Run TIGRFAM and Pfam HMM searches for one genome's proteins, writing
/// GTDB-Tk's `{id}_tigrfam.tsv`, `{id}_pfam.tsv` and top-hit files.
pub(crate) fn search_markers(
    hmmsearch: &Path,
    data_root: &Path,
    pfam_info: &PfamInfo,
    genome_id: &str,
    protein_faa: &Path,
    genome_dir: &Path,
) -> Result<(TopHits, TopHits), Box<dyn Error>> {
    let path = |suffix: &str| genome_dir.join(format!("{genome_id}{suffix}"));
    run(Command::new(hmmsearch)
        .arg("-o")
        .arg(path("_tigrfam.out"))
        .arg("--tblout")
        .arg(path("_tigrfam.tsv"))
        .args(["--noali", "--notextw", "--cut_nc", "--cpu", "1"])
        .arg(data_root.join("markers/tigrfam/tigrfam.hmm"))
        .arg(protein_faa))?;
    let mut tigrfam = TopHits::default();
    for line in BufReader::new(File::open(path("_tigrfam.tsv"))?).lines() {
        let line = line?;
        if line.starts_with('#') {
            continue;
        }
        let fields = line.split_whitespace().collect::<Vec<_>>();
        if fields.len() < 6 {
            continue;
        }
        tigrfam.add_single_family(Hit {
            gene_id: fields[0].to_owned(),
            hmm_id: fields[3].to_owned(),
            evalue: fields[4].parse()?,
            bitscore: fields[5].parse()?,
        });
    }
    tigrfam.write(&path("_tigrfam_tophit.tsv"))?;

    let domtbl = path("_pfam.domtblout");
    let pfam_hmm = data_root.join("markers/pfam/Pfam-A.hmm");
    run(Command::new(hmmsearch)
        .args(["--notextw", "--cut_ga", "--cpu", "1", "-o", "/dev/null"])
        .arg("--domtblout")
        .arg(&domtbl)
        .arg(&pfam_hmm)
        .arg(protein_faa))?;
    let domains = resolve_clan_overlaps(read_pfam_domains(&domtbl)?, pfam_info);
    fs::remove_file(&domtbl)?;
    write_pfam_table(
        &path("_pfam.tsv"),
        &domains,
        pfam_info,
        protein_faa,
        &pfam_hmm,
        longest_sequence_id(protein_faa)?,
    )?;
    let mut pfam = TopHits::default();
    for domain in domains {
        pfam.add(Hit {
            gene_id: domain.gene_id,
            hmm_id: domain.accession,
            evalue: domain.evalue,
            bitscore: domain.bitscore,
        });
    }
    pfam.write(&path("_pfam_tophit.tsv"))?;
    Ok((pfam, tigrfam))
}

fn longest_sequence_id(fasta: &Path) -> Result<usize, Box<dyn Error>> {
    let mut longest = 0;
    for line in BufReader::new(File::open(fasta)?).lines() {
        if let Some(header) = line?.strip_prefix('>') {
            longest = longest.max(header.split_whitespace().next().unwrap_or_default().len());
        }
    }
    Ok(longest)
}

/// GTDB-Tk's per-marker copy-number classes for one genome and marker set.
#[derive(Debug, Default, Clone)]
pub(crate) struct CopyNumber {
    /// Unique and multiple-but-identical markers: marker -> (gene, sequence).
    pub(crate) unique: BTreeMap<String, (String, String)>,
    pub(crate) multiple_unique: BTreeMap<String, (String, String)>,
    pub(crate) multiple: BTreeSet<String>,
    pub(crate) missing: BTreeSet<String>,
}

impl CopyNumber {
    pub(crate) fn single_copy_count(&self) -> usize {
        self.unique.len() + self.multiple_unique.len()
    }
}

pub(crate) fn copy_number(
    marker_set: &str,
    proteins: &HashMap<String, String>,
    pfam: &TopHits,
    tigrfam: &TopHits,
) -> CopyNumber {
    let mut by_marker: HashMap<&str, Vec<&Hit>> = HashMap::new();
    for hit in pfam.iter().chain(tigrfam.iter()) {
        by_marker.entry(hit.hmm_id.as_str()).or_default().push(hit);
    }
    let mut result = CopyNumber::default();
    for marker in marker_ids(marker_set) {
        let sequence = |hit: &Hit| proteins.get(&hit.gene_id).cloned().unwrap_or_default();
        match by_marker.get(marker).map(Vec::as_slice) {
            None | Some([]) => {
                result.missing.insert(marker.to_owned());
            }
            Some([hit]) => {
                result
                    .unique
                    .insert(marker.to_owned(), (hit.gene_id.clone(), sequence(hit)));
            }
            Some(hits) => {
                let sequences = hits
                    .iter()
                    .map(|hit| sequence(hit))
                    .collect::<BTreeSet<_>>();
                if sequences.len() == 1 {
                    let top = hits
                        .iter()
                        .max_by(|left, right| left.significance(right))
                        .expect("non-empty");
                    result
                        .multiple_unique
                        .insert(marker.to_owned(), (top.gene_id.clone(), sequence(top)));
                } else {
                    result.multiple.insert(marker.to_owned());
                }
            }
        }
    }
    result
}

/// Proteins keyed by the first word of each header, trailing `*` removed.
pub(crate) fn read_proteins(path: &Path) -> Result<HashMap<String, String>, Box<dyn Error>> {
    let mut reader = parse_fastx_file(path)?;
    let mut proteins = HashMap::new();
    while let Some(record) = reader.next() {
        let record = record?;
        let id = String::from_utf8_lossy(record.id());
        let id = id.split_whitespace().next().unwrap_or_default().to_owned();
        let mut sequence = String::from_utf8_lossy(&record.seq()).into_owned();
        if sequence.ends_with('*') {
            sequence.pop();
        }
        proteins.insert(id, sequence);
    }
    Ok(proteins)
}

/// `identify/{prefix}.{marker_set}.markers_summary.tsv`.
pub(crate) fn write_markers_summary(
    path: &Path,
    genomes: &BTreeMap<String, CopyNumber>,
) -> std::io::Result<()> {
    let mut output = BufWriter::new(File::create(path)?);
    writeln!(
        output,
        "name\tnumber_unique_genes\tnumber_multiple_genes\tnumber_multiple_unique_genes\tnumber_missing_genes\tlist_unique_genes\tlist_multiple_genes\tlist_multiple_unique_genes\tlist_missing_genes"
    )?;
    let join = |names: Vec<&String>| {
        names
            .into_iter()
            .map(String::as_str)
            .collect::<Vec<_>>()
            .join(",")
    };
    for (genome, copies) in genomes {
        writeln!(
            output,
            "{genome}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            copies.unique.len(),
            copies.multiple.len(),
            copies.multiple_unique.len(),
            copies.missing.len(),
            join(copies.unique.keys().collect()),
            join(copies.multiple.iter().collect()),
            join(copies.multiple_unique.keys().collect()),
            join(copies.missing.iter().collect()),
        )?;
    }
    output.flush()
}

/// Domain call from single-copy marker percentages, as GTDB-Tk's
/// `genome_domain`; returns the marker set and, when the call is within 10
/// points, the rounded (bacterial, archaeal) percentages for the warning.
pub(crate) fn genome_domain(
    bac_single_copy: usize,
    ar_single_copy: usize,
) -> (&'static str, Option<(f64, f64)>) {
    let bacterial = bac_single_copy as f64 * 100.0 / 120.0;
    let archaeal = ar_single_copy as f64 * 100.0 / 53.0;
    let marker_set = if bacterial >= archaeal {
        "bac120"
    } else {
        "ar53"
    };
    let questionable = ((bacterial - archaeal).abs() <= 10.0).then(|| {
        (
            py_round(bacterial, 1).parse().unwrap(),
            py_round(archaeal, 1).parse().unwrap(),
        )
    });
    (marker_set, questionable)
}

/// GTDB-Tk's multiple-hit warning for the genome's assigned marker set.
pub(crate) fn multihit_warning(copies: &CopyNumber, marker_set: &str) -> Option<String> {
    let percent = copies.multiple.len() as f64 * 100.0 / marker_ids(marker_set).len() as f64;
    (percent >= DEFAULT_MULTIHIT_THRESHOLD).then(|| {
        format!(
            "Genome has more than {}% of markers with multiple hits",
            py_round(percent, 1)
        )
    })
}

/// Align the genomes' single-copy markers with one `hmmalign` call per
/// marker (in parallel) and concatenate match columns in marker order;
/// missing markers become gaps of the HMM length. Genomes without any
/// single-copy marker are left out, as in GTDB-Tk.
pub(crate) fn align_markers(
    hmmalign: &Path,
    hmms: &[MarkerHmm],
    genomes: &BTreeMap<String, CopyNumber>,
    work_dir: &Path,
) -> Result<BTreeMap<String, String>, Box<dyn Error>> {
    fs::create_dir_all(work_dir)?;
    let ids = genomes.keys().collect::<Vec<_>>();
    let per_marker = hmms
        .par_iter()
        .map(|hmm| -> Result<HashMap<String, String>, String> {
            let input = work_dir.join(format!("{}.faa", hmm.id));
            let mut present = Vec::new();
            {
                let mut writer = BufWriter::new(File::create(&input).map_err(|e| e.to_string())?);
                for (index, id) in ids.iter().enumerate() {
                    let copies = &genomes[*id];
                    if let Some((_, sequence)) = copies
                        .unique
                        .get(hmm.id)
                        .or_else(|| copies.multiple_unique.get(hmm.id))
                    {
                        writeln!(writer, ">g{index}\n{sequence}").map_err(|e| e.to_string())?;
                        present.push(index);
                    }
                }
                writer.flush().map_err(|e| e.to_string())?;
            }
            if present.is_empty() {
                fs::remove_file(&input).map_err(|e| e.to_string())?;
                return Ok(HashMap::new());
            }
            let output = Command::new(hmmalign)
                .args(["--outformat", "Pfam"])
                .arg(&hmm.path)
                .arg(&input)
                .output()
                .map_err(|error| format!("failed to start hmmalign: {error}"))?;
            fs::remove_file(&input).map_err(|e| e.to_string())?;
            if !output.status.success() {
                return Err(format!(
                    "hmmalign failed for {}: {}",
                    hmm.id,
                    String::from_utf8_lossy(&output.stderr).trim()
                ));
            }
            let aligned = match_columns(&String::from_utf8_lossy(&output.stdout))
                .map_err(|error| format!("{}: {error}", hmm.id))?;
            if aligned.len() != present.len() {
                return Err(format!("Not all genomes could be aligned to {}", hmm.id));
            }
            for row in aligned.values() {
                if row.len() != hmm.length {
                    return Err(format!(
                        "{} alignment has {} columns; expected {}",
                        hmm.id,
                        row.len(),
                        hmm.length
                    ));
                }
            }
            Ok(aligned)
        })
        .collect::<Result<Vec<_>, _>>()?;

    let mut alignments = BTreeMap::new();
    for (index, id) in ids.iter().enumerate() {
        if genomes[*id].single_copy_count() == 0 {
            continue;
        }
        let tag = format!("g{index}");
        let mut alignment = String::new();
        for (hmm, aligned) in hmms.iter().zip(&per_marker) {
            match aligned.get(&tag) {
                Some(row) => alignment.push_str(row),
                None => alignment.push_str(&"-".repeat(hmm.length)),
            }
        }
        alignments.insert((*id).clone(), alignment);
    }
    Ok(alignments)
}
/// Match-state columns (`#=GC RF` == 'x') of each sequence in Pfam-format output.
fn match_columns(stockholm: &str) -> Result<HashMap<String, String>, Box<dyn Error>> {
    let mut rows = Vec::new();
    let mut reference = None;
    for line in stockholm.lines() {
        if line.starts_with("#=GC RF") {
            reference = line.rsplit(' ').next().map(str::to_owned);
        } else if !line.starts_with('#') && !line.starts_with("//") && !line.trim().is_empty() {
            let (name, row) = line.split_once(' ').ok_or("malformed hmmalign row")?;
            rows.push((name.to_owned(), row.trim().to_owned()));
        }
    }
    let reference = reference.ok_or("hmmalign output has no #=GC RF line")?;
    Ok(rows
        .into_iter()
        .map(|(name, row)| {
            let kept = row
                .chars()
                .zip(reference.chars())
                .filter(|(_, mask)| *mask == 'x')
                .map(|(residue, _)| residue)
                .collect();
            (name, kept)
        })
        .collect())
}

pub(crate) fn read_mask(
    data: &GtdbtkData,
    root: &Path,
    marker_set: &str,
) -> Result<Vec<bool>, Box<dyn Error>> {
    let path = root
        .join("masks")
        .join(format!("gtdb_{}_{marker_set}.mask", data.version));
    Ok(fs::read_to_string(&path)
        .map_err(|error| format!("cannot read {}: {error}", path.display()))?
        .trim()
        .bytes()
        .map(|byte| byte == b'1')
        .collect())
}

/// Masked alignment, its share of non-gap characters (0-1), and GTDB-Tk's
/// `aa_percent_msa` (alphabetic characters, percent).
pub(crate) fn apply_mask(
    alignment: &str,
    mask: &[bool],
) -> Result<(String, f64, f64), Box<dyn Error>> {
    if alignment.len() != mask.len() {
        return Err(format!(
            "Mask ({}) and alignment ({}) length do not match.",
            mask.len(),
            alignment.len()
        )
        .into());
    }
    let masked = alignment
        .chars()
        .zip(mask)
        .filter(|(_, keep)| **keep)
        .map(|(residue, _)| residue)
        .collect::<String>();
    if masked.is_empty() {
        return Ok((masked, 0.0, 0.0));
    }
    let valid = masked.chars().filter(|&c| c != '.' && c != '-').count();
    let alphabetic = masked.chars().filter(char::is_ascii_alphabetic).count();
    let length = masked.len() as f64;
    Ok((
        masked,
        valid as f64 / length,
        alphabetic as f64 / length * 100.0,
    ))
}

/// GTDB's concatenated reference MSA (`msa/gtdb_{release}_{marker_set}.faa`).
pub(crate) fn read_reference_msa(
    root: &Path,
    release: &str,
    marker_set: &str,
) -> Result<BTreeMap<String, String>, Box<dyn Error>> {
    read_msa(
        &root
            .join("msa")
            .join(format!("gtdb_{release}_{marker_set}.faa")),
    )
}

/// An aligned FASTA file, plain or gzipped, keyed by sequence ID.
pub(crate) fn read_msa(path: &Path) -> Result<BTreeMap<String, String>, Box<dyn Error>> {
    let mut msa = BTreeMap::new();
    // Test packages may ship an empty MSA, which needletail rejects.
    if fs::metadata(path).is_ok_and(|metadata| metadata.len() == 0) {
        return Ok(msa);
    }
    let mut reader = parse_fastx_file(path)
        .map_err(|error| format!("cannot read {}: {error}", path.display()))?;
    while let Some(record) = reader.next() {
        let record = record?;
        let id = String::from_utf8_lossy(record.id());
        let id = id.split_whitespace().next().unwrap_or_default().to_owned();
        msa.insert(id, String::from_utf8_lossy(&record.seq()).into_owned());
    }
    Ok(msa)
}

/// Write sorted records gzip-compressed, as GTDB-Tk's `_write_msa`: GTDB
/// genomes carry their taxonomy in the header, user genomes only their ID.
pub(crate) fn write_msa_gz(
    path: &Path,
    sequences: &BTreeMap<String, String>,
    data: &GtdbtkData,
) -> std::io::Result<()> {
    let mut output = GzEncoder::new(BufWriter::new(File::create(path)?), Compression::default());
    for (id, sequence) in sequences {
        match data.taxonomy_of(id) {
            Some(taxonomy) if id.starts_with("GB_") || id.starts_with("RS_") => {
                writeln!(output, ">{id} {}\n{sequence}", taxonomy.join(";"))?
            }
            _ => writeln!(output, ">{id}\n{sequence}")?,
        }
    }
    output.finish()?.flush()
}

/// GTDB-Tk's reason for pruning a masked alignment from the MSA.
pub(crate) fn filtered_reason(masked: &str) -> String {
    let percent = if masked.is_empty() {
        0.0
    } else {
        masked.chars().filter(char::is_ascii_alphabetic).count() as f64 * 100.0
            / masked.len() as f64
    };
    format!("Insufficient number of amino acids in MSA ({percent:.1}%)")
}

/// `align/{prefix}.{marker_set}.filtered.tsv` lines for pruned user genomes.
pub(crate) fn write_filtered(
    path: &Path,
    pruned: &BTreeMap<String, String>,
) -> std::io::Result<()> {
    let mut output = BufWriter::new(File::create(path)?);
    for (id, masked) in pruned {
        writeln!(output, "{id}\t{}", filtered_reason(masked))?;
    }
    output.flush()
}
#[cfg(test)]
mod tests {
    use super::*;

    fn hit(gene: &str, hmm: &str, evalue: f64, bitscore: f64) -> Hit {
        Hit {
            gene_id: gene.into(),
            hmm_id: hmm.into(),
            evalue,
            bitscore,
        }
    }

    #[test]
    fn marker_sets_have_gtdb_sizes() {
        assert_eq!(marker_ids("bac120").len(), 120);
        assert_eq!(marker_ids("ar53").len(), 53);
        assert_eq!(marker_ids("ar53")[0], "PF00410.20");
    }

    #[test]
    fn tigrfam_keeps_only_the_most_significant_family_per_gene() {
        let mut hits = TopHits::default();
        hits.add_single_family(hit("g1", "TIGR00001", 1e-10, 50.0));
        hits.add_single_family(hit("g1", "TIGR00002", 1e-20, 80.0));
        hits.add_single_family(hit("g1", "TIGR00003", 1e-5, 20.0));
        assert_eq!(
            hits.iter()
                .map(|hit| hit.hmm_id.as_str())
                .collect::<Vec<_>>(),
            ["TIGR00002"]
        );
    }

    #[test]
    fn copy_number_classes_follow_gtdbtk() {
        let mut tigrfam = TopHits::default();
        tigrfam.add_single_family(hit("g1", "TIGR00037", 1e-40, 140.0));
        tigrfam.add_single_family(hit("g2", "TIGR00064", 1e-40, 140.0));
        tigrfam.add_single_family(hit("g3", "TIGR00064", 1e-30, 100.0));
        tigrfam.add_single_family(hit("g4", "TIGR00111", 1e-40, 140.0));
        tigrfam.add_single_family(hit("g5", "TIGR00111", 1e-30, 100.0));
        let proteins = [
            ("g1", "MA"),
            ("g2", "MB"),
            ("g3", "MC"),
            ("g4", "MD"),
            ("g5", "MD"),
        ]
        .into_iter()
        .map(|(id, seq)| (id.to_owned(), seq.to_owned()))
        .collect();
        let copies = copy_number("ar53", &proteins, &TopHits::default(), &tigrfam);
        assert!(copies.unique.contains_key("TIGR00037"));
        assert!(copies.multiple.contains("TIGR00064"));
        assert_eq!(copies.multiple_unique["TIGR00111"].0, "g4");
        assert_eq!(copies.missing.len(), 50);
    }

    #[test]
    fn domain_call_and_questionable_warning() {
        assert_eq!(genome_domain(110, 5), ("bac120", None));
        assert_eq!(genome_domain(10, 50).0, "ar53");
        assert_eq!(genome_domain(60, 25), ("bac120", Some((50.0, 47.2))));
    }

    #[test]
    fn match_columns_drop_insert_states() {
        let stockholm = "# STOCKHOLM 1.0\n\ng0           MK-ae.LV\n#=GC RF      xxx...xx\n//\n";
        assert_eq!(match_columns(stockholm).unwrap()["g0"], "MK-LV");
    }

    #[test]
    fn mask_counts_valid_and_alphabetic_columns() {
        let (masked, valid, percent) =
            apply_mask("MK-.AB", &[true, true, true, true, false, true]).unwrap();
        assert_eq!(masked, "MK-.B");
        assert!((valid - 0.6).abs() < 1e-12);
        assert!((percent - 60.0).abs() < 1e-12);
    }
}

#[cfg(test)]
mod gtdbtk_comparison {
    use super::*;

    /// AEON_MARKER_DATA=<gtdbtk data> AEON_MARKER_IN=<marker_genes dir> AEON_MARKER_OUT=<dir>
    #[test]
    #[ignore]
    fn identify_like_gtdbtk() {
        let data = PathBuf::from(std::env::var("AEON_MARKER_DATA").unwrap());
        let input = PathBuf::from(std::env::var("AEON_MARKER_IN").unwrap());
        let out = PathBuf::from(std::env::var("AEON_MARKER_OUT").unwrap());
        let info = PfamInfo::load(&data).unwrap();
        let hmmsearch = crate::classify::external::find_program("hmmsearch").unwrap();
        let mut summaries: BTreeMap<&str, BTreeMap<String, CopyNumber>> = BTreeMap::new();
        let mut ids = fs::read_dir(&input)
            .unwrap()
            .map(|e| e.unwrap().file_name().into_string().unwrap())
            .collect::<Vec<_>>();
        ids.sort();
        for id in ids {
            let dir = out.join(&id);
            fs::create_dir_all(&dir).unwrap();
            let faa = input.join(&id).join(format!("{id}_protein.faa"));
            let (pfam, tigrfam) =
                search_markers(&hmmsearch, &data, &info, &id, &faa, &dir).unwrap();
            let proteins = read_proteins(&faa).unwrap();
            for set in MARKER_SETS {
                summaries
                    .entry(set)
                    .or_default()
                    .insert(id.clone(), copy_number(set, &proteins, &pfam, &tigrfam));
            }
        }
        for (set, genomes) in &summaries {
            write_markers_summary(
                &out.join(format!("gtdbtk.{set}.markers_summary.tsv")),
                genomes,
            )
            .unwrap();
        }
        let gtdb = GtdbtkData::load(data.to_str().unwrap()).unwrap();
        let hmmalign = crate::classify::external::find_program("hmmalign").unwrap();
        for set in MARKER_SETS {
            let domain_genomes = summaries[set]
                .iter()
                .filter(|(id, _)| {
                    let (bac, ar) = (&summaries["bac120"][*id], &summaries["ar53"][*id]);
                    genome_domain(bac.single_copy_count(), ar.single_copy_count()).0 == set
                })
                .map(|(id, copies)| (id.clone(), copies.clone()))
                .collect::<BTreeMap<_, _>>();
            if domain_genomes.is_empty() {
                continue;
            }
            let hmms = marker_hmms(&data, set).unwrap();
            write_marker_info(&out.join(format!("gtdbtk.{set}.marker_info.tsv")), &hmms).unwrap();
            let user =
                align_markers(&hmmalign, &hmms, &domain_genomes, &out.join("align_tmp")).unwrap();
            let mask = read_mask(&gtdb, &data, set).unwrap();
            let mut masked = BTreeMap::new();
            for (id, alignment) in &user {
                let (sequence, _, percent) = apply_mask(alignment, &mask).unwrap();
                eprintln!("{set} {id} msa_percent={}", py_round(percent, 2));
                masked.insert(id.clone(), sequence);
            }
            write_msa_gz(
                &out.join(format!("gtdbtk.{set}.user_msa.fasta.gz")),
                &masked,
                &gtdb,
            )
            .unwrap();
            let mut all = masked.clone();
            for (id, alignment) in read_reference_msa(&data, &gtdb.version, set).unwrap() {
                all.insert(id, apply_mask(&alignment, &mask).unwrap().0);
            }
            write_msa_gz(&out.join(format!("gtdbtk.{set}.msa.fasta.gz")), &all, &gtdb).unwrap();
        }
    }
}
