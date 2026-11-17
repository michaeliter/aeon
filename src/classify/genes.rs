//! Gene calling with FragGeneScanRs's complete-genome model, run as a
//! separate program (FragGeneScanRs is GPL-3.0; aeon is MIT).

use needletail::parse_fastx_file;
use std::env;
use std::error::Error;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::Command;

use crate::classify::external::{find_program, run_with_input};

/// FragGeneScanRs translates with the standard bacterial code.
pub(crate) const TRANSLATION_TABLE: u8 = 11;
const TRAINING_MODEL: &str = "complete";
/// FragGeneScanRs reports no gene shorter than this in complete genomes, and
/// fails on records of one or two bases.
const MIN_STRETCH: usize = 120;

/// FragGeneScanRs and the directory holding its training models.
pub(crate) struct GeneCaller {
    program: PathBuf,
    train_dir: PathBuf,
}

impl GeneCaller {
    /// Find FragGeneScanRs and its `train/` directory: `train_dir` if given,
    /// else a `train/` next to the program or in its `share/FragGeneScanRs`,
    /// or the sources `cargo install` left behind.
    pub(crate) fn locate(program: &str, train_dir: Option<&str>) -> Result<Self, Box<dyn Error>> {
        let program = find_program(program)?;
        let train_dir = match train_dir.map(PathBuf::from) {
            Some(dir) => dir,
            None => default_train_dirs(&program)
                .into_iter()
                .find(|dir| is_train_dir(dir))
                .ok_or(
                    "FragGeneScanRs training models were not found; pass --fgs-train-dir (the train/ directory of the FragGeneScanRs sources)",
                )?,
        };
        if !is_train_dir(&train_dir) {
            return Err(format!(
                "{} does not contain the FragGeneScanRs '{TRAINING_MODEL}' model",
                train_dir.display()
            )
            .into());
        }
        Ok(Self { program, train_dir })
    }

    /// Call genes for one genome into `out_dir/{genome_id}/`, writing
    /// `{id}_protein.faa` and, with `gene_files`, `{id}_protein.fna` and
    /// `{id}_protein.gff` like GTDB-Tk's `identify` step. Returns the protein
    /// FASTA, or `None` for a genome without sequence or genes.
    /// FragGeneScanRs reads the genome from its standard input and splits
    /// the contigs over `threads`. It stops calling genes on a record at its
    /// first non-ACGT base, so each stretch between ambiguous bases is fed as
    /// its own record and the gene coordinates are mapped back to the contigs.
    pub(crate) fn call(
        &self,
        genome_path: &str,
        genome_id: &str,
        out_dir: &Path,
        threads: usize,
        gene_files: bool,
    ) -> Result<Option<PathBuf>, Box<dyn Error>> {
        let genome_dir = out_dir.join(genome_id);
        fs::create_dir_all(&genome_dir)?;
        // Genome IDs such as GCF_000005845.2 contain dots, so names are built
        // explicitly rather than with `Path::with_extension`.
        let output = |suffix: &str| genome_dir.join(format!("{genome_id}{suffix}"));
        let protein_faa = output("_protein.faa");
        let mut command = Command::new(&self.program);
        command
            .args(["-w", "1", "-t", TRAINING_MODEL])
            .args(["-p", &threads.max(1).to_string()])
            .arg("-r")
            .arg(&self.train_dir)
            .arg("-a")
            .arg(&protein_faa);
        if gene_files {
            command
                .arg("-n")
                .arg(output("_protein.fna"))
                .arg("-g")
                .arg(output("_protein.gff"));
        }
        let mut stretches = Vec::new();
        run_with_input(&mut command, |stdin| {
            stretches = write_stretches(genome_path, stdin)?;
            Ok(())
        })?;
        if stretches.is_empty() {
            return Ok(None);
        }
        restore_coordinates(&protein_faa, &stretches)?;
        if gene_files {
            restore_coordinates(&output("_protein.fna"), &stretches)?;
            restore_coordinates(&output("_protein.gff"), &stretches)?;
        }
        Ok((count_records(&protein_faa)? > 0).then_some(protein_faa))
    }
}

/// A run of unambiguous bases within a contig, fed to FragGeneScanRs as the
/// record named by its index.
struct Stretch {
    contig: String,
    offset: usize,
}

/// Write the ACGT stretches long enough to hold a gene of every contig in
/// `input` as records named by their index, and return them.
fn write_stretches(input: &str, output: &mut dyn Write) -> Result<Vec<Stretch>, Box<dyn Error>> {
    let mut reader =
        parse_fastx_file(input).map_err(|error| format!("invalid FASTA {input}: {error}"))?;
    let mut stretches = Vec::new();
    while let Some(record) = reader.next() {
        let record = record.map_err(|error| format!("error reading {input}: {error}"))?;
        // FragGeneScanRs names genes after the first word of the header.
        let id = record.id();
        let contig = id.split(u8::is_ascii_whitespace).next().unwrap_or(id);
        let contig = String::from_utf8_lossy(contig).into_owned();
        let sequence = record.seq();
        for (offset, bases) in
            acgt_stretches(&sequence).filter(|(_, bases)| bases.len() >= MIN_STRETCH)
        {
            write_fasta_record(output, stretches.len().to_string().as_bytes(), bases)?;
            stretches.push(Stretch {
                contig: contig.clone(),
                offset,
            });
        }
    }
    output.flush()?;
    Ok(stretches)
}

/// The non-empty maximal runs of A, C, G and T (either case) in `sequence`,
/// with their offsets.
fn acgt_stretches(sequence: &[u8]) -> impl Iterator<Item = (usize, &[u8])> {
    sequence
        .split(|base| !matches!(base.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T'))
        .scan(0, |start, bases| {
            let offset = *start;
            *start += bases.len() + 1;
            Some((offset, bases))
        })
        .filter(|(_, bases)| !bases.is_empty())
}

/// Rewrite FragGeneScanRs's FASTA headers and GFF rows from stretch records
/// to contig names and coordinates.
fn restore_coordinates(path: &Path, stretches: &[Stretch]) -> Result<(), Box<dyn Error>> {
    let text = fs::read_to_string(path)?;
    let mut out = BufWriter::new(File::create(path)?);
    for line in text.lines() {
        if let Some(gene) = line.strip_prefix('>') {
            writeln!(out, ">{}", contig_gene_id(gene, stretches)?)?;
        } else if line.contains('\t') {
            // record, source, type, start, end, score, strand, phase, `ID={gene};...`
            let mut fields = line.split('\t').map(str::to_owned).collect::<Vec<_>>();
            let stretch = stretch(&fields[0], stretches)?;
            fields[0] = stretch.contig.clone();
            for position in &mut fields[3..5] {
                *position = (position.parse::<usize>()? + stretch.offset).to_string();
            }
            if let Some((gene, rest)) = fields[8]
                .strip_prefix("ID=")
                .and_then(|a| a.split_once(';'))
            {
                fields[8] = format!("ID={};{rest}", contig_gene_id(gene, stretches)?);
            }
            writeln!(out, "{}", fields.join("\t"))?;
        } else {
            writeln!(out, "{line}")?;
        }
    }
    Ok(out.flush()?)
}

/// `{record}_{start}_{end}_{strand}` from FragGeneScanRs, with the record a
/// stretch index, as `{contig}_{start}_{end}_{strand}` in contig coordinates.
fn contig_gene_id(gene: &str, stretches: &[Stretch]) -> Result<String, Box<dyn Error>> {
    let fields = gene.splitn(4, '_').collect::<Vec<_>>();
    let [record, start, end, strand] = fields[..] else {
        return Err(format!("unexpected FragGeneScanRs gene ID {gene}").into());
    };
    let stretch = stretch(record, stretches)?;
    let start = start.parse::<usize>()? + stretch.offset;
    let end = end.parse::<usize>()? + stretch.offset;
    Ok(format!("{}_{start}_{end}_{strand}", stretch.contig))
}

fn stretch<'a>(record: &str, stretches: &'a [Stretch]) -> Result<&'a Stretch, Box<dyn Error>> {
    record
        .parse::<usize>()
        .ok()
        .and_then(|index| stretches.get(index))
        .ok_or_else(|| format!("unexpected FragGeneScanRs record {record}").into())
}

fn is_train_dir(dir: &Path) -> bool {
    dir.join(TRAINING_MODEL).is_file()
}

fn default_train_dirs(program: &Path) -> Vec<PathBuf> {
    let mut dirs = Vec::new();
    if let Some(bin) = program.parent() {
        dirs.push(bin.join("train"));
        dirs.push(bin.join("../share/FragGeneScanRs/train"));
    }
    // `cargo install frag_gene_scan_rs` keeps the crate sources, which
    // include the models, in the registry.
    let cargo_home = env::var_os("CARGO_HOME")
        .map(PathBuf::from)
        .or_else(|| env::var_os("HOME").map(|home| Path::new(&home).join(".cargo")));
    if let Some(registry) = cargo_home.map(|home| home.join("registry").join("src")) {
        let mut crates = fs::read_dir(registry)
            .into_iter()
            .flatten()
            .flatten()
            .flat_map(|index| fs::read_dir(index.path()).into_iter().flatten().flatten())
            .map(|entry| entry.path())
            .filter(|path| {
                path.file_name()
                    .and_then(|name| name.to_str())
                    .is_some_and(|name| name.starts_with("frag_gene_scan_rs-"))
            })
            .collect::<Vec<_>>();
        // Newest version last, so prefer it.
        crates.sort();
        dirs.extend(crates.into_iter().rev().map(|path| path.join("train")));
    }
    dirs
}

/// Decompress and normalise a FASTA/FASTQ file to plain FASTA; returns its base count.
pub(crate) fn write_plain_fasta(
    input: &str,
    output: &mut dyn Write,
) -> Result<usize, Box<dyn Error>> {
    let mut reader =
        parse_fastx_file(input).map_err(|error| format!("invalid FASTA {input}: {error}"))?;
    let mut total = 0;
    while let Some(record) = reader.next() {
        let record = record.map_err(|error| format!("error reading {input}: {error}"))?;
        let sequence = record.seq();
        total += sequence.len();
        write_fasta_record(output, record.id(), &sequence)?;
    }
    output.flush()?;
    Ok(total)
}

/// `write_plain_fasta` into the file `output`.
pub(crate) fn write_plain_fasta_file(input: &str, output: &Path) -> Result<usize, Box<dyn Error>> {
    write_plain_fasta(input, &mut BufWriter::new(File::create(output)?))
}

/// One FASTA record: `>id`, then the sequence on one line.
pub(crate) fn write_fasta_record(
    writer: &mut (impl Write + ?Sized),
    id: &[u8],
    sequence: &[u8],
) -> std::io::Result<()> {
    writer.write_all(b">")?;
    writer.write_all(id)?;
    writer.write_all(b"\n")?;
    writer.write_all(sequence)?;
    writer.write_all(b"\n")
}

pub(crate) fn count_records(fasta: &Path) -> Result<usize, Box<dyn Error>> {
    let mut count = 0;
    for line in BufReader::new(File::open(fasta)?).lines() {
        if line?.starts_with('>') {
            count += 1;
        }
    }
    Ok(count)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn stretches_skip_ambiguous_bases() {
        let stretches = acgt_stretches(b"NNacGTNRTTANN").collect::<Vec<_>>();
        assert_eq!(stretches, [(2, &b"acGT"[..]), (8, &b"TTA"[..])]);
    }

    #[test]
    fn gene_ids_are_restored_to_contig_coordinates() {
        let stretches = [
            Stretch {
                contig: "contig_1".into(),
                offset: 0,
            },
            Stretch {
                contig: "contig_1".into(),
                offset: 120,
            },
        ];
        assert_eq!(
            contig_gene_id("1_5_310_-", &stretches).unwrap(),
            "contig_1_125_430_-"
        );
        assert!(contig_gene_id("2_5_310_-", &stretches).is_err());
    }
}
