//! Gene calling with FragGeneScanRs's complete-genome model, run as a
//! separate program (FragGeneScanRs is GPL-3.0; aeon is MIT).

use needletail::parse_fastx_file;
use std::env;
use std::error::Error;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::Command;

use crate::classify::external::{find_program, run};

/// FragGeneScanRs translates with the standard bacterial code.
pub(crate) const TRANSLATION_TABLE: u8 = 11;
const TRAINING_MODEL: &str = "complete";

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct GeneCalls {
    pub(crate) protein_faa: PathBuf,
    pub(crate) proteins: usize,
}

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
    /// `{id}_protein.faa`, `{id}_protein.fna` and `{id}_protein.gff` like
    /// GTDB-Tk's `identify` step. Returns `None` for a genome without
    /// sequence or genes. FragGeneScanRs splits the contigs over `threads`.
    pub(crate) fn call(
        &self,
        genome_path: &str,
        genome_id: &str,
        out_dir: &Path,
        threads: usize,
    ) -> Result<Option<GeneCalls>, Box<dyn Error>> {
        let genome_dir = out_dir.join(genome_id);
        fs::create_dir_all(&genome_dir)?;
        let input = genome_dir.join(format!("{genome_id}.input.fna"));
        let total_bases = write_plain_fasta(genome_path, &input)?;
        if total_bases == 0 {
            fs::remove_file(&input)?;
            return Ok(None);
        }
        // Genome IDs such as GCF_000005845.2 contain dots, so names are built
        // explicitly rather than with `Path::with_extension`.
        let output = |suffix: &str| genome_dir.join(format!("{genome_id}{suffix}"));
        let protein_faa = output("_protein.faa");
        run(Command::new(&self.program)
            .args(["-w", "1", "-t", TRAINING_MODEL])
            .args(["-p", &threads.max(1).to_string()])
            .arg("-r")
            .arg(&self.train_dir)
            .arg("-s")
            .arg(&input)
            .arg("-a")
            .arg(&protein_faa)
            .arg("-n")
            .arg(output("_protein.fna"))
            .arg("-g")
            .arg(output("_protein.gff")))?;
        fs::remove_file(&input)?;
        let proteins = count_records(&protein_faa)?;
        Ok((proteins > 0).then_some(GeneCalls {
            protein_faa,
            proteins,
        }))
    }
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
pub(crate) fn write_plain_fasta(input: &str, output: &Path) -> Result<usize, Box<dyn Error>> {
    let mut reader =
        parse_fastx_file(input).map_err(|error| format!("invalid FASTA {input}: {error}"))?;
    let mut writer = BufWriter::new(File::create(output)?);
    let mut total = 0;
    while let Some(record) = reader.next() {
        let record = record.map_err(|error| format!("error reading {input}: {error}"))?;
        let sequence = record.seq();
        total += sequence.len();
        writer.write_all(b">")?;
        writer.write_all(record.id())?;
        writer.write_all(b"\n")?;
        writer.write_all(&sequence)?;
        writer.write_all(b"\n")?;
    }
    writer.flush()?;
    Ok(total)
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
