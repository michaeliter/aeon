//! Gene calling with pyrodigal, choosing between translation tables 4 and 11
//! exactly as GTDB-Tk's `Prodigal._producer` does.

use needletail::parse_fastx_file;
use std::collections::HashMap;
use std::error::Error;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::Command;

use crate::classify::external::run;
use crate::classify::summary::py_float;

/// Genomes shorter than this use Prodigal's metagenomic mode.
const MIN_SINGLE_MODE_BASES: usize = 100_000;

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct GeneCalls {
    pub(crate) protein_faa: PathBuf,
    pub(crate) translation_table: u8,
    pub(crate) coding_density_4: f64,
    pub(crate) coding_density_11: f64,
    pub(crate) proteins: usize,
}

/// GTDB-Tk's table choice: table 4 only when it codes markedly more of the genome.
pub(crate) fn best_translation_table(coding_density_4: f64, coding_density_11: f64) -> u8 {
    if coding_density_4 - coding_density_11 > 0.05 && coding_density_4 > 0.7 {
        4
    } else {
        11
    }
}

/// Call genes for one genome into `out_dir/{genome_id}/`, writing
/// `{id}_protein.faa`, `{id}_protein.fna`, `{id}_protein.gff`, and
/// `{id}_translation_table.tsv` like GTDB-Tk's `identify` step.
/// Returns `None` for a genome without sequence.
pub(crate) fn call_genes(
    pyrodigal: &Path,
    genome_path: &str,
    genome_id: &str,
    out_dir: &Path,
) -> Result<Option<GeneCalls>, Box<dyn Error>> {
    let genome_dir = out_dir.join(genome_id);
    fs::create_dir_all(&genome_dir)?;
    let input = genome_dir.join(format!("{genome_id}.input.fna"));
    let total_bases = write_plain_fasta(genome_path, &input)?;
    if total_bases == 0 {
        fs::remove_file(&input)?;
        return Ok(None);
    }
    let mode = if total_bases < MIN_SINGLE_MODE_BASES {
        "meta"
    } else {
        "single"
    };

    // Genome IDs such as GCF_000005845.2 contain dots, so names are built
    // explicitly rather than with `Path::with_extension`.
    let table_path = |table: u8, extension: &str| {
        genome_dir.join(format!("{genome_id}.table{table}.{extension}"))
    };
    let mut density = HashMap::new();
    for table in [4_u8, 11] {
        let (faa, fna, gff) = (
            table_path(table, "faa"),
            table_path(table, "fna"),
            table_path(table, "gff"),
        );
        run(Command::new(pyrodigal)
            .args(["-m", "-p", mode, "-f", "gff", "-g", &table.to_string()])
            .arg("-a")
            .arg(&faa)
            .arg("-d")
            .arg(&fna)
            .arg("-i")
            .arg(&input)
            .arg("-o")
            .arg(&gff))?;
        density.insert(table, coding_bases(&gff)? as f64 / total_bases as f64);
    }
    let (coding_density_4, coding_density_11) = (density[&4], density[&11]);
    let table = best_translation_table(coding_density_4, coding_density_11);

    let final_path = |suffix: &str| genome_dir.join(format!("{genome_id}{suffix}"));
    for table_used in [4_u8, 11] {
        for (extension, suffix) in [
            ("faa", "_protein.faa"),
            ("fna", "_protein.fna"),
            ("gff", "_protein.gff"),
        ] {
            let path = table_path(table_used, extension);
            if table_used == table {
                fs::rename(&path, final_path(suffix))?;
            } else {
                fs::remove_file(&path)?;
            }
        }
    }
    fs::remove_file(&input)?;
    fs::write(
        final_path("_translation_table.tsv"),
        format!(
            "best_translation_table\t{table}\ncoding_density_4\t{}\ncoding_density_11\t{}\n",
            py_float(coding_density_4),
            py_float(coding_density_11)
        ),
    )?;
    let protein_faa = final_path("_protein.faa");
    Ok(Some(GeneCalls {
        proteins: count_records(&protein_faa)?,
        protein_faa,
        translation_table: table,
        coding_density_4,
        coding_density_11,
    }))
}

/// Decompress and normalise a FASTA/FASTQ file to plain FASTA; returns its base count.
fn write_plain_fasta(input: &str, output: &Path) -> Result<usize, Box<dyn Error>> {
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

/// Coding bases in a Prodigal GFF, reproducing GTDB-Tk's
/// `ProdigalGeneFeatureParser`: per sequence, a mask of length equal to the
/// last coding base, with each gene marking indices `start..=end` (1-based
/// coordinates used as 0-based indices).
pub(crate) fn coding_bases(gff: &Path) -> Result<usize, Box<dyn Error>> {
    let mut genes: HashMap<String, Vec<(usize, usize)>> = HashMap::new();
    for line in BufReader::new(File::open(gff)?).lines() {
        let line = line?;
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let fields = line.split('\t').collect::<Vec<_>>();
        if fields.len() < 5 {
            continue;
        }
        genes
            .entry(fields[0].to_owned())
            .or_default()
            .push((fields[3].parse()?, fields[4].parse()?));
    }
    Ok(genes.values().map(|genes| masked_length(genes)).sum())
}

fn masked_length(genes: &[(usize, usize)]) -> usize {
    let length = genes.iter().map(|&(_, end)| end).max().unwrap_or(0);
    let mut mask = vec![false; length];
    for &(start, end) in genes {
        let stop = (end + 1).min(length);
        if start < stop {
            mask[start..stop].fill(true);
        }
    }
    mask.into_iter().filter(|&coding| coding).count()
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
    fn table_4_needs_a_clear_coding_density_gain() {
        assert_eq!(best_translation_table(0.90, 0.80), 4);
        assert_eq!(best_translation_table(0.90, 0.86), 11);
        assert_eq!(best_translation_table(0.69, 0.50), 11);
    }

    #[test]
    fn coding_mask_matches_gtdbtk_indexing() {
        // numpy: mask = zeros(10); mask[2:6] = True; mask[5:11] = True -> indices 2..=9
        assert_eq!(masked_length(&[(2, 5), (5, 10)]), 8);
        // A single gene loses its last base to the mask length.
        assert_eq!(masked_length(&[(1, 100)]), 99);
        assert_eq!(masked_length(&[]), 0);
    }
}
