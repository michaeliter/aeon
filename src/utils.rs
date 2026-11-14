use rayon::ThreadPoolBuilder;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::sync::Once;

static INIT_RAYON: Once = Once::new();

pub(crate) fn init_rayon_global(threads: usize) {
    INIT_RAYON.call_once(|| {
        ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()
            .expect("failed to initialize the global Rayon thread pool");
    });
    log::debug!("Rayon thread pool ready threads={threads}");
}

pub(crate) fn requested_threads(matches: &clap::ArgMatches) -> usize {
    matches
        .get_one::<usize>("threads")
        .copied()
        .unwrap_or_else(num_cpus::get)
}

pub(crate) fn read_list_file(filepath: &str) -> io::Result<Vec<String>> {
    let reader = BufReader::new(File::open(filepath)?);
    reader
        .lines()
        .map(|line| line.map(|value| value.trim().to_owned()))
        .filter(|line| line.as_ref().map_or(true, |value| !value.is_empty()))
        .collect()
}

pub(crate) fn write_genome_list(filepath: &str, genomes: &[String]) -> io::Result<()> {
    let mut writer = BufWriter::new(File::create(filepath)?);
    for genome in genomes {
        writeln!(writer, "{genome}")?;
    }
    writer.flush()
}

/// Write the durable genome metadata in current vector-ID order.
///
/// Older Aeon indexes have two fields (`id<TAB>path`).  New indexes add the
/// normalized compressed genome-content SHA-256 as a third field.  Keeping the
/// digest beside its path avoids a second sidecar while preserving the old
/// static DiskANN index format unchanged.
pub(crate) fn write_idmap_tsv(
    filepath: &str,
    genomes: &[String],
    hashes: Option<&[[u8; 32]]>,
) -> io::Result<()> {
    if let Some(hashes) = hashes {
        if hashes.len() != genomes.len() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "cannot write id map with {} genomes and {} content hashes",
                    genomes.len(),
                    hashes.len()
                ),
            ));
        }
    }

    let mut writer = BufWriter::new(File::create(filepath)?);
    for (id, genome) in genomes.iter().enumerate() {
        if let Some(hashes) = hashes {
            writeln!(writer, "{id}\t{genome}\t{}", hash_to_hex(&hashes[id]))?;
        } else {
            writeln!(writer, "{id}\t{genome}")?;
        }
    }
    writer.flush()
}

/// Read content hashes carried by the third field of `PREFIX.idmap.tsv`.
///
/// Returns `Ok(None)` for a valid legacy two-column map.  All rows must use
/// the same format and must exactly match the current `.genomes.txt` ordering;
/// this prevents a stale mapping from silently assigning a digest to the wrong
/// vector after a compacting update.
pub(crate) fn read_idmap_content_hashes(
    filepath: &str,
    genomes: &[String],
) -> io::Result<Option<Vec<[u8; 32]>>> {
    let reader = BufReader::new(File::open(filepath)?);
    let mut hashes = Vec::with_capacity(genomes.len());
    let mut contains_hashes = None;
    let mut rows = 0_usize;

    for (row, line) in reader.lines().enumerate() {
        let line = line?;
        let mut fields = line.split('\t');
        let id = fields.next().unwrap_or_default();
        let genome = fields.next().ok_or_else(|| {
            invalid_idmap(filepath, row + 1, "expected at least an ID and genome path")
        })?;
        let hash = fields.next();
        if fields.next().is_some() {
            return Err(invalid_idmap(
                filepath,
                row + 1,
                "expected two or three tab-separated fields",
            ));
        }

        let parsed_id = id.parse::<usize>().map_err(|_| {
            invalid_idmap(filepath, row + 1, "vector ID is not a non-negative integer")
        })?;
        if parsed_id != row {
            return Err(invalid_idmap(
                filepath,
                row + 1,
                &format!("expected vector ID {row}, found {parsed_id}"),
            ));
        }
        let expected_genome = genomes.get(row).ok_or_else(|| {
            invalid_idmap(
                filepath,
                row + 1,
                "contains more rows than PREFIX.genomes.txt",
            )
        })?;
        if genome != expected_genome {
            return Err(invalid_idmap(
                filepath,
                row + 1,
                "genome path does not match PREFIX.genomes.txt",
            ));
        }

        let row_has_hash = hash.is_some();
        match contains_hashes {
            Some(expected) if expected != row_has_hash => {
                return Err(invalid_idmap(
                    filepath,
                    row + 1,
                    "mixes two-column and three-column rows",
                ));
            }
            None => contains_hashes = Some(row_has_hash),
            _ => {}
        }
        if let Some(hash) = hash {
            hashes.push(parse_hash_hex(hash, filepath, row + 1)?);
        }
        rows += 1;
    }

    if rows != genomes.len() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "id map has {rows} rows but PREFIX.genomes.txt has {}: {filepath}",
                genomes.len()
            ),
        ));
    }

    match contains_hashes {
        Some(true) => Ok(Some(hashes)),
        Some(false) | None => Ok(None),
    }
}

fn invalid_idmap(filepath: &str, line: usize, message: &str) -> io::Error {
    io::Error::new(
        io::ErrorKind::InvalidData,
        format!("invalid id map {filepath} at line {line}: {message}"),
    )
}

fn hash_to_hex(hash: &[u8; 32]) -> String {
    const HEX: &[u8; 16] = b"0123456789abcdef";
    let mut output = String::with_capacity(hash.len() * 2);
    for byte in hash {
        output.push(HEX[(byte >> 4) as usize] as char);
        output.push(HEX[(byte & 0x0f) as usize] as char);
    }
    output
}

fn parse_hash_hex(value: &str, filepath: &str, line: usize) -> io::Result<[u8; 32]> {
    if value.len() != 64 {
        return Err(invalid_idmap(
            filepath,
            line,
            "content hash must contain exactly 64 hexadecimal characters",
        ));
    }
    let bytes = value.as_bytes();
    let mut hash = [0_u8; 32];
    for (index, target) in hash.iter_mut().enumerate() {
        let high = hex_nibble(bytes[index * 2]).ok_or_else(|| {
            invalid_idmap(
                filepath,
                line,
                "content hash contains a non-hexadecimal character",
            )
        })?;
        let low = hex_nibble(bytes[index * 2 + 1]).ok_or_else(|| {
            invalid_idmap(
                filepath,
                line,
                "content hash contains a non-hexadecimal character",
            )
        })?;
        *target = (high << 4) | low;
    }
    Ok(hash)
}

fn hex_nibble(byte: u8) -> Option<u8> {
    match byte {
        b'0'..=b'9' => Some(byte - b'0'),
        b'a'..=b'f' => Some(byte - b'a' + 10),
        b'A'..=b'F' => Some(byte - b'A' + 10),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::{read_idmap_content_hashes, write_idmap_tsv};
    use std::fs;

    fn temporary_path(name: &str) -> String {
        std::env::temp_dir()
            .join(format!("aeon-{name}-{}", std::process::id()))
            .to_string_lossy()
            .into_owned()
    }

    #[test]
    fn idmap_round_trips_content_hashes() {
        let path = temporary_path("idmap-content-hash.tsv");
        let genomes = vec!["first.fna".to_owned(), "second.fna".to_owned()];
        let hashes = [[0x12; 32], [0xab; 32]];

        write_idmap_tsv(&path, &genomes, Some(&hashes)).unwrap();
        assert_eq!(
            read_idmap_content_hashes(&path, &genomes).unwrap(),
            Some(hashes.to_vec())
        );

        fs::remove_file(path).unwrap();
    }

    #[test]
    fn legacy_two_column_idmap_has_no_hashes() {
        let path = temporary_path("idmap-legacy.tsv");
        let genomes = vec!["first.fna".to_owned(), "second.fna".to_owned()];

        write_idmap_tsv(&path, &genomes, None).unwrap();
        assert_eq!(read_idmap_content_hashes(&path, &genomes).unwrap(), None);

        fs::remove_file(path).unwrap();
    }
}
