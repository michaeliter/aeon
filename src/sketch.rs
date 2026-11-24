use needletail::errors::{ParseError, ParseErrorKind};
use needletail::{FastxReader, Sequence, parse_fastx_file};
use num_traits::{NumCast, PrimInt, ToPrimitive};
use rayon::prelude::*;
use sha2::{Digest, Sha256};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::thread;
use std::time::{Duration, Instant};
use xxhash_rust::xxh3::xxh3_64_with_seed;

// DNA / base k-mer machinery
use kmerutils::base::{
    CompressedKmerT, Kmer16b32bit, Kmer32bit, Kmer64bit, KmerBuilder,
    alphabet::Alphabet2b,
    kmergenerator::{KmerGenerationPattern, KmerGenerator},
    sequence::Sequence as SequenceStruct,
};
use kmerutils::sketcharg::{DataType, SeqSketcherParams, SketchAlgo};
use kmerutils::sketching::setsketchert::{OptDensHashSketch, RevOptDensHashSketch, SeqSketcherT};

// Amino-acid k-mers and AA sketchers (AA-specific generator + pattern)
use kmerutils::aautils::kmeraa::{
    Alphabet as AAAlphabet, KmerAA32bit, KmerAA64bit,
    KmerGenerationPattern as AAKmerGenerationPattern, KmerGenerator as AAKmerGenerator, SequenceAA,
};

use kmerutils::aautils::setsketchert as aasketch;
use kmerutils::aautils::setsketchert::SeqSketcherAAT;

use crate::metadata::PrefixParams;

/// SHA-256 of the complete, normalized content that Aeon passes to its
/// sketcher. The digest includes record boundaries and therefore represents a
/// genome rather than an unordered collection of bases.
pub(crate) type GenomeContentHash = [u8; 32];

#[derive(Debug)]
pub(crate) struct SketchedGenome {
    pub(crate) vector: Vec<u16>,
    pub(crate) content_hash: GenomeContentHash,
}

/// Converts ASCII-encoded bases (from Needletail) into our `SequenceStruct`.
fn ascii_to_seq(bases: &[u8]) -> Result<SequenceStruct, ()> {
    let alphabet = Alphabet2b::new();
    let mut seq = SequenceStruct::with_capacity(2, bases.len());
    seq.encode_and_add(bases, &alphabet);
    Ok(seq)
}

/// Converts ASCII-encoded amino-acid sequence into `SequenceAA`.
/// - uppercases all residues
/// - drops any residue not in the 20-AA alphabet (ACDEFGHIKLMNPQRSTVWY),
///   including '*', 'X', gaps, etc.
fn ascii_to_seq_aa(bases: &[u8]) -> Result<SequenceAA, ()> {
    // Normalize to uppercase to match Alphabet's bases
    let mut upper = Vec::with_capacity(bases.len());
    for &b in bases {
        upper.push(b.to_ascii_uppercase());
    }

    let alphabet = AAAlphabet::new();
    Ok(SequenceAA::new_filtered(&upper, &alphabet))
}

/// Attempts at reading a sequence file before giving up: shared filesystems
/// such as NFS can fail reads transiently.
const READ_ATTEMPTS: u32 = 5;
/// Wait before the first retry; it doubles before each further retry.
const FIRST_RETRY_DELAY: Duration = Duration::from_secs(1);

/// Open `path` and collect its records with `read`, retrying the whole file
/// after an I/O error with exponential backoff. Format errors panic at once.
fn read_with_retries<T>(
    path: &str,
    first_delay: Duration,
    read: impl Fn(&mut dyn FastxReader) -> Result<Vec<T>, ParseError>,
) -> Vec<T> {
    let mut delay = first_delay;
    let mut attempt = 1;
    loop {
        let (result, opened) = match parse_fastx_file(path) {
            Ok(mut reader) => (read(reader.as_mut()), true),
            Err(error) => (Err(error), false),
        };
        match result {
            Ok(records) => return records,
            Err(error) if error.kind == ParseErrorKind::Io && attempt < READ_ATTEMPTS => {
                eprintln!(
                    "WARNING: I/O error reading {path} (attempt {attempt}/{READ_ATTEMPTS}), retrying in {delay:?}: {error}"
                );
                thread::sleep(delay);
                delay *= 2;
                attempt += 1;
            }
            Err(error) if opened => panic!("Error reading record in {path}: {error}"),
            Err(error) => panic!("Invalid FASTA/Q {path}: {error}"),
        }
    }
}

fn read_dna_sequences(path: &str) -> Vec<SequenceStruct> {
    read_with_retries(path, FIRST_RETRY_DELAY, |reader| {
        let mut sequences = Vec::new();
        while let Some(record) = reader.next() {
            let normalized = record?.normalize(false).into_owned();
            sequences.push(ascii_to_seq(&normalized).expect("DNA compression must succeed"));
        }
        Ok(sequences)
    })
}

fn read_aa_sequences(path: &str) -> Vec<SequenceAA> {
    read_with_retries(path, FIRST_RETRY_DELAY, |reader| {
        let mut sequences = Vec::new();
        let mut record_num = 0_u64;
        while let Some(record) = reader.next() {
            let record = record?;
            let bases = record.seq();
            let sequence_id = String::from_utf8_lossy(record.id());
            if bases.is_empty() {
                eprintln!(
                    "ERROR: sequence of null length, file: {path:?}, record num: {record_num}, sequence id: {sequence_id}"
                );
                record_num += 1;
                continue;
            }
            let sequence =
                ascii_to_seq_aa(&bases).unwrap_or_else(|_| panic!("AA parse error in {path}"));
            if sequence.is_empty() {
                eprintln!(
                    "ERROR: null encoded sequence, file: {path:?}, record num: {record_num}, sequence id: {sequence_id}"
                );
                record_num += 1;
                continue;
            }
            sequences.push(sequence);
            record_num += 1;
        }
        Ok(sequences)
    })
}

fn hash_dna_genome(sequences: &[SequenceStruct]) -> GenomeContentHash {
    let mut hasher = Sha256::new();
    hasher.update(b"aeon/genome-content/dna/v1\0");
    hasher.update((sequences.len() as u64).to_le_bytes());
    for sequence in sequences {
        // `decompress` reconstructs the exact canonical A/C/G/T stream from
        // kmerutils' two-bit representation, independently of FASTA headers,
        // case, wrapping, and invalid-base filtering.
        let bases = sequence.decompress();
        hasher.update((bases.len() as u64).to_le_bytes());
        hasher.update(&bases);
    }
    hasher.finalize().into()
}

fn hash_aa_genome(sequences: &[SequenceAA]) -> GenomeContentHash {
    let mut hasher = Sha256::new();
    hasher.update(b"aeon/genome-content/aa/v1\0");
    hasher.update((sequences.len() as u64).to_le_bytes());
    let mut buffer = [0_u8; 8192];
    for sequence in sequences {
        let length = sequence.size();
        hasher.update((length as u64).to_le_bytes());
        for start in (0..length).step_by(buffer.len()) {
            let width = (length - start).min(buffer.len());
            for (offset, byte) in buffer[..width].iter_mut().enumerate() {
                *byte = sequence.get_base(start + offset);
            }
            hasher.update(&buffer[..width]);
        }
    }
    hasher.finalize().into()
}

fn hash_dna_files_ordered(paths: &[String]) -> Vec<GenomeContentHash> {
    let completed = AtomicUsize::new(0);
    let total = paths.len();
    let mut indexed = paths
        .par_iter()
        .enumerate()
        .map(|(index, path)| {
            let hash = hash_dna_genome(&read_dna_sequences(path));
            let done = completed.fetch_add(1, Ordering::Relaxed) + 1;
            if done % 1000 == 0 {
                eprintln!("Hashed {done}/{total} genomes");
            }
            (index, hash)
        })
        .collect::<Vec<_>>();
    indexed.sort_by_key(|(index, _)| *index);
    indexed.into_iter().map(|(_, hash)| hash).collect()
}

fn hash_aa_files_ordered(paths: &[String]) -> Vec<GenomeContentHash> {
    let completed = AtomicUsize::new(0);
    let total = paths.len();
    let mut indexed = paths
        .par_iter()
        .enumerate()
        .map(|(index, path)| {
            let hash = hash_aa_genome(&read_aa_sequences(path));
            let done = completed.fetch_add(1, Ordering::Relaxed) + 1;
            if done % 1000 == 0 {
                eprintln!("Hashed {done}/{total} genomes");
            }
            (index, hash)
        })
        .collect::<Vec<_>>();
    indexed.sort_by_key(|(index, _)| *index);
    indexed.into_iter().map(|(_, hash)| hash).collect()
}

/// Strict, order-preserving parallel DNA sketching. Content hashes are computed
/// in the same parse/compression pass when requested by build or insertion.
fn sketch_files_ordered_u16<Kmer, Sketcher, F>(
    file_paths: &[String],
    sketcher: &Sketcher,
    kmer_hash_fn: F,
    include_content_hash: bool,
) -> Vec<SketchedGenome>
where
    Kmer: CompressedKmerT + KmerBuilder<Kmer> + Send + Sync,
    <Kmer as CompressedKmerT>::Val: num::PrimInt + Send + Sync + std::fmt::Debug,
    KmerGenerator<Kmer>: KmerGenerationPattern<Kmer>,
    Sketcher: SeqSketcherT<Kmer, Sig = u16> + Sync,
    F: Fn(&Kmer) -> <Kmer as CompressedKmerT>::Val + Send + Sync + Copy,
{
    let completed = AtomicUsize::new(0);
    let total = file_paths.len();
    // We keep a local Vec<(original_index, result)> and then sort by index
    // to guarantee that final vectors follow file_paths order exactly.
    let mut indexed: Vec<(usize, SketchedGenome)> = file_paths
        .par_iter()
        .enumerate()
        .map(|(i, path)| {
            let sequences = read_dna_sequences(path);
            let content_hash = include_content_hash
                .then(|| hash_dna_genome(&sequences))
                .unwrap_or([0_u8; 32]);

            let sequences_ref: Vec<&SequenceStruct> = sequences.iter().collect();

            // returns Vec<Vec<u16>>, and for "seqs" interface inner vec has size 1
            let signature = sketcher.sketch_compressedkmer_seqs(&sequences_ref, kmer_hash_fn);
            let sig_u16 = signature
                .first()
                .cloned()
                .unwrap_or_else(|| panic!("sketcher returned empty signature for {path}"));

            let done = completed.fetch_add(1, Ordering::Relaxed) + 1;
            if done % 1000 == 0 {
                eprintln!("Sketched {done}/{total} genomes");
            }

            (
                i,
                SketchedGenome {
                    vector: sig_u16,
                    content_hash,
                },
            )
        })
        .collect();

    indexed.sort_by_key(|(i, _)| *i);
    indexed.into_iter().map(|(_, record)| record).collect()
}

/// Strict, order-preserving parallel AA sketching (Sig = u16).
fn sketch_files_ordered_aa_u16<Kmer, Sketcher, F>(
    file_paths: &[String],
    sketcher: &Sketcher,
    kmer_hash_fn: F,
    include_content_hash: bool,
) -> Vec<SketchedGenome>
where
    Kmer: CompressedKmerT + KmerBuilder<Kmer> + Send + Sync,
    <Kmer as CompressedKmerT>::Val: num::PrimInt + Send + Sync + std::fmt::Debug,
    // *** AA-specific generator + pattern bound (this fixes E0277) ***
    AAKmerGenerator<Kmer>: AAKmerGenerationPattern<Kmer>,
    Sketcher: SeqSketcherAAT<Kmer, Sig = u16> + Sync,
    F: Fn(&Kmer) -> <Kmer as CompressedKmerT>::Val + Send + Sync + Copy,
{
    let completed = AtomicUsize::new(0);
    let total = file_paths.len();
    let mut indexed: Vec<(usize, SketchedGenome)> = file_paths
        .par_iter()
        .enumerate()
        .map(|(i, path)| {
            let sequences = read_aa_sequences(path);
            let content_hash = include_content_hash
                .then(|| hash_aa_genome(&sequences))
                .unwrap_or([0_u8; 32]);

            let sequences_ref: Vec<&SequenceAA> = sequences.iter().collect();

            // AA sketcher interface: sketch_compressedkmeraa_seqs
            let signature = sketcher.sketch_compressedkmeraa_seqs(&sequences_ref, kmer_hash_fn);
            let sig_u16 = signature
                .first()
                .cloned()
                .unwrap_or_else(|| panic!("sketcher returned empty signature for {path}"));

            let done = completed.fetch_add(1, Ordering::Relaxed) + 1;
            if done % 1000 == 0 {
                eprintln!("Sketched {done}/{total} genomes");
            }

            (
                i,
                SketchedGenome {
                    vector: sig_u16,
                    content_hash,
                },
            )
        })
        .collect();

    indexed.sort_by_key(|(i, _)| *i);
    indexed.into_iter().map(|(_, record)| record).collect()
}

/// Canonicalize (min(kmer, revcomp)) -> pack -> XXH3 -> return Kmer::Val (u32/u64).
/// Safe for k<=32 (never shifts by 64).
pub fn make_xxh3_canonical_kmer_hash_fn<Kmer>(
    seed: u64,
) -> impl Fn(&Kmer) -> Kmer::Val + Copy + Send + Sync
where
    Kmer: CompressedKmerT + KmerBuilder<Kmer> + Copy,
    Kmer::Val: PrimInt + ToPrimitive + NumCast,
{
    move |kmer: &Kmer| -> Kmer::Val {
        // canonicalize
        let rc = kmer.reverse_complement();
        let canonical = if rc < *kmer { rc } else { *kmer };

        // DNA packing: 2 bits per base
        let k: usize = canonical.get_nb_base() as usize;
        let bits: usize = 2usize * k;

        // safe mask (avoid 1<<64)
        let mask_u64: u64 = if bits >= 64 {
            u64::MAX
        } else {
            (1u64 << bits) - 1
        };

        // convert packed bits to u64 safely (works for u32/u64)
        let packed_u64: u64 = canonical
            .get_compressed_value()
            .to_u64()
            .expect("Kmer::Val must fit in u64")
            & mask_u64;

        // hash packed bytes
        let h64: u64 = xxh3_64_with_seed(&packed_u64.to_le_bytes(), seed);

        // output type:
        // - if Val is u32-ish => use low 32 bits
        // - if Val is u64-ish => keep full 64
        if std::mem::size_of::<Kmer::Val>() <= 4 {
            let low32 = h64 as u32 as u64; // intentional truncation
            NumCast::from(low32).expect("casting hash->Val failed")
        } else {
            NumCast::from(h64).expect("casting hash->Val failed")
        }
    }
}

/// Amino-acid packed k-mer hashing (no reverse-complement).
/// AA alphabet needs 5 bits (2^5 = 32) for 20 residues.
pub fn make_xxh3_aa_kmer_hash_fn<Kmer>(
    seed: u64,
) -> impl Fn(&Kmer) -> Kmer::Val + Copy + Send + Sync
where
    Kmer: CompressedKmerT + KmerBuilder<Kmer>,
    Kmer::Val: PrimInt + NumCast + ToPrimitive,
{
    const AA_BITS_PER_RESIDUE: usize = 5;

    move |kmer: &Kmer| -> Kmer::Val {
        let k: usize = kmer.get_nb_base() as usize;
        let bits: usize = AA_BITS_PER_RESIDUE * k;

        let mask_u64: u64 = if bits >= 64 {
            u64::MAX
        } else {
            (1u64 << bits) - 1
        };

        let packed_u64: u64 = kmer
            .get_compressed_value()
            .to_u64()
            .expect("Kmer::Val must be convertible to u64")
            & mask_u64;

        let h64: u64 = xxh3_64_with_seed(&packed_u64.to_le_bytes(), seed);

        let out_u64 = if std::mem::size_of::<Kmer::Val>() <= 4 {
            (h64 as u32) as u64
        } else {
            h64
        };

        NumCast::from(out_u64).expect("cast to Kmer::Val must succeed")
    }
}

/// Dispatch sketching by k-mer size.
/// Fixed Sig=u16 for the sketcher output (DNA).
fn sketch_with_kmer_dispatch_u16(
    paths: &[String],
    kmer_size: usize,
    sketch_size: usize,
    densification: usize, // 0 optdens, 1 revoptdens
    hash_seed: u64,
    include_content_hash: bool,
) -> Vec<SketchedGenome> {
    let sketch_args = SeqSketcherParams::new(
        kmer_size,
        sketch_size,
        SketchAlgo::OPTDENS, // label; actual is controlled by sketcher type
        DataType::DNA,
    );

    if kmer_size <= 14 {
        let kmer_hash_fn = make_xxh3_canonical_kmer_hash_fn::<Kmer32bit>(hash_seed);

        match densification {
            0 => {
                let sketcher = OptDensHashSketch::<Kmer32bit, f32>::new(&sketch_args);
                sketch_files_ordered_u16::<Kmer32bit, _, _>(
                    paths,
                    &sketcher,
                    kmer_hash_fn,
                    include_content_hash,
                )
            }
            1 => {
                let sketcher = RevOptDensHashSketch::<Kmer32bit, f32>::new(&sketch_args);
                sketch_files_ordered_u16::<Kmer32bit, _, _>(
                    paths,
                    &sketcher,
                    kmer_hash_fn,
                    include_content_hash,
                )
            }
            _ => panic!("densification must be 0 or 1"),
        }
    } else if kmer_size == 16 {
        let kmer_hash_fn = make_xxh3_canonical_kmer_hash_fn::<Kmer16b32bit>(hash_seed);

        match densification {
            0 => {
                let sketcher = OptDensHashSketch::<Kmer16b32bit, f32>::new(&sketch_args);
                sketch_files_ordered_u16::<Kmer16b32bit, _, _>(
                    paths,
                    &sketcher,
                    kmer_hash_fn,
                    include_content_hash,
                )
            }
            1 => {
                let sketcher = RevOptDensHashSketch::<Kmer16b32bit, f32>::new(&sketch_args);
                sketch_files_ordered_u16::<Kmer16b32bit, _, _>(
                    paths,
                    &sketcher,
                    kmer_hash_fn,
                    include_content_hash,
                )
            }
            _ => panic!("densification must be 0 or 1"),
        }
    } else if kmer_size <= 32 {
        let kmer_hash_fn = make_xxh3_canonical_kmer_hash_fn::<Kmer64bit>(hash_seed);

        match densification {
            0 => {
                let sketcher = OptDensHashSketch::<Kmer64bit, f32>::new(&sketch_args);
                sketch_files_ordered_u16::<Kmer64bit, _, _>(
                    paths,
                    &sketcher,
                    kmer_hash_fn,
                    include_content_hash,
                )
            }
            1 => {
                let sketcher = RevOptDensHashSketch::<Kmer64bit, f32>::new(&sketch_args);
                sketch_files_ordered_u16::<Kmer64bit, _, _>(
                    paths,
                    &sketcher,
                    kmer_hash_fn,
                    include_content_hash,
                )
            }
            _ => panic!("densification must be 0 or 1"),
        }
    } else {
        panic!("kmer_size cannot exceed 32 and must not be 15!");
    }
}

/// Amino-acid k-mer dispatcher (Sig = u16).
/// AA k must be <= 12 (KmerAA32bit for k<=6, KmerAA64bit for 7..=12).
fn sketch_aa_with_kmer_dispatch_u16(
    paths: &[String],
    kmer_size: usize,
    sketch_size: usize,
    densification: usize, // 0 optdens, 1 revoptdens
    hash_seed: u64,
    include_content_hash: bool,
) -> Vec<SketchedGenome> {
    let sketch_args = SeqSketcherParams::new(
        kmer_size,
        sketch_size,
        SketchAlgo::OPTDENS, // label; actual controlled by sketcher type
        DataType::AA,
    );

    if kmer_size <= 6 {
        let kmer_hash_fn = make_xxh3_aa_kmer_hash_fn::<KmerAA32bit>(hash_seed);

        match densification {
            0 => {
                let sketcher = aasketch::OptDensHashSketch::<KmerAA32bit, f32>::new(&sketch_args);
                sketch_files_ordered_aa_u16::<KmerAA32bit, _, _>(
                    paths,
                    &sketcher,
                    kmer_hash_fn,
                    include_content_hash,
                )
            }
            1 => {
                let sketcher =
                    aasketch::RevOptDensHashSketch::<KmerAA32bit, f32>::new(&sketch_args);
                sketch_files_ordered_aa_u16::<KmerAA32bit, _, _>(
                    paths,
                    &sketcher,
                    kmer_hash_fn,
                    include_content_hash,
                )
            }
            _ => panic!("densification must be 0 or 1"),
        }
    } else if kmer_size <= 12 {
        let kmer_hash_fn = make_xxh3_aa_kmer_hash_fn::<KmerAA64bit>(hash_seed);

        match densification {
            0 => {
                let sketcher = aasketch::OptDensHashSketch::<KmerAA64bit, f32>::new(&sketch_args);
                sketch_files_ordered_aa_u16::<KmerAA64bit, _, _>(
                    paths,
                    &sketcher,
                    kmer_hash_fn,
                    include_content_hash,
                )
            }
            1 => {
                let sketcher =
                    aasketch::RevOptDensHashSketch::<KmerAA64bit, f32>::new(&sketch_args);
                sketch_files_ordered_aa_u16::<KmerAA64bit, _, _>(
                    paths,
                    &sketcher,
                    kmer_hash_fn,
                    include_content_hash,
                )
            }
            _ => panic!("densification must be 0 or 1"),
        }
    } else {
        panic!("kmer_size for amino acids must be <= 12");
    }
}

fn sketch_records_from_params(
    paths: &[String],
    params: &PrefixParams,
    include_content_hash: bool,
) -> Vec<SketchedGenome> {
    let started = Instant::now();
    log::debug!(
        "sketching files={} seq_type={} kmer_size={} sketch_size={} densification={}",
        paths.len(),
        params.seq_type,
        params.kmer_size,
        params.sketch_size,
        params.densification
    );
    let sketches = match params.seq_type.as_str() {
        "dna" => sketch_with_kmer_dispatch_u16(
            paths,
            params.kmer_size,
            params.sketch_size,
            params.densification,
            params.hash_seed,
            include_content_hash,
        ),
        "aa" => sketch_aa_with_kmer_dispatch_u16(
            paths,
            params.kmer_size,
            params.sketch_size,
            params.densification,
            params.hash_seed,
            include_content_hash,
        ),
        other => panic!("unknown sequence type in parameters: {other}"),
    };
    eprintln!(
        "Sketched {}/{} genomes in {:.3}s",
        paths.len(),
        paths.len(),
        started.elapsed().as_secs_f64()
    );
    sketches
}

/// Sketch query genomes without the additional content-hash work needed by
/// build and update operations.
pub(crate) fn sketch_from_params(paths: &[String], params: &PrefixParams) -> Vec<Vec<u16>> {
    sketch_records_from_params(paths, params, false)
        .into_iter()
        .map(|record| record.vector)
        .collect()
}

/// Parse, kmerutils-compress, content-hash, and sketch genomes in one pass.
/// Build and insertion use this so content deduplication happens before graph
/// construction or Vamana insertion.
pub(crate) fn sketch_and_hash_from_params(
    paths: &[String],
    params: &PrefixParams,
) -> Vec<SketchedGenome> {
    sketch_records_from_params(paths, params, true)
}

/// Recomputes content hashes for a legacy index whose `idmap.tsv` has only the
/// original two fields. It deliberately does not sketch or allocate graph
/// vectors.
pub(crate) fn genome_content_hashes_from_params(
    paths: &[String],
    params: &PrefixParams,
) -> Vec<GenomeContentHash> {
    let started = Instant::now();
    let hashes = match params.seq_type.as_str() {
        "dna" => hash_dna_files_ordered(paths),
        "aa" => hash_aa_files_ordered(paths),
        other => panic!("unknown sequence type in parameters: {other}"),
    };
    eprintln!(
        "Hashed {}/{} genomes in {:.3}s",
        paths.len(),
        paths.len(),
        started.elapsed().as_secs_f64()
    );
    hashes
}

// Database operations, CLI parsing, and process orchestration live in their
// own modules. This module only owns sequence parsing and sketch generation.

#[cfg(test)]
mod tests {
    use super::{ascii_to_seq, ascii_to_seq_aa, hash_dna_genome, read_with_retries};
    use needletail::errors::{ParseError, ParseErrorKind};
    use std::cell::Cell;
    use std::time::Duration;
    use std::{fs, io};

    /// A one-record FASTA file in the temporary directory.
    fn fasta_file(name: &str) -> String {
        let path = std::env::temp_dir().join(format!("aeon-{}-{name}.fa", std::process::id()));
        fs::write(&path, ">r1\nACGT\n").unwrap();
        path.to_str().unwrap().to_owned()
    }

    #[test]
    fn reads_are_retried_after_io_errors() {
        let path = fasta_file("retry");
        let attempts = Cell::new(0);
        let records = read_with_retries(&path, Duration::ZERO, |reader| {
            attempts.set(attempts.get() + 1);
            if attempts.get() < 3 {
                return Err(ParseError::from(io::Error::other("transient")));
            }
            let mut ids = Vec::new();
            while let Some(record) = reader.next() {
                ids.push(record?.id().to_vec());
            }
            Ok(ids)
        });
        assert_eq!((records, attempts.get()), (vec![b"r1".to_vec()], 3));
        fs::remove_file(path).unwrap();
    }

    #[test]
    fn format_errors_are_not_retried() {
        let path = fasta_file("format");
        let attempts = Cell::new(0);
        let outcome = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            read_with_retries::<()>(&path, Duration::ZERO, |_| {
                attempts.set(attempts.get() + 1);
                let mut error = ParseError::from(io::Error::other("bad record"));
                error.kind = ParseErrorKind::InvalidStart;
                Err(error)
            })
        }));
        assert!(outcome.is_err());
        assert_eq!(attempts.get(), 1);
        fs::remove_file(path).unwrap();
    }

    #[test]
    fn aa_conversion_exposes_empty_and_fully_filtered_records() {
        assert!(ascii_to_seq_aa(b"").unwrap().is_empty());
        assert!(ascii_to_seq_aa(b"***XXX---").unwrap().is_empty());
        assert!(!ascii_to_seq_aa(b"MTEQIELIK").unwrap().is_empty());
    }

    #[test]
    fn dna_content_hash_tracks_compressed_content_and_record_boundaries() {
        let canonical = ascii_to_seq(b"ACGTN").unwrap();
        let case_variant = ascii_to_seq(b"acgt").unwrap();
        assert_eq!(
            hash_dna_genome(&[canonical]),
            hash_dna_genome(&[case_variant])
        );

        let first = ascii_to_seq(b"AC").unwrap();
        let second = ascii_to_seq(b"GT").unwrap();
        let single = ascii_to_seq(b"ACGT").unwrap();
        assert_ne!(
            hash_dna_genome(&[first, second]),
            hash_dna_genome(&[single])
        );
    }
}
