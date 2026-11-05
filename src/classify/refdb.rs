//! The `aeon classify-db` reference database: three aeon indexes over GTDB
//! species representatives (genome nucleotides, proteomes, and concatenated
//! marker genes) plus the taxonomy and radii needed to report GTDB-Tk-style
//! classifications.

use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, HashMap, HashSet};
use std::error::Error;
use std::fs;
use std::path::{Path, PathBuf};
use std::time::Instant;

use needletail::parse_fastx_file;

use crate::classify::external::find_tool;
use crate::classify::genes::GeneCaller;
use crate::classify::markers::{
    MARKER_SETS, MarkerTools, PfamInfo, apply_mask, extract_masked_markers, marker_sequence,
    read_mask, reference_msa_path,
};
use crate::classify::refdata::{
    GtdbtkData, accession_from_path, canonical_gid, marker_set_for_taxonomy, strip_db_prefix,
};
use crate::database::build_database;
use crate::metadata::{PrefixParams, genomes_path, load_params};
use crate::utils::read_list_file;
use rayon::prelude::*;
use std::sync::atomic::{AtomicUsize, Ordering};

const MANIFEST_FILE: &str = "manifest.json";
const NT_PREFIX: &str = "nt";
const AA_PREFIX: &str = "aa";
const MARKER_PREFIX: &str = "markers";
const PROTEIN_DIR: &str = "proteins";
const MARKER_DIR: &str = "marker_genes";

/// Where reference proteomes for the AAI and marker indexes come from.
pub(crate) enum ProteinSource {
    /// Precomputed proteomes (e.g. GTDB `protein_faa_reps`), named by accession.
    List(String),
    /// Call genes for every reference genome with FragGeneScanRs.
    CallGenes(GeneCaller),
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub(crate) struct Manifest {
    pub(crate) aeon_version: String,
    pub(crate) gtdb_release: String,
    pub(crate) gtdbtk_data: String,
    pub(crate) references: usize,
    /// References with markers in the marker index.
    pub(crate) marker_references: usize,
}

/// One aeon index of the database.
pub(crate) struct Index {
    pub(crate) prefix: String,
    pub(crate) params: PrefixParams,
    /// Indexed files in vector-ID order.
    pub(crate) paths: Vec<String>,
    /// GTDB accessions (`GCF_000005845.2`) aligned with `paths`.
    pub(crate) accessions: Vec<String>,
}

impl Index {
    fn open(dir: &Path, name: &str) -> Result<Self, Box<dyn Error>> {
        let prefix = dir.join(name).to_string_lossy().into_owned();
        let paths = read_list_file(&genomes_path(&prefix))?;
        Ok(Self {
            params: load_params(&prefix)?,
            accessions: paths.iter().map(|path| accession_from_path(path)).collect(),
            paths,
            prefix,
        })
    }

    pub(crate) fn len(&self) -> usize {
        self.paths.len()
    }
}

/// An opened classify database.
pub(crate) struct ClassifyDb {
    pub(crate) manifest: Manifest,
    pub(crate) data: GtdbtkData,
    /// Reference genomes.
    pub(crate) nt: Index,
    /// Reference proteomes, aligned with `nt` vector IDs.
    pub(crate) aa: Index,
    /// Concatenated marker genes of the references that have them.
    pub(crate) markers: Index,
}

impl ClassifyDb {
    pub(crate) fn open(dir: &str) -> Result<Self, Box<dyn Error>> {
        let dir = PathBuf::from(dir);
        let manifest: Manifest = serde_json::from_str(
            &fs::read_to_string(dir.join(MANIFEST_FILE)).map_err(|error| {
                format!(
                    "{} is not an aeon classify database: {error}",
                    dir.display()
                )
            })?,
        )
        .map_err(|error| {
            format!(
                "{} was built by an older aeon classify-db; rebuild it: {error}",
                dir.display()
            )
        })?;
        let (nt, aa) = (Index::open(&dir, NT_PREFIX)?, Index::open(&dir, AA_PREFIX)?);
        if nt.len() != aa.len() {
            return Err("amino-acid and nucleotide indexes have different sizes".into());
        }
        Ok(Self {
            data: GtdbtkData::load(&manifest.gtdbtk_data)?,
            markers: Index::open(&dir, MARKER_PREFIX)?,
            nt,
            aa,
            manifest,
        })
    }

    pub(crate) fn root(&self) -> &Path {
        Path::new(&self.manifest.gtdbtk_data)
    }
}

/// Index parameters of the three indexes.
pub(crate) struct IndexParams {
    pub(crate) nt: PrefixParams,
    pub(crate) aa: PrefixParams,
    pub(crate) markers: PrefixParams,
}

pub(crate) fn build_classify_db(
    db_dir: &str,
    genome_list: &str,
    gtdbtk_data: &str,
    proteins: ProteinSource,
    bin_dir: Option<&str>,
    params: IndexParams,
) -> Result<(), Box<dyn Error>> {
    let started = Instant::now();
    let genomes = absolute_paths(read_list_file(genome_list)?)?;
    if genomes.is_empty() {
        return Err("reference genome list is empty".into());
    }
    let root = std::path::absolute(gtdbtk_data)?;
    let data = GtdbtkData::load(&root.to_string_lossy())?;
    eprintln!(
        "Loaded GTDB {} reference data: {} taxonomy records, {} species radii",
        data.version,
        data.taxonomy.len(),
        data.radii.len()
    );

    let accessions = genomes
        .iter()
        .map(|path| accession_from_path(path))
        .collect::<Vec<_>>();
    check_references(&accessions, &data)?;

    let db = &std::path::absolute(db_dir)?;
    fs::create_dir_all(db)?;
    // The nucleotide index needs only the genomes, so it builds while the
    // proteomes, protein index and marker index are prepared.
    let (nt, rest) = rayon::join(
        || {
            build_database(&db.join(NT_PREFIX).to_string_lossy(), &genomes, params.nt)
                .map_err(|error| error.to_string())
        },
        || -> Result<_, String> {
            let proteomes = match proteins {
                ProteinSource::List(list) => read_list_file(&list)
                    .map_err(Into::into)
                    .and_then(absolute_paths)
                    .and_then(|paths| match_proteomes(&accessions, &paths)),
                ProteinSource::CallGenes(caller) => {
                    call_reference_genes(&caller, &genomes, &accessions, db)
                }
            }
            .map_err(|error| error.to_string())?;
            let (aa, markers) = rayon::join(
                || {
                    build_database(&db.join(AA_PREFIX).to_string_lossy(), &proteomes, params.aa)
                        .map_err(|error| error.to_string())
                },
                || {
                    let markers =
                        write_reference_markers(&root, &data, &accessions, &proteomes, bin_dir, db)
                            .map_err(|error| error.to_string())?;
                    build_database(
                        &db.join(MARKER_PREFIX).to_string_lossy(),
                        &markers,
                        params.markers,
                    )
                    .map_err(|error| error.to_string())?;
                    Ok::<_, String>(markers)
                },
            );
            aa?;
            markers
        },
    );
    nt?;
    let markers = rest?;
    let manifest = Manifest {
        aeon_version: env!("CARGO_PKG_VERSION").into(),
        gtdb_release: data.version.clone(),
        gtdbtk_data: root.to_string_lossy().into_owned(),
        references: genomes.len(),
        marker_references: markers.len(),
    };
    fs::write(
        db.join(MANIFEST_FILE),
        serde_json::to_string_pretty(&manifest)?,
    )?;
    eprintln!(
        "Built classify database: dir={}, references={}, with markers={}, elapsed={:.3}s",
        db_dir,
        genomes.len(),
        markers.len(),
        started.elapsed().as_secs_f64()
    );
    Ok(())
}

/// Write each reference's masked marker residues to `marker_genes/{accession}.faa`
/// and return the files in reference order. Rows come from the GTDB-Tk
/// package's reference MSA; references missing from it (custom or test
/// packages) get theirs from hmmsearch and hmmalign on their proteomes.
fn write_reference_markers(
    root: &Path,
    data: &GtdbtkData,
    accessions: &[String],
    proteomes: &[String],
    bin_dir: Option<&str>,
    db: &Path,
) -> Result<Vec<String>, Box<dyn Error>> {
    let started = Instant::now();
    let marker_dir = db.join(MARKER_DIR);
    fs::create_dir_all(&marker_dir)?;
    let file_of = |accession: &str| marker_dir.join(format!("{accession}.faa"));
    let write = |accession: &str, masked: &str| -> std::io::Result<bool> {
        let sequence = marker_sequence(masked);
        if sequence.is_empty() {
            return Ok(false);
        }
        fs::write(file_of(accession), format!(">{accession}\n{sequence}\n"))?;
        Ok(true)
    };
    let mut with_markers = HashSet::new();
    let mut from_package = 0;
    for set in MARKER_SETS {
        let mut pending = accessions
            .iter()
            .zip(proteomes)
            .filter(|(accession, _)| {
                data.taxonomy_of(accession)
                    .is_some_and(|taxonomy| marker_set_for_taxonomy(taxonomy) == set)
            })
            .map(|(accession, proteome)| (canonical_gid(accession), (accession, proteome)))
            .collect::<HashMap<_, _>>();
        if pending.is_empty() {
            continue;
        }
        let mask = read_mask(data, root, set)?;
        // The R232 bac120 MSA is several GB, so rows are streamed.
        let msa = reference_msa_path(root, &data.version, set);
        if fs::metadata(&msa).is_ok_and(|metadata| metadata.len() > 0) {
            let mut reader = parse_fastx_file(&msa)
                .map_err(|error| format!("cannot read {}: {error}", msa.display()))?;
            while let Some(record) = reader.next() {
                let record = record?;
                let id = String::from_utf8_lossy(record.id());
                let id = strip_db_prefix(id.split_whitespace().next().unwrap_or_default());
                let Some((accession, _)) = pending.remove(&canonical_gid(id)) else {
                    continue;
                };
                let (masked, _, _) = apply_mask(&String::from_utf8_lossy(&record.seq()), &mask)?;
                if write(accession, &masked)? {
                    with_markers.insert(accession.clone());
                    from_package += 1;
                }
                if pending.is_empty() {
                    break;
                }
            }
        }
        if pending.is_empty() {
            continue;
        }
        eprintln!(
            "Extracting {set} markers with HMMER for {} reference(s) missing from {}",
            pending.len(),
            msa.display()
        );
        let pfam_info = PfamInfo::load(root)?;
        let (hmmsearch, hmmalign) = (
            find_tool(bin_dir, "hmmsearch")?,
            find_tool(bin_dir, "hmmalign")?,
        );
        let tools = MarkerTools {
            hmmsearch: &hmmsearch,
            hmmalign: &hmmalign,
            root,
            pfam_info: &pfam_info,
        };
        let inputs = pending
            .values()
            .map(|(accession, proteome)| ((*accession).clone(), PathBuf::from(proteome)))
            .collect::<Vec<_>>();
        let work_dir = db.join(format!("{MARKER_DIR}_tmp"));
        let masked: BTreeMap<String, String> = extract_masked_markers(
            &tools,
            &mask,
            set,
            &inputs,
            &work_dir,
            rayon::current_num_threads(),
        )?;
        fs::remove_dir_all(&work_dir)?;
        for (accession, masked) in masked {
            if write(&accession, &masked)? {
                with_markers.insert(accession);
            }
        }
    }
    eprintln!(
        "Wrote marker genes for {}/{} references ({from_package} from the package MSA) in {:.3}s",
        with_markers.len(),
        accessions.len(),
        started.elapsed().as_secs_f64()
    );
    Ok(accessions
        .iter()
        .filter(|accession| with_markers.contains(*accession))
        .map(|accession| file_of(accession).to_string_lossy().into_owned())
        .collect())
}

/// Paths stored in the database must not depend on the directory classify
/// runs from.
fn absolute_paths(paths: Vec<String>) -> Result<Vec<String>, Box<dyn Error>> {
    paths
        .into_iter()
        .map(|path| Ok(std::path::absolute(&path)?.to_string_lossy().into_owned()))
        .collect()
}

/// Order proteome files to match the reference genomes by accession.
fn match_proteomes(
    accessions: &[String],
    proteomes: &[String],
) -> Result<Vec<String>, Box<dyn Error>> {
    let by_accession = proteomes
        .iter()
        .map(|path| (canonical_gid(&accession_from_path(path)), path))
        .collect::<HashMap<_, _>>();
    accessions
        .iter()
        .map(|accession| {
            by_accession
                .get(&canonical_gid(accession))
                .map(|path| (*path).clone())
                .ok_or_else(|| format!("no proteome in the protein list for {accession}").into())
        })
        .collect()
}

/// Call genes for every reference and keep only the proteomes.
fn call_reference_genes(
    caller: &GeneCaller,
    genomes: &[String],
    accessions: &[String],
    db: &Path,
) -> Result<Vec<String>, Box<dyn Error>> {
    let started = Instant::now();
    let protein_dir = db.join(PROTEIN_DIR);
    let done = AtomicUsize::new(0);
    let threads = (rayon::current_num_threads() / genomes.len()).max(1);
    let proteomes = genomes
        .par_iter()
        .zip(accessions)
        .map(|(genome, accession)| -> Result<String, String> {
            let calls = caller
                .call(genome, accession, &protein_dir, threads)
                .map_err(|error| format!("gene calling failed for {genome}: {error}"))?
                .ok_or_else(|| format!("no genes were called for reference genome {genome}"))?;
            let genome_dir = protein_dir.join(accession);
            for suffix in ["_protein.fna", "_protein.gff"] {
                let _ = fs::remove_file(genome_dir.join(format!("{accession}{suffix}")));
            }
            let finished = done.fetch_add(1, Ordering::Relaxed) + 1;
            if finished % 100 == 0 {
                eprintln!("Called genes for {finished}/{} references", genomes.len());
            }
            Ok(calls.protein_faa.to_string_lossy().into_owned())
        })
        .collect::<Result<Vec<_>, _>>()?;
    eprintln!(
        "Called genes for {} references in {:.3}s",
        genomes.len(),
        started.elapsed().as_secs_f64()
    );
    Ok(proteomes)
}

/// Every reference must be a GTDB species representative with a taxonomy
/// and an ANI radius, and accessions must be unique.
fn check_references(accessions: &[String], data: &GtdbtkData) -> Result<(), Box<dyn Error>> {
    let mut seen = HashSet::with_capacity(accessions.len());
    let mut missing = Vec::new();
    for accession in accessions {
        if !seen.insert(canonical_gid(accession)) {
            return Err(format!("reference list contains accession {accession} twice").into());
        }
        if data.taxonomy_of(accession).is_none() || data.radius_of(accession).is_none() {
            missing.push(accession.as_str());
        }
    }
    if missing.is_empty() {
        return Ok(());
    }
    Err(format!(
        "{} reference genome(s) are not GTDB {} species representatives with a taxonomy and ANI radius (e.g. {}). Reference files must be named by accession, like GCF_000005845.2_genomic.fna.gz",
        missing.len(),
        data.version,
        missing[..missing.len().min(5)].join(", ")
    )
    .into())
}
