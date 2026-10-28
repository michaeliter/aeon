//! The `aeon classify-db` reference database: a nucleotide aeon index over
//! GTDB species representatives plus the taxonomy and radii needed to report
//! GTDB-Tk-style classifications.

use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};
use std::error::Error;
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

use crate::classify::genes::call_genes;
use crate::classify::refdata::{
    GtdbtkData, accession_from_path, canonical_gid, read_sp_clusters_radii, read_taxonomy,
};
use crate::database::build_database;
use crate::metadata::{PrefixParams, genomes_path, load_params};
use crate::utils::read_list_file;
use rayon::prelude::*;
use std::sync::atomic::{AtomicUsize, Ordering};

const MANIFEST_FILE: &str = "manifest.json";
const TAXONOMY_FILE: &str = "taxonomy.tsv";
const RADII_FILE: &str = "radii.tsv";
const NT_PREFIX: &str = "nt";
const AA_PREFIX: &str = "aa";
const PROTEIN_DIR: &str = "proteins";

/// Where reference taxonomy and species radii come from.
pub(crate) enum TaxonomySource {
    /// A GTDB-Tk reference data package (GTDBTK_DATA_PATH).
    GtdbtkPackage(String),
    /// GTDB release files: `*_taxonomy_r*.tsv[.gz]` and `sp_clusters_r*.tsv`.
    ReleaseFiles {
        release: String,
        taxonomy: Vec<String>,
        sp_clusters: String,
    },
}

/// Where reference proteomes for the AAI tier come from.
pub(crate) enum ProteinSource {
    /// No amino-acid index; classify stops after the ANI tier.
    None,
    /// Precomputed proteomes (e.g. GTDB `protein_faa_reps`), named by accession.
    List(String),
    /// Call genes for every reference genome with pyrodigal.
    CallGenes(PathBuf),
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub(crate) struct Manifest {
    pub(crate) aeon_version: String,
    pub(crate) gtdb_release: String,
    pub(crate) gtdbtk_data: Option<String>,
    pub(crate) references: usize,
    /// Whether `aa.*` holds an amino-acid index aligned with `nt.*` vector IDs.
    #[serde(default)]
    pub(crate) proteins: bool,
}

/// An opened classify database.
pub(crate) struct ClassifyDb {
    pub(crate) dir: PathBuf,
    pub(crate) manifest: Manifest,
    pub(crate) data: GtdbtkData,
    pub(crate) nt_params: PrefixParams,
    /// Reference genome paths in nucleotide index vector-ID order.
    pub(crate) nt_paths: Vec<String>,
    /// GTDB accessions (`GCF_000005845.2`) aligned with `nt_paths`.
    pub(crate) nt_accessions: Vec<String>,
    /// Amino-acid sketch parameters and proteome paths, aligned with `nt_paths`.
    pub(crate) aa: Option<(PrefixParams, Vec<String>)>,
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
        )?;
        // A GTDB-Tk package's tables cover every genome in its reference
        // trees, which tier 3 needs; the database's own copies cover only
        // the indexed references.
        let data = match &manifest.gtdbtk_data {
            Some(root) => GtdbtkData::load(root)?,
            None => GtdbtkData::from_files(
                &manifest.gtdb_release,
                &[dir.join(TAXONOMY_FILE)],
                &dir.join(RADII_FILE),
            )?,
        };
        let nt_prefix = nt_prefix(&dir);
        let nt_paths = read_list_file(&genomes_path(&nt_prefix))?;
        let nt_accessions = nt_paths
            .iter()
            .map(|path| accession_from_path(path))
            .collect();
        let aa = if manifest.proteins {
            let aa_prefix = dir.join(AA_PREFIX).to_string_lossy().into_owned();
            let paths = read_list_file(&genomes_path(&aa_prefix))?;
            if paths.len() != nt_paths.len() {
                return Err("amino-acid and nucleotide indexes have different sizes".into());
            }
            Some((load_params(&aa_prefix)?, paths))
        } else {
            None
        };
        Ok(Self {
            nt_params: load_params(&nt_prefix)?,
            aa,
            dir,
            manifest,
            data,
            nt_paths,
            nt_accessions,
        })
    }

    pub(crate) fn nt_prefix(&self) -> String {
        nt_prefix(&self.dir)
    }

    pub(crate) fn aa_prefix(&self) -> String {
        self.dir.join(AA_PREFIX).to_string_lossy().into_owned()
    }
}

fn nt_prefix(dir: &Path) -> String {
    dir.join(NT_PREFIX).to_string_lossy().into_owned()
}

pub(crate) fn build_classify_db(
    db_dir: &str,
    genome_list: &str,
    source: TaxonomySource,
    proteins: ProteinSource,
    nt_params: PrefixParams,
    aa_params: PrefixParams,
) -> Result<(), Box<dyn Error>> {
    let started = Instant::now();
    let genomes = absolute_paths(read_list_file(genome_list)?)?;
    if genomes.is_empty() {
        return Err("reference genome list is empty".into());
    }
    let (data, gtdbtk_data) = match source {
        TaxonomySource::GtdbtkPackage(root) => {
            let root = std::path::absolute(root)?.to_string_lossy().into_owned();
            (GtdbtkData::load(&root)?, Some(root))
        }
        TaxonomySource::ReleaseFiles {
            release,
            taxonomy,
            sp_clusters,
        } => {
            let mut ranks = HashMap::new();
            for path in &taxonomy {
                ranks.extend(read_taxonomy(Path::new(path))?);
            }
            let data = GtdbtkData {
                version: release,
                taxonomy: ranks,
                radii: read_sp_clusters_radii(Path::new(&sp_clusters))?,
                ..GtdbtkData::default()
            };
            (data, None)
        }
    };
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
    write_reference_tables(db, &accessions, &data)?;
    build_database(&nt_prefix(db), &genomes, nt_params)?;
    let proteomes = match proteins {
        ProteinSource::None => None,
        ProteinSource::List(list) => Some(match_proteomes(
            &accessions,
            &absolute_paths(read_list_file(&list)?)?,
        )?),
        ProteinSource::CallGenes(pyrodigal) => {
            Some(call_reference_genes(&pyrodigal, &genomes, &accessions, db)?)
        }
    };
    if let Some(proteomes) = &proteomes {
        build_database(&db.join(AA_PREFIX).to_string_lossy(), proteomes, aa_params)?;
    }
    let manifest = Manifest {
        aeon_version: env!("CARGO_PKG_VERSION").into(),
        gtdb_release: data.version.clone(),
        gtdbtk_data,
        references: genomes.len(),
        proteins: proteomes.is_some(),
    };
    fs::write(
        db.join(MANIFEST_FILE),
        serde_json::to_string_pretty(&manifest)?,
    )?;
    eprintln!(
        "Built classify database: dir={}, references={}, elapsed={:.3}s",
        db_dir,
        genomes.len(),
        started.elapsed().as_secs_f64()
    );
    Ok(())
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
    pyrodigal: &Path,
    genomes: &[String],
    accessions: &[String],
    db: &Path,
) -> Result<Vec<String>, Box<dyn Error>> {
    let started = Instant::now();
    let protein_dir = db.join(PROTEIN_DIR);
    let done = AtomicUsize::new(0);
    let proteomes = genomes
        .par_iter()
        .zip(accessions)
        .map(|(genome, accession)| -> Result<String, String> {
            let calls = call_genes(pyrodigal, genome, accession, &protein_dir)
                .map_err(|error| format!("gene calling failed for {genome}: {error}"))?
                .ok_or_else(|| format!("reference genome has no sequence: {genome}"))?;
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

fn write_reference_tables(
    db: &Path,
    accessions: &[String],
    data: &GtdbtkData,
) -> Result<(), Box<dyn Error>> {
    let mut taxonomy = BufWriter::new(File::create(db.join(TAXONOMY_FILE))?);
    let mut radii = BufWriter::new(File::create(db.join(RADII_FILE))?);
    for accession in accessions {
        let ranks = data.taxonomy_of(accession).expect("checked above");
        let radius = &data.radii[&canonical_gid(accession)];
        writeln!(taxonomy, "{}\t{}", accession, ranks.join(";"))?;
        writeln!(radii, "{}\t{}\t{}", radius.species, accession, radius.ani)?;
    }
    taxonomy.flush()?;
    radii.flush()?;
    Ok(())
}
