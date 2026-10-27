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

use crate::classify::refdata::{
    GtdbtkData, accession_from_path, canonical_gid, read_sp_clusters_radii, read_taxonomy,
};
use crate::database::build_database;
use crate::metadata::{PrefixParams, genomes_path, load_params};
use crate::utils::read_list_file;

const MANIFEST_FILE: &str = "manifest.json";
const TAXONOMY_FILE: &str = "taxonomy.tsv";
const RADII_FILE: &str = "radii.tsv";
const NT_PREFIX: &str = "nt";

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

#[derive(Debug, Clone, Serialize, Deserialize)]
pub(crate) struct Manifest {
    pub(crate) aeon_version: String,
    pub(crate) gtdb_release: String,
    pub(crate) gtdbtk_data: Option<String>,
    pub(crate) references: usize,
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
        let mut data = GtdbtkData::from_files(
            &manifest.gtdb_release,
            &[dir.join(TAXONOMY_FILE)],
            &dir.join(RADII_FILE),
        )?;
        if let Some(root) = &manifest.gtdbtk_data {
            let package = GtdbtkData::load(root)?;
            data.root = package.root;
            data.red_dist_bac = package.red_dist_bac;
            data.red_dist_arc = package.red_dist_arc;
        }
        let nt_prefix = nt_prefix(&dir);
        let nt_paths = read_list_file(&genomes_path(&nt_prefix))?;
        let nt_accessions = nt_paths
            .iter()
            .map(|path| accession_from_path(path))
            .collect();
        Ok(Self {
            nt_params: load_params(&nt_prefix)?,
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
}

fn nt_prefix(dir: &Path) -> String {
    dir.join(NT_PREFIX).to_string_lossy().into_owned()
}

pub(crate) fn build_classify_db(
    db_dir: &str,
    genome_list: &str,
    source: TaxonomySource,
    nt_params: PrefixParams,
) -> Result<(), Box<dyn Error>> {
    let started = Instant::now();
    let genomes = read_list_file(genome_list)?;
    if genomes.is_empty() {
        return Err("reference genome list is empty".into());
    }
    let (data, gtdbtk_data) = match source {
        TaxonomySource::GtdbtkPackage(root) => (GtdbtkData::load(&root)?, Some(root)),
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

    let db = Path::new(db_dir);
    fs::create_dir_all(db)?;
    write_reference_tables(db, &accessions, &data)?;
    build_database(&nt_prefix(db), genome_list, nt_params)?;
    let manifest = Manifest {
        aeon_version: env!("CARGO_PKG_VERSION").into(),
        gtdb_release: data.version.clone(),
        gtdbtk_data,
        references: genomes.len(),
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
