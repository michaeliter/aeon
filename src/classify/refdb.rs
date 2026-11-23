//! The `aeon classify-db` reference database: two aeon indexes over the
//! reference genomes (e.g. GTDB species representatives or WoL2; genome
//! nucleotides and proteomes), plus a `reference_data/` directory with the
//! references' taxonomy, species radii, masked marker alignment and taxon
//! sibling radii, and the marker HMMs (and masks, when built from a GTDB-Tk
//! package). Classify needs nothing else.

use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, HashMap, HashSet};
use std::error::Error;
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

use needletail::parse_fastx_file;

use crate::classify::external::find_tool;
use crate::classify::genes::GeneCaller;
use crate::classify::marker_rank::{
    encode_row, read_marker_radii, resolve_radii, taxon_radii, write_marker_alignment,
    write_marker_radii,
};
use crate::classify::markers::{
    MARKER_SETS, MarkerTools, PfamInfo, apply_mask, extract_masked_markers, read_mask,
    reference_msa_path, write_marker_data,
};
use crate::classify::refdata::{
    DEFAULT_RADIUS, GtdbtkData, Radius, accession_from_path, canonical_gid,
    marker_set_for_taxonomy, strip_db_prefix,
};
use crate::database::build_database;
use crate::metadata::{PrefixParams, genomes_path, load_params};
use crate::utils::{read_list_file, write_genome_list};
use rayon::prelude::*;
use std::sync::atomic::{AtomicUsize, Ordering};

const MANIFEST_FILE: &str = "manifest.json";
const NT_PREFIX: &str = "nt";
const AA_PREFIX: &str = "aa";
const PROTEIN_DIR: &str = "proteins";
const MARKER_TMP_DIR: &str = "marker_genes_tmp";
const REFERENCE_DIR: &str = "reference_data";
const TAXONOMY_FILE: &str = "taxonomy.tsv";
const RADII_FILE: &str = "radii.tsv";

/// Where reference proteomes for the AAI index and marker alignment come from.
pub(crate) enum ProteinSource {
    /// Precomputed proteomes (e.g. GTDB `protein_faa_reps`), named by genome ID.
    List(String),
    /// Call genes for every reference genome with FragGeneScanRs.
    CallGenes(GeneCaller),
}

/// Where reference taxonomy, radii and marker data come from.
pub(crate) struct ReferenceSources {
    /// A GTDB-Tk reference data package: taxonomy, radii, marker HMMs,
    /// masks and the reference MSA.
    pub(crate) gtdbtk_data: Option<String>,
    /// Taxonomy files, used instead of the package's.
    pub(crate) taxonomy: Vec<String>,
    /// A radii file, used instead of the package's.
    pub(crate) radii: Option<String>,
    /// Release label, by default the package's.
    pub(crate) release: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub(crate) struct Manifest {
    pub(crate) aeon_version: String,
    pub(crate) gtdb_release: String,
    /// The GTDB-Tk package the database was built from, for the record.
    pub(crate) gtdbtk_data: Option<String>,
    pub(crate) references: usize,
    /// References with rows in the marker alignment.
    pub(crate) marker_references: usize,
}

/// One aeon index of the database.
pub(crate) struct Index {
    pub(crate) prefix: String,
    pub(crate) params: PrefixParams,
    /// Indexed files in vector-ID order.
    pub(crate) paths: Vec<String>,
    /// Genome IDs (`GCF_000005845.2`, `G000005825`, ...) aligned with `paths`.
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
    /// Taxonomy, radii and marker data (`reference_data/`).
    reference_dir: PathBuf,
    /// Reference genomes.
    pub(crate) nt: Index,
    /// Reference proteomes, named by their genomes' accessions.
    pub(crate) aa: Index,
    /// The radius of each reference taxon from phylum to genus for the
    /// marker rank call.
    pub(crate) marker_radii: HashMap<String, f64>,
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
        let nt = Index::open(&dir, NT_PREFIX)?;
        let mut aa = Index::open(&dir, AA_PREFIX)?;
        // Proteome file names may write an accession differently from the
        // genome's (`RS_` prefix, GenBank or RefSeq form), so name each
        // proteome by its genome's accession.
        let genome_accessions = nt
            .accessions
            .iter()
            .map(|accession| (canonical_gid(accession), accession))
            .collect::<HashMap<_, _>>();
        for accession in &mut aa.accessions {
            if let Some(genome) = genome_accessions.get(&canonical_gid(accession)) {
                accession.clone_from(genome);
            }
        }
        let reference_dir = dir.join(REFERENCE_DIR);
        let file = |name: &str| reference_dir.join(name).to_string_lossy().into_owned();
        Ok(Self {
            marker_radii: read_marker_radii(&reference_dir)?,
            data: GtdbtkData::from_files(
                &manifest.gtdb_release,
                &[file(TAXONOMY_FILE)],
                Some(&file(RADII_FILE)),
            )?,
            nt,
            aa,
            reference_dir,
            manifest,
        })
    }

    /// The marker-data directory (GTDB-Tk package layout).
    pub(crate) fn root(&self) -> &Path {
        &self.reference_dir
    }
}

/// Index parameters of the two indexes.
pub(crate) struct IndexParams {
    pub(crate) nt: PrefixParams,
    pub(crate) aa: PrefixParams,
}

pub(crate) fn build_classify_db(
    db_dir: &str,
    genome_list: &str,
    sources: ReferenceSources,
    proteins: ProteinSource,
    bin_dir: Option<&str>,
    params: IndexParams,
) -> Result<(), Box<dyn Error>> {
    let started = Instant::now();
    let genomes = absolute_paths(read_list_file(genome_list)?)?;
    if genomes.is_empty() {
        return Err("reference genome list is empty".into());
    }
    let package = sources
        .gtdbtk_data
        .as_deref()
        .map(std::path::absolute)
        .transpose()?;
    let mut data = match &package {
        Some(root) => GtdbtkData::load(&root.to_string_lossy())?,
        None => GtdbtkData::default(),
    };
    let package_release = package.as_ref().map(|_| data.version.clone());
    let files = GtdbtkData::from_files("", &sources.taxonomy, sources.radii.as_deref())?;
    if !sources.taxonomy.is_empty() {
        data.taxonomy = files.taxonomy;
    }
    if sources.radii.is_some() {
        data.radii = files.radii;
    }
    data.version = sources
        .release
        .or_else(|| package_release.clone())
        .unwrap_or_else(|| "custom".into());
    eprintln!(
        "Loaded {} reference data: {} taxonomy records, {} species radii",
        data.version,
        data.taxonomy.len(),
        data.radii.len()
    );

    let accessions = genomes
        .iter()
        .map(|path| accession_from_path(path))
        .collect::<Vec<_>>();
    check_references(&accessions, &data)?;
    let data = reference_subset(&data, &accessions);

    let db = &std::path::absolute(db_dir)?;
    let reference_dir = db.join(REFERENCE_DIR);
    fs::create_dir_all(&reference_dir)?;
    write_reference_data(&reference_dir, &data, &accessions)?;
    let package = package.as_deref().zip(package_release.as_deref());
    write_marker_data(&reference_dir, package)?;
    eprintln!(
        "Wrote marker HMMs and alignment masks from {}",
        package.map_or("aeon's bundled GTDB R232 copy".into(), |(root, _)| root
            .display()
            .to_string())
    );
    // The nucleotide index needs only the genomes, so it builds while the
    // proteomes, protein index and marker alignment are prepared.
    let (nt, rest) = rayon::join(
        || build_index(&db.join(NT_PREFIX), &genomes, params.nt),
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
                || build_index(&db.join(AA_PREFIX), &proteomes, params.aa),
                || {
                    write_reference_markers(
                        &reference_dir,
                        package,
                        &data,
                        &accessions,
                        &proteomes,
                        bin_dir,
                        db,
                    )
                    .map_err(|error| error.to_string())
                },
            );
            aa?;
            markers
        },
    );
    nt?;
    let marker_references = rest?;
    // The nt index build drops exact-content duplicates, so count what was
    // indexed. The marker alignment keeps every reference with markers.
    let references = Index::open(db, NT_PREFIX)?.len();
    let manifest = Manifest {
        aeon_version: env!("CARGO_PKG_VERSION").into(),
        gtdb_release: data.version.clone(),
        gtdbtk_data: package.map(|(root, _)| root.to_string_lossy().into_owned()),
        references,
        marker_references,
    };
    fs::write(
        db.join(MANIFEST_FILE),
        serde_json::to_string_pretty(&manifest)?,
    )?;
    eprintln!(
        "Built classify database: dir={}, references={}, with markers={}, elapsed={:.3}s",
        db_dir,
        references,
        marker_references,
        started.elapsed().as_secs_f64()
    );
    Ok(())
}

/// Write the references' masked marker alignment rows to `root`, one file
/// per marker set, and the radii of their taxa from phylum to genus, and
/// return how many references have rows. With a GTDB-Tk package, rows come
/// from its reference MSA; other references get theirs from hmmsearch and
/// hmmalign on their proteomes, using the marker data in `root`.
fn write_reference_markers(
    root: &Path,
    package: Option<(&Path, &str)>,
    data: &GtdbtkData,
    accessions: &[String],
    proteomes: &[String],
    bin_dir: Option<&str>,
    db: &Path,
) -> Result<usize, Box<dyn Error>> {
    let started = Instant::now();
    let (mut with_markers, mut from_package) = (0, 0);
    let mut radii = BTreeMap::new();
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
            write_marker_alignment(root, set, &[])?;
            continue;
        }
        let mut rows = Vec::new();
        let mask = read_mask(root, set)?;
        // The R232 bac120 MSA is several GB, so rows are streamed.
        let msa = package.map(|(package, release)| reference_msa_path(package, release, set));
        if let Some(msa) = msa
            .as_ref()
            .filter(|msa| fs::metadata(msa).is_ok_and(|metadata| metadata.len() > 0))
        {
            let mut reader = parse_fastx_file(&msa)
                .map_err(|error| format!("cannot read {}: {error}", msa.display()))?;
            while let Some(record) = reader.next() {
                let record = record?;
                let id = String::from_utf8_lossy(record.id());
                let id = strip_db_prefix(id.split_whitespace().next().unwrap_or_default());
                let Some((accession, _)) = pending.remove(&canonical_gid(id)) else {
                    continue;
                };
                let (masked, _, _) =
                    apply_mask(&String::from_utf8_lossy(&record.seq()), mask.as_deref())?;
                push_marker_row(&mut rows, accession, &masked);
                if pending.is_empty() {
                    break;
                }
            }
            from_package += rows.len();
        }
        if !pending.is_empty() {
            eprintln!(
                "Extracting {set} markers with HMMER for {} reference(s){}",
                pending.len(),
                msa.map_or(String::new(), |msa| format!(
                    " missing from {}",
                    msa.display()
                ))
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
            let work_dir = db.join(MARKER_TMP_DIR);
            let masked: BTreeMap<String, String> = extract_masked_markers(
                &tools,
                mask.as_deref(),
                set,
                &inputs,
                &work_dir,
                rayon::current_num_threads(),
            )?;
            fs::remove_dir_all(&work_dir)?;
            for (accession, masked) in masked {
                push_marker_row(&mut rows, &accession, &masked);
            }
        }
        write_marker_alignment(root, set, &rows)?;
        with_markers += rows.len();
        let radii_started = Instant::now();
        let references = rows
            .iter()
            .filter_map(|(accession, row)| Some((data.taxonomy_of(accession)?, row.as_slice())))
            .collect::<Vec<_>>();
        let set_radii = taxon_radii(&references);
        eprintln!(
            "Computed {set} marker radii of {} taxa from {} references in {:.3}s",
            set_radii.len(),
            references.len(),
            radii_started.elapsed().as_secs_f64()
        );
        radii.extend(set_radii);
    }
    // Rank medians pool both marker sets.
    write_marker_radii(root, &radii, &resolve_radii(&radii))?;
    eprintln!(
        "Wrote marker alignment rows for {with_markers}/{} references ({from_package} from the package MSA) in {:.3}s",
        accessions.len(),
        started.elapsed().as_secs_f64()
    );
    Ok(with_markers)
}

/// Add a reference's masked alignment row, unless it has no residues.
fn push_marker_row(rows: &mut Vec<(String, Vec<u8>)>, accession: &str, masked: &str) {
    let row = encode_row(masked.as_bytes());
    if row.iter().any(|&residue| residue != 0) {
        rows.push((accession.to_owned(), row));
    }
}

/// Build the aeon index at `prefix` over `paths`. `build_database` reads its
/// inputs from a list file, so the paths are written to one beside the index
/// for the build.
fn build_index(prefix: &Path, paths: &[String], params: PrefixParams) -> Result<(), String> {
    let prefix = prefix.to_string_lossy();
    let list = format!("{prefix}.inputs.txt");
    write_genome_list(&list, paths).map_err(|error| error.to_string())?;
    let built = build_database(&prefix, &list, params).map_err(|error| error.to_string());
    fs::remove_file(&list).map_err(|error| error.to_string())?;
    built
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
            let protein_faa = caller
                .call(genome, accession, &protein_dir, threads, false)
                .map_err(|error| format!("gene calling failed for {genome}: {error}"))?
                .ok_or_else(|| format!("no genes were called for reference genome {genome}"))?;
            let finished = done.fetch_add(1, Ordering::Relaxed) + 1;
            if finished % 100 == 0 {
                eprintln!("Called genes for {finished}/{} references", genomes.len());
            }
            Ok(protein_faa.to_string_lossy().into_owned())
        })
        .collect::<Result<Vec<_>, _>>()?;
    eprintln!(
        "Called genes for {} references in {:.3}s",
        genomes.len(),
        started.elapsed().as_secs_f64()
    );
    Ok(proteomes)
}

/// Every reference needs a seven-rank taxonomy, and accessions must be unique.
fn check_references(accessions: &[String], data: &GtdbtkData) -> Result<(), Box<dyn Error>> {
    let mut seen = HashSet::with_capacity(accessions.len());
    let mut missing = Vec::new();
    for accession in accessions {
        if !seen.insert(canonical_gid(accession)) {
            return Err(format!("reference list contains accession {accession} twice").into());
        }
        if data.taxonomy_of(accession).is_none() {
            missing.push(accession.as_str());
        }
    }
    if missing.is_empty() {
        return Ok(());
    }
    Err(format!(
        "{} reference genome(s) have no taxonomy (e.g. {}). Reference files must be named by accession, like GCF_000005845.2_genomic.fna.gz, and listed in the taxonomy file",
        missing.len(),
        missing[..missing.len().min(5)].join(", ")
    )
    .into())
}

/// The references' taxonomy and radii; species without a radius get
/// GTDB-Tk's default of 95%.
fn reference_subset(data: &GtdbtkData, accessions: &[String]) -> GtdbtkData {
    let mut subset = GtdbtkData {
        version: data.version.clone(),
        ..GtdbtkData::default()
    };
    let mut defaulted = 0;
    for accession in accessions {
        let gid = canonical_gid(accession);
        let taxonomy = data.taxonomy[&gid].clone();
        let radius = data.radii.get(&gid).cloned().unwrap_or_else(|| {
            defaulted += 1;
            Radius {
                species: taxonomy[6].clone(),
                ani: DEFAULT_RADIUS,
            }
        });
        subset.radii.insert(gid.clone(), radius);
        subset.taxonomy.insert(gid, taxonomy);
    }
    if defaulted > 0 {
        eprintln!(
            "{defaulted} reference(s) have no species radius; using {DEFAULT_RADIUS}% for them"
        );
    }
    subset
}

/// `taxonomy.tsv` and `radii.tsv` in GTDB-Tk's package formats.
fn write_reference_data(
    dir: &Path,
    data: &GtdbtkData,
    accessions: &[String],
) -> Result<(), Box<dyn Error>> {
    let mut taxonomy = BufWriter::new(File::create(dir.join(TAXONOMY_FILE))?);
    let mut radii = BufWriter::new(File::create(dir.join(RADII_FILE))?);
    for accession in accessions {
        let gid = canonical_gid(accession);
        writeln!(taxonomy, "{accession}\t{}", data.taxonomy[&gid].join(";"))?;
        let radius = &data.radii[&gid];
        writeln!(radii, "{}\t{accession}\t{}", radius.species, radius.ani)?;
    }
    taxonomy.flush()?;
    radii.flush()?;
    Ok(())
}
