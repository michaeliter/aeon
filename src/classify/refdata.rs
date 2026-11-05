//! Readers for the GTDB-Tk reference data package (GTDBTK_DATA_PATH layout).

use flate2::read::MultiGzDecoder;
use std::collections::HashMap;
use std::error::Error;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

pub(crate) const RANK_PREFIXES: [&str; 7] = ["d__", "p__", "c__", "o__", "f__", "g__", "s__"];

/// Canonical form of an NCBI accession, mirroring GTDB-Tk's `canonical_gid`:
/// `RS_GCF_005435135.1` -> `G005435135`. UBA identifiers are returned unchanged.
pub(crate) fn canonical_gid(gid: &str) -> String {
    if gid.starts_with('U') {
        return gid.to_owned();
    }
    let gid = gid
        .replace("RS_", "")
        .replace("GB_", "")
        .replace("GCA_", "G")
        .replace("GCF_", "G");
    match gid.find('.') {
        Some(dot) => gid[..dot].to_owned(),
        None => gid,
    }
}

/// Strip the GTDB `RS_`/`GB_` database prefix, keeping the versioned accession.
pub(crate) fn strip_db_prefix(gid: &str) -> &str {
    gid.strip_prefix("RS_")
        .or_else(|| gid.strip_prefix("GB_"))
        .unwrap_or(gid)
}

/// Derive a reference accession from a GTDB file path such as
/// `.../GCF_000005845.2_genomic.fna.gz` or `.../GCA_000006155.2_protein.faa`.
pub(crate) fn accession_from_path(path: &str) -> String {
    let name = Path::new(path)
        .file_name()
        .and_then(|name| name.to_str())
        .unwrap_or(path);
    let mut stem = name;
    for suffix in [
        ".gz", ".fna", ".faa", ".fa", ".fasta", "_genomic", "_protein",
    ] {
        stem = stem.strip_suffix(suffix).unwrap_or(stem);
    }
    let stem = strip_db_prefix(stem);
    // GCF_000005845.2_ASM584v2 -> GCF_000005845.2
    if stem.starts_with("GCA_") || stem.starts_with("GCF_") {
        if let Some(dot) = stem.find('.') {
            let version_end = stem[dot + 1..]
                .find(|c: char| !c.is_ascii_digit())
                .map_or(stem.len(), |offset| dot + 1 + offset);
            return stem[..version_end].to_owned();
        }
    }
    stem.to_owned()
}

/// Pad an incomplete taxonomy string to seven ranks, like GTDB-Tk's
/// `standardise_taxonomy`.
pub(crate) fn standardise_taxonomy(taxonomy: &str, marker_set: Option<&str>) -> String {
    let mut ranks = taxonomy
        .split(';')
        .filter(|rank| !rank.is_empty())
        .map(str::to_owned)
        .collect::<Vec<_>>();
    let domain = match marker_set {
        Some("bac120") => Some("d__Bacteria"),
        Some("ar53") => Some("d__Archaea"),
        _ => None,
    };
    if let Some(domain) = domain {
        if ranks.first().map(String::as_str) != Some(domain) {
            ranks.insert(0, domain.to_owned());
        }
    }
    ranks.extend(
        RANK_PREFIXES[ranks.len().min(RANK_PREFIXES.len())..]
            .iter()
            .map(|prefix| (*prefix).to_owned()),
    );
    ranks.join(";")
}

/// The marker set GTDB-Tk reports a reference-derived taxonomy under.
pub(crate) fn marker_set_for_taxonomy(taxonomy: &[String]) -> &'static str {
    match taxonomy.first().map(String::as_str) {
        Some("d__Archaea") => "ar53",
        _ => "bac120",
    }
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Radius {
    pub(crate) species: String,
    pub(crate) ani: f64,
}

/// Reference taxonomy and species radii from a GTDB-Tk package.
#[derive(Debug, Default)]
pub(crate) struct GtdbtkData {
    pub(crate) version: String,
    /// Seven-rank taxonomy keyed by canonical genome ID.
    pub(crate) taxonomy: HashMap<String, Vec<String>>,
    /// Species ANI radius keyed by canonical representative genome ID.
    pub(crate) radii: HashMap<String, Radius>,
}

impl GtdbtkData {
    pub(crate) fn load(root: &str) -> Result<Self, Box<dyn Error>> {
        let root = Path::new(root);
        Ok(Self {
            version: read_version(&root.join("metadata").join("metadata.txt"))?,
            taxonomy: read_taxonomy(&root.join("taxonomy").join("gtdb_taxonomy.tsv"))?,
            radii: read_radii(&root.join("radii").join("gtdb_radii.tsv"))?,
        })
    }

    pub(crate) fn taxonomy_of(&self, gid: &str) -> Option<&[String]> {
        self.taxonomy.get(&canonical_gid(gid)).map(Vec::as_slice)
    }

    pub(crate) fn radius_of(&self, gid: &str) -> Option<f64> {
        self.radii.get(&canonical_gid(gid)).map(|radius| radius.ani)
    }
}

/// Lines of a plain or gzip-compressed text file.
pub(crate) fn read_lines(
    path: &Path,
) -> Result<impl Iterator<Item = std::io::Result<String>>, Box<dyn Error>> {
    let file =
        File::open(path).map_err(|error| format!("cannot open {}: {error}", path.display()))?;
    let reader: Box<dyn Read> = if path.extension().is_some_and(|extension| extension == "gz") {
        Box::new(MultiGzDecoder::new(file))
    } else {
        Box::new(file)
    };
    Ok(BufReader::new(reader).lines())
}

pub(crate) fn read_taxonomy(path: &Path) -> Result<HashMap<String, Vec<String>>, Box<dyn Error>> {
    let mut taxonomy = HashMap::new();
    for line in read_lines(path)? {
        let line = line?;
        let line = line.trim_end();
        if line.is_empty() {
            continue;
        }
        let (gid, taxa) = line
            .split_once('\t')
            .ok_or_else(|| format!("malformed taxonomy line in {}: {line}", path.display()))?;
        let ranks = taxa.split(';').map(str::to_owned).collect::<Vec<_>>();
        if ranks.len() != RANK_PREFIXES.len() {
            return Err(format!("expected a 7-rank taxonomy for {gid}: {taxa}").into());
        }
        taxonomy.insert(canonical_gid(gid), ranks);
    }
    Ok(taxonomy)
}

pub(crate) fn read_radii(path: &Path) -> Result<HashMap<String, Radius>, Box<dyn Error>> {
    let mut radii = HashMap::new();
    for line in read_lines(path)? {
        let line = line?;
        let fields = line.trim().split('\t').collect::<Vec<_>>();
        if fields.len() == 1 && fields[0].is_empty() {
            continue;
        }
        let [species, genome, ani] = fields[..] else {
            return Err(format!("malformed radii line in {}: {line}", path.display()).into());
        };
        radii.insert(
            canonical_gid(genome),
            Radius {
                species: species.to_owned(),
                ani: ani.parse()?,
            },
        );
    }
    Ok(radii)
}

/// `VERSION_DATA` (e.g. `r232`) from the package's `metadata/metadata.txt`.
pub(crate) fn read_version(path: &Path) -> Result<String, Box<dyn Error>> {
    let text = fs::read_to_string(path)
        .map_err(|error| format!("cannot open GTDB-Tk file {}: {error}", path.display()))?;
    text.lines()
        .find_map(|line| line.trim().strip_prefix("VERSION_DATA="))
        .map(str::to_owned)
        .ok_or_else(|| format!("{} is missing VERSION_DATA", path.display()).into())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn canonical_gid_matches_gtdbtk_examples() {
        assert_eq!(canonical_gid("G005435135"), "G005435135");
        assert_eq!(canonical_gid("GCF_005435135.1"), "G005435135");
        assert_eq!(
            canonical_gid("GCF_005435135.1_ASM543513v1_genomic"),
            "G005435135"
        );
        assert_eq!(canonical_gid("RS_GCF_005435135.1"), "G005435135");
        assert_eq!(canonical_gid("GB_GCA_005435135.1"), "G005435135");
        assert_eq!(canonical_gid("UBA1234"), "UBA1234");
    }

    #[test]
    fn accessions_come_from_gtdb_file_names() {
        assert_eq!(
            accession_from_path("/db/GCF_000005845.2_genomic.fna.gz"),
            "GCF_000005845.2"
        );
        assert_eq!(
            accession_from_path("GCA_000006155.2_protein.faa"),
            "GCA_000006155.2"
        );
        assert_eq!(
            accession_from_path("RS_GCF_000005845.2_ASM584v2_genomic.fna"),
            "GCF_000005845.2"
        );
        assert_eq!(accession_from_path("/mags/bin.12.fa"), "bin.12");
    }

    #[test]
    fn taxonomy_is_padded_to_seven_ranks() {
        assert_eq!(
            standardise_taxonomy("d__Bacteria;p__Bacillota;c__Bacilli", None),
            "d__Bacteria;p__Bacillota;c__Bacilli;o__;f__;g__;s__"
        );
        assert_eq!(
            standardise_taxonomy("p__Halobacteriota", Some("ar53")),
            "d__Archaea;p__Halobacteriota;c__;o__;f__;g__;s__"
        );
        assert_eq!(
            standardise_taxonomy("", Some("bac120")),
            "d__Bacteria;p__;c__;o__;f__;g__;s__"
        );
    }

    #[test]
    fn reads_package_version() {
        let dir = std::env::temp_dir().join(format!("aeon-refdata-{}", std::process::id()));
        fs::create_dir_all(&dir).unwrap();
        let path = dir.join("metadata.txt");
        fs::write(
            &path,
            "VERSION_DATA=r232\nRED_DIST_BAC_DICT={\"p__\": 0.35}\n",
        )
        .unwrap();
        let version = read_version(&path).unwrap();
        fs::remove_dir_all(&dir).unwrap();
        assert_eq!(version, "r232");
    }
}
