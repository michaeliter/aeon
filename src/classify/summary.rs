//! Writers for GTDB-Tk-compatible classification output files.

use std::fs::{self, File};
use std::io::{self, BufWriter, Write};
use std::path::Path;

use crate::classify::refdata::GtdbtkData;

pub(crate) const NONE_VALUE: &str = "N/A";

pub(crate) const SUMMARY_COLUMNS: [&str; 20] = [
    "user_genome",
    "classification",
    "closest_genome_reference",
    "closest_genome_reference_radius",
    "closest_genome_taxonomy",
    "closest_genome_ani",
    "closest_genome_af",
    "closest_placement_reference",
    "closest_placement_radius",
    "closest_placement_taxonomy",
    "closest_placement_ani",
    "closest_placement_af",
    "pplacer_taxonomy",
    "classification_method",
    "note",
    "other_related_references(genome_id,species_name,radius,ANI,AF)",
    "msa_percent",
    "translation_table",
    "red_value",
    "warnings",
];

/// Python `repr(float)`: shortest round-trip digits, always with a decimal
/// point or exponent, switching to exponent notation outside [1e-4, 1e16).
pub(crate) fn py_float(value: f64) -> String {
    if value.is_nan() {
        return "nan".into();
    }
    if value.is_infinite() {
        return if value > 0.0 { "inf" } else { "-inf" }.into();
    }
    let magnitude = value.abs();
    if magnitude != 0.0 && !(1e-4..1e16).contains(&magnitude) {
        let formatted = format!("{value:e}");
        let (mantissa, exponent) = formatted.split_once('e').expect("exponent format");
        let (sign, digits) = match exponent.strip_prefix('-') {
            Some(digits) => ('-', digits),
            None => ('+', exponent),
        };
        return format!("{mantissa}e{sign}{digits:0>2}");
    }
    let formatted = value.to_string();
    if formatted.contains('.') {
        formatted
    } else {
        format!("{formatted}.0")
    }
}

/// Python `str(round(value, digits))`.
pub(crate) fn py_round(value: f64, digits: usize) -> String {
    let rounded: f64 = format!("{value:.digits$}")
        .parse()
        .expect("formatted float parses");
    py_float(rounded)
}

/// One `other_related_references` entry, as GTDB-Tk's `Classify.formatnote`.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct RelatedReference {
    pub(crate) genome_id: String,
    pub(crate) ani: f64,
    pub(crate) af: f64,
}

/// Format closest-first references, skipping `exclude`, capped at 50 entries.
pub(crate) fn format_related_references(
    references: &[RelatedReference],
    exclude: &[&str],
    data: &GtdbtkData,
) -> Option<String> {
    let entries = references
        .iter()
        .filter(|reference| !exclude.contains(&reference.genome_id.as_str()))
        .take(50)
        .map(|reference| {
            let species = data
                .taxonomy_of(&reference.genome_id)
                .map_or("None", |taxonomy| taxonomy[6].as_str());
            let radius = data
                .radius_of(&reference.genome_id)
                .map_or_else(|| "None".to_owned(), py_float);
            format!(
                "{}, {}, {}, {}, {}",
                reference.genome_id,
                species,
                radius,
                py_round(reference.ani, 2),
                py_round(reference.af, 3)
            )
        })
        .collect::<Vec<_>>();
    (!entries.is_empty()).then(|| entries.join("; "))
}

#[derive(Debug, Clone, Default, PartialEq)]
pub(crate) struct SummaryRow {
    pub(crate) user_genome: String,
    pub(crate) classification: Option<String>,
    pub(crate) closest_genome_reference: Option<String>,
    pub(crate) closest_genome_reference_radius: Option<f64>,
    pub(crate) closest_genome_taxonomy: Option<String>,
    pub(crate) closest_genome_ani: Option<f64>,
    pub(crate) closest_genome_af: Option<f64>,
    pub(crate) closest_placement_reference: Option<String>,
    pub(crate) closest_placement_radius: Option<f64>,
    pub(crate) closest_placement_taxonomy: Option<String>,
    pub(crate) closest_placement_ani: Option<f64>,
    pub(crate) closest_placement_af: Option<f64>,
    pub(crate) pplacer_taxonomy: Option<String>,
    pub(crate) classification_method: Option<String>,
    pub(crate) note: Option<String>,
    pub(crate) other_related_references: Option<String>,
    pub(crate) msa_percent: Option<f64>,
    pub(crate) translation_table: Option<u8>,
    pub(crate) red_value: Option<f64>,
    /// Joined with `;` in sorted order, so output is deterministic.
    pub(crate) warnings: Vec<String>,
}

impl SummaryRow {
    pub(crate) fn new(user_genome: impl Into<String>) -> Self {
        Self {
            user_genome: user_genome.into(),
            ..Self::default()
        }
    }

    fn fields(&self) -> [String; 20] {
        let text = |value: &Option<String>| value.clone().unwrap_or_else(|| NONE_VALUE.into());
        let float = |value: Option<f64>| value.map_or_else(|| NONE_VALUE.into(), py_float);
        let rounded = |value: Option<f64>, digits| {
            value.map_or_else(|| NONE_VALUE.into(), |value| py_round(value, digits))
        };
        let mut warnings = self.warnings.clone();
        warnings.sort();
        warnings.dedup();
        [
            self.user_genome.clone(),
            text(&self.classification),
            text(&self.closest_genome_reference),
            float(self.closest_genome_reference_radius),
            text(&self.closest_genome_taxonomy),
            rounded(self.closest_genome_ani, 2),
            rounded(self.closest_genome_af, 3),
            text(&self.closest_placement_reference),
            float(self.closest_placement_radius),
            text(&self.closest_placement_taxonomy),
            rounded(self.closest_placement_ani, 2),
            rounded(self.closest_placement_af, 3),
            text(&self.pplacer_taxonomy),
            text(&self.classification_method),
            text(&self.note),
            text(&self.other_related_references),
            rounded(self.msa_percent, 2),
            self.translation_table
                .map_or_else(|| NONE_VALUE.into(), |table| table.to_string()),
            rounded(self.red_value, 5),
            if warnings.is_empty() {
                NONE_VALUE.into()
            } else {
                warnings.join(";")
            },
        ]
    }
}

/// Write `{prefix}.{marker_set}.summary.tsv` rows sorted by `user_genome`.
pub(crate) fn write_summary(path: &Path, rows: &[SummaryRow]) -> io::Result<()> {
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent)?;
    }
    let mut output = BufWriter::new(File::create(path)?);
    write_summary_to(&mut output, rows)?;
    output.flush()
}

pub(crate) fn write_summary_to(output: &mut impl Write, rows: &[SummaryRow]) -> io::Result<()> {
    writeln!(output, "{}", SUMMARY_COLUMNS.join("\t"))?;
    let mut sorted = rows.iter().collect::<Vec<_>>();
    sorted.sort_by(|left, right| left.user_genome.cmp(&right.user_genome));
    for row in sorted {
        writeln!(output, "{}", row.fields().join("\t"))?;
    }
    Ok(())
}

/// One genome assigned by the ANI screen, for `classify/ani_screen/{p}.{ms}.ani_summary.tsv`.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct AniScreenRow {
    pub(crate) user_genome: String,
    pub(crate) reference_genome: String,
    pub(crate) ani: f64,
    pub(crate) af: f64,
    pub(crate) reference_taxonomy: String,
    pub(crate) other_related_references: Option<String>,
}

pub(crate) fn write_ani_summary(path: &Path, rows: &[AniScreenRow]) -> io::Result<()> {
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent)?;
    }
    let mut output = BufWriter::new(File::create(path)?);
    writeln!(
        output,
        "user_genome\treference_genome\tskani_ani\tskani_af\treference_taxonomy\tother_related_references(genome_id,species_name,radius,ANI,AF)"
    )?;
    let mut sorted = rows.iter().collect::<Vec<_>>();
    sorted.sort_by(|left, right| left.user_genome.cmp(&right.user_genome));
    for row in sorted {
        writeln!(
            output,
            "{}\t{}\t{}\t{}\t{}\t{}",
            row.user_genome,
            row.reference_genome,
            py_float(row.ani),
            py_round(row.af, 3),
            row.reference_taxonomy,
            row.other_related_references
                .as_deref()
                .unwrap_or(NONE_VALUE)
        )?;
    }
    output.flush()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::classify::refdata::Radius;

    #[test]
    fn floats_match_python_repr() {
        assert_eq!(py_float(95.0), "95.0");
        assert_eq!(py_float(0.94), "0.94");
        assert_eq!(py_float(100.0), "100.0");
        assert_eq!(py_float(0.0), "0.0");
        assert_eq!(py_float(0.0005), "0.0005");
        assert_eq!(py_float(0.00005), "5e-05");
        assert_eq!(py_float(1e-7), "1e-07");
        assert_eq!(py_float(0.1 + 0.2), "0.30000000000000004");
    }

    #[test]
    fn rounding_matches_python_round() {
        assert_eq!(py_round(0.125, 2), "0.12");
        assert_eq!(py_round(0.375, 2), "0.38");
        assert_eq!(py_round(2.675, 2), "2.67");
        assert_eq!(py_round(97.123456, 2), "97.12");
        assert_eq!(py_round(0.95, 3), "0.95");
        assert_eq!(py_round(0.912345678, 5), "0.91235");
    }

    #[test]
    fn summary_rows_are_sorted_and_use_na() {
        let mut ani = SummaryRow::new("query_b");
        ani.classification = Some("d__Bacteria;p__P;c__C;o__O;f__F;g__G;s__G s".into());
        ani.closest_genome_reference = Some("GCF_000005845.2".into());
        ani.closest_genome_reference_radius = Some(95.0);
        ani.closest_genome_ani = Some(98.7654);
        ani.closest_genome_af = Some(0.91234);
        ani.classification_method = Some("ani_screen".into());
        ani.note = Some("classification based on ANI only".into());
        ani.translation_table = Some(11);
        ani.warnings = vec!["b warning".into(), "a warning".into()];
        let unclassified = SummaryRow {
            classification: Some("Unclassified Bacteria".into()),
            ..SummaryRow::new("query_a")
        };

        let mut buffer = Vec::new();
        write_summary_to(&mut buffer, &[ani, unclassified]).unwrap();
        let text = String::from_utf8(buffer).unwrap();
        let lines = text.lines().collect::<Vec<_>>();
        assert_eq!(lines[0].split('\t').count(), 20);
        assert!(lines[1].starts_with("query_a\tUnclassified Bacteria\tN/A\t"));
        let fields = lines[2].split('\t').collect::<Vec<_>>();
        assert_eq!(fields[2], "GCF_000005845.2");
        assert_eq!(fields[3], "95.0");
        assert_eq!(fields[5], "98.77");
        assert_eq!(fields[6], "0.912");
        assert_eq!(fields[17], "11");
        assert_eq!(fields[18], "N/A");
        assert_eq!(fields[19], "a warning;b warning");
    }

    #[test]
    fn reproduces_gtdbtk_documented_summary_example() {
        let taxonomy = |genus_species: &str| {
            format!(
                "d__Archaea;p__Thermoplasmatota;c__Thermoplasmata;o__Methanomassiliicoccales;f__Methanomethylophilaceae;{genus_species}"
            )
        };
        let ruminantium = "d__Archaea;p__Methanobacteriota;c__Methanobacteria;o__Methanobacteriales;f__Methanobacteriaceae;g__Methanobrevibacter;s__Methanobrevibacter ruminantium";
        let nitrosopumilaceae = "d__Archaea;p__Thermoproteota;c__Nitrososphaeria;o__Nitrososphaerales;f__Nitrosopumilaceae;g__;s__";
        let rows = [
            SummaryRow {
                classification: Some(ruminantium.into()),
                closest_genome_reference: Some("GCF_000024185.1".into()),
                closest_genome_reference_radius: Some(95.0),
                closest_genome_taxonomy: Some(ruminantium.into()),
                closest_genome_ani: Some(100.0),
                closest_genome_af: Some(1.0),
                classification_method: Some("ani_screen".into()),
                note: Some("classification based on ANI only".into()),
                other_related_references: Some("GCA_900321995.1, s__Methanobrevibacter sp900321995, 95.0, 80.9, 0.7; GCF_900114585.1, s__Methanobrevibacter olleyae, 95.0, 79.96, 0.55".into()),
                ..SummaryRow::new("genome_1")
            },
            SummaryRow {
                classification: Some(taxonomy("g__VadinCA11;s__VadinCA11 sp002498365")),
                closest_genome_reference: Some("GCA_002498365.1".into()),
                closest_genome_reference_radius: Some(95.0),
                closest_genome_taxonomy: Some(taxonomy("g__VadinCA11;s__VadinCA11 sp002498365")),
                closest_genome_ani: Some(99.16),
                closest_genome_af: Some(0.94),
                closest_placement_reference: Some("GCA_002498365.1".into()),
                closest_placement_radius: Some(95.0),
                closest_placement_taxonomy: Some(taxonomy("g__VadinCA11;s__VadinCA11 sp002498365")),
                closest_placement_ani: Some(99.16),
                closest_placement_af: Some(0.94),
                pplacer_taxonomy: Some(taxonomy("g__VadinCA11;s__")),
                classification_method: Some("taxonomic classification defined by topology and ANI".into()),
                note: Some("topological placement and ANI have congruent species assignments".into()),
                other_related_references: Some("GCA_002505345.1, s__VadinCA11 sp002505345, 95.0, 89.92, 0.89; GCA_002509405.1, s__VadinCA11 sp002509405, 95.0, 88.13, 0.89".into()),
                msa_percent: Some(87.1),
                translation_table: Some(11),
                ..SummaryRow::new("genome_2")
            },
            SummaryRow {
                classification: Some(nitrosopumilaceae.into()),
                pplacer_taxonomy: Some(nitrosopumilaceae.into()),
                classification_method: Some("taxonomic novelty determined using RED".into()),
                msa_percent: Some(85.47),
                translation_table: Some(11),
                red_value: Some(0.895231),
                ..SummaryRow::new("genome_3")
            },
            SummaryRow {
                classification: Some("Unclassified".into()),
                warnings: vec!["No bacterial or archaeal marker".into()],
                ..SummaryRow::new("genome_4")
            },
        ];
        let mut buffer = Vec::new();
        write_summary_to(&mut buffer, &rows).unwrap();
        let expected = include_str!("../../tests/data/gtdbtk_2.8.0_example.summary.tsv");
        assert_eq!(String::from_utf8(buffer).unwrap(), expected);
    }

    #[test]
    fn related_references_follow_formatnote() {
        let mut data = GtdbtkData::default();
        data.taxonomy.insert(
            "G000005845".into(),
            "d__Bacteria;p__P;c__C;o__O;f__F;g__Escherichia;s__Escherichia coli"
                .split(';')
                .map(str::to_owned)
                .collect(),
        );
        data.radii.insert(
            "G000005845".into(),
            Radius {
                species: "s__Escherichia coli".into(),
                ani: 95.0,
            },
        );
        let references = vec![
            RelatedReference {
                genome_id: "GCF_000001.1".into(),
                ani: 99.0,
                af: 0.9,
            },
            RelatedReference {
                genome_id: "GCF_000005845.2".into(),
                ani: 96.123,
                af: 0.8756,
            },
        ];
        assert_eq!(
            format_related_references(&references, &["GCF_000001.1"], &data).unwrap(),
            "GCF_000005845.2, s__Escherichia coli, 95.0, 96.12, 0.876"
        );
        assert_eq!(
            format_related_references(&references[..1], &["GCF_000001.1"], &data),
            None
        );
    }
}
