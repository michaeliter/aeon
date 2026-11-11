//! GTDB-Tk's `{prefix}.json` run record: the command, the database and, per
//! step, its timing, options and output files.
//!
//! GTDB-Tk reads this file back to resume a run and rejects step names and
//! fields it does not know, so aeon's tiers are filed under GTDB-Tk's four
//! `classify_wf` steps: gene calling and marker search under `identify`,
//! the AAI tier and the marker-index search under `classify`.

use chrono::{DateTime, Local};
use serde::Serialize;
use std::collections::BTreeMap;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};

pub(crate) const ANI_SCREEN: &str = "ANI screen";
pub(crate) const IDENTIFY: &str = "identify";
pub(crate) const ALIGN: &str = "align";
pub(crate) const CLASSIFY: &str = "classify";
const STEPS: [&str; 4] = [ANI_SCREEN, IDENTIFY, ALIGN, CLASSIFY];

/// Output files keyed by marker set, or `all`.
pub(crate) type OutputFiles = BTreeMap<String, Vec<String>>;

#[derive(Default)]
struct StepLog {
    spans: Vec<(DateTime<Local>, DateTime<Local>)>,
    output_files: OutputFiles,
}

/// Collects the wall-clock spans and output files behind each step.
pub(crate) struct RunRecord {
    last_mark: DateTime<Local>,
    steps: BTreeMap<&'static str, StepLog>,
}

impl RunRecord {
    pub(crate) fn start() -> Self {
        Self {
            last_mark: Local::now(),
            steps: BTreeMap::new(),
        }
    }

    /// Attribute the time since the previous mark to `step`.
    pub(crate) fn mark(&mut self, step: &'static str) {
        let now = Local::now();
        self.steps
            .entry(step)
            .or_default()
            .spans
            .push((self.last_mark, now));
        self.last_mark = now;
    }

    /// Record `path` as an output of `step` for `marker_set` (or `all`).
    pub(crate) fn output(&mut self, step: &'static str, marker_set: &str, path: &Path) {
        self.steps
            .entry(step)
            .or_default()
            .output_files
            .entry(marker_set.to_owned())
            .or_default()
            .push(absolute(path));
    }

    /// Write the record to `{out_dir}/{prefix}.json`. A step spans its
    /// first start to its last end; its duration is the time spent in it.
    pub(crate) fn write(self, run: &RunInfo) -> io::Result<()> {
        let path = run.out_dir.join(format!("{}.json", run.prefix));
        let output_dir = absolute(&run.out_dir);
        let mut steps = Vec::new();
        for name in STEPS {
            let Some(log) = self.steps.get(name) else {
                continue;
            };
            let (Some(first), Some(last)) = (log.spans.first(), log.spans.last()) else {
                continue;
            };
            let seconds = log
                .spans
                .iter()
                .map(|(start, end)| (*end - *start).num_microseconds().unwrap_or(0))
                .sum::<i64>()
                / 1_000_000;
            let mut output_files = log.output_files.clone();
            let details = match name {
                ANI_SCREEN => StepDetails::AniScreen {
                    genome_dir: run.genome_dir.clone(),
                    batchfile: run.batchfile.clone(),
                    min_af: run.min_af,
                    // GTDB-Tk maps each marker set to its one ANI summary.
                    unassigned_hits_file: output_files
                        .remove("all")
                        .and_then(|mut files| files.pop()),
                    output_files: output_files
                        .into_iter()
                        .filter_map(|(set, mut files)| Some((set, files.pop()?)))
                        .collect(),
                },
                IDENTIFY => StepDetails::Identify {
                    genes: false,
                    extension: run.extension.clone(),
                    write_single_copy_genes: false,
                    genome_dir: run.genome_dir.clone(),
                    batchfile: run.batchfile.clone(),
                    output_files,
                },
                ALIGN => StepDetails::Align {
                    identify_dir: output_dir.clone(),
                    // The MSA holds the user genomes only.
                    skip_gtdb_refs: true,
                    taxa_filter: None,
                    min_perc_aa: run.min_perc_aa,
                    custom_msa_filters: false,
                    skip_trimming: false,
                    rnd_seed: None,
                    cols_per_gene: None,
                    min_consensus: None,
                    max_consensus: None,
                    min_perc_taxa: None,
                    outgroup_taxon: None,
                    output_files,
                },
                _ => StepDetails::Classify {
                    align_dir: output_dir.clone(),
                    genome_dir: run.genome_dir.clone(),
                    batchfile: run.batchfile.clone(),
                    scratch_dir: None,
                    debug_option: false,
                    full_tree: false,
                    place_species: false,
                    output_files,
                },
            };
            steps.push(Step {
                name,
                output_dir: output_dir.clone(),
                starts_at: isoformat(first.0),
                ends_at: isoformat(last.1),
                duration: py_timedelta(seconds),
                status: "completed",
                details,
            });
        }
        let record = Record {
            version: env!("CARGO_PKG_VERSION"),
            command_line: &run.command_line,
            database_version: &run.database_version,
            database_path: absolute(&run.database_path),
            steps,
            output_dir,
            path: absolute(&path),
        };
        let mut output = BufWriter::new(File::create(&path)?);
        // Python's `json.dump(..., indent=4)`.
        let formatter = serde_json::ser::PrettyFormatter::with_indent(b"    ");
        let mut serializer = serde_json::Serializer::with_formatter(&mut output, formatter);
        record
            .serialize(&mut serializer)
            .map_err(io::Error::other)?;
        output.flush()
    }
}

/// The run-wide fields of the record.
pub(crate) struct RunInfo {
    pub(crate) command_line: String,
    pub(crate) database_version: String,
    pub(crate) database_path: PathBuf,
    pub(crate) out_dir: PathBuf,
    pub(crate) prefix: String,
    pub(crate) genome_dir: Option<String>,
    pub(crate) batchfile: Option<String>,
    pub(crate) extension: Option<String>,
    pub(crate) min_af: f64,
    pub(crate) min_perc_aa: f64,
}

#[derive(Serialize)]
struct Record<'a> {
    version: &'static str,
    command_line: &'a str,
    database_version: &'a str,
    database_path: String,
    steps: Vec<Step>,
    output_dir: String,
    path: String,
}

#[derive(Serialize)]
struct Step {
    name: &'static str,
    output_dir: String,
    starts_at: String,
    ends_at: String,
    duration: String,
    status: &'static str,
    #[serde(flatten)]
    details: StepDetails,
}

/// The fields of GTDB-Tk's step dataclasses, in their order.
#[derive(Serialize)]
#[serde(untagged)]
enum StepDetails {
    AniScreen {
        genome_dir: Option<String>,
        batchfile: Option<String>,
        min_af: f64,
        output_files: BTreeMap<String, String>,
        unassigned_hits_file: Option<String>,
    },
    Identify {
        genes: bool,
        extension: Option<String>,
        write_single_copy_genes: bool,
        genome_dir: Option<String>,
        batchfile: Option<String>,
        output_files: OutputFiles,
    },
    Align {
        identify_dir: String,
        skip_gtdb_refs: bool,
        taxa_filter: Option<String>,
        min_perc_aa: f64,
        custom_msa_filters: bool,
        skip_trimming: bool,
        rnd_seed: Option<u64>,
        cols_per_gene: Option<u64>,
        min_consensus: Option<f64>,
        max_consensus: Option<f64>,
        min_perc_taxa: Option<f64>,
        outgroup_taxon: Option<String>,
        output_files: OutputFiles,
    },
    Classify {
        align_dir: String,
        genome_dir: Option<String>,
        batchfile: Option<String>,
        scratch_dir: Option<String>,
        debug_option: bool,
        full_tree: bool,
        place_species: bool,
        output_files: OutputFiles,
    },
}

fn absolute(path: &Path) -> String {
    std::path::absolute(path)
        .unwrap_or_else(|_| path.to_path_buf())
        .to_string_lossy()
        .into_owned()
}

/// Python's `datetime.isoformat()` for a local time.
fn isoformat(time: DateTime<Local>) -> String {
    time.naive_local()
        .format("%Y-%m-%dT%H:%M:%S%.6f")
        .to_string()
}

/// Python's `str(timedelta(seconds=...))`: `H:MM:SS`, with days before.
fn py_timedelta(seconds: i64) -> String {
    let (days, rest) = (seconds / 86_400, seconds % 86_400);
    let clock = format!("{}:{:02}:{:02}", rest / 3600, rest % 3600 / 60, rest % 60);
    match days {
        0 => clock,
        1 => format!("1 day, {clock}"),
        days => format!("{days} days, {clock}"),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn durations_format_like_python_timedelta() {
        assert_eq!(py_timedelta(9), "0:00:09");
        assert_eq!(py_timedelta(3 * 3600 + 62), "3:01:02");
        assert_eq!(py_timedelta(86_400 + 5), "1 day, 0:00:05");
        assert_eq!(py_timedelta(2 * 86_400), "2 days, 0:00:00");
    }

    #[test]
    fn steps_follow_gtdbtk_order_and_fields() {
        let dir = std::env::temp_dir().join(format!("aeon-record-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let mut record = RunRecord::start();
        record.mark(ANI_SCREEN);
        record.output(ANI_SCREEN, "all", &dir.join("hits.tsv.gz"));
        record.mark(IDENTIFY);
        record.mark(CLASSIFY);
        record.output(CLASSIFY, "bac120", &dir.join("summary.tsv"));
        record
            .write(&RunInfo {
                command_line: "aeon classify".into(),
                database_version: "r232".into(),
                database_path: dir.clone(),
                out_dir: dir.clone(),
                prefix: "gtdbtk".into(),
                genome_dir: Some("queries".into()),
                batchfile: None,
                extension: Some("fna.gz".into()),
                min_af: 0.5,
                min_perc_aa: 10.0,
            })
            .unwrap();
        let text = std::fs::read_to_string(dir.join("gtdbtk.json")).unwrap();
        std::fs::remove_dir_all(&dir).unwrap();
        let json: serde_json::Value = serde_json::from_str(&text).unwrap();
        let names = json["steps"]
            .as_array()
            .unwrap()
            .iter()
            .map(|step| step["name"].as_str().unwrap())
            .collect::<Vec<_>>();
        assert_eq!(names, [ANI_SCREEN, IDENTIFY, CLASSIFY]);
        assert_eq!(json["steps"][0]["output_files"], serde_json::json!({}));
        assert!(
            json["steps"][0]["unassigned_hits_file"]
                .as_str()
                .unwrap()
                .ends_with("hits.tsv.gz")
        );
        assert!(text.contains("\n    \"version\": "));
    }
}
