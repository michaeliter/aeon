mod classify;
mod database;
mod metadata;
mod sketch;
mod utils;

use std::error::Error;

use crate::classify::refdb::{ProteinSource, TaxonomySource, build_classify_db};
use crate::classify::{ClassifyOptions, QueryInput, classify_genomes};
use crate::database::{build_database, search_database, update_database};
use crate::metadata::PrefixParams;
use crate::utils::{init_rayon_global, read_list_file, requested_threads};

use clap::{Arg, ArgAction, Command};

fn prefix_arg() -> Arg {
    Arg::new("prefix")
        .long("prefix")
        .short('p')
        .required(true)
        .action(ArgAction::Set)
        .value_name("PREFIX")
        .help("Database prefix shared by the .diskann, .genomes.txt, .idmap.tsv, and .params.json files")
}

fn threads_arg() -> Arg {
    Arg::new("threads")
        .long("threads")
        .short('t')
        .value_name("N")
        .value_parser(clap::value_parser!(usize))
        .visible_alias("cpus")
        .help("Worker threads [default: all logical CPUs]")
}

fn pyrodigal_arg() -> Arg {
    Arg::new("pyrodigal")
        .long("pyrodigal")
        .default_value("pyrodigal")
        .value_name("PATH")
        .help("pyrodigal executable used for gene calling")
}

fn classify_index_params(args: &clap::ArgMatches, seq_type: &str) -> PrefixParams {
    let (kmer_size, sketch_size) = match seq_type {
        "aa" => ("aa_kmer_size", "aa_sketch_size"),
        _ => ("kmer_size", "sketch_size"),
    };
    PrefixParams {
        kmer_size: *args.get_one(kmer_size).unwrap(),
        sketch_size: *args.get_one(sketch_size).unwrap(),
        densification: 0,
        seq_type: seq_type.into(),
        max_degree: *args.get_one("max_degree").unwrap(),
        build_beam_width: *args.get_one("build_beam_width").unwrap(),
        alpha: 1.2,
        extra_seeds: 1,
        sketch_elem_size: 2,
        sketch_elem_type: "u16".into(),
        distance: "anndists::dist::DistHamming".into(),
        hash_seed: 1337,
    }
}

fn update_beam_arg() -> Arg {
    Arg::new("beam_width")
        .long("beam-width")
        .value_name("N")
        .value_parser(clap::value_parser!(usize))
        .help("Candidate beam for inserted genomes [default: original build beam]")
        .long_help("Candidate beam used to connect inserted genomes to the existing Vamana graph. Larger values generally improve graph quality but increase update time. When omitted, Aeon reuses the build beam stored in PREFIX.params.json.")
}

fn build_cli() -> Command {
    Command::new("aeon")
        .version(env!("CARGO_PKG_VERSION"))
        .about("Efficient and robust genome similarity search with evolving disk-based proximity graphs")
        .long_about("Build and search mmap-backed genome sketch indexes. Insertions and MERIT deletions run in a temporary dynamic workspace and are committed back to the ordinary static DiskANN format.")
        .subcommand_required(true)
        .arg_required_else_help(true)
        .subcommand(
            Command::new("todisk")
                .about("Build a new static diskANN index using the Vamana graph construction algorithm")
                .long_about("Sketch every sequence file in the reference list, construct a Vamana graph, and write the static mmap-searchable index plus its genome-name mapping and parameter metadata.")
                .arg(prefix_arg())
                .arg(
                    Arg::new("reference_list")
                        .long("reference-list")
                        .short('r')
                        .required(true)
                        .value_name("FILE")
                        .help("Text file containing one reference FASTA/FASTQ path per line")
                        .long_help("Text file containing one reference FASTA/FASTQ path per line. Needletail-supported compressed inputs are accepted. Paths are retained as the genome names used by search and deletion."),
                )
                .arg(
                    Arg::new("kmer_size")
                        .long("kmer-size")
                        .short('k')
                        .default_value("16")
                        .value_name("K")
                        .value_parser(clap::value_parser!(usize))
                        .help("k-mer length; DNA supports <=32 except 15, amino acid supports <=12"),
                )
                .arg(
                    Arg::new("sketch_size")
                        .long("sketch-size")
                        .short('s')
                        .default_value("8192")
                        .value_name("N")
                        .value_parser(clap::value_parser!(usize))
                        .help("Number of u16 values in each MinHash sketch"),
                )
                .arg(
                    Arg::new("densification")
                        .long("densification")
                        .short('d')
                        .default_value("0")
                        .value_name("MODE")
                        .value_parser(["0", "1"])
                        .help("Densification algorithm: 0=OPTDENS, 1=RevOPTDENS"),
                )
                .arg(
                    Arg::new("seq_type")
                        .long("seq-type")
                        .default_value("dna")
                        .value_name("TYPE")
                        .value_parser(["dna", "aa"])
                        .help("Input alphabet: dna or aa (amino acid)"),
                )
                .arg(
                    Arg::new("max_degree")
                        .long("max-degree")
                        .default_value("64")
                        .value_name("R")
                        .value_parser(clap::value_parser!(usize))
                        .help("Maximum Vamana out-degree; larger values use more space and usually improve recall"),
                )
                .arg(
                    Arg::new("build_beam_width")
                        .long("build-beam-width")
                        .default_value("256")
                        .value_name("L")
                        .value_parser(clap::value_parser!(usize))
                        .help("Construction beam; larger values improve graph quality at higher build cost"),
                )
                .arg(
                    Arg::new("alpha")
                        .long("alpha")
                        .default_value("1.2")
                        .value_name("FLOAT")
                        .value_parser(clap::value_parser!(f32))
                        .help("Vamana robust-pruning alpha parameter"),
                )
                .arg(
                    Arg::new("extra_seeds")
                        .long("extra-seeds")
                        .default_value("1")
                        .value_name("N")
                        .value_parser(clap::value_parser!(usize))
                        .help("Additional random construction-search entry points per node"),
                )
                .arg(
                    Arg::new("hash_seed")
                        .long("hash-seed")
                        .default_value("1337")
                        .value_name("INTEGER")
                        .value_parser(clap::value_parser!(u64))
                        .help("XXH3 seed stored in metadata and reused for every future sketch"),
                )
                .arg(threads_arg()),
        )
        .subcommand(
            Command::new("insert")
                .about("Insert sketch vectors into an existing diskANN index/Vamana graph after sketching new genomes")
                .long_about("Sketch new genome files with the original database parameters, connect them to a temporary dynamic Vamana graph, and commit a compact static index.")
                .arg(prefix_arg())
                .arg(
                    Arg::new("genome_list")
                        .long("genome-list")
                        .short('i')
                        .required(true)
                        .value_name("FILE")
                        .help("Text file containing new genome FASTA/FASTQ file paths to sketch and then insert, one per line"),
                )
                .arg(update_beam_arg())
                .arg(threads_arg()),
        )
        .subcommand(
            Command::new("delete")
                .about("Delete sketch vectors from Vamana graph with routability-aware MERIT graph repair")
                .long_about("The temporary DiskANN session starts with its default routability guard enabled. Before deleting stored genome paths or unique genome filenames, Aeon virtually masks the requested batch and verifies that a small set of dispersed graph landmarks remains reachable. An admitted batch is repaired with MERIT, compacted, and committed as an ordinary static index. A deferred batch does not alter the database. Sequence files are not read during deletion.")
                .arg(prefix_arg())
                .arg(
                    Arg::new("name_list")
                        .long("name-list")
                        .short('d')
                        .required(true)
                        .value_name("FILE")
                        .help("Text file containing stored paths or unique genome filenames from PREFIX.genomes.txt")
                        .long_help("Text file containing one genome name per line. A complete stored path always matches. A filename without directories is accepted when it uniquely identifies one stored genome; use the complete path when the same filename occurs more than once."),
                )
                .arg(threads_arg()),
        )
        .subcommand(
            Command::new("update")
                .about("Routability-guarded Vamana graph deletion and MERIT graph repair followed by insertion in one transaction")
                .long_about("The temporary DiskANN session starts with its default routability guard enabled. Virtually validate the requested deletion batch against dispersed graph landmarks, apply MERIT repair only when it is admitted, then insert new sketches in one temporary update session and perform one compact static commit. A deferred deletion leaves the database unchanged, so no insertion is applied. A deleted path may be reinserted in the same command.")
                .arg(prefix_arg())
                .arg(
                    Arg::new("delete_list")
                        .long("delete-list")
                        .short('d')
                        .required(true)
                        .value_name("FILE")
                        .help("Text file containing stored paths or unique genome filenames to delete")
                        .long_help("Text file containing one genome name per line. A complete stored path always matches. A filename without directories is accepted when it uniquely identifies one stored genome; use the complete path when the same filename occurs more than once."),
                )
                .arg(
                    Arg::new("insert_list")
                        .long("insert-list")
                        .short('i')
                        .required(true)
                        .value_name("FILE")
                        .help("Text file containing genome FASTA/FASTQ paths to sketch and then insert"),
                )
                .arg(update_beam_arg())
                .arg(threads_arg()),
        )
        .subcommand(
            Command::new("search")
                .about("Search query genomes against a static diskANN index")
                .long_about("Sketch query genomes using parameters recovered from database metadata, search the mmap-backed static DiskANN graph in parallel, and report ranked neighbors.")
                .arg(prefix_arg())
                .arg(
                    Arg::new("query_list")
                        .long("query-list")
                        .short('q')
                        .required(true)
                        .value_name("FILE")
                        .help("Text file containing one query FASTA/FASTQ path per line"),
                )
                .arg(
                    Arg::new("k")
                        .long("k")
                        .default_value("10")
                        .value_name("N")
                        .value_parser(clap::value_parser!(usize))
                        .help("Number of nearest neighbors reported per query"),
                )
                .arg(
                    Arg::new("beam_width")
                        .long("beam-width")
                        .default_value("512")
                        .value_name("N")
                        .value_parser(clap::value_parser!(usize))
                        .help("Search beam; larger values usually improve recall but reduce throughput"),
                )
                .arg(
                    Arg::new("output")
                        .long("output")
                        .short('o')
                        .value_name("FILE")
                        .help("Write TSV results to FILE instead of standard output"),
                )
                .arg(threads_arg()),
        )
        .subcommand(
            Command::new("classify-db")
                .about("Build a GTDB reference database for aeon classify")
                .long_about("Sketch GTDB species-representative genomes into a static DiskANN index and store the GTDB taxonomy and species ANI radii needed by aeon classify. Taxonomy comes from a GTDB-Tk reference data package, or from GTDB release taxonomy and sp_clusters files. Reference files must be named by accession (e.g. GCF_000005845.2_genomic.fna.gz) and must stay at their listed paths, because classify computes exact ANI against them.")
                .arg(
                    Arg::new("db")
                        .long("db")
                        .required(true)
                        .value_name("DIR")
                        .help("Output database directory"),
                )
                .arg(
                    Arg::new("genome_list")
                        .long("genome-list")
                        .short('r')
                        .required(true)
                        .value_name("FILE")
                        .help("Text file with one GTDB representative genome FASTA path per line"),
                )
                .arg(
                    Arg::new("gtdbtk_data")
                        .long("gtdbtk-data")
                        .value_name("DIR")
                        .env("GTDBTK_DATA_PATH")
                        .conflicts_with_all(["taxonomy", "sp_clusters"])
                        .help("GTDB-Tk reference data package; required later for marker-based placement"),
                )
                .arg(
                    Arg::new("taxonomy")
                        .long("taxonomy")
                        .value_name("FILE")
                        .action(ArgAction::Append)
                        .requires_all(["sp_clusters", "release"])
                        .help("GTDB release taxonomy file(s), e.g. bac120_taxonomy_r232.tsv.gz (repeatable)"),
                )
                .arg(
                    Arg::new("sp_clusters")
                        .long("sp-clusters")
                        .value_name("FILE")
                        .requires("taxonomy")
                        .help("GTDB release sp_clusters_r*.tsv providing species ANI radii"),
                )
                .arg(
                    Arg::new("release")
                        .long("release")
                        .value_name("NAME")
                        .help("GTDB release name recorded with --taxonomy, e.g. r232"),
                )
                .arg(
                    Arg::new("protein_list")
                        .long("protein-list")
                        .value_name("FILE")
                        .conflicts_with("call_genes")
                        .help("Text file with one reference proteome path per line, named by accession (e.g. GTDB protein_faa_reps)"),
                )
                .arg(
                    Arg::new("call_genes")
                        .long("call-genes")
                        .action(ArgAction::SetTrue)
                        .help("Call reference genes with pyrodigal to build the AAI protein index"),
                )
                .arg(pyrodigal_arg())
                .arg(
                    Arg::new("kmer_size")
                        .long("kmer-size")
                        .short('k')
                        .default_value("16")
                        .value_name("K")
                        .value_parser(clap::value_parser!(usize))
                        .help("Nucleotide k-mer length"),
                )
                .arg(
                    Arg::new("aa_kmer_size")
                        .long("aa-kmer-size")
                        .default_value("7")
                        .value_name("K")
                        .value_parser(clap::value_parser!(usize))
                        .help("Amino-acid k-mer length for the protein index"),
                )
                .arg(
                    Arg::new("aa_sketch_size")
                        .long("aa-sketch-size")
                        .default_value("12000")
                        .value_name("N")
                        .value_parser(clap::value_parser!(usize))
                        .help("Number of u16 values in each protein MinHash sketch"),
                )
                .arg(
                    Arg::new("sketch_size")
                        .long("sketch-size")
                        .short('s')
                        .default_value("12000")
                        .value_name("N")
                        .value_parser(clap::value_parser!(usize))
                        .help("Number of u16 values in each MinHash sketch"),
                )
                .arg(
                    Arg::new("max_degree")
                        .long("max-degree")
                        .default_value("64")
                        .value_name("R")
                        .value_parser(clap::value_parser!(usize))
                        .help("Maximum Vamana out-degree"),
                )
                .arg(
                    Arg::new("build_beam_width")
                        .long("build-beam-width")
                        .default_value("256")
                        .value_name("L")
                        .value_parser(clap::value_parser!(usize))
                        .help("Construction beam"),
                )
                .arg(threads_arg()),
        )
        .subcommand(
            Command::new("classify")
                .about("Classify genomes against GTDB with GTDB-Tk-compatible output")
                .long_about("Classify genomes in a cascade. DiskANN sketch search proposes the closest GTDB representatives and TurboANI computes their exact ANI. A genome within the closest reference's ANI radius (AF >= --min_af) receives that species, as in the GTDB-Tk ANI screen; one at or above --genus-ani to a strict majority of one genus receives that genus. Remaining genomes have genes called with pyrodigal; protein sketch search proposes candidate reference proteomes and DIAMOND reciprocal best hits give AAI, assigning a genus (--genus-aai) or family (--family-aai) by majority among the --aai-top closest references. Genomes still unresolved go through GTDB-Tk's marker pipeline: hmmsearch finds the bac120/ar53 markers, the copy-number profile picks the domain, hmmalign builds the masked concatenated alignment, and pplacer places it in the GTDB reference tree (bacteria in the backbone tree, then a class-level tree, unless --full_tree), where RED and topology decide the rank. This tier requires a classify-db built from a GTDB-Tk data package. Results are written in the GTDB-Tk classify_wf layout: {out_dir}/{prefix}.{bac120,ar53}.summary.tsv, classify/, and {prefix}.log.")
                .arg(
                    Arg::new("db")
                        .long("db")
                        .required(true)
                        .value_name("DIR")
                        .help("Database directory built by aeon classify-db"),
                )
                .arg(
                    Arg::new("genome_dir")
                        .long("genome_dir")
                        .visible_alias("genome-dir")
                        .value_name("DIR")
                        .required_unless_present("batchfile")
                        .conflicts_with("batchfile")
                        .help("Directory containing genome FASTA files"),
                )
                .arg(
                    Arg::new("extension")
                        .long("extension")
                        .short('x')
                        .default_value("fna")
                        .value_name("EXT")
                        .help("Extension of genome files in --genome_dir (e.g. fna, fna.gz)"),
                )
                .arg(
                    Arg::new("batchfile")
                        .long("batchfile")
                        .value_name("FILE")
                        .help("Tab-separated file of genome paths with optional genome IDs"),
                )
                .arg(
                    Arg::new("out_dir")
                        .long("out_dir")
                        .visible_alias("out-dir")
                        .required(true)
                        .value_name("DIR")
                        .help("Output directory"),
                )
                .arg(
                    Arg::new("prefix")
                        .long("prefix")
                        .default_value("gtdbtk")
                        .value_name("PREFIX")
                        .help("Prefix of all output files"),
                )
                .arg(
                    Arg::new("min_af")
                        .long("min_af")
                        .visible_alias("min-af")
                        .default_value("0.5")
                        .value_name("FLOAT")
                        .value_parser(clap::value_parser!(f64))
                        .help("Minimum alignment fraction for a species assignment"),
                )
                .arg(
                    Arg::new("genus_ani")
                        .long("genus-ani")
                        .default_value("80")
                        .value_name("PERCENT")
                        .value_parser(clap::value_parser!(f64))
                        .help("ANI needed for genus consensus; more distant genomes are left for the AAI tier"),
                )
                .arg(
                    Arg::new("genus_aai")
                        .long("genus-aai")
                        .default_value("75")
                        .value_name("PERCENT")
                        .value_parser(clap::value_parser!(f64))
                        .help("AAI needed for genus consensus"),
                )
                .arg(
                    Arg::new("family_aai")
                        .long("family-aai")
                        .default_value("65")
                        .value_name("PERCENT")
                        .value_parser(clap::value_parser!(f64))
                        .help("AAI needed for family consensus"),
                )
                .arg(
                    Arg::new("aai_top")
                        .long("aai-top")
                        .default_value("10")
                        .value_name("N")
                        .value_parser(clap::value_parser!(usize))
                        .help("Closest references by AAI that vote on genus and family"),
                )
                .arg(pyrodigal_arg())
                .arg(
                    Arg::new("diamond")
                        .long("diamond")
                        .default_value("diamond")
                        .value_name("PATH")
                        .help("DIAMOND executable used for AAI"),
                )
                .arg(
                    Arg::new("min_perc_aa")
                        .long("min_perc_aa")
                        .visible_alias("min-perc-aa")
                        .default_value("10")
                        .value_name("PERCENT")
                        .value_parser(clap::value_parser!(f64))
                        .help("Leave out genomes with amino acids in fewer than this percentage of masked MSA columns"),
                )
                .arg(
                    Arg::new("full_tree")
                        .short('f')
                        .long("full_tree")
                        .visible_alias("full-tree")
                        .action(ArgAction::SetTrue)
                        .help("Place bacterial genomes in the full bac120 tree instead of the backbone and class-level trees"),
                )
                .arg(
                    Arg::new("bin_dir")
                        .long("bin-dir")
                        .value_name("DIR")
                        .help("Directory holding hmmsearch, hmmalign, pplacer and guppy (default: search PATH)"),
                )
                .arg(
                    Arg::new("candidates")
                        .long("candidates")
                        .default_value("10")
                        .value_name("N")
                        .value_parser(clap::value_parser!(usize))
                        .help("Sketch-search candidates re-scored with exact ANI per genome"),
                )
                .arg(
                    Arg::new("beam_width")
                        .long("beam-width")
                        .default_value("512")
                        .value_name("N")
                        .value_parser(clap::value_parser!(usize))
                        .help("DiskANN search beam"),
                )
                .arg(threads_arg()),
        )
}

fn main() -> Result<(), Box<dyn Error>> {
    println!("\n ************** initializing logger *****************\n");
    let _ = env_logger::Builder::from_default_env().try_init();
    log::info!("Logger initialized from default environment");
    let matches = build_cli().get_matches();
    let (command, args) = matches.subcommand().expect("subcommand is required");
    init_rayon_global(requested_threads(args));
    eprintln!("Running aeon {command}");

    match command {
        "todisk" => build_database(
            args.get_one::<String>("prefix").unwrap(),
            args.get_one::<String>("reference_list").unwrap(),
            PrefixParams {
                kmer_size: *args.get_one("kmer_size").unwrap(),
                sketch_size: *args.get_one("sketch_size").unwrap(),
                densification: args
                    .get_one::<String>("densification")
                    .unwrap()
                    .parse()
                    .unwrap(),
                seq_type: args.get_one::<String>("seq_type").unwrap().clone(),
                max_degree: *args.get_one("max_degree").unwrap(),
                build_beam_width: *args.get_one("build_beam_width").unwrap(),
                alpha: *args.get_one("alpha").unwrap(),
                extra_seeds: *args.get_one("extra_seeds").unwrap(),
                sketch_elem_size: 2,
                sketch_elem_type: "u16".into(),
                distance: "anndists::dist::DistHamming".into(),
                hash_seed: *args.get_one("hash_seed").unwrap(),
            },
        )?,
        "insert" => update_database(
            args.get_one::<String>("prefix").unwrap(),
            read_list_file(args.get_one::<String>("genome_list").unwrap())?,
            Vec::new(),
            args.get_one::<usize>("beam_width").copied(),
        )?,
        "delete" => update_database(
            args.get_one::<String>("prefix").unwrap(),
            Vec::new(),
            read_list_file(args.get_one::<String>("name_list").unwrap())?,
            None,
        )?,
        "update" => update_database(
            args.get_one::<String>("prefix").unwrap(),
            read_list_file(args.get_one::<String>("insert_list").unwrap())?,
            read_list_file(args.get_one::<String>("delete_list").unwrap())?,
            args.get_one::<usize>("beam_width").copied(),
        )?,
        "search" => search_database(
            args.get_one::<String>("prefix").unwrap(),
            args.get_one::<String>("query_list").unwrap(),
            *args.get_one::<usize>("k").unwrap(),
            *args.get_one::<usize>("beam_width").unwrap(),
            args.get_one::<String>("output").map(String::as_str),
        )?,
        "classify-db" => build_classify_db(
            args.get_one::<String>("db").unwrap(),
            args.get_one::<String>("genome_list").unwrap(),
            match args.get_one::<String>("gtdbtk_data") {
                Some(root) if !args.contains_id("taxonomy") => {
                    TaxonomySource::GtdbtkPackage(root.clone())
                }
                _ => TaxonomySource::ReleaseFiles {
                    taxonomy: args
                        .get_many::<String>("taxonomy")
                        .ok_or("provide --gtdbtk-data, or --taxonomy with --sp-clusters")?
                        .cloned()
                        .collect(),
                    release: args.get_one::<String>("release").unwrap().clone(),
                    sp_clusters: args.get_one::<String>("sp_clusters").unwrap().clone(),
                },
            },
            if args.get_flag("call_genes") {
                ProteinSource::CallGenes(classify::external::find_program(
                    args.get_one::<String>("pyrodigal").unwrap(),
                )?)
            } else if let Some(list) = args.get_one::<String>("protein_list") {
                ProteinSource::List(list.clone())
            } else {
                ProteinSource::None
            },
            classify_index_params(args, "dna"),
            classify_index_params(args, "aa"),
        )?,
        "classify" => classify_genomes(&ClassifyOptions {
            db: args.get_one::<String>("db").unwrap().clone(),
            input: match args.get_one::<String>("batchfile") {
                Some(batchfile) => QueryInput::Batchfile(batchfile.clone()),
                None => QueryInput::GenomeDir {
                    dir: args.get_one::<String>("genome_dir").unwrap().clone(),
                    extension: args.get_one::<String>("extension").unwrap().clone(),
                },
            },
            out_dir: args.get_one::<String>("out_dir").unwrap().clone(),
            prefix: args.get_one::<String>("prefix").unwrap().clone(),
            min_af: *args.get_one("min_af").unwrap(),
            candidates: *args.get_one("candidates").unwrap(),
            beam_width: *args.get_one("beam_width").unwrap(),
            genus_ani: *args.get_one("genus_ani").unwrap(),
            genus_aai: *args.get_one("genus_aai").unwrap(),
            family_aai: *args.get_one("family_aai").unwrap(),
            aai_top: *args.get_one("aai_top").unwrap(),
            pyrodigal: args.get_one::<String>("pyrodigal").unwrap().clone(),
            diamond: args.get_one::<String>("diamond").unwrap().clone(),
            min_perc_aa: *args.get_one("min_perc_aa").unwrap(),
            bin_dir: args.get_one::<String>("bin_dir").cloned(),
            full_tree: args.get_flag("full_tree"),
            threads: requested_threads(args),
        })?,
        _ => unreachable!(),
    }
    Ok(())
}
