<div align="center">
  <img width="25%" src ="aeon-logo.png">
</div>


# [Aeon](https://dictionary.cambridge.org/us/dictionary/english/aeon) (/ˈiːɒn/) : efficient and robust genome similarity search with evolving disk-based proximity graphs 🦀

***"Built once, updated forever. Nodes and edges evolve like cells, while Aeon remains eternal."***


Aeon is a dynamic genome similarity search platform built on mmap-backed DiskANN index, MERIT deletion repair, and b-bit one-permutation MinHash sketches. Search indexes remain in the ordinary static DiskANN format: updates use a temporary dynamic workspace and commit a compact static index.

## Install

Aeon requires nightly Rust because the `stdsimd` features of its dependencies use `portable_simd`. The included `rust-toolchain.toml` selects nightly automatically (install it with `rustup toolchain install nightly`).

```bash
git clone https://github.com/jianshu93/aeon
cd aeon
cargo build --release
./target/release/aeon --help
```

Input lists contain one genome path or stored genome name per line. FASTA/FASTQ and compressed inputs supported by Needletail are accepted. DNA and amino-acid databases are supported.

Run `aeon COMMAND --help` for the full description of any command and its options. `--threads` defaults to all logical CPUs for every command.

## Logging

Aeon writes a small number of plain progress messages directly to stderr, including how many genomes were sketched, inserted, deleted, or searched. These messages do not use the logging framework and never mix with TSV search results written to stdout.

During DNA or amino-acid sketching, Aeon prints progress after every 1,000 completed genomes and once at completion.

```bash
RUST_LOG=debug aeon delete --prefix streptomyces --name-list delete.txt
RUST_LOG=rust_diskann=debug aeon delete --prefix streptomyces --name-list delete.txt
```

Optional `RUST_LOG` settings add debugging details; they do not control the plain progress messages.

## Commands

### Build

```bash
aeon todisk \
  --prefix index_out \
  --reference-list genomes.txt \
  --kmer-size 16 \
  --sketch-size 12000 \
  --max-degree 64 \
  --build-beam-width 256
```

This creates:

- `index_out.diskann`: static mmap-searchable DiskANN index.
- `index_out.genomes.txt`: genome names in current vector-ID order.
- `index_out.idmap.tsv`: vector ID to genome name mapping.
- `index_out.params.json`: sketch and graph parameters.

Important build options:

| Option | Default | Meaning |
|---|---:|---|
| `--seq-type <dna\|aa>` | `dna` | Select the DNA or amino-acid sketcher. |
| `--kmer-size <K>` | `16` | k-mer length. DNA supports up to 32 except 15; amino acid supports up to 12. |
| `--sketch-size <N>` | `8192` | Number of `u16` values per b-bit MinHash sketch. Larger sketches improve resolution but increase index size and distance cost. |
| `--densification <0\|1>` | `0` | `0` selects OPTDENS; `1` selects RevOPTDENS. |
| `--max-degree <R>` | `64` | Maximum Vamana out-degree. Larger values usually improve recall while increasing storage and traversal work. |
| `--build-beam-width <L>` | `256` | Construction beam. Larger values generally improve graph quality at higher build cost. |
| `--alpha <FLOAT>` | `1.2` | Vamana robust-pruning alpha. |
| `--extra-seeds <N>` | `1` | Additional random construction-search entry points per node. |
| `--hash-seed <INTEGER>` | `1337` | XXH3 seed saved in metadata and reused by insert and search. |

### Insert

The insertion list must contain paths to genome files. Aeon sketches them with the database's original parameters.

```bash
aeon insert --prefix index_out --genome-list insert.txt
```

`--beam-width` controls candidate exploration while connecting inserted genomes. It defaults to the original build beam stored in `PREFIX.params.json`; increasing it may improve graph quality at the cost of update time.

Aeon uses rust-diskann's parallel 256-node Vamana micro-batch insertion path. Outgoing graph searches and grouped reverse-edge RobustPrune work run through the configured Rayon pool; only the final deterministic mmap row writes are serialized. This is especially important for large b-bit MinHash sketches, where individual Hamming comparisons are memory-bandwidth intensive.

### Delete

Deletion does not read or sketch genomes. Each line may be either a complete stored path from `PREFIX.genomes.txt` or a unique genome filename such as `GCA_965228215.1_genomic.fna.gz`.

```bash
aeon delete --prefix index_out --name-list delete.txt
```

Aeon's `rust-diskann` update session enables routability-aware deletion admission control by default. Aeon reads that built-in guard rather than initializing a second one. It virtually masks the proposed batch, checks whether the same residual entry point can still reach a small dispersed landmark set, and applies MERIT versioned-edge invalidation and local repair only when no previously reachable landmark route is lost. A deferred batch leaves the existing static database unchanged and reports the affected landmark IDs; split that deletion list or schedule a rebuild. This is a lightweight navigability safeguard, not a formal Recall@k guarantee.

For an admitted batch, rust-diskann groups MERIT's ordered directed repair operations by mutable Vamana adjacency row. Operations for one row remain sequential, while independent rows run through Aeon's Rayon pool; immutable Hamming distances are cached within each row queue. This preserves the exact per-edge RobustPrune decisions while avoiding conflict waves that can collapse delete parallelism on closely clustered genomes.

After an admitted deletion, Aeon removes dynamic version state while committing the updated ordinary static index.

Complete stored paths always match. A basename without directories is accepted only when it identifies exactly one stored genome. If the database contains the same filename under multiple paths, Aeon rejects the ambiguous basename and asks for the complete stored path.

### Combined Update

Delete and insert in one update session and one static commit:

```bash
aeon update \
  --prefix index_out \
  --delete-list delete.txt \
  --insert-list insert.txt
```

A path may be deleted and reinserted in the same command. Names that remain in the database cannot be inserted again.

`--delete-list` follows the same complete-path-or-unique-basename rule as `--name-list`.

`--beam-width` has the same meaning as for `insert`. Before either `delete` or `update` mutates graph data, Aeon uses rust-diskann's default routability admission check. When it defers a batch, neither the deletion nor the pending insertion is committed. Admitted deletions use the rust-diskann MERIT defaults (`repair beam = 2R`, `k_r = 2`).

When an update contains insertions, it uses the same parallel Vamana micro-batch insertion path as `insert` after an admitted deletion repair.

### Search

```bash
aeon search \
  --prefix index_out \
  --query-list queries.txt \
  --k 10 \
  --beam-width 512 \
  --output neighbors.tsv
```

Search always opens the static index. Results contain query path, matched genome name, current vector ID, raw sketch Jaccard estimate, normalized Hamming distance (or hash collision probability), and `ANI` as the final column. For DNA sketches, Aeon calculates `ANI = 1 + ln(2J / (1 + J)) / k`, where `J` is the sketch collision-derived Jaccard estimate and `k` is the stored nucleotide k-mer length. `ANI` is written as a fraction from 0 to 1; a value above `0.95` is strong evidence for an already represented species-level genome. For incomplete MAGs, a value below `0.95` does not exclude a same-species match because reduced alignment fraction also lowers Jaccard. Use the value to prioritize candidates and confirm species calls with an alignment-based ANI and alignment fraction. Amino-acid searches write `NA` because the nucleotide Mash model does not apply.

Search options:

| Option | Default | Meaning |
|---|---:|---|
| `--k <N>` | `10` | Nearest neighbors reported for each query. |
| `--beam-width <N>` | `512` | Graph-search beam. Larger values generally improve recall but reduce throughput. |
| `--output <FILE>` | stdout | Write tab-separated results to a file. |

### Classify

`aeon classify` assigns GTDB taxonomy and writes GTDB-Tk 2.8 `classify_wf` outputs: `{prefix}.bac120.summary.tsv`, `{prefix}.ar53.summary.tsv`, and the `identify/`, `align/` and `classify/` directories with GTDB-Tk's file names and columns. Results should be close to GTDB-Tk's, but the method differs. aeon has no pplacer, RED or reference-tree step. Below species it classifies by AAI and by a search of marker genes.

**Three indexes.** `aeon classify-db` builds three aeon indexes over the GTDB species representatives, and every classify run uses all three:

| Index | Sequences | k-mers | Sketch size |
|---|---|---|---|
| `nt` | reference genomes | nucleotide, k=16 (`-k`) | 15,000 (`-s`) |
| `aa` | reference proteomes | amino acid, k=7 (`--aa-kmer-size`) | 10,000 (`--aa-sketch-size`) |
| `markers` | concatenated bac120/ar53 marker genes (masked GTDB-Tk alignment, gaps removed) | amino acid, k=5 (`--marker-kmer-size`) | 3,000 (`--marker-sketch-size`) |

All three graphs are built with max degree 128 (`--max-degree`), construction beam 512 (`--build-beam-width`) and alpha 1.2 (`--alpha`).

Marker genes for references come from the GTDB-Tk package's reference MSA. References missing from it get theirs from hmmsearch and hmmalign on their proteomes.

**Cascade.** Each genome stops at the first step that resolves it:

1. **nt search, then TurboANI when candidates span species.** The genome sketch search proposes `--candidates` (10) species representatives with a sketch-estimated ANI.
   - **One candidate at or above `--min-sketch-ani` (90%):** if its estimate clears that species' radius by at least 1 point, the genome gets the species without TurboANI.
   - **Several candidates at or above 90%, or one without that margin:** TurboANI computes their exact ANI and AF, and GTDB-Tk's ANI screen rule decides. The closest reference with AF ≥ `--min_af` (0.5) gives its species if the ANI is within that species' GTDB radius (95% for most species).
   - **No candidate at or above 90%:** the genome goes to step 2.
2. **FragGeneScanRs, aa search, DIAMOND AAI.**
   - Genes are called with FragGeneScanRs's complete-genome model (`-t complete -w 1`).
   - The proteome sketch search proposes candidate reference proteomes.
   - DIAMOND reciprocal best hits give AAI. The implementation is inspired by the enveomics collection (Rodriguez-R LM & Konstantinidis KT (2016). The enveomics collection: a toolbox for specialized analyses of microbial genomes and metagenomes. *PeerJ Preprints* 4:e1900v1, [doi:10.7287/peerj.preprints.1900v1](https://doi.org/10.7287/peerj.preprints.1900v1)):
     - `diamond blastp` runs in both directions in DIAMOND's default (fast) mode, with its default e-value (0.001).
     - Each protein's best hit is its top hit with ≥20% identity. There is no alignment-length or coverage filter.
     - AAI is the unweighted mean identity of the reciprocal best hits. A reference needs at least 50 of them to get an AAI.
     - A query is searched against all its candidates in one database. Each hit's e-value is rescaled to its own reference proteome's size, and each protein keeps its best hit per reference, so the result matches separate searches of each genome pair.
   - The closest reference by AAI gives the call:
     - **Genus:** its genus is assigned if its AAI is ≥ `--genus-aai` (65%) and beats every reference of another genus by at least `--aai-margin` (5) points.
     - **Family:** failing that, its family is assigned under the same rule with ≥ `--family-aai` (45%).
     - **Otherwise** the genome goes to step 3.
3. **hmmer marker extraction, marker search.**
   - hmmsearch and hmmalign extract the bac120/ar53 markers exactly as GTDB-Tk's `identify` and `align` steps do, with the same copy-number rules, domain call, mask and `--min_perc_aa` filter.
   - The masked concatenated markers are searched against the marker index.
   - Only hits from the genome's own domain (bac120 or ar53) count, and only at or above `--min-marker-identity` (25% estimated identity).
   - The genome gets the closest such reference's lineage down to order.
   - A genome with no such hit is left unclassified below its domain.

**Design decisions:**
- **pplacer is replaced by AAI and the marker search.** This removes the reference tree and its large memory requirement.
- **FragGeneScanRs replaces Prodigal/pyrodigal.** It is faster, and accurate enough for AAI and marker extraction. It runs as a separate program because it is GPL-3.0 and aeon is MIT. It translates with table 11 only, so `translation_table` is always 11.
- **Species calls use GTDB's species radius**, so they match GTDB-Tk's ANI screen. 95% applies to most species.
- **TurboANI runs only when the sketch search cannot settle the species on its own.**
  - A species call made from the sketch estimate alone reports `closest_genome_af` as `N/A`. Its Jaccard-based estimate already drops with the shared fraction.
  - The 1-point margin sends borderline genomes to TurboANI. These include chimeras whose AF would fail, because the sketch cannot check AF.
- **Genus 65% AAI and family 45% AAI** are the defaults, the genus and family boundaries proposed by Konstantinidis, Rosselló-Móra & Amann ([ISME J 2017](https://doi.org/10.1038/ismej.2017.113)). With the AAI margin below, they made 3 wrong calls in 104 holdout queries.
- **AAI calls use a margin, not a vote.** A novel genus sits about equally far from several genera in its family, so its closest reference leads the next genus by only a point or two. A genus that is in the database stands out from its neighbours. Unlike a majority vote among the closest references, the margin doesn't penalise genera with few representatives; many GTDB genera have only one.
- **Marker-search calls stop at order**, using the closest marker reference's lineage.
  - The 25% floor only rejects chance k-mer matches. Random marker sequences share almost no 5-mers (J ≈ 0.001).
  - It does not separate distant phyla from real but fast-evolving relatives. In the holdouts, other phyla reached up to 56% estimated identity, and Mycoplasmatales matched their own order at 42–50%.
  - So a genome from a phylum missing from the database still gets the order of its closest phylum.

**Output differences from GTDB-Tk:**
- `pplacer_taxonomy` and `red_value` are always `N/A`.
- For marker-search calls, `closest_placement_reference`, `closest_placement_radius` and `closest_placement_taxonomy` describe the closest marker reference.
- `classification_method` is one of:
  - `ani_screen`;
  - `taxonomic classification defined by AAI`;
  - `taxonomic classification defined by marker gene search`.
- Extra files:
  - `classify/aai/{prefix}.aai_summary.tsv`: every AAI hit, with the SD of the reciprocal-best-hit identities, their count, and the smaller proteome's protein count;
  - `classify/markers/{prefix}.marker_summary.tsv`: every marker-index hit, with an estimated marker identity.
- `align/` holds the user MSA only, not the full reference MSA.

```bash
# External tools: DIAMOND and HMMER (environment.yml), and FragGeneScanRs.
conda env create -f environment.yml && conda activate aeon-classify
cargo install frag_gene_scan_rs

# Reference database. Genomes and proteomes are named by GTDB accession.
aeon classify-db --db gtdb_db --gtdbtk-data $GTDBTK_DATA_PATH \
  -r gtdb_reps.txt --protein-list gtdb_protein_reps.txt   # or --call-genes

aeon classify --db gtdb_db --genome_dir genomes -x fna.gz --out_dir out --cpus 8
```

`--call-genes` calls reference genes with FragGeneScanRs instead of reading `--protein-list` (e.g. GTDB `protein_faa_reps`). aeon finds FragGeneScanRs's `train/` models next to the executable or in the sources `cargo install` keeps. Otherwise, pass `--fgs-train-dir`.

**Holdout results.** Holdouts of GTDB R232 species representatives were used, each with the query's species, whole genus or whole family withheld from the database.

| Set | Orders | References | Queries (species / genus / family withheld) |
|---|---|---|---|
| Training | Lactobacillales, Campylobacterales, Bacteroidales, Flavobacteriales, Oscillospirales | 272 | 55 (20 / 15 / 20) |
| Validation | Enterobacterales, Rhodobacterales, Lachnospirales, Mycoplasmatales, Actinomycetales | 221 | 49 (20 / 13 / 16) |

A genus or family holdout has no reference of its true taxon, so the right answer is to stop one rank higher. A "wrong" call assigns a neighbouring genus or family; in every case the ranks above it were right.

All rules below use 65%/45% AAI. Results are over both sets (104 queries):

| AAI rule | Wrong calls | Species holdouts given their genus (of 40) | Genus holdouts given their family (of 28) | Family holdouts stopped at order (of 36) |
|---|---|---|---|---|
| Margin 5 points (default) | 3 | 35 | 20 | 35 |
| Margin 3 points | 5 | 36 | 23 | 33 |
| Margin 7 points | 1 | 32 | 20 | 36 |
| Margin 10 points | 1 | 30 | 17 | 36 |
| Strict majority of the 10 closest references | 18 | 22 | 20 | 23 |
| Closest reference only, no margin | 46 | 37 | 18 | 0 |

The default row comes from full classify runs. The other rows replay the rule over the same AAI hits. The margin was chosen on these same queries.

- **Marker search:** all 42 order-level calls were correct.
- **Why a fixed AAI cutoff is not enough.** GTDB draws genus and family boundaries by relative evolutionary divergence, not AAI, so a fixed AAI cutoff cuts across them:

  | Comparison (per query) | Min | Median | Max |
  |---|---|---|---|
  | Same genus (species holdouts vs references of their own genus) | 53.0 | 74.0 | 93.9 |
  | Closest other genus in the same family (genus holdouts) | 54.0 | 61.6 | 73.0 |
  | Closest other family in the same order (family holdouts) | 48.8 | 55.3 | 62.5 |

  - 10 of 28 genus holdouts had a sister genus at ≥65% AAI, while 10% of same-genus pairs fell below 64.3%.
  - Every family holdout (36 of 36) had another family at ≥45% AAI, so the 45% cutoff alone stops none of them.
- **Why the margin works.** When the true taxon is absent, the closest reference led the next taxon by less than 3 points in 41 of 46 cases (the largest lead was 13.0). When the closest reference was right, it usually led by 6 points or more; 8 of 65 led by less than 5.
- **Why a margin, not a majority vote.** The vote blocked 18 correct calls, 15 of them genera with only 2 references in the database that were outvoted by a larger neighbour.
- **DIAMOND mode:** the default mode gives AAI about 2 points higher than `--sensitive` on the same pairs (median +2.2 and +2.3 on the two sets). It made about the same number of wrong calls under the majority vote (19 against 21) in about a third of the time.
- **Species calls from the sketch alone:** reference genomes cut into 10 kb pieces, keeping 50–90% of them, were all assigned their own species, 24 of 24, with no TurboANI call.

**Runtime.** On the 55-genome training set with 8 threads, classify took 80 s wall time and 570 CPU-seconds, with a peak RSS of 0.95 GB. Building the 272-reference database with `--call-genes` took 126 s, with a peak RSS of 1.7 GB. Index scaling at full GTDB size has not been measured yet.

## Update Semantics

- Updates are written to temporary files before replacing database files.
- Vector IDs may change after deletion because live records are compacted.
- Aeon rebuilds `.genomes.txt` and `.idmap.tsv` from rust-diskann's old-to-new ID mapping after every update.
- `insert`, `delete`, and `update` finish with the same static format produced by `todisk`.
- `delete` and deletion-containing `update` commands are admitted only when the temporary landmark probe preserves every route that was reachable before the proposed batch; a deferred batch writes no replacement index.

## Benchmark

Apple Silicon M4 Max benchmark on 2,799 RefSeq Streptomycetaceae genomes (6.4 GB input), using 8 threads, DNA k=16, OPTDENS sketch size 8,192, graph degree 64, and construction beam/ef 256:

| Operation | aeon | GSearch | Comparison |
|---|---:|---:|---:|
| Build wall time | 50.38 s | 42.22 s | aeon 19% slower |
| Build peak RSS | 675 MiB | 991 MiB | aeon 32% lower |
| Database size, including mappings | 48.23 MB | 47.43 MB | aeon 1.7% larger |
| Search 100 genomes | 1.63 s | 2.40 s | aeon 1.47x faster |
| Search peak RSS | 394 MiB | 571 MiB | aeon 31% lower |
| Top-1 zero-distance matches | 100/100 | 100/100 | equal |

For the larger GTDB search benchmark, the database contained 199,923 representative genomes. Both indexes used DNA k=16, OPTDENS sketches of 12,000 `u16` values, graph degree 128, and 8 search threads. Aeon used search beam 512; GSearch used its stored HNSW `ef=1600` configuration. Values below are medians from three warm mmap runs of the same 100 query genomes on an Apple Silicon M4 Max:

| Search 100 genomes against GTDB-NT | aeon | GSearch | Comparison |
|---|---:|---:|---:|
| Wall time | 0.89 s | 6.92 s | aeon 7.78x faster |
| Peak RSS | 4.44 GiB | 5.17 GiB | aeon 14% lower |

For large number of reference genome databases, e.g., millions of genomes, aeon will be orders of magnitude faster than GSearch


A combined 1% replacement (28 MERIT deletions followed by 28 insertions and one static commit) took 4.53 seconds with 222 MiB peak RSS. This was 11.1x faster than rebuilding the aeon database, and all 28 replacement queries returned a zero-distance top hit.

## References

- Subramanya et al., 2019 *DiskANN: Fast Accurate Billion-point Nearest Neighbor Search on a Single Node*. NeurIPS 2019.
- Wu et.al., 2026. *MERIT: Efficient In-Place Deletion for Dynamic Graph-Based Approximate Nearest Neighbor Indexes*. arXiv preprint arXiv:2607.29173.
- Li et.al., 2012. *One permutation hashing*. Advances in Neural Information Processing Systems 25.
- Li et al., 2010. *b-Bit Minwise Hashing*. WWW 2010.
- Ondov et al., 2016. *Mash: Fast genome and metagenome distance estimation using MinHash*. Genome Biology 17:132.
