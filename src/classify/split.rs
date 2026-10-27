//! Split-tree placement of bacterial genomes, ported from GTDB-Tk 2.8's
//! `split.py` and the class-level branch of `_classify_red_topology`:
//! pplacer places genomes into a backbone tree, the class each genome gets
//! there picks a class-level tree, and rules 1-6 reconcile the two
//! classifications.

use std::collections::{BTreeMap, HashMap, HashSet};
use std::error::Error;
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::Path;

use crate::classify::placement::{
    PlacedTree, Placement, rank_index, rank_prefix, split_taxa, taxon_of,
};
use crate::classify::refdata::{GtdbtkData, RANK_PREFIXES, standardise_taxonomy};
use crate::classify::summary::py_float;

const PHYLUM: usize = 1;
const CLASS: usize = 2;
const MARKER_SET: &str = "bac120";

const BACKBONE_NOTE: &str = "classification based on placement in backbone tree";
const CLASS_LEVEL_NOTE: &str = "classification based on placement in class-level tree";
const CONSENSUS_NOTE: &str =
    "classification based on consensus between backbone and class-level tree";

/// A genome's classification in the backbone tree
/// (`Split.get_high_pplacer_taxonomy`).
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct BackboneCall {
    pub(crate) taxonomy_red: String,
    pub(crate) taxonomy_terminal: String,
    pub(crate) pplacer_taxonomy: String,
    pub(crate) is_terminal: bool,
    pub(crate) red: f64,
}

/// Classify every user genome in a backbone tree, in leaf order.
pub(crate) fn classify_backbone(
    tree_path: &Path,
    mrca_red: &Path,
    users: &HashSet<String>,
    data: &GtdbtkData,
    red_dict: &HashMap<String, f64>,
) -> Result<Vec<(String, BackboneCall)>, Box<dyn Error>> {
    let placed = PlacedTree::load(tree_path, mrca_red, users, data)?;
    let backbone = Backbone {
        placed: &placed,
        users,
        data,
        red_dict,
        lineages: lineage_index(data),
    };
    placed
        .user_leaves(users)
        .map(|(leaf, label)| Ok((label.to_owned(), backbone.classify(leaf)?)))
        .collect()
}

/// For each taxon, the lineage of a genome that has it. GTDB lineages are
/// consistent, so any genome gives the same ranks above the taxon.
fn lineage_index(data: &GtdbtkData) -> HashMap<&str, &[String]> {
    let mut lineages = HashMap::new();
    for lineage in data.taxonomy.values() {
        for taxon in lineage {
            lineages.entry(taxon.as_str()).or_insert(lineage.as_slice());
        }
    }
    lineages
}

/// Python's `items[start:end]` for non-negative bounds.
fn py_slice(items: &[String], start: usize, end: usize) -> &[String] {
    let end = end.min(items.len());
    if start >= end {
        &[]
    } else {
        &items[start..end]
    }
}

struct Backbone<'a> {
    placed: &'a PlacedTree,
    users: &'a HashSet<String>,
    data: &'a GtdbtkData,
    red_dict: &'a HashMap<String, f64>,
    lineages: HashMap<&'a str, &'a [String]>,
}

impl<'a> Backbone<'a> {
    fn red_for(&self, taxon: &str) -> Result<f64, String> {
        self.red_dict
            .get(rank_prefix(taxon))
            .copied()
            .ok_or_else(|| format!("no median RED for rank {}", rank_prefix(taxon)))
    }

    fn taxonomy(&self, node: usize) -> Result<&'a [String], Box<dyn Error>> {
        let label = self
            .placed
            .tree
            .label(node)
            .unwrap_or_default()
            .replace('\'', "");
        self.data
            .taxonomy_of(&label)
            .ok_or_else(|| format!("{label} has no GTDB taxonomy").into())
    }

    fn classify(&self, leaf: usize) -> Result<BackboneCall, Box<dyn Error>> {
        let tree = &self.placed.tree;
        // Walk to the root collecting decorations, the first RED below one,
        // and whether the genome sits on a terminal branch (its first
        // internal ancestor holds a single other genome).
        let mut taxa = Vec::new();
        let mut current_red = 1.0;
        let (mut terminal, mut tested) = (None, false);
        let mut node = leaf;
        while let Some(parent) = tree.parent(node) {
            if let Some(red) = self.placed.red[node] {
                if current_red == 1.0 && red < 1.0 {
                    current_red = red;
                }
            }
            if !tree.is_leaf(node) && !tested {
                let others = tree
                    .leaves_under(node)
                    .filter(|&n| !tree.label(n).is_some_and(|l| self.users.contains(l)))
                    .collect::<Vec<_>>();
                if others.len() == 1 {
                    terminal = Some(self.taxonomy(others[0])?);
                }
                tested = !others.is_empty();
            }
            if let Some(taxon) = taxon_of(tree, node) {
                taxa.extend(split_taxa(&taxon).into_iter().rev());
            }
            node = parent;
        }
        taxa.reverse();
        let pplacer_taxonomy = taxa.join(";");
        let mut taxa_str = pplacer_taxonomy.clone();

        let mut terminal_str = String::new();
        if let Some(term) = terminal {
            // The reference's ranks below the lowest decoration, down to genus.
            let leaf_ranks = if taxa.len() > 1 {
                let lowest = taxa.last().expect("non-empty");
                let position = term.iter().position(|t| t == lowest).ok_or_else(|| {
                    format!("{lowest} is not in the taxonomy of the neighbouring reference")
                })?;
                py_slice(term, position + 1, term.len() - 1)
            } else {
                taxa_str = "d__Bacteria".into();
                py_slice(term, 1, term.len() - 1)
            };
            let parent_rank = rank_prefix(taxa_str.rsplit(';').next().unwrap_or_default());
            terminal_str =
                self.classify_on_terminal_branch(leaf_ranks, current_red, parent_rank, term)?;
        }

        let mut parent_node = tree.parent(leaf).ok_or("placed genome is the root")?;
        let parent_taxon = loop {
            if let Some(taxon) = taxon_of(tree, parent_node) {
                break split_taxa(&taxon).pop().expect("split yields a taxon");
            }
            parent_node = tree
                .parent(parent_node)
                .ok_or("no decorated ancestor; the tree root needs a domain label")?;
        };
        let mut red_str = String::new();
        if rank_prefix(&parent_taxon) != "g__" {
            let mut node = leaf;
            while self.placed.references_under(node).is_empty() {
                node = tree
                    .parent(node)
                    .ok_or("no reference genome above placement")?;
            }
            let child_index = rank_index(rank_prefix(&parent_taxon))
                .ok_or_else(|| format!("unknown rank {parent_taxon}"))?
                + 1;
            let child_rank = *RANK_PREFIXES
                .get(child_index)
                .ok_or("species-level decoration")?;
            let references = self.placed.references_under(node);
            let names = references
                .iter()
                .map(|&reference| self.taxonomy(reference).map(|t| t[child_index].as_str()))
                .collect::<Result<HashSet<_>, _>>()?;
            if names.len() == 1 {
                let mut child_taxons = Vec::new();
                let mut child_red = 0.0;
                for &sub in tree.subtree(node) {
                    if tree.is_leaf(sub) {
                        continue;
                    }
                    if let Some(taxon) =
                        taxon_of(tree, sub).filter(|t| t.trim().starts_with(child_rank))
                    {
                        child_taxons = split_taxa(&taxon);
                        child_red = self.placed.red[sub].ok_or_else(|| {
                            format!("{taxon} has no RED value; is the RED file complete?")
                        })?;
                        break;
                    }
                }
                (red_str, terminal_str) = self.classify_on_internal_branch(
                    &child_taxons,
                    current_red,
                    child_red,
                    &references,
                    &parent_taxon,
                    child_index,
                    taxa_str,
                    terminal_str,
                    terminal.is_some(),
                )?;
            } else {
                red_str = taxa_str;
            }
        }

        Ok(BackboneCall {
            taxonomy_red: standardise_taxonomy(&red_str, Some(MARKER_SET)),
            taxonomy_terminal: standardise_taxonomy(&terminal_str, Some(MARKER_SET)),
            pplacer_taxonomy: standardise_taxonomy(&pplacer_taxonomy, Some(MARKER_SET)),
            is_terminal: terminal.is_some(),
            red: current_red,
        })
    }

    /// `_classify_on_terminal_branch`: guess the lower ranks from the
    /// neighbouring reference by RED, never above class.
    fn classify_on_terminal_branch(
        &self,
        leaf_ranks: &[String],
        current_red: f64,
        parent_rank: &str,
        term: &[String],
    ) -> Result<String, Box<dyn Error>> {
        let mut closest = None;
        for (index, taxon) in leaf_ranks.iter().enumerate().rev() {
            let other = if index == 0 {
                parent_rank
            } else {
                &leaf_ranks[index - 1]
            };
            if (current_red - self.red_for(taxon)?).abs()
                < (current_red - self.red_for(other)?).abs()
            {
                closest = Some(rank_prefix(taxon));
                break;
            }
        }
        let closest = closest.unwrap_or(parent_rank);
        let index = rank_index(closest).ok_or_else(|| format!("unknown rank {closest}"))?;
        Ok(py_slice(term, 1, index.max(CLASS) + 1).join(";"))
    }

    /// `_classify_on_internal_branch`: the RED classification of a genome
    /// whose references all share one child rank, and the class implied for
    /// genomes off terminal branches.
    #[allow(clippy::too_many_arguments)]
    fn classify_on_internal_branch(
        &self,
        child_taxons: &[String],
        current_red: f64,
        child_red: f64,
        references: &[usize],
        parent_taxon: &str,
        child_index: usize,
        mut taxa_str: String,
        mut terminal_str: String,
        is_terminal: bool,
    ) -> Result<(String, String), Box<dyn Error>> {
        let mut closest: Option<String> = None;
        if child_taxons.is_empty() {
            if references.len() != 1 {
                // GTDB-Tk looks up the full taxon in its rank list here and fails.
                return Err(format!(
                    "cannot classify a genome above {} references of one {} without a decorated {} node (GTDB-Tk fails here too)",
                    references.len(),
                    RANK_PREFIXES[child_index],
                    RANK_PREFIXES[child_index]
                )
                .into());
            }
            let leaf_ranks = py_slice(self.taxonomy(references[0])?, child_index, 6);
            for (index, taxon) in leaf_ranks.iter().enumerate().rev() {
                let other = if index == 0 {
                    parent_taxon
                } else {
                    &leaf_ranks[index - 1]
                };
                if (current_red - self.red_for(taxon)?).abs()
                    < (current_red - self.red_for(other)?).abs()
                {
                    closest = Some(taxon.clone());
                    break;
                }
            }
            closest.get_or_insert_with(|| parent_taxon.to_owned());
        }
        for (index, taxon) in child_taxons.iter().enumerate().rev() {
            let distance = (current_red - self.red_for(taxon)?).abs();
            let child_distance = (child_red - self.red_for(taxon)?).abs();
            if index == 0 {
                if distance < child_distance
                    && distance < (current_red - self.red_for(parent_taxon)?).abs()
                {
                    closest = Some(taxon.clone());
                } else if closest.is_none() {
                    closest = Some(parent_taxon.to_owned());
                }
            } else if distance < (current_red - self.red_for(&child_taxons[index - 1])?).abs()
                && distance < child_distance
            {
                closest = Some(taxon.clone());
                break;
            }
        }
        if let Some(lineage) = closest.and_then(|taxon| {
            let index = rank_index(rank_prefix(&taxon))?;
            self.lineages
                .get(taxon.as_str())
                .map(|lineage| (*lineage, index))
        }) {
            let (lineage, index) = lineage;
            taxa_str = py_slice(lineage, 1, index + 1).join(";");
            // Genomes placed on an internal branch above a single class are
            // classified at least to class.
            if child_taxons.iter().any(|t| t.starts_with("c__"))
                && index < CLASS
                && (terminal_str.split(';').any(|t| t == "c__") || !is_terminal)
            {
                terminal_str = py_slice(lineage, 1, CLASS + 1).join(";");
            }
        }
        Ok((taxa_str, terminal_str))
    }
}

/// The class-level tree a backbone call maps to (`Split.map_high_taxonomy`):
/// by its RED class, else its terminal-branch class, else `None` (the
/// backbone classification stands).
pub(crate) fn class_tree<'a>(
    call: &BackboneCall,
    mapping: &'a HashMap<String, String>,
) -> Result<Option<&'a str>, Box<dyn Error>> {
    for taxonomy in [&call.taxonomy_red, &call.taxonomy_terminal] {
        let class = taxonomy.split(';').nth(CLASS).unwrap_or_default();
        if class.len() > 3 {
            return mapping
                .get(class)
                .map(|tree| Some(tree.as_str()))
                .ok_or_else(|| format!("{class} is not in the class-level tree mapping").into());
        }
    }
    Ok(None)
}

/// A genome the backbone tree alone classifies (rule 1 of `map_high_taxonomy`).
pub(crate) fn backbone_placement(call: &BackboneCall) -> Placement {
    Placement {
        pplacer_taxonomy: call.pplacer_taxonomy.clone(),
        classification: call.taxonomy_red.clone(),
        red_value: call.red,
        method: None,
        reference: None,
        note: Some(BACKBONE_NOTE),
        warning: None,
    }
}

/// Where rules 1-6 left a genome placed in a class-level tree.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Reconciled {
    pub(crate) placement: Placement,
    pub(crate) rule: &'static str,
    /// Rule 1 falls back to the backbone tree.
    pub(crate) on_backbone: bool,
}

/// Rules 1-6 of GTDB-Tk's `_classify_red_topology` for a genome classified
/// in both trees. GTDB-Tk reports the class-level RED taxonomy, not the
/// pplacer decorations, as `pplacer_taxonomy` here.
pub(crate) fn reconcile(
    backbone: &BackboneCall,
    mut placement: Placement,
    valid_classes: &HashSet<String>,
    valid_phyla: &HashSet<String>,
) -> Reconciled {
    let backbone_tax = backbone.taxonomy_red.split(';').collect::<Vec<_>>();
    let class_tax = placement.classification.split(';').collect::<Vec<_>>();
    fn rank<'t>(taxonomy: &[&'t str], index: usize) -> &'t str {
        taxonomy.get(index).copied().unwrap_or_default()
    }
    // `limit_rank`: keep ranks down to `index`, empty below.
    fn limit<'t>(taxonomy: &[&'t str], index: usize) -> Vec<&'t str> {
        taxonomy
            .iter()
            .take(index + 1)
            .copied()
            .chain(RANK_PREFIXES.iter().skip(index + 1).copied())
            .collect()
    }
    let (backbone_class, class) = (rank(&backbone_tax, CLASS), rank(&class_tax, CLASS));
    let (backbone_phylum, phylum) = (rank(&backbone_tax, PHYLUM), rank(&class_tax, PHYLUM));
    let (classification, rule, note, warning) = if backbone_class == "c__" && class == "c__" {
        let classification = backbone.taxonomy_red.clone();
        placement.pplacer_taxonomy = classification.clone();
        placement.red_value = backbone.red;
        placement.classification = classification;
        placement.note = Some(BACKBONE_NOTE);
        return Reconciled {
            placement,
            rule: "Rule 1",
            on_backbone: true,
        };
    } else if valid_classes.contains(class) && backbone_class == class {
        (class_tax.clone(), "Rule 2", CLASS_LEVEL_NOTE, None)
    } else if valid_classes.contains(class) {
        (
            class_tax.clone(),
            "Rule 3",
            CLASS_LEVEL_NOTE,
            Some("backbone tree has incongruent class classification"),
        )
    } else if valid_phyla.contains(phylum) && backbone_phylum == phylum {
        (
            limit(&class_tax, PHYLUM),
            "Rule 4",
            CLASS_LEVEL_NOTE,
            Some("classification restricted to valid phylum"),
        )
    } else if valid_phyla.contains(phylum) {
        (
            limit(&class_tax, PHYLUM),
            "Rule 5",
            CLASS_LEVEL_NOTE,
            Some("backbone tree has incongruent phylum classification"),
        )
    } else {
        let consensus = backbone_tax
            .iter()
            .zip(&class_tax)
            .enumerate()
            .map(|(index, (backbone_taxon, class_taxon))| {
                if backbone_taxon == class_taxon {
                    *backbone_taxon
                } else {
                    RANK_PREFIXES[index]
                }
            })
            .collect();
        (consensus, "Rule 6", CONSENSUS_NOTE, None)
    };
    placement.pplacer_taxonomy = placement.classification.clone();
    placement.classification = classification.join(";");
    placement.note = Some(note);
    placement.warning = warning;
    Reconciled {
        placement,
        rule,
        on_backbone: false,
    }
}

/// `split/class_level/tree_mapping.tsv`: class name to class-level tree.
pub(crate) fn read_tree_mapping(path: &Path) -> Result<HashMap<String, String>, Box<dyn Error>> {
    let text = fs::read_to_string(path)
        .map_err(|error| format!("cannot read {}: {error}", path.display()))?;
    text.lines()
        .filter(|line| !line.trim().is_empty())
        .map(|line| {
            let mut fields = line.split_whitespace();
            match (fields.next(), fields.next()) {
                (Some(class), Some(tree)) => Ok((class.to_owned(), tree.to_owned())),
                _ => Err(format!("malformed tree mapping line: {line}").into()),
            }
        })
        .collect()
}

/// Phyla of every genome in the given classes (`get_authorised_rank`).
pub(crate) fn phyla_of_classes(classes: &HashSet<String>, data: &GtdbtkData) -> HashSet<String> {
    data.taxonomy
        .values()
        .filter(|lineage| {
            lineage
                .get(CLASS)
                .is_some_and(|class| classes.contains(class))
        })
        .filter_map(|lineage| lineage.get(PHYLUM).cloned())
        .collect()
}

/// `classify/intermediate_results/{prefix}.bac120.backbone.classification_pplacer.tsv`.
pub(crate) fn write_backbone_classification(
    path: &Path,
    calls: &[(String, BackboneCall)],
) -> std::io::Result<()> {
    let mut output = BufWriter::new(File::create(path)?);
    writeln!(
        output,
        "user_genome\tgtdb_taxonomy_red\tgtdb_taxonomy_terminal\tpplacer_taxonomy\tis_terminal\tred"
    )?;
    let sorted = calls
        .iter()
        .map(|(g, c)| (g, c))
        .collect::<BTreeMap<_, _>>();
    for (genome, call) in sorted {
        writeln!(
            output,
            "{genome}\t{}\t{}\t{}\t{}\t{}",
            call.taxonomy_red,
            call.taxonomy_terminal,
            call.pplacer_taxonomy,
            if call.is_terminal { "True" } else { "False" },
            py_float(call.red)
        )?;
    }
    output.flush()
}

/// `classify/{prefix}.bac120.tree.mapping.tsv`: the tree and rule behind
/// each genome's classification.
pub(crate) fn write_tree_mapping(
    path: &Path,
    rows: &BTreeMap<String, (String, &'static str)>,
) -> std::io::Result<()> {
    let mut output = BufWriter::new(File::create(path)?);
    writeln!(
        output,
        "user_genome\tis_ani_classification\tclass_tree_mapped\tclassification_rule"
    )?;
    for (genome, (tree, rule)) in rows {
        writeln!(output, "{genome}\tFalse\t{tree}\t{rule}")?;
    }
    output.flush()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn call(red: &str) -> BackboneCall {
        BackboneCall {
            taxonomy_red: red.into(),
            taxonomy_terminal: "d__Bacteria;p__;c__;o__;f__;g__;s__".into(),
            pplacer_taxonomy: "d__Bacteria;p__;c__;o__;f__;g__;s__".into(),
            is_terminal: false,
            red: 0.4,
        }
    }

    fn placement(classification: &str) -> Placement {
        Placement {
            pplacer_taxonomy: "d__Bacteria;p__P;c__;o__;f__;g__;s__".into(),
            classification: classification.into(),
            red_value: 0.9,
            method: Some("taxonomic novelty determined using RED"),
            reference: None,
            note: None,
            warning: None,
        }
    }

    fn set(items: &[&str]) -> HashSet<String> {
        items.iter().map(|s| (*s).to_owned()).collect()
    }

    #[test]
    fn rules_follow_gtdbtk_precedence() {
        let classes = set(&["c__A", "c__B"]);
        let phyla = set(&["p__P"]);
        let low = "d__Bacteria;p__P;c__A;o__O;f__;g__;s__";
        let cases = [
            (
                "d__Bacteria;p__P;c__;o__;f__;g__;s__",
                "d__Bacteria;p__P;c__;o__;f__;g__;s__",
                "Rule 1",
            ),
            ("d__Bacteria;p__P;c__A;o__;f__;g__;s__", low, "Rule 2"),
            ("d__Bacteria;p__P;c__B;o__;f__;g__;s__", low, "Rule 3"),
        ];
        for (backbone, class_level, rule) in cases {
            let result = reconcile(&call(backbone), placement(class_level), &classes, &phyla);
            assert_eq!(result.rule, rule);
        }
        let rule1 = reconcile(
            &call("d__Bacteria;p__P;c__;o__;f__;g__;s__"),
            placement("d__Bacteria;p__P;c__;o__;f__;g__;s__"),
            &classes,
            &phyla,
        );
        assert!(rule1.on_backbone);
        assert_eq!(rule1.placement.red_value, 0.4);
        assert_eq!(rule1.placement.note, Some(BACKBONE_NOTE));

        // A class outside the tree's classes but a valid phylum.
        let off = "d__Bacteria;p__P;c__Z;o__O;f__F;g__;s__";
        let rule4 = reconcile(
            &call("d__Bacteria;p__P;c__A;o__;f__;g__;s__"),
            placement(off),
            &classes,
            &phyla,
        );
        assert_eq!(rule4.rule, "Rule 4");
        assert_eq!(
            rule4.placement.classification,
            "d__Bacteria;p__P;c__;o__;f__;g__;s__"
        );
        assert_eq!(rule4.placement.pplacer_taxonomy, off);
        assert_eq!(
            rule4.placement.warning,
            Some("classification restricted to valid phylum")
        );
        let rule5 = reconcile(
            &call("d__Bacteria;p__Q;c__X;o__;f__;g__;s__"),
            placement(off),
            &classes,
            &phyla,
        );
        assert_eq!(rule5.rule, "Rule 5");

        let rule6 = reconcile(
            &call("d__Bacteria;p__Q;c__X;o__O;f__;g__;s__"),
            placement("d__Bacteria;p__R;c__Y;o__O;f__;g__;s__"),
            &classes,
            &phyla,
        );
        assert_eq!(rule6.rule, "Rule 6");
        assert_eq!(
            rule6.placement.classification,
            "d__Bacteria;p__;c__;o__O;f__;g__;s__"
        );
        assert_eq!(rule6.placement.note, Some(CONSENSUS_NOTE));
    }

    #[test]
    fn backbone_terminal_branch_takes_ranks_from_the_neighbour() {
        let dir = std::env::temp_dir().join(format!("aeon-split-{}", std::process::id()));
        fs::create_dir_all(&dir).unwrap();
        let tree = dir.join("backbone.tree");
        // Q was placed beside R3, 0.2 below p__A on a 0.3 edge.
        fs::write(
            &tree,
            "(((R1:0.1,R2:0.1)'100:c__X':0.2,(R3:0.1,Q:0.05):0.2)'100:p__A':0.3,(R4:0.5,R5:0.5)'100:p__B':0.3)d__Bacteria;",
        )
        .unwrap();
        let red = dir.join("red.tsv");
        fs::write(
            &red,
            "R1|R5\t0.0\nR1|R3\t0.3\nR1|R2\t0.6\nR4|R5\t0.3\nR1\t1.0\nR2\t1.0\nR3\t1.0\nR4\t1.0\nR5\t1.0\n",
        )
        .unwrap();
        let lineage = |ranks: &str| {
            format!("d__Bacteria;{ranks}")
                .split(';')
                .map(str::to_owned)
                .collect::<Vec<_>>()
        };
        let x = lineage("p__A;c__X;o__O1;f__F1;g__G1;s__G1 a");
        let data = GtdbtkData {
            taxonomy: HashMap::from([
                ("R1".to_owned(), x.clone()),
                ("R2".to_owned(), x),
                (
                    "R3".to_owned(),
                    lineage("p__A;c__Y;o__O2;f__F2;g__G2;s__G2 a"),
                ),
                (
                    "R4".to_owned(),
                    lineage("p__B;c__Z;o__O3;f__F3;g__G3;s__G3 a"),
                ),
                (
                    "R5".to_owned(),
                    lineage("p__B;c__Z;o__O3;f__F3;g__G3;s__G3 b"),
                ),
            ]),
            ..GtdbtkData::default()
        };
        let red_dict = [
            ("d__", 0.0),
            ("p__", 0.34),
            ("c__", 0.476),
            ("o__", 0.647),
            ("f__", 0.797),
            ("g__", 0.944),
        ]
        .into_iter()
        .map(|(rank, red)| (rank.to_owned(), red))
        .collect::<HashMap<_, _>>();
        let users = set(&["Q"]);
        let calls = classify_backbone(&tree, &red, &users, &data, &red_dict).unwrap();
        fs::remove_dir_all(dir).unwrap();
        let (genome, call) = &calls[0];
        assert_eq!(genome, "Q");
        assert!(call.is_terminal);
        assert!((call.red - (0.3 + 0.7 * 0.2 / 0.3)).abs() < 1e-12);
        // RED 0.77 is nearest the family median, so both calls stop at f__F2.
        let expected = "d__Bacteria;p__A;c__Y;o__O2;f__F2;g__;s__";
        assert_eq!(call.taxonomy_red, expected);
        assert_eq!(call.taxonomy_terminal, expected);
        assert_eq!(
            call.pplacer_taxonomy,
            "d__Bacteria;p__A;c__;o__;f__;g__;s__"
        );
    }

    #[test]
    fn class_tree_prefers_the_red_class() {
        let mapping = HashMap::from([
            ("c__A".to_owned(), "1".to_owned()),
            ("c__T".to_owned(), "2".to_owned()),
        ]);
        let mut backbone = call("d__Bacteria;p__P;c__A;o__;f__;g__;s__");
        backbone.taxonomy_terminal = "d__Bacteria;p__P;c__T;o__;f__;g__;s__".into();
        assert_eq!(class_tree(&backbone, &mapping).unwrap(), Some("1"));
        backbone.taxonomy_red = "d__Bacteria;p__P;c__;o__;f__;g__;s__".into();
        assert_eq!(class_tree(&backbone, &mapping).unwrap(), Some("2"));
        backbone.taxonomy_terminal = backbone.taxonomy_red.clone();
        assert_eq!(class_tree(&backbone, &mapping).unwrap(), None);
    }
}
