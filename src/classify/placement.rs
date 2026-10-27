//! Tier 3 placement: pplacer puts genomes into the GTDB reference tree and
//! relative evolutionary divergence (RED) turns the placement into a
//! taxonomy, ported from GTDB-Tk 2.8's `classify.py` (`_assign_mrca_red`,
//! `_get_pplacer_taxonomy`, `_get_pplacer_reference`,
//! `_classify_red_topology` and `_get_redtax`). Split-tree reconciliation
//! lives in `split.rs`.

use std::collections::{BTreeMap, HashMap, HashSet, VecDeque};
use std::error::Error;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;
use std::process::Command;

use crate::classify::refdata::{GtdbtkData, RANK_PREFIXES, standardise_taxonomy};
use crate::classify::summary::py_float;

/// A rooted tree; a leaf's `label` is its taxon name, an internal node's
/// label is GTDB's `support:taxa` decoration.
#[derive(Debug, Default)]
pub(crate) struct Tree {
    nodes: Vec<Node>,
    /// Node indices in preorder, and each node's position in it.
    preorder: Vec<usize>,
    position: Vec<usize>,
    /// Exclusive end of each node's subtree in `preorder`.
    subtree_end: Vec<usize>,
}

#[derive(Debug, Default, Clone)]
struct Node {
    parent: Option<usize>,
    children: Vec<usize>,
    label: Option<String>,
    length: f64,
}

const ROOT: usize = 0;

impl Tree {
    /// Parse Newick as dendropy does with `preserve_underscores=True`.
    pub(crate) fn parse(text: &str) -> Result<Self, Box<dyn Error>> {
        let mut tree = Self {
            nodes: vec![Node::default()],
            ..Self::default()
        };
        let chars = text.trim().chars().collect::<Vec<_>>();
        let mut current = ROOT;
        let mut index = 0;
        while index < chars.len() {
            match chars[index] {
                '(' => current = tree.add_child(current),
                ',' => {
                    let parent = tree.nodes[current].parent.ok_or("unbalanced Newick ','")?;
                    current = tree.add_child(parent);
                }
                ')' => current = tree.nodes[current].parent.ok_or("unbalanced Newick ')'")?,
                ';' => break,
                ':' => {
                    let start = index + 1;
                    while index + 1 < chars.len() && !"(),;[".contains(chars[index + 1]) {
                        index += 1;
                    }
                    let number = chars[start..=index].iter().collect::<String>();
                    tree.nodes[current].length = number
                        .trim()
                        .parse()
                        .map_err(|_| format!("invalid Newick branch length {number:?}"))?;
                }
                '[' => {
                    while index < chars.len() && chars[index] != ']' {
                        index += 1;
                    }
                }
                '\'' => {
                    let mut label = String::new();
                    index += 1;
                    while index < chars.len() {
                        if chars[index] == '\'' {
                            if chars.get(index + 1) == Some(&'\'') {
                                label.push('\'');
                                index += 1;
                            } else {
                                break;
                            }
                        } else {
                            label.push(chars[index]);
                        }
                        index += 1;
                    }
                    tree.nodes[current].label = Some(label);
                }
                c if c.is_whitespace() => {}
                _ => {
                    let start = index;
                    while index + 1 < chars.len()
                        && !"():,;[".contains(chars[index + 1])
                        && !chars[index + 1].is_whitespace()
                    {
                        index += 1;
                    }
                    tree.nodes[current].label = Some(chars[start..=index].iter().collect());
                }
            }
            index += 1;
        }
        tree.index();
        Ok(tree)
    }

    fn add_child(&mut self, parent: usize) -> usize {
        let child = self.nodes.len();
        self.nodes.push(Node {
            parent: Some(parent),
            ..Node::default()
        });
        self.nodes[parent].children.push(child);
        child
    }

    fn index(&mut self) {
        let count = self.nodes.len();
        self.preorder = Vec::with_capacity(count);
        self.position = vec![0; count];
        self.subtree_end = vec![0; count];
        let mut stack = vec![ROOT];
        while let Some(node) = stack.pop() {
            self.position[node] = self.preorder.len();
            self.preorder.push(node);
            stack.extend(self.nodes[node].children.iter().rev());
        }
        for &node in self.preorder.iter().rev() {
            self.subtree_end[node] = self.nodes[node]
                .children
                .last()
                .map_or(self.position[node] + 1, |&child| self.subtree_end[child]);
        }
    }

    pub(super) fn is_leaf(&self, node: usize) -> bool {
        self.nodes[node].children.is_empty()
    }

    pub(super) fn label(&self, node: usize) -> Option<&str> {
        self.nodes[node].label.as_deref()
    }

    pub(super) fn parent(&self, node: usize) -> Option<usize> {
        self.nodes[node].parent
    }

    /// Subtree of `node` in preorder (including `node`).
    pub(super) fn subtree(&self, node: usize) -> &[usize] {
        &self.preorder[self.position[node]..self.subtree_end[node]]
    }

    /// Leaves under `node`, left to right.
    pub(super) fn leaves_under(&self, node: usize) -> impl Iterator<Item = usize> + '_ {
        self.subtree(node)
            .iter()
            .copied()
            .filter(|&n| self.is_leaf(n))
    }

    fn leaves(&self) -> impl Iterator<Item = usize> + '_ {
        self.leaves_under(ROOT)
    }

    fn is_ancestor(&self, ancestor: usize, node: usize) -> bool {
        self.position[ancestor] <= self.position[node]
            && self.position[node] < self.subtree_end[ancestor]
    }

    fn mrca(&self, left: usize, right: usize) -> usize {
        let mut node = left;
        while !self.is_ancestor(node, right) {
            node = self.parent(node).expect("nodes share the root");
        }
        node
    }

    fn levelorder(&self, node: usize) -> impl Iterator<Item = usize> + '_ {
        let mut queue = VecDeque::from([node]);
        std::iter::from_fn(move || {
            let next = queue.pop_front()?;
            queue.extend(self.nodes[next].children.iter().copied());
            Some(next)
        })
    }
}

/// GTDB-Tk's `parse_label`: `support:taxon|aux`, a bare support value, or a taxon.
fn parse_label(label: Option<&str>) -> (Option<f64>, Option<String>) {
    let Some(label) = label.map(str::trim) else {
        return (None, None);
    };
    let label = label.split('|').next().unwrap_or(label);
    if let Some((support, taxon)) = label.split_once(':') {
        return (support.parse().ok(), Some(taxon.to_owned()));
    }
    match label.parse::<f64>() {
        Ok(support) => (Some(support), None),
        Err(_) if !label.is_empty() => (None, Some(label.to_owned())),
        Err(_) => (None, None),
    }
}

/// An internal node's taxon decoration; dendropy keeps leaf names apart
/// from labels, so leaves have none.
pub(super) fn taxon_of(tree: &Tree, node: usize) -> Option<String> {
    if tree.is_leaf(node) {
        return None;
    }
    parse_label(tree.label(node)).1
}

/// Ranks of a decoration such as `p__X;c__Y`, trimmed.
pub(super) fn split_taxa(taxon: &str) -> Vec<String> {
    taxon.split(';').map(|t| t.trim().to_owned()).collect()
}

pub(super) fn rank_prefix(taxon: &str) -> &str {
    taxon.get(..3).unwrap_or(taxon)
}

pub(super) fn rank_index(prefix: &str) -> Option<usize> {
    RANK_PREFIXES.iter().position(|rank| *rank == prefix)
}

/// Run pplacer on a user MSA against a reference package and convert its
/// placements to a tree with `guppy tog`.
#[allow(clippy::too_many_arguments)]
pub(crate) fn run_pplacer(
    pplacer: &Path,
    guppy: &Path,
    refpkg: &Path,
    user_msa: &Path,
    json_out: &Path,
    pplacer_out: &Path,
    tree_out: &Path,
    threads: usize,
) -> Result<(), Box<dyn Error>> {
    if let Some(dir) = json_out.parent() {
        fs::create_dir_all(dir)?;
    }
    let output = Command::new(pplacer)
        .args(["-m", "wag", "-j", &threads.clamp(1, 64).to_string(), "-c"])
        .arg(refpkg)
        .arg("-o")
        .arg(json_out)
        .arg(user_msa)
        .output()
        .map_err(|error| format!("failed to start pplacer: {error}"))?;
    fs::write(pplacer_out, &output.stdout)?;
    if !output.status.success() || !json_out.is_file() {
        return Err(format!(
            "An error was encountered while running pplacer, check the log file: {} ({})",
            pplacer_out.display(),
            String::from_utf8_lossy(&output.stderr).trim()
        )
        .into());
    }
    let output = Command::new(guppy)
        .arg("tog")
        .arg("-o")
        .arg(tree_out)
        .arg(json_out)
        .output()
        .map_err(|error| format!("failed to start guppy: {error}"))?;
    if !output.status.success() || !tree_out.is_file() {
        return Err(format!(
            "An error was encountered while running tog: {}",
            String::from_utf8_lossy(&output.stderr).trim()
        )
        .into());
    }
    Ok(())
}

pub(crate) fn program_version(program: &Path) -> String {
    Command::new(program)
        .arg("--version")
        .output()
        .ok()
        .map(|output| String::from_utf8_lossy(&output.stdout).trim().to_owned())
        .filter(|version| !version.is_empty())
        .unwrap_or_else(|| "unknown".to_owned())
}

/// RED of reference nodes from GTDB-Tk's `mrca_red` file, interpolated onto
/// the nodes pplacer inserted (`_assign_mrca_red`).
fn assign_mrca_red(tree: &Tree, mrca_red: &Path) -> Result<Vec<Option<f64>>, Box<dyn Error>> {
    let leaf_map = tree
        .leaves()
        .filter_map(|leaf| tree.label(leaf).map(|label| (label.to_owned(), leaf)))
        .collect::<HashMap<_, _>>();
    let mut red = vec![None; tree.nodes.len()];
    let mut reference = vec![false; tree.nodes.len()];
    let lookup = |label: &str| {
        leaf_map
            .get(label)
            .copied()
            .ok_or_else(|| format!("{label} from the RED file is not in the tree"))
    };
    for line in BufReader::new(File::open(mrca_red)?).lines() {
        let line = line?;
        let Some((labels, value)) = line.trim().split_once('\t') else {
            continue;
        };
        let node = match labels.split_once('|') {
            Some((left, right)) => tree.mrca(lookup(left)?, lookup(right)?),
            None => lookup(labels)?,
        };
        red[node] = Some(value.parse()?);
        reference[node] = true;
    }

    let leaves = tree.leaves().collect::<Vec<_>>();
    for &leaf in &leaves {
        if reference[leaf] {
            continue;
        }
        red[leaf] = Some(1.0);
        let mut placed = leaf;
        let mut parent = tree
            .parent(placed)
            .ok_or("user genome placed at the root")?;
        while !tree.leaves_under(placed).any(|n| reference[n]) {
            placed = parent;
            parent = tree.parent(placed).ok_or("no reference above placement")?;
        }
        let child = tree
            .levelorder(placed)
            .find(|&n| reference[n])
            .ok_or("no reference below placement")?;
        while !reference[parent] {
            parent = tree.parent(parent).ok_or("no reference ancestor")?;
        }
        let mut current = tree.parent(child).ok_or("reference child at root")?;
        let mut edge_length = tree.nodes[child].length.max(1e-9);
        let (mut on_branch, mut placed_length) = (false, 0.0);
        while current != parent {
            if on_branch || current == placed {
                on_branch = true;
                placed_length += tree.nodes[current].length;
            }
            edge_length += tree.nodes[current].length;
            current = tree.parent(current).ok_or("walked past the root")?;
        }
        let (child_red, parent_red) = (red[child].unwrap_or(1.0), red[parent].unwrap_or(0.0));
        red[placed] = Some(parent_red + (child_red - parent_red) * (placed_length / edge_length));
    }
    Ok(red)
}

/// A user genome's placement and its RED-based classification.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Placement {
    pub(crate) pplacer_taxonomy: String,
    pub(crate) classification: String,
    pub(crate) red_value: f64,
    /// `None` for genomes classified from the backbone tree alone, which
    /// GTDB-Tk reports without a method.
    pub(crate) method: Option<&'static str>,
    /// The single reference genome nearest the placement, as labelled in the tree.
    pub(crate) reference: Option<String>,
    /// Split-tree note and warning (`classification based on placement in ...`).
    pub(crate) note: Option<&'static str>,
    pub(crate) warning: Option<&'static str>,
}

/// A pplacer tree with RED on its reference and placement nodes.
pub(super) struct PlacedTree {
    pub(super) tree: Tree,
    pub(super) red: Vec<Option<f64>>,
    /// Leaves that are GTDB references rather than user genomes.
    pub(super) is_reference: Vec<bool>,
}

impl PlacedTree {
    pub(super) fn load(
        tree_path: &Path,
        mrca_red: &Path,
        users: &HashSet<String>,
        data: &GtdbtkData,
    ) -> Result<Self, Box<dyn Error>> {
        let tree = Tree::parse(&fs::read_to_string(tree_path)?)?;
        let red = assign_mrca_red(&tree, mrca_red)?;
        let is_reference = (0..tree.nodes.len())
            .map(|node| {
                tree.is_leaf(node)
                    && tree.label(node).is_some_and(|label| {
                        !users.contains(label)
                            && data.taxonomy_of(&label.replace('\'', "")).is_some()
                    })
            })
            .collect();
        Ok(Self {
            tree,
            red,
            is_reference,
        })
    }

    /// Placed user genomes, left to right.
    pub(super) fn user_leaves<'a>(
        &'a self,
        users: &'a HashSet<String>,
    ) -> impl Iterator<Item = (usize, &'a str)> + 'a {
        self.tree.leaves().filter_map(move |leaf| {
            self.tree
                .label(leaf)
                .filter(|label| users.contains(*label))
                .map(|label| (leaf, label))
        })
    }

    /// Reference leaves under `node`.
    pub(super) fn references_under(&self, node: usize) -> Vec<usize> {
        self.tree
            .leaves_under(node)
            .filter(|&n| self.is_reference[n])
            .collect()
    }
}

const NOVELTY: &str = "taxonomic novelty determined using RED";
const TOPOLOGY: &str = "taxonomic classification fully defined by topology";

/// Classify every user genome placed in `tree_path` (GTDB-Tk's full-tree path).
pub(crate) fn classify_tree(
    tree_path: &Path,
    mrca_red: &Path,
    marker_set: &str,
    users: &HashSet<String>,
    data: &GtdbtkData,
    red_dict: &HashMap<String, f64>,
) -> Result<BTreeMap<String, Placement>, Box<dyn Error>> {
    let placed = PlacedTree::load(tree_path, mrca_red, users, data)?;
    let tree = &placed.tree;
    let references = pplacer_references(tree, &placed.is_reference);
    let mut placements = BTreeMap::new();
    for (leaf, label) in placed.user_leaves(users) {
        let (classification, red_value, method) = classify_red_topology(
            tree,
            leaf,
            &placed.red,
            &placed.is_reference,
            data,
            red_dict,
        )?;
        placements.insert(
            label.to_owned(),
            Placement {
                pplacer_taxonomy: pplacer_taxonomy(tree, leaf, marker_set),
                classification,
                red_value,
                method: Some(method),
                reference: references
                    .get(&leaf)
                    .and_then(|&node| tree.label(node).map(str::to_owned)),
                note: None,
                warning: None,
            },
        );
    }
    Ok(placements)
}

/// Decorations above a placed genome (`_get_pplacer_taxonomy`).
fn pplacer_taxonomy(tree: &Tree, leaf: usize, marker_set: &str) -> String {
    let mut taxa = Vec::new();
    let mut node = leaf;
    while let Some(parent) = tree.parent(node) {
        if let Some(taxon) = taxon_of(tree, node) {
            taxa.extend(split_taxa(&taxon).into_iter().rev());
        }
        node = parent;
    }
    taxa.reverse();
    standardise_taxonomy(&taxa.join(";"), Some(marker_set))
}

/// For each placed leaf, the reference genome under the first ancestor that
/// holds exactly one reference (`_get_pplacer_reference`).
fn pplacer_references(tree: &Tree, is_reference: &[bool]) -> HashMap<usize, usize> {
    let count = tree.nodes.len();
    let (mut references, mut single) = (vec![0usize; count], vec![None; count]);
    for &node in tree.preorder.iter().rev() {
        if tree.is_leaf(node) {
            if is_reference[node] {
                references[node] = 1;
                single[node] = Some(node);
            }
            continue;
        }
        for &child in &tree.nodes[node].children {
            references[node] += references[child];
            if references[child] == 1 {
                single[node] = single[child];
            }
        }
        if references[node] > 1 {
            single[node] = None;
        }
    }
    let mut placed = HashMap::new();
    for leaf in tree.leaves().filter(|&leaf| !is_reference[leaf]) {
        let mut node = tree.parent(leaf);
        while let Some(current) = node {
            if references[current] == 1 {
                if let Some(reference) = single[current] {
                    placed.insert(leaf, reference);
                }
                break;
            } else if references[current] > 1 {
                break;
            }
            node = tree.parent(current);
        }
    }
    placed
}

/// `_classify_red_topology` for one placed genome: its taxonomy, RED value
/// and classification method.
fn classify_red_topology(
    tree: &Tree,
    leaf: usize,
    red: &[Option<f64>],
    is_reference: &[bool],
    data: &GtdbtkData,
    red_dict: &HashMap<String, f64>,
) -> Result<(String, f64, &'static str), Box<dyn Error>> {
    let reference_leaves = |node: usize| {
        tree.leaves_under(node)
            .filter(|&n| is_reference[n])
            .collect::<Vec<_>>()
    };
    let taxonomy = |node: usize| -> Result<Vec<String>, Box<dyn Error>> {
        let label = tree.label(node).unwrap_or_default().replace('\'', "");
        data.taxonomy_of(&label)
            .map(<[String]>::to_vec)
            .ok_or_else(|| format!("{label} has no GTDB taxonomy").into())
    };
    let red_of = |node: usize| {
        red[node]
            .ok_or_else(|| format!("node {node} has no RED value; is the mrca_red file complete?"))
    };
    let red_for = |taxon: &str| -> Result<f64, String> {
        red_dict
            .get(rank_prefix(taxon))
            .copied()
            .ok_or_else(|| format!("no median RED for rank {}", rank_prefix(taxon)))
    };

    let mut node = tree.parent(leaf).ok_or("placed genome is the root")?;
    while reference_leaves(node).is_empty() {
        node = tree
            .parent(node)
            .ok_or("no reference genome above placement")?;
    }
    let current_red = red_of(node)?;
    let mut parent_taxon_node = tree
        .parent(node)
        .ok_or("placement directly below the root")?;
    let parent_taxon = loop {
        if let Some(taxon) = taxon_of(tree, parent_taxon_node) {
            break taxon;
        }
        parent_taxon_node = tree
            .parent(parent_taxon_node)
            .ok_or("no decorated ancestor; the tree root needs a domain label")?;
    };
    let parent_rank = rank_prefix(split_taxa(&parent_taxon).last().expect("split")).to_owned();
    let parent_index =
        rank_index(&parent_rank).ok_or_else(|| format!("unknown rank {parent_rank}"))?;

    let mut child_taxons: Vec<String> = Vec::new();
    let mut child_red = 0.0;
    let mut closest_rank: Option<String> = None;
    let mut method = NOVELTY;
    let child_rank = RANK_PREFIXES.get(parent_index + 1).copied();
    if parent_rank != "g__" {
        let child_rank = child_rank.ok_or("species-level decoration")?;
        let child_index = parent_index + 1;
        let ranks = reference_leaves(node)
            .into_iter()
            .map(|reference| taxonomy(reference).map(|t| t[child_index].clone()))
            .collect::<Result<HashSet<_>, _>>()?;
        if ranks.len() == 1 {
            for &sub in tree.subtree(node) {
                if tree.is_leaf(sub) {
                    continue;
                }
                if let Some(taxon) =
                    taxon_of(tree, sub).filter(|t| t.trim().starts_with(child_rank))
                {
                    child_taxons = split_taxa(&taxon);
                    child_red = red_of(sub)?;
                    break;
                }
            }
        } else {
            closest_rank = Some(parent_rank.clone());
            method = TOPOLOGY;
        }
    } else {
        closest_rank = Some(parent_rank.clone());
        method = TOPOLOGY;
    }

    if child_taxons.is_empty() && closest_rank.is_none() {
        let child_index = parent_index + 1;
        let leaves = reference_leaves(node);
        if leaves.len() != 1 {
            let subranks = leaves
                .iter()
                .map(|&reference| taxonomy(reference).map(|t| t[child_index].clone()))
                .collect::<Result<HashSet<_>, _>>()?;
            if subranks.len() == 1 {
                return Err("There should be only one leaf.".into());
            }
            closest_rank = Some(parent_rank.clone());
            method = TOPOLOGY;
        }
        // GTDB-Tk compares RED even after a topology call above.
        let leaf_ranks = taxonomy(leaves[0])?[child_index..6].to_vec();
        for (index, leaf_taxon) in leaf_ranks.iter().enumerate().rev() {
            let distance = (current_red - red_for(leaf_taxon)?).abs();
            if *leaf_taxon == leaf_ranks[0] {
                if distance < (current_red - red_for(&parent_rank)?).abs() {
                    closest_rank = Some(rank_prefix(leaf_taxon).to_owned());
                    break;
                }
            } else if distance < (current_red - red_for(&leaf_ranks[index - 1])?).abs() {
                closest_rank = Some(rank_prefix(leaf_taxon).to_owned());
                break;
            }
        }
        if closest_rank.is_none() {
            closest_rank = Some(parent_rank.clone());
        }
    }

    for (index, child_taxon) in child_taxons.iter().enumerate().rev() {
        let distance = (current_red - red_for(child_taxon)?).abs();
        let child_distance = (child_red - red_for(child_taxon)?).abs();
        if *child_taxon == child_taxons[0] {
            if distance < child_distance && distance < (current_red - red_for(&parent_rank)?).abs()
            {
                closest_rank = Some(rank_prefix(child_taxon).to_owned());
            } else if closest_rank.is_none() {
                closest_rank = Some(parent_rank.clone());
            }
        } else if distance < (current_red - red_for(&child_taxons[index - 1])?).abs()
            && distance < child_distance
        {
            // GTDB-Tk keeps the full taxon here, not just its rank prefix.
            closest_rank = Some(child_taxon.clone());
            break;
        }
    }
    let closest_rank = closest_rank.ok_or("closest rank is None")?;

    let subtree = tree.leaves_under(node).collect::<Vec<_>>();
    let red_taxonomy = red_taxonomy(&subtree, &closest_rank, is_reference, &taxonomy)?;
    Ok((
        standardise_taxonomy(&red_taxonomy.join(";"), None),
        current_red,
        method,
    ))
}

/// `_get_redtax`: the first reference's ranks above `closest_rank`, plus the
/// closest-rank taxon when every reference agrees on it.
fn red_taxonomy(
    leaves: &[usize],
    closest_rank: &str,
    is_reference: &[bool],
    taxonomy: &dyn Fn(usize) -> Result<Vec<String>, Box<dyn Error>>,
) -> Result<Vec<String>, Box<dyn Error>> {
    let (mut subtax, mut multirefrank) = (Vec::new(), Vec::new());
    let mut initial = true;
    for &leaf in leaves.iter().filter(|&&leaf| is_reference[leaf]) {
        for rank in taxonomy(leaf)? {
            if rank.starts_with(closest_rank) {
                multirefrank.push(rank);
                initial = false;
                break;
            } else if initial {
                subtax.push(rank);
            }
        }
    }
    let distinct = multirefrank.iter().collect::<HashSet<_>>();
    if distinct.len() == 1 {
        subtax.push(multirefrank[0].clone());
    } else {
        subtax.push(format!("{closest_rank}undefined"));
    }
    Ok(subtax)
}

/// `classify/intermediate_results/{prefix}.{marker_set}.classification_pplacer.tsv`.
pub(crate) fn write_pplacer_classification(
    path: &Path,
    placements: &BTreeMap<String, Placement>,
) -> std::io::Result<()> {
    let mut output = BufWriter::new(File::create(path)?);
    for (genome, placement) in placements {
        writeln!(output, "{genome}\t{}", placement.pplacer_taxonomy)?;
    }
    output.flush()
}

/// `classify/intermediate_results/{prefix}.{marker_set}.red_dictionary.tsv`.
pub(crate) fn write_red_dictionary(
    path: &Path,
    red_dict: &HashMap<String, f64>,
) -> std::io::Result<()> {
    let mut output = BufWriter::new(File::create(path)?);
    for (name, prefix) in [
        ("Phylum", "p__"),
        ("Class", "c__"),
        ("Order", "o__"),
        ("Family", "f__"),
        ("Genus", "g__"),
    ] {
        let value = red_dict
            .get(prefix)
            .map_or_else(|| "None".to_owned(), |v| py_float(*v));
        writeln!(output, "{name}\t{value}")?;
    }
    output.flush()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_gtdb_newick_labels() {
        let tree = Tree::parse("((A:1,B:2)'100.0:g__X':0.5,C:1e-3)d__Bacteria;").unwrap();
        let labels = tree
            .leaves()
            .map(|leaf| tree.label(leaf).unwrap().to_owned())
            .collect::<Vec<_>>();
        assert_eq!(labels, ["A", "B", "C"]);
        let inner = tree.parent(tree.leaves().next().unwrap()).unwrap();
        assert_eq!(
            parse_label(tree.label(inner)),
            (Some(100.0), Some("g__X".into()))
        );
        assert_eq!(tree.nodes[inner].length, 0.5);
        assert_eq!(tree.label(ROOT), Some("d__Bacteria"));
        assert_eq!(tree.mrca(2, 4), ROOT);
        assert_eq!(tree.mrca(2, 3), inner);
    }

    #[test]
    fn parse_label_matches_gtdbtk() {
        assert_eq!(parse_label(Some("82.0")), (Some(82.0), None));
        assert_eq!(
            parse_label(Some("p__A;c__B")),
            (None, Some("p__A;c__B".into()))
        );
        assert_eq!(
            parse_label(Some("99:o__O|x")),
            (Some(99.0), Some("o__O".into()))
        );
        assert_eq!(parse_label(None), (None, None));
    }

    #[test]
    fn interpolates_red_on_the_placement_edge() {
        // Reference edge root(0.0) -> AB(0.8) of length 1.0; pplacer split it
        // 0.25 below the root and hung Q off the new node.
        let dir = std::env::temp_dir().join(format!("aeon-red-{}", std::process::id()));
        fs::create_dir_all(&dir).unwrap();
        let red_file = dir.join("red.tsv");
        fs::write(&red_file, "A|C\t0.0\nA|B\t0.8\nA\t1.0\nB\t1.0\nC\t1.0\n").unwrap();
        let tree = Tree::parse("(((A:1,B:1):0.75,Q:0.1):0.25,C:2);").unwrap();
        let red = assign_mrca_red(&tree, &red_file).unwrap();
        let q = tree.leaves().find(|&n| tree.label(n) == Some("Q")).unwrap();
        let placed = tree.parent(q).unwrap();
        assert!((red[placed].unwrap() - 0.2).abs() < 1e-12);
        assert_eq!(red[q], Some(1.0));
        fs::remove_dir_all(dir).unwrap();
    }
}
