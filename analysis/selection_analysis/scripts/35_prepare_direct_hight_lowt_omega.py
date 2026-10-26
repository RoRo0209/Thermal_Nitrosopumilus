#!/usr/bin/env python3
"""Prepare direct HighT-versus-LowT three-category codeml branch models.

The input foreground trees are reused from the established Figure 7 pipeline.
For each gene:

* Other crown-external branches are class #0.
* HighT within-crown branches are class #1.
* LowT within-crown branches are class #2 in the alternative model.
* HighT and LowT within-crown branches share class #1 in the null model.

The stem entering each crown is excluded, matching the scope of the existing
separate-foreground branch models. Both terminal and internal branches inside
each crown are labelled; QC explicitly rejects tip-only foreground labelling.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import shutil
from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class Node:
    name: str = ""
    children: list["Node"] = field(default_factory=list)
    branch_class: int = 0

    @property
    def is_leaf(self) -> bool:
        return not self.children


class NewickParser:
    def __init__(self, text: str):
        self.text = text.strip()
        self.i = 0

    def parse(self) -> Node:
        node = self._subtree()
        self._space()
        if self.i < len(self.text) and self.text[self.i] == ";":
            self.i += 1
        self._space()
        if self.i != len(self.text):
            raise ValueError(
                f"Unexpected Newick text at position {self.i}: "
                f"{self.text[self.i:self.i + 40]}"
            )
        return node

    def _space(self) -> None:
        while self.i < len(self.text) and self.text[self.i].isspace():
            self.i += 1

    def _token(self) -> str:
        self._space()
        start = self.i
        while self.i < len(self.text) and self.text[self.i] not in ":,();":
            self.i += 1
        return self.text[start:self.i].strip()

    def _skip_length(self) -> None:
        self._space()
        if self.i < len(self.text) and self.text[self.i] == ":":
            self.i += 1
            while self.i < len(self.text) and self.text[self.i] not in ",();":
                self.i += 1

    @staticmethod
    def _split_class(token: str) -> tuple[str, int]:
        for marker in ("#2", "#1", "#0"):
            if token.endswith(marker):
                return token[: -len(marker)], int(marker[1:])
        return token, 0

    def _subtree(self) -> Node:
        self._space()
        if self.i >= len(self.text):
            raise ValueError("Unexpected end of Newick")
        if self.text[self.i] == "(":
            self.i += 1
            children = [self._subtree()]
            while True:
                self._space()
                if self.text[self.i] == ",":
                    self.i += 1
                    children.append(self._subtree())
                elif self.text[self.i] == ")":
                    self.i += 1
                    break
                else:
                    raise ValueError(f"Expected ',' or ')' at position {self.i}")
            token = self._token()
            name, branch_class = self._split_class(token)
            self._skip_length()
            return Node(name=name, children=children, branch_class=branch_class)
        token = self._token()
        if not token:
            raise ValueError(f"Missing leaf name at position {self.i}")
        name, branch_class = self._split_class(token)
        self._skip_length()
        return Node(name=name, branch_class=branch_class)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-job-dir", type=Path, required=True)
    parser.add_argument("--tree-dir", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument(
        "--starts", default="0.05,0.2,1.0",
        help="Comma-separated starting omega values for both models."
    )
    return parser.parse_args()


def leaf_set(node: Node) -> set[str]:
    if node.is_leaf:
        return {node.name}
    result: set[str] = set()
    for child in node.children:
        result.update(leaf_set(child))
    return result


def find_lca(node: Node, targets: set[str]) -> Node:
    for child in node.children:
        if targets <= leaf_set(child):
            return find_lca(child, targets)
    return node


def marked_leaf_set(node: Node) -> set[str]:
    if node.is_leaf:
        return {node.name} if node.branch_class > 0 else set()
    result: set[str] = set()
    for child in node.children:
        result.update(marked_leaf_set(child))
    return result


def clone(node: Node) -> Node:
    return Node(
        name=node.name,
        children=[clone(child) for child in node.children],
        branch_class=node.branch_class,
    )


def clear_classes(node: Node) -> None:
    node.branch_class = 0
    for child in node.children:
        clear_classes(child)


def mark_subtree_branches(node: Node, branch_class: int) -> None:
    node.branch_class = branch_class
    for child in node.children:
        mark_subtree_branches(child, branch_class)


def mark_crown(root: Node, targets: set[str], branch_class: int) -> Node:
    if len(targets) < 2:
        raise ValueError("Crown has fewer than two retained targets")
    lca = find_lca(root, targets)
    descendants = leaf_set(lca)
    if descendants != targets:
        raise ValueError(
            "Crown targets are not monophyletic: "
            f"targets={len(targets)}, descendants={len(descendants)}"
        )
    if lca is root:
        raise ValueError("Crown LCA is the root")
    # Exclude the stem entering the LCA. Label every descendant branch.
    for child in lca.children:
        mark_subtree_branches(child, branch_class)
    return lca


def branch_counts(node: Node, branch_class: int) -> tuple[int, int]:
    terminal = 0
    internal = 0
    if node.branch_class == branch_class:
        if node.is_leaf:
            terminal += 1
        else:
            internal += 1
    for child in node.children:
        child_terminal, child_internal = branch_counts(child, branch_class)
        terminal += child_terminal
        internal += child_internal
    return terminal, internal


def serialize(node: Node) -> str:
    if node.is_leaf:
        text = node.name
    else:
        text = "(" + ",".join(serialize(child) for child in node.children) + ")"
    if node.branch_class > 0:
        text += f"#{node.branch_class}"
    return text


def read_phylip_header(path: Path) -> tuple[int, int]:
    first = path.read_text(encoding="ascii", errors="strict").splitlines()[0].split()
    if len(first) < 2:
        raise ValueError(f"Invalid PHYLIP header: {path}")
    return int(first[0]), int(first[1])


def ctl_text(treefile: str, omega: float) -> str:
    return f"""seqfile = ../alignment.phy
treefile = ../{treefile}
outfile = results.txt

noisy = 3
verbose = 1
runmode = 0
seqtype = 1
CodonFreq = 2
clock = 0
aaDist = 0
model = 2
NSsites = 0
icode = 0
Mgene = 0
fix_kappa = 0
kappa = 2
fix_omega = 0
omega = {omega}
fix_alpha = 1
alpha = 0
Malpha = 0
ncatG = 10
getSE = 0
RateAncestor = 1
Small_Diff = 1e-6
cleandata = 1
method = 0
"""


def write_tree(path: Path, tree: Node, n_taxa: int) -> str:
    newick = serialize(tree) + ";"
    path.write_text(f"{n_taxa} 1\n{newick}\n", encoding="ascii", newline="\n")
    return hashlib.sha256(newick.encode("ascii")).hexdigest()


def write_runs(model_dir: Path, model_name: str, treefile: str, starts: list[float]) -> None:
    for index, omega in enumerate(starts, start=1):
        run_dir = model_dir / f"{model_name}_rep{index}"
        run_dir.mkdir(parents=True, exist_ok=True)
        (run_dir / "codeml.ctl").write_text(
            ctl_text(treefile, omega), encoding="ascii", newline="\n"
        )


def main() -> None:
    args = parse_args()
    starts = [float(value) for value in args.starts.split(",") if value.strip()]
    if len(starts) < 2:
        raise ValueError("At least two starting omega values are required")
    args.outdir.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, object]] = []
    source_families = sorted(path for path in args.source_job_dir.iterdir() if path.is_dir())
    for source_family in source_families:
        family = source_family.name
        alignment_source = source_family / "alignment.phy"
        plain_path = args.tree_dir / "species_tree" / f"{family}.nwk"
        high_path = args.tree_dir / "foreground_HighT_clade" / f"{family}.nwk"
        low_path = args.tree_dir / "foreground_LowT_clade" / f"{family}.nwk"
        required = (alignment_source, plain_path, high_path, low_path)
        missing = [str(path) for path in required if not path.is_file()]
        if missing:
            raise FileNotFoundError(f"{family}: missing inputs: {missing}")

        n_taxa, alignment_nt = read_phylip_header(alignment_source)
        plain = NewickParser(plain_path.read_text(encoding="utf-8-sig")).parse()
        high_existing = NewickParser(high_path.read_text(encoding="utf-8-sig")).parse()
        low_existing = NewickParser(low_path.read_text(encoding="utf-8-sig")).parse()
        high_targets = marked_leaf_set(high_existing)
        low_targets = marked_leaf_set(low_existing)
        all_tips = leaf_set(plain)

        if len(all_tips) != n_taxa:
            raise ValueError(
                f"{family}: tree/alignment taxon mismatch: {len(all_tips)} versus {n_taxa}"
            )
        if high_targets & low_targets:
            raise ValueError(f"{family}: HighT and LowT tip sets overlap")
        if not high_targets <= all_tips or not low_targets <= all_tips:
            raise ValueError(f"{family}: focal targets absent from plain family tree")

        null_tree = clone(plain)
        clear_classes(null_tree)
        mark_crown(null_tree, high_targets, 1)
        mark_crown(null_tree, low_targets, 1)

        alt_tree = clone(plain)
        clear_classes(alt_tree)
        mark_crown(alt_tree, high_targets, 1)
        mark_crown(alt_tree, low_targets, 2)

        high_terminal, high_internal = branch_counts(alt_tree, 1)
        low_terminal, low_internal = branch_counts(alt_tree, 2)
        null_terminal, null_internal = branch_counts(null_tree, 1)
        if min(high_internal, low_internal) < 1:
            raise ValueError(f"{family}: foreground labelling is tip-only")
        if high_terminal != len(high_targets) or low_terminal != len(low_targets):
            raise ValueError(f"{family}: not all focal terminal branches were labelled")
        if null_terminal != high_terminal + low_terminal:
            raise ValueError(f"{family}: null foreground does not pool both crowns")

        family_dir = args.outdir / family
        family_dir.mkdir(exist_ok=True)
        shutil.copy2(alignment_source, family_dir / "alignment.phy")
        null_hash = write_tree(family_dir / "null_H_equals_L.tree", null_tree, n_taxa)
        alt_hash = write_tree(family_dir / "alternative_H_L_O.tree", alt_tree, n_taxa)
        write_runs(family_dir, "null", "null_H_equals_L.tree", starts)
        write_runs(family_dir, "alternative", "alternative_H_L_O.tree", starts)

        rows.append({
            "split_family": family,
            "n_taxa": n_taxa,
            "alignment_nt": alignment_nt,
            "alignment_codons": alignment_nt // 3,
            "n_HighT_tips": len(high_targets),
            "n_LowT_tips": len(low_targets),
            "HighT_terminal_branches_labelled": high_terminal,
            "HighT_internal_branches_labelled": high_internal,
            "LowT_terminal_branches_labelled": low_terminal,
            "LowT_internal_branches_labelled": low_internal,
            "null_pooled_terminal_branches_labelled": null_terminal,
            "null_pooled_internal_branches_labelled": null_internal,
            "foreground_tip_sets_overlap": 0,
            "HighT_monophyletic": 1,
            "LowT_monophyletic": 1,
            "stem_branches_excluded": 1,
            "tip_only_labelling": 0,
            "null_tree_sha256": null_hash,
            "alternative_tree_sha256": alt_hash,
            "n_start_values_per_model": len(starts),
            "start_values": ",".join(str(value) for value in starts),
            "tree_qc_pass": 1,
        })

    if len(rows) != 110:
        raise ValueError(f"Expected 110 codeml-ready families; observed {len(rows)}")

    fields = list(rows[0])
    with (args.outdir / "direct_model_tree_qc.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    n_runs = len(rows) * 2 * len(starts)
    print(f"Families prepared: {len(rows)}")
    print(f"Models per family: 2; starting values per model: {len(starts)}")
    print(f"codeml runs prepared: {n_runs}")
    print("Tree QC: all HighT/LowT crowns monophyletic, non-overlapping,")
    print("         internal + terminal crown branches labelled, stems excluded")


if __name__ == "__main__":
    main()
