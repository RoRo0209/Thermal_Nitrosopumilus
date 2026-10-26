#!/usr/bin/env python3
"""Prune the master species tree per family and label HighT/LowT foreground stems."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class Node:
    name: str = ""
    children: list["Node"] = field(default_factory=list)
    foreground: bool = False

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
            raise ValueError(f"Unexpected Newick text at position {self.i}: {self.text[self.i:self.i+30]}")
        return node

    def _space(self):
        while self.i < len(self.text) and self.text[self.i].isspace():
            self.i += 1

    def _token(self) -> str:
        self._space()
        start = self.i
        while self.i < len(self.text) and self.text[self.i] not in ":,();":
            self.i += 1
        return self.text[start:self.i].strip()

    def _skip_length(self):
        self._space()
        if self.i < len(self.text) and self.text[self.i] == ":":
            self.i += 1
            while self.i < len(self.text) and self.text[self.i] not in ",();":
                self.i += 1

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
            name = self._token()  # Internal support/name is not needed by codeml.
            self._skip_length()
            return Node(name=name, children=children)
        name = self._token()
        if not name:
            raise ValueError(f"Missing leaf name at position {self.i}")
        self._skip_length()
        return Node(name=name)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--master-tree", type=Path, required=True)
    p.add_argument("--protein-dir", type=Path, required=True)
    p.add_argument("--genome-tree-map", type=Path, required=True)
    p.add_argument("--hight-list", type=Path, required=True)
    p.add_argument("--lowt-list", type=Path, required=True)
    p.add_argument("--outdir", type=Path, required=True)
    return p.parse_args()


def read_fasta_ids(path: Path) -> list[str]:
    ids = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                ids.append(line[1:].split()[0])
    return ids


def clone(node: Node) -> Node:
    return Node(node.name, [clone(child) for child in node.children], node.foreground)


def prune(node: Node, keep: set[str]) -> Node | None:
    if node.is_leaf:
        return Node(name=node.name) if node.name in keep else None
    children = [x for child in node.children if (x := prune(child, keep)) is not None]
    if not children:
        return None
    if len(children) == 1:
        return children[0]
    return Node(children=children)


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


def serialize(node: Node) -> str:
    if node.is_leaf:
        text = node.name
    else:
        text = "(" + ",".join(serialize(child) for child in node.children) + ")"
    if node.foreground:
        text += "#1"
    return text


def mark_stem_foreground(root: Node, targets: set[str]):
    if len(targets) < 2:
        return False, set(), "fewer_than_two_targets"
    lca = find_lca(root, targets)
    descendants = leaf_set(lca)
    monophyletic = descendants == targets
    if not monophyletic:
        return False, descendants, "not_monophyletic"
    if lca is root:
        return False, descendants, "foreground_lca_is_root"
    lca.foreground = True
    return True, descendants, "ok"


def mark_all_branches(node: Node) -> None:
    node.foreground = True
    for child in node.children:
        mark_all_branches(child)


def mark_crown_foreground(root: Node, targets: set[str]):
    if len(targets) < 2:
        return False, set(), "fewer_than_two_targets"
    lca = find_lca(root, targets)
    descendants = leaf_set(lca)
    if descendants != targets:
        return False, descendants, "not_monophyletic"
    # Label every branch inside the crown but not the stem entering the crown LCA.
    for child in lca.children:
        mark_all_branches(child)
    return True, descendants, "ok"


def read_names(path: Path) -> list[str]:
    return [line.strip() for line in path.read_text(encoding="utf-8-sig").splitlines() if line.strip()]


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    plain_dir = args.outdir / "species_tree"
    high_stem_dir = args.outdir / "foreground_HighT_stem"
    low_stem_dir = args.outdir / "foreground_LowT_stem"
    high_clade_dir = args.outdir / "foreground_HighT_clade"
    low_clade_dir = args.outdir / "foreground_LowT_clade"
    for directory in (plain_dir, high_stem_dir, low_stem_dir, high_clade_dir, low_clade_dir):
        directory.mkdir(exist_ok=True)

    mapping = {}
    with args.genome_tree_map.open(encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            mapping[row["genome"]] = row["tree_id"]
    high_all = {mapping[name] for name in read_names(args.hight_list)}
    low_all = {mapping[name] for name in read_names(args.lowt_list)}
    master = NewickParser(args.master_tree.read_text(encoding="utf-8-sig")).parse()
    master_tips = leaf_set(master)
    missing_map = (high_all | low_all) - master_tips
    if missing_map:
        raise ValueError("Focal IDs absent from master tree: " + ", ".join(sorted(missing_map)))

    rows = []
    for protein_path in sorted(args.protein_dir.glob("*.faa")):
        family = protein_path.stem
        ids = read_fasta_ids(protein_path)
        if len(ids) != len(set(ids)):
            raise ValueError(f"Duplicate FASTA IDs in {protein_path}")
        absent = set(ids) - master_tips
        if absent:
            raise ValueError(f"Family {family} IDs absent from master tree: {sorted(absent)}")
        family_tree = prune(master, set(ids))
        if family_tree is None:
            raise ValueError(f"No retained tree for {family}")
        high_targets = set(ids) & high_all
        low_targets = set(ids) & low_all

        high_stem_tree = clone(family_tree)
        high_stem_ok, high_desc, high_reason = mark_stem_foreground(high_stem_tree, high_targets)
        low_stem_tree = clone(family_tree)
        low_stem_ok, low_desc, low_reason = mark_stem_foreground(low_stem_tree, low_targets)
        high_clade_tree = clone(family_tree)
        high_clade_ok, _, high_clade_reason = mark_crown_foreground(high_clade_tree, high_targets)
        low_clade_tree = clone(family_tree)
        low_clade_ok, _, low_clade_reason = mark_crown_foreground(low_clade_tree, low_targets)

        (plain_dir / f"{family}.nwk").write_text(serialize(family_tree) + ";\n", encoding="utf-8")
        if high_stem_ok:
            (high_stem_dir / f"{family}.nwk").write_text(serialize(high_stem_tree) + ";\n", encoding="utf-8")
        if low_stem_ok:
            (low_stem_dir / f"{family}.nwk").write_text(serialize(low_stem_tree) + ";\n", encoding="utf-8")
        if high_clade_ok:
            (high_clade_dir / f"{family}.nwk").write_text(serialize(high_clade_tree) + ";\n", encoding="utf-8")
        if low_clade_ok:
            (low_clade_dir / f"{family}.nwk").write_text(serialize(low_clade_tree) + ";\n", encoding="utf-8")
        rows.append({
            "split_family": family,
            "n_tree_tips": len(ids),
            "n_HighT": len(high_targets),
            "n_LowT": len(low_targets),
            "HighT_monophyletic": int(high_desc == high_targets),
            "LowT_monophyletic": int(low_desc == low_targets),
            "HighT_stem_tree_ready": int(high_stem_ok),
            "LowT_stem_tree_ready": int(low_stem_ok),
            "HighT_clade_tree_ready": int(high_clade_ok),
            "LowT_clade_tree_ready": int(low_clade_ok),
            "HighT_stem_status": high_reason,
            "LowT_stem_status": low_reason,
            "HighT_clade_status": high_clade_reason,
            "LowT_clade_status": low_clade_reason,
        })

    fields = list(rows[0]) if rows else []
    with (args.outdir / "tree_qc.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Families processed: {len(rows)}")
    print(f"HighT stem/clade trees ready: {sum(row['HighT_stem_tree_ready'] for row in rows)}/{sum(row['HighT_clade_tree_ready'] for row in rows)}")
    print(f"LowT stem/clade trees ready: {sum(row['LowT_stem_tree_ready'] for row in rows)}/{sum(row['LowT_clade_tree_ready'] for row in rows)}")


if __name__ == "__main__":
    main()
