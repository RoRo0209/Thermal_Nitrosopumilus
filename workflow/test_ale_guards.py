#!/usr/bin/env python3
"""Synthetic regression tests only; these are not study observations."""
import argparse
import csv
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
RSCRIPT = os.environ.get("RSCRIPT", "Rscript")


class AleGuards(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix="synthetic_ale_")
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.uml = self.root / "raw"
        self.uml.mkdir()
        self.out = self.root / "verified"
        self.tree = self.root / "species.nwk"
        self.tree.write_text("((D:1,C:1):1,(B:1,A:1):1);\n")
        self.groups = self.root / "groups.tsv"
        self.groups.write_text("genome_id\tgroup\nA\tHighT\nB\tHighT\nC\tLowT\nD\tLowT\n")
        self.raw = "#ALEml_undated synthetic test only\nS: ((A:1,B:1)40:1,(C:1,D:1)50:1)60;\n"
        for tip, aid in zip("ABCD", range(10, 14)):
            self.raw += f"S_terminal_branch {tip}({aid}) 0 0 0 0\n"
        for branch in (40, 50, 60):
            self.raw += f"S_internal_branch {branch} 1.25 2.5 3.75 0\n"
        for name in ("synthetic_one", "synthetic_two"):
            (self.uml / (name + ".ale.uml_rec")).write_text(self.raw)

    def run_verify(self, success=True):
        result = subprocess.run([RSCRIPT, "--vanilla", str(ROOT / "analysis/ALE_2/verify_ale_mapping.R"),
                                 str(self.tree), str(self.uml / "synthetic_one.ale.uml_rec"),
                                 str(self.groups), str(self.out)], capture_output=True, text=True)
        self.assertEqual(result.returncode == 0, success, result.stdout + result.stderr)

    def extract(self, success=True):
        result = subprocess.run([sys.executable, str(ROOT / "analysis/ALE_2/extract_ale_events.py"),
                                 "--uml-dir", str(self.uml), "--verified-dir", str(self.out),
                                 "--expected-families", "2"], capture_output=True, text=True)
        self.assertEqual(result.returncode == 0, success, result.stdout + result.stderr)

    def test_reordered_tips_keep_actual_ale_ids_and_exact_sums(self):
        self.run_verify()
        with (self.out / "crowns.tsv").open() as h:
            rows = list(csv.DictReader(h, delimiter="\t"))
        self.assertEqual({r["group"]: r["ALE_ID"] for r in rows}, {"HighT": "40", "LowT": "50"})
        self.assertTrue(all(r["R_node"] != r["ALE_ID"] for r in rows))
        self.extract()
        with (self.out / "DTLO_table_rgy.txt").open() as h:
            values = {r["Branch"]: r for r in csv.DictReader(h, delimiter="\t")}
        self.assertEqual(values["40"]["D"], "2.50")
        self.assertEqual(len(values), 7)

    def test_different_clades_rejected(self):
        self.tree.write_text("((A:1,C:1):1,(B:1,D:1):1);\n")
        self.run_verify(False)

    def test_nonclade_group_rejected(self):
        self.groups.write_text("genome_id\tgroup\nA\tHighT\nC\tHighT\nB\tLowT\nD\tLowT\n")
        self.run_verify(False)

    def test_different_run_numbering_rejected(self):
        self.run_verify()
        p = self.uml / "synthetic_two.ale.uml_rec"
        p.write_text(self.raw.replace(")40:", ")41:").replace("branch 40 ", "branch 41 "))
        self.extract(False)

    def test_missing_branch_in_all_families_rejected(self):
        self.run_verify()
        for p in self.uml.glob("*.uml_rec"):
            p.write_text(self.raw.replace("S_internal_branch 60 1.25 2.5 3.75 0\n", ""))
        self.extract(False)

    def test_duplicate_record_rejected(self):
        self.run_verify()
        (self.uml / "synthetic_two.ale.uml_rec").write_text(self.raw + "S_internal_branch 40 1 1 1 1\n")
        self.extract(False)

    def test_plot_rejects_changed_event_table(self):
        self.run_verify()
        self.extract()
        with (self.out / "branch40_gene_DTLO.tsv").open("a") as h:
            h.write("tampered\t9\t9\t9\t9\n")
        r = subprocess.run([RSCRIPT, "--vanilla", str(ROOT / "analysis/ALE_2/plot_functional_events.R"),
                            str(self.out), "absent.xlsx", "absent.json", str(self.root / "figures")],
                           capture_output=True, text=True)
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("Verified event files changed", r.stdout + r.stderr)

    def test_plot_rejects_other_numbered_source(self):
        self.run_verify()
        self.extract()
        p = self.root / "other.ale.uml_rec"
        p.write_text(self.raw.replace(")40:", ")41:"))
        r = subprocess.run([RSCRIPT, "--vanilla", str(ROOT / "analysis/ALE_2/plot_branch_events.R"),
                            str(self.tree), str(p), "absent.xlsx", str(self.out), str(self.root / "figures")],
                           capture_output=True, text=True)
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("Plot ALE source differs", r.stdout + r.stderr)


if __name__ == "__main__":
    unittest.main()
