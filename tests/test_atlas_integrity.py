#!/usr/bin/env python3
"""
Automated Reproducibility, Schema, Portability & Mathematical Integrity Suite
for miRBioGenAtlas v2.0.0 (runs with zero third-party dependencies via `unittest`).

Usage:
    python3 -m unittest discover -s tests -v
"""

import csv
import gzip
import hashlib
import json
import pathlib
import subprocess
import sys
import unittest

REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent
DATA_DIR = REPO_ROOT / "data"


class TestMiRBioGenAtlasIntegrity(unittest.TestCase):
    def test_01_snapshot_manifest_sha256_all_45_archives(self):
        manifest_path = DATA_DIR / "raw_snapshot" / "snapshot_manifest.json"
        self.assertTrue(manifest_path.exists(), "Missing snapshot_manifest.json")
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        self.assertEqual(len(manifest["files"]), 45)
        for item in manifest["files"]:
            fp = REPO_ROOT / item["path"]
            self.assertTrue(fp.exists(), f"Missing frozen snapshot: {item['path']}")
            digest = hashlib.sha256(fp.read_bytes()).hexdigest()
            self.assertEqual(digest, item["sha256"], f"SHA-256 mismatch for {item['path']}")

    def test_02_core_table_dimensions_and_union_identity(self):
        with gzip.open(DATA_DIR / "tcga_pancan_miRNA_biogenesis_sample_level_alterations.csv.gz", "rt", encoding="utf-8") as f:
            n_sl = sum(1 for _ in csv.DictReader(f))
        self.assertEqual(n_sl, 152925)

        with gzip.open(DATA_DIR / "miRNA_biogenesis_saturation_mutagenesis_213940_atlas.csv.gz", "rt", encoding="utf-8") as f:
            n_sat = sum(1 for _ in csv.DictReader(f))
        self.assertEqual(n_sat, 213940)

        with open(DATA_DIR / "tcga_pancan_miRNA_biogenesis_patient_level_burden.csv", "r", encoding="utf-8") as f:
            pb_rows = list(csv.DictReader(f))
        self.assertEqual(len(pb_rows), 10195)
        self.assertIn("osMonths", pb_rows[0])
        self.assertIn("hrdScore", pb_rows[0])

        with open(DATA_DIR / "tcga_pancan_miRNA_biogenesis_alteration_freqs_long_FINAL.csv", "r", encoding="utf-8") as f:
            long_rows = list(csv.DictReader(f))
        self.assertEqual(len(long_rows), 480)
        for r in long_rows:
            self.assertEqual(
                int(r["any_alter_n"]),
                int(r["mut_only_n"]) + int(r["cna_only_n"]) + int(r["mut_and_cna_n"]),
            )

    def test_03_cbioportal_20_spotchecks_exact_zero_discrepancy(self):
        with open(DATA_DIR / "validation_20_cbioportal_spotchecks.csv", "r", encoding="utf-8") as f:
            rows = list(csv.DictReader(f))
        self.assertEqual(len(rows), 20)
        for r in rows:
            self.assertAlmostEqual(float(r["abs_freq_discrepancy"]), 0.0, places=4)
            self.assertIn(str(r["union_identity_holds"]).strip().lower(), ("1", "true"))

    def test_04_clinvar_roc_auc_and_multivariable_cox_benchmarks(self):
        with open(DATA_DIR / "validation_clinvar_roc_auc_benchmarks.csv", "r", encoding="utf-8") as f:
            roc_rows = list(csv.DictReader(f))
        am_all = [
            r for r in roc_rows
            if r["benchmark_subset"] == "ClinVar_All15_Missense_AllStars" and r["predictor"] == "AlphaMissense"
        ]
        self.assertEqual(len(am_all), 1)
        self.assertGreater(float(am_all[0]["roc_auc"]), 0.97)

        with open(DATA_DIR / "validation_clinical_survival_benchmarks.csv", "r", encoding="utf-8") as f:
            surv_rows = list(csv.DictReader(f))
        self.assertGreaterEqual(len(surv_rows), 12)
        self.assertIn("multivar_cox_hr", surv_rows[0])
        lin28_row = [r for r in surv_rows if r["test_id"] == "S6_LIN28AB_Oncofetal_Reactivation_Stratified_OS"][0]
        self.assertGreater(float(lin28_row["multivar_cox_hr"]), 1.2)

    def test_05_no_hardcoded_user_paths_in_portable_code(self):
        check_files = list((REPO_ROOT / "scripts").glob("*.py")) + [
            REPO_ROOT / "notebooks" / "notebook.ipynb",
            REPO_ROOT / "miRBioGenAtlas_figures.R",
        ]
        for fp in check_files:
            txt = fp.read_text(encoding="utf-8")
            self.assertNotIn("/Users/", txt, f"Hard-coded macOS path in {fp.name}")
            self.assertNotIn("/home/", txt, f"Hard-coded Linux path in {fp.name}")

    def test_06_web_explorer_and_15_saturation_chunks_present(self):
        idx = REPO_ROOT / "docs" / "index.html"
        self.assertTrue(idx.exists())
        self.assertGreater(idx.stat().st_size, 500_000)
        sat_js = list((REPO_ROOT / "docs" / "sat").glob("*.js"))
        self.assertEqual(len(sat_js), 15)

    def test_07_cli_variant_query_utility(self):
        res = subprocess.run(
            [sys.executable, str(REPO_ROOT / "scripts" / "query_variant.py"), "DICER1", "E1705K", "--json"],
            capture_output=True,
            text=True,
            check=True,
        )
        payload = json.loads(res.stdout)
        self.assertEqual(payload["geneSymbol"], "DICER1")
        self.assertEqual(payload["variant"], "p.E1705K")
        self.assertGreater(payload["tcgaPanCanSomaticCount_N10195"], 0)
        self.assertGreater(payload["alphamissenseScore"], 0.9)


if __name__ == "__main__":
    unittest.main()
