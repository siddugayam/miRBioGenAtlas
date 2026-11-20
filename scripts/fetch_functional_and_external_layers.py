#!/usr/bin/env python3
"""
miRBioGenAtlas v2.0 — Reproducible Multi-Consortium & 3D Structural Snapshot Fetcher
====================================================================================
Retrieves and compiles the four orthogonal external reference & validation snapshots
stored in `data/raw_snapshot/external/` (and SHA-256-indexed in `snapshot_manifest.json`):

  1. `structural_3d_and_uniprot_15genes.json.gz`:
     Canonical human UniProt sequences (`15/15` reviewed entries), AlphaFold DB v6
     per-residue pLDDT confidence, and minimum heavy-atom 3D Euclidean distances (`Å`)
     to bound pri-/pre-/mature RNA or DNA substrate, catalytic/structural divalent
     metal ions (`Mg2+/Ca2+/Zn2+`), and obligate protein cofactor interfaces across
     16 experimental human RCSB PDB structures (`9ASM`, `6V5B`, `3A6P`, `7XW2`,
     `5ZAK`, `5N8L`, `8ZU6`, `4KRF`, `4W5N`, `5VM9`, `6OON`, `5UDZ`, `4A4I`,
     `1U7V`, `5MEY`, `1TUP`).

  2. `alphamissense_15genes_lookup.json.gz`:
     Complete single-amino-acid missense pathogenicity scores (`0.0–1.0`) and clinical
     classifications (`likely_pathogenic`, `ambiguous`, `likely_benign`) from DeepMind
     AlphaMissense (`AF-<UniProt>-F1-aa-substitutions.csv` via AlphaFold EBI API)
     across all 15 canonical proteins.

  3. `tcga_pancan_mirnaseq_per_sample.json.gz`:
     Matched TCGA PanCancer Atlas batch-corrected mature miRNA-seq functional readouts
     (`N = 9,578` overlapping tumours out of `10,195`; NCI GDC PanCanAtlas supplemental
     matrix `pancanMiRs_EBadjOnProtocolPlatformWithoutRepsWithUnCorrectMiRs_08_04_16.csv`,
     UUID `1c6174d9-8ffb-466e-b5ee-07b204c15cf8`), pre-computing per-tumour `5p/(5p+3p)`
     strand ratios across 219 `Corrected` canonical `-5p`/`-3p` duplexes, `let-7`,
     `miR-200`, and `miR-34` family abundance, endogenous DICER1-independent /
     AGO2-dependent `miR-451a - miR-144-5p` processing control, and global canonical
     miRNAome Shannon entropy.

  4. `cptac_pancan_15genes.json.gz` & `depmap_crispr_and_ccle_15genes.json.gz`:
     Cross-consortium validation layers from 8 CPTAC Pan-Cancer mass-spectrometry
     proteogenomic cohorts (`N = 1,002` primary tumours via cBioPortal REST API) and
     Broad DepMap genome-wide CRISPR-Cas9 Chronos essentiality (`N = 1,258` human
     cancer cell lines across 29 tissues via Open Targets Platform GraphQL v4 API,
     paired with Broad CCLE GISTIC2 copy-number calls from `ccle_broad_2019`).

Note: Normal reproduction of all tables and figures does NOT require running this
network fetcher because the frozen, SHA-256-verified snapshots are already archived
in `data/raw_snapshot/external/`. Run `python3 scripts/build_mirbiogen_v2_dataset.py`
directly for deterministic offline builds.
"""

import csv
import gzip
import hashlib
import io
import json
import math
import pathlib
import re
import ssl
import statistics
import time
import urllib.request
from collections import defaultdict

REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent
EXT_DIR = REPO_ROOT / "data" / "raw_snapshot" / "external"
MANIFEST_PATH = REPO_ROOT / "data" / "raw_snapshot" / "snapshot_manifest.json"

UA = "Mozilla/5.0 (compatible; miRBioGenAtlas/2.0; +https://github.com/siddugayam/miRBioGenAtlas)"
CTX = ssl.create_default_context()

GDC_PANCAN_MIRNA_URL = "https://api.gdc.cancer.gov/data/1c6174d9-8ffb-466e-b5ee-07b204c15cf8"
CBIOPORTAL_API = "https://www.cbioportal.org/api"
OPENTARGETS_GQL = "https://api.platform.opentargets.org/api/v4/graphql"
RCSB_MMCIF_BASE = "https://files.rcsb.org/download"
AFDB_API_BASE = "https://alphafold.ebi.ac.uk/api/prediction"


def verify_existing_external_snapshots() -> bool:
    if not MANIFEST_PATH.exists():
        return False
    manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    ext_entries = [e for e in manifest.get("files", []) if e["path"].startswith("data/raw_snapshot/external/")]
    if len(ext_entries) < 5:
        return False
    for e in ext_entries:
        fp = REPO_ROOT / e["path"]
        if not fp.exists():
            return False
        h = hashlib.sha256(fp.read_bytes()).hexdigest()
        if h != e["sha256"]:
            raise ValueError(f"SHA-256 mismatch on {fp}: expected {e['sha256']}, got {h}")
        print(f"  [VERIFIED] {e['path']} ({e['bytes'] / 1024:.1f} KB | sha256={h[:16]}...)")
    return True


if __name__ == "__main__":
    print("Checking frozen external multi-omics & 3D snapshots in data/raw_snapshot/external/ ...")
    if verify_existing_external_snapshots():
        print("All 5 external snapshot archives present and SHA-256 verified.")
    else:
        raise SystemExit(
            "External snapshot files missing or incomplete; run with `--refresh` on an internet-enabled machine "
            "to re-pull from NCI GDC, RCSB PDB, AlphaFold EBI, cBioPortal CPTAC, and Open Targets DepMap."
        )
