#!/usr/bin/env python3
"""Verify SHA-256 checksums of external reference snapshots in data/raw_snapshot/external/."""

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
        print(f"  [VERIFIED] {e['path']} ({e.get('size_bytes', e.get('bytes', 0)) / 1024:.1f} KB | sha256={h[:16]}...)")
    return True


if __name__ == "__main__":
    print("Checking frozen external multi-omics & 3D snapshots in data/raw_snapshot/external/ ...")
    if verify_existing_external_snapshots():
        print("All external snapshot archives present and SHA-256 verified.")
    else:
        raise SystemExit(
            "External snapshot files missing or incomplete; run with `--refresh` on an internet-enabled machine "
            "to re-pull from NCI GDC, RCSB PDB, AlphaFold EBI, cBioPortal CPTAC, and Open Targets DepMap."
        )
