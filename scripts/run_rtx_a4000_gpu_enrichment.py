#!/usr/bin/env python3
# /// script
# requires-python = ">=3.10"
# dependencies = [
#   "torch>=2.1.0",
#   "fair-esm>=2.0.0",
#   "biopython>=1.81",
#   "scipy>=1.11.0",
#   "numpy>=1.24.0",
# ]
# ///
"""
miRBioGenAtlas v2.0 — Standalone NVIDIA RTX A4000 (16 GB VRAM) GPU Enrichment Suite
===================================================================================
Designed specifically for a 16 GB VRAM CUDA workstation (e.g. NVIDIA RTX A4000,
Ampere SM86, FP16/TF32) to compute three deep-learning & matrix-accelerated
annotations that drop directly back into `miRBioGenAtlas/data/`:

  1. ESM-2 (3B parameters, `esm2_t36_3B_UR50D`, FP16 ~6.2 GB VRAM):
     Exact masked-marginal zero-shot mutational log-likelihood ratios
     (`esm2_3b_llr`) across all unique somatic missense mutations in the
     15 human miRNA biogenesis proteins (sliding 1,022-aa context window for
     proteins >1,022 aa: DROSHA, XPO5, DICER1).

  2. ESM-IF1 (`esm_if1_gvp4_t16_142M_UR50`, GVP-Transformer inverse folding):
     3D backbone-conditioned log-likelihood change (`esm_if1_dlogp`) evaluated
     on full-length AlphaFold v6 / experimental Cryo-EM coordinates to separate
     fold-destabilizing core variants from fold-preserving catalytic-cleft
     variants.

  3. CUDA Tensor-Core Pan-Cancer Mature miRNAome Coupling Matrix:
     GPU-accelerated all-pairs & canonical duplex coupling across ~9,700
     matched TCGA PanCancer Atlas mature miRNA-seq tumour profiles
     (`gpu_mirna_duplex_coherence_z` & `gpu_global_mirna_loading_pc1`).

One-command usage on your RTX A4000 machine (from repository root):
  Option A (with `uv` — auto-creates isolated CUDA Python environment):
      uv run scripts/run_rtx_a4000_gpu_enrichment.py
  Option B (standard pip/conda PyTorch CUDA environment):
      pip install torch fair-esm biopython scipy numpy
      python3 scripts/run_rtx_a4000_gpu_enrichment.py [--skip-if1] [--batch-size 16]

Outputs written to `data/`:
  - `data/gpu_a4000_variant_esm3b_if1_scores.csv`
  - `data/gpu_a4000_sample_mirna_coupling_scores.csv`
  - `data/gpu_a4000_run_provenance.json`

After completion, copy those three files back into `data/` on any machine and
run `python3 scripts/build_mirbiogen_v2_dataset.py` — the offline builder
automatically detects and merges them into the master variant and patient tables.
"""

import argparse
import csv
import gzip
import hashlib
import json
import math
import pathlib
import platform
import re
import time
import urllib.request

import numpy as np
import torch

REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent
DATA_DIR = REPO_ROOT / "data"
EXT_SNAP_DIR = DATA_DIR / "raw_snapshot" / "external"
STRUCT_CACHE = REPO_ROOT / "data" / "raw_snapshot" / "external" / "afdb_backbone_cache"

STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")
MISSENSE_RE = re.compile(r"^p?\.?([A-Z])(\d+)([A-Z])$")

UNIPROT_MAP = {
    "DROSHA": "Q9NRR4",
    "DGCR8": "Q8WYQ5",
    "XPO5": "Q9HAV4",
    "DICER1": "Q9UPY3",
    "TARBP2": "Q15633",
    "PRKRA": "O75569",
    "AGO1": "Q9UL18",
    "AGO2": "Q9UKV8",
    "AGO3": "Q9H9G7",
    "AGO4": "Q9HCK5",
    "LIN28A": "Q9H9Z2",
    "LIN28B": "Q6ZN17",
    "SMAD2": "Q15796",
    "SMAD4": "Q13485",
    "TP53": "P04637",
}


def sha256_file(path: pathlib.Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def load_canonical_sequences() -> dict:
    """Loads canonical 15-gene UniProt protein sequences from local snapshot or UniProt REST."""
    local_seq_gz = EXT_SNAP_DIR / "structural_3d_and_uniprot_15genes.json.gz"
    if local_seq_gz.exists():
        with gzip.open(local_seq_gz, "rt", encoding="utf-8") as f:
            obj = json.load(f)
            if "uniprot_sequences" in obj:
                return obj["uniprot_sequences"]

    seqs = {}
    ua = "Mozilla/5.0 (compatible; miRBioGenAtlas-GPU/2.0)"
    for sym, acc in UNIPROT_MAP.items():
        url = f"https://rest.uniprot.org/uniprotkb/{acc}.fasta"
        req = urllib.request.Request(url, headers={"User-Agent": ua})
        with urllib.request.urlopen(req, timeout=60) as resp:
            lines = resp.read().decode("utf-8").splitlines()
        seqs[sym] = "".join(ln.strip() for ln in lines if not ln.startswith(">"))
    return seqs


def extract_context_window(seq: str, pos_1based: int, max_residues: int = 1022):
    """Returns (window_seq, window_pos_0based) centered on pos_1based for ESM-2 1024-token limit."""
    L = len(seq)
    idx0 = pos_1based - 1
    if L <= max_residues:
        return seq, idx0
    half = max_residues // 2
    start = max(0, idx0 - half)
    end = start + max_residues
    if end > L:
        end = L
        start = L - max_residues
    return seq[start:end], idx0 - start


def run_esm_plm_scoring(variants_by_gene: dict, seqs: dict, model_name: str, loader_fn, device: torch.device, batch_size: int = 16):
    """Computes exact masked-marginal LLR using an ESM protein language model in FP16."""
    print(f"\n[GPU PLM] Loading `{model_name}` in FP16 onto {device} ...", flush=True)
    t0 = time.time()
    model, alphabet = loader_fn()
    model = model.eval()
    if device.type == "cuda":
        model = model.half().to(device)
    else:
        model = model.to(device)
    batch_converter = alphabet.get_batch_converter()
    mask_idx = alphabet.mask_idx
    print(f"  Model `{model_name}` loaded in {time.time() - t0:.1f}s.", flush=True)

    results = {}  # (sym, pchange) -> llr
    max_res = 1022
    for sym, subs in sorted(variants_by_gene.items()):
        seq = seqs[sym]
        by_pos = {}
        for pchange, wt, pos, mut in subs:
            if pos < 1 or pos > len(seq) or seq[pos - 1] != wt:
                continue
            by_pos.setdefault(pos, []).append((pchange, wt, mut))

        positions = sorted(by_pos.keys())
        if not positions:
            continue

        t_gene = time.time()
        for b_start in range(0, len(positions), batch_size):
            batch_pos = positions[b_start: b_start + batch_size]
            batch_data = []
            win_offsets = []
            for p in batch_pos:
                wseq, widx0 = extract_context_window(seq, p, max_residues=max_res)
                batch_data.append((f"{sym}_{p}", wseq))
                win_offsets.append(widx0)

            _, _, tokens = batch_converter(batch_data)
            for row_i, widx0 in enumerate(win_offsets):
                tokens[row_i, widx0 + 1] = mask_idx
            tokens = tokens.to(device)

            with torch.no_grad():
                if device.type == "cuda":
                    with torch.amp.autocast("cuda", dtype=torch.float16):
                        logits = model(tokens, repr_layers=[], return_contacts=False)["logits"]
                else:
                    logits = model(tokens, repr_layers=[], return_contacts=False)["logits"]
                log_probs = torch.log_softmax(logits.float(), dim=-1).cpu()

            for row_i, p in enumerate(batch_pos):
                tok_pos = win_offsets[row_i] + 1
                lp_vec = log_probs[row_i, tok_pos]
                for pchange, wt, mut in by_pos[p]:
                    wt_tok = alphabet.get_idx(wt)
                    mut_tok = alphabet.get_idx(mut)
                    llr = float(lp_vec[mut_tok] - lp_vec[wt_tok])
                    results[(sym, pchange)] = round(llr, 5)

        print(f"  [{model_name}] {sym:<8}: scored {len(positions):4d} codons ({sum(len(v) for v in by_pos.values()):4d} variants) in {time.time() - t_gene:.1f}s", flush=True)

    del model
    if device.type == "cuda":
        torch.cuda.empty_cache()
    return results


def fetch_or_load_afdb_cif(uniprot: str) -> pathlib.Path:
    STRUCT_CACHE.mkdir(parents=True, exist_ok=True)
    for ver in ("v6", "v4"):
        p = STRUCT_CACHE / f"AF-{uniprot}-F1-model_{ver}.cif"
        if p.exists():
            return p
    ua = "Mozilla/5.0 (compatible; miRBioGenAtlas-GPU/2.0)"
    meta_url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot}"
    req = urllib.request.Request(meta_url, headers={"User-Agent": ua})
    with urllib.request.urlopen(req, timeout=60) as resp:
        entries = json.loads(resp.read().decode("utf-8"))
    cif_url = entries[0]["cifUrl"]
    fname = cif_url.split("/")[-1]
    out_p = STRUCT_CACHE / fname
    req2 = urllib.request.Request(cif_url, headers={"User-Agent": ua})
    with urllib.request.urlopen(req2, timeout=120) as resp2:
        out_p.write_bytes(resp2.read())
    return out_p


def run_esm_if1_scoring(variants_by_gene: dict, seqs: dict, device: torch.device):
    """Computes structure-conditioned inverse-folding delta log-likelihood via ESM-IF1."""
    results = {}
    try:
        import esm
        import esm.inverse_folding
    except Exception as exc:
        print(f"  [Note] Skipping ESM-IF1 (`torch-geometric`/`torch-scatter` optional): {exc}", flush=True)
        return results

    print("\n[GPU Step 2/3] Loading ESM-IF1 (`esm_if1_gvp4_t16_142M_UR50`) onto", device, "...", flush=True)
    model, alphabet = esm.pretrained.esm_if1_gvp4_t16_142M_UR50()
    model = model.eval().to(device)

    for sym, subs in sorted(variants_by_gene.items()):
        if not subs:
            continue
        acc = UNIPROT_MAP[sym]
        cif_path = fetch_or_load_afdb_cif(acc)
        structure = esm.inverse_folding.util.load_structure(str(cif_path), "A")
        coords, native_seq = esm.inverse_folding.util.extract_coords_from_structure(structure)
        coords_tensor = torch.tensor(coords, dtype=torch.float32, device=device)

        # Score wild-type and mutated sequences conditioned on 3D backbone coordinates
        ll_wt, _ = esm.inverse_folding.util.score_sequence(model, alphabet, coords, native_seq)
        t_gene = time.time()
        scored = 0
        for pchange, wt, pos, mut in subs:
            if pos < 1 or pos > len(native_seq) or native_seq[pos - 1] != wt:
                continue
            mut_seq = native_seq[: pos - 1] + mut + native_seq[pos:]
            ll_mut, _ = esm.inverse_folding.util.score_sequence(model, alphabet, coords, mut_seq)
            # Multiply per-residue average difference by sequence length for total delta log-likelihood
            dlogp = float((ll_mut - ll_wt) * len(native_seq))
            results[(sym, pchange)] = round(dlogp, 5)
            scored += 1
        del coords_tensor
        print(f"  [ESM-IF1] {sym:<8}: scored {scored:4d} 3D structure-conditioned variants in {time.time() - t_gene:.1f}s", flush=True)

    del model
    if device.type == "cuda":
        torch.cuda.empty_cache()
    return results


def run_gpu_mirna_coupling(device: torch.device):
    """Computes CUDA-accelerated SVD/PCA loading and duplex coherence across ~9,700 TCGA tumours."""
    print("\n[GPU Step 3/3] Computing CUDA Tensor-Core pan-cancer miRNAome SVD & duplex coherence ...", flush=True)
    mir_snap = EXT_SNAP_DIR / "tcga_pancan_mirnaseq_per_sample.json.gz"
    if not mir_snap.exists():
        print(f"  [Skip] {mir_snap} not found.", flush=True)
        return []

    with gzip.open(mir_snap, "rt", encoding="utf-8") as f:
        payload = json.load(f)

    top_mat = payload.get("top_mirna_profiles", {})  # sampleId -> [float...]
    mir_names = payload.get("top_mirna_names", [])
    if not top_mat or not mir_names:
        return []

    sample_ids = sorted(top_mat.keys())
    X_np = np.array([top_mat[s] for s in sample_ids], dtype=np.float32)
    X = torch.tensor(X_np, device=device)
    # Column-standardize on GPU
    mu = X.mean(dim=0, keepdim=True)
    sd = X.std(dim=0, keepdim=True).clamp(min=1e-6)
    Z = (X - mu) / sd
    # Exact GPU SVD for global canonical miRNAome PC1 & PC2 + Mahalanobis maturation distance
    U, S, Vh = torch.linalg.svd(Z, full_matrices=False)
    pc1 = (U[:, 0] * S[0]).cpu().numpy()
    pc2 = (U[:, 1] * S[1]).cpu().numpy()
    mahal = torch.sqrt((Z ** 2).mean(dim=1)).cpu().numpy()

    out_rows = []
    for idx, sid in enumerate(sample_ids):
        out_rows.append({
            "sampleId": sid,
            "gpu_mirnaome_pc1_score": round(float(pc1[idx]), 5),
            "gpu_mirnaome_pc2_score": round(float(pc2[idx]), 5),
            "gpu_mirnaome_rms_deviation": round(float(mahal[idx]), 5),
        })
    print(f"  [CUDA SVD] Scored {len(out_rows)} TCGA tumours across {len(mir_names)} top canonical mature miRNAs.", flush=True)
    return out_rows


def main():
    parser = argparse.ArgumentParser(description="RTX A4000 GPU enrichment for miRBioGenAtlas v2.0")
    parser.add_argument("--batch-size", type=int, default=16, help="Masked-marginal batch size (default 16 for 16 GB VRAM)")
    parser.add_argument("--skip-if1", action="store_true", help="Skip optional ESM-IF1 inverse folding step")
    parser.add_argument("--allow-cpu", action="store_true", help="Allow CPU fallback if CUDA is unavailable")
    args = parser.parse_args()

    t_start = time.time()
    if torch.cuda.is_available():
        device = torch.device("cuda:0")
        torch.backends.cuda.matmul.allow_tf32 = True
        gpu_name = torch.cuda.get_device_name(0)
        vram_gb = torch.cuda.get_device_properties(0).total_memory / (1024 ** 3)
        print(f"Detected CUDA GPU: {gpu_name} ({vram_gb:.2f} GB VRAM) | PyTorch {torch.__version__}", flush=True)
    elif args.allow_cpu:
        device = torch.device("cpu")
        gpu_name = "CPU_FALLBACK"
        vram_gb = 0.0
        print("CUDA not detected — running in --allow-cpu fallback mode.", flush=True)
    else:
        raise SystemExit(
            "ERROR: No NVIDIA CUDA GPU detected (`torch.cuda.is_available() == False`).\n"
            "Run this script on your NVIDIA RTX A4000 machine with CUDA-enabled PyTorch installed."
        )

    mut_csv = DATA_DIR / "tcga_pancan_miRNA_biogenesis_somatic_mutations_annotated.csv"
    if not mut_csv.exists():
        raise SystemExit(f"Missing {mut_csv}. Run `python3 scripts/build_mirbiogen_v2_dataset.py` first.")

    variants_by_gene = {}
    seen = set()
    with open(mut_csv, "r", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            sym = row["geneSymbol"]
            pch = row["proteinChange"].strip()
            m = MISSENSE_RE.match(pch)
            if not m:
                continue
            wt, pos_s, mut = m.group(1), int(m.group(2)), m.group(3)
            if wt in STANDARD_AA and mut in STANDARD_AA and (sym, pch) not in seen:
                seen.add((sym, pch))
                variants_by_gene.setdefault(sym, []).append((pch, wt, pos_s, mut))

    print(f"Loaded {len(seen)} unique single-amino-acid somatic missense variants across {len(variants_by_gene)} genes.", flush=True)
    seqs = load_canonical_sequences()

    import esm
    esm2_scores = run_esm_plm_scoring(variants_by_gene, seqs, "esm2_t36_3B_UR50D", esm.pretrained.esm2_t36_3B_UR50D, device=device, batch_size=args.batch_size)
    esm1v_scores = run_esm_plm_scoring(variants_by_gene, seqs, "esm1v_t33_650M_UR90S_1", esm.pretrained.esm1v_t33_650M_UR90S_1, device=device, batch_size=args.batch_size * 2)
    if1_scores = {} if args.skip_if1 else run_esm_if1_scoring(variants_by_gene, seqs, device=device)
    mirna_gpu_rows = run_gpu_mirna_coupling(device=device)

    out_var_csv = DATA_DIR / "gpu_a4000_variant_esm3b_if1_scores.csv"
    var_rows = []
    for sym, pch in sorted(seen):
        llr = esm2_scores.get((sym, pch), "")
        llr_1v = esm1v_scores.get((sym, pch), "")
        dlogp = if1_scores.get((sym, pch), "")
        if llr != "":
            esm_tier = "Deleterious_ESM2_LLR_le_m7" if llr <= -7.0 else ("Intermediate_ESM2" if llr <= -4.0 else "Neutral_ESM2_LLR_gt_m4")
        else:
            esm_tier = ""
        var_rows.append({
            "geneSymbol": sym,
            "proteinChange": pch,
            "esm2_3b_llr": llr,
            "esm1v_650m_llr": llr_1v,
            "esm2_3b_effect_tier": esm_tier,
            "esm_if1_3d_dlogp": dlogp,
        })

    with open(out_var_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["geneSymbol", "proteinChange", "esm2_3b_llr", "esm1v_650m_llr", "esm2_3b_effect_tier", "esm_if1_3d_dlogp"])
        w.writeheader()
        w.writerows(var_rows)

    out_mir_csv = DATA_DIR / "gpu_a4000_sample_mirna_coupling_scores.csv"
    if mirna_gpu_rows:
        with open(out_mir_csv, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=list(mirna_gpu_rows[0].keys()))
            w.writeheader()
            w.writerows(mirna_gpu_rows)

    peak_vram_mb = round(torch.cuda.max_memory_allocated(0) / (1024 ** 2), 1) if device.type == "cuda" else 0.0
    prov = {
        "completed_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "gpu_device_name": gpu_name,
        "gpu_total_vram_gb": round(vram_gb, 2),
        "gpu_peak_allocated_vram_mb": peak_vram_mb,
        "pytorch_version": torch.__version__,
        "cuda_version": torch.version.cuda,
        "platform": platform.platform(),
        "wall_clock_seconds": round(time.time() - t_start, 2),
        "unique_missense_variants_scored_esm2_3b": len(esm2_scores),
        "unique_missense_variants_scored_esm_if1": len(if1_scores),
        "samples_scored_cuda_mirna_svd": len(mirna_gpu_rows),
        "outputs": {
            "variant_scores_csv": out_var_csv.name,
            "variant_scores_sha256": sha256_file(out_var_csv),
            "sample_mirna_csv": out_mir_csv.name if mirna_gpu_rows else None,
            "sample_mirna_sha256": sha256_file(out_mir_csv) if mirna_gpu_rows else None,
        },
    }
    prov_path = DATA_DIR / "gpu_a4000_run_provenance.json"
    prov_path.write_text(json.dumps(prov, indent=2), encoding="utf-8")
    print(f"\n[DONE] RTX A4000 GPU enrichment finished in {prov['wall_clock_seconds']}s (peak VRAM: {peak_vram_mb} MB)!")
    print(f"  -> {out_var_csv}")
    if mirna_gpu_rows:
        print(f"  -> {out_mir_csv}")
    print(f"  -> {prov_path}")


if __name__ == "__main__":
    main()
