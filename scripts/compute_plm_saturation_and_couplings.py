#!/usr/bin/env python3
"""
Optional deep-learning and matrix regeneration script for miRBioGenAtlas.
Recomputes ESM-IF1 3D inverse-folding scores, ESM-2 (3B) and ESM-1v saturation
mutagenesis log-likelihood ratios across 213,940 substitutions, and pan-cancer
mature miRNA x mRNA correlations across dual-profiled TCGA tumours.

Normal dataset builds do not require running this script because pre-computed,
SHA-256-verified snapshots are archived in `data/raw_snapshot/`.
"""

import argparse
import csv
import gzip
import hashlib
import json
import math
import pathlib
import re
import ssl
import time
import urllib.request
from collections import Counter, defaultdict

import numpy as np
import torch

if hasattr(torch.serialization, "add_safe_globals"):
    torch.serialization.add_safe_globals([argparse.Namespace])

WORK_ROOT = pathlib.Path(__file__).resolve().parent.parent
DATA_DIR = WORK_ROOT / "data"
EXT_DIR = DATA_DIR / "raw_snapshot" / "external"
AFDB_CACHE = EXT_DIR / "afdb_backbone_cache"

STANDARD_AA = "ACDEFGHIKLMNPQRSTVWY"
STANDARD_AA_SET = set(STANDARD_AA)
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
GENE_ORDER = list(UNIPROT_MAP.keys())

CTX = ssl.create_default_context()
UA = "Mozilla/5.0 (compatible; miRBioGenAtlas/2.0; +https://github.com/siddugayam/miRBioGenAtlas)"


def select_device(requested: str) -> torch.device:
    if requested != "auto":
        return torch.device(requested)
    if torch.cuda.is_available():
        return torch.device("cuda")
    if hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
        return torch.device("mps")
    return torch.device("cpu")


def sha256_file(path: pathlib.Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_cif_backbone_coords(cif_path: pathlib.Path, expected_len: int) -> np.ndarray:
    coords = np.full((expected_len, 3, 3), np.nan, dtype=np.float32)
    atom_idx = {"N": 0, "CA": 1, "C": 2}
    cols = []
    in_atom = False
    with open(cif_path, "r", encoding="utf-8", errors="replace") as f:
        for ln in f:
            if ln.startswith("_atom_site."):
                in_atom = True
                cols.append(ln.strip().split()[0])
                continue
            if in_atom:
                if not ln.startswith("ATOM "):
                    if cols and ln.strip() in ("#", "loop_"):
                        break
                    continue
                parts = ln.split()
                if len(parts) < len(cols):
                    continue
                rec = dict(zip(cols, parts))
                aname = rec.get("_atom_site.label_atom_id", "")
                if aname in atom_idx:
                    pos = int(rec["_atom_site.label_seq_id"]) - 1
                    if 0 <= pos < expected_len:
                        coords[pos, atom_idx[aname], 0] = float(rec["_atom_site.Cartn_x"])
                        coords[pos, atom_idx[aname], 1] = float(rec["_atom_site.Cartn_y"])
                        coords[pos, atom_idx[aname], 2] = float(rec["_atom_site.Cartn_z"])
    return coords


def extract_context_window(seq: str, pos_1based: int, max_residues: int = 1022):
    seq_len = len(seq)
    idx0 = pos_1based - 1
    if seq_len <= max_residues:
        return seq, idx0
    half = max_residues // 2
    start = max(0, idx0 - half)
    end = start + max_residues
    if end > seq_len:
        end = seq_len
        start = seq_len - max_residues
    return seq[start:end], idx0 - start


def run_esm_if1(seqs: dict, tcga_variants_by_gene: dict, device: torch.device, cache_dir: pathlib.Path):
    ckpt_npz = cache_dir / "ckpt_if1.npz"
    ckpt_json = cache_dir / "ckpt_if1_tcga_full.json"
    if ckpt_npz.exists() and ckpt_json.exists():
        z = np.load(ckpt_npz)
        if1_saturation = {sym: z[sym] for sym in GENE_ORDER}
        raw_j = json.loads(ckpt_json.read_text(encoding="utf-8"))
        if1_tcga_full = {tuple(k.split("|", 1)): float(v) for k, v in raw_j.items()}
        return if1_saturation, if1_tcga_full

    import esm
    import esm.inverse_folding

    model, alphabet = esm.pretrained.esm_if1_gvp4_t16_142M_UR50()
    model = model.eval().to(device)
    coord_converter = esm.inverse_folding.util.CoordBatchConverter(alphabet)

    if1_saturation = {}
    if1_tcga_full = {}
    aa_to_alpha_idx = {aa: alphabet.get_idx(aa) for aa in STANDARD_AA}

    for sym in GENE_ORDER:
        seq = seqs[sym]
        seq_len = len(seq)
        acc = UNIPROT_MAP[sym]
        cif_candidates = sorted(AFDB_CACHE.glob(f"AF-{acc}-F1-model_*.cif"))
        if not cif_candidates:
            raise FileNotFoundError(f"Missing AlphaFold mmCIF for {sym} ({acc}) in {AFDB_CACHE}")
        coords = parse_cif_backbone_coords(cif_candidates[0], seq_len)

        coords_t, confidence_t, _, tokens_t, padding_mask_t = coord_converter([(coords, None, seq)], device=device)
        prev_output_tokens = tokens_t[:, :-1]
        valid_mask_1d = torch.isfinite(coords_t[:, 1: seq_len + 1, 0, 0]).float()

        with torch.no_grad():
            enc_out = model.encoder(coords_t, padding_mask_t, confidence_t)
            logits, _ = model.decoder(prev_output_tokens, encoder_out=enc_out)
            if logits.shape[0] == len(alphabet):
                lp_t = torch.log_softmax(logits[:, 0, :seq_len].float(), dim=0).transpose(0, 1)
            elif logits.shape[1] == len(alphabet):
                lp_t = torch.log_softmax(logits[0, :, :seq_len].float(), dim=0).transpose(0, 1)
            else:
                lp_t = torch.log_softmax(logits[0, :seq_len, :].float(), dim=-1)
            lp = lp_t.cpu().numpy()

        mat = np.zeros((seq_len, 20), dtype=np.float32)
        for pos0 in range(seq_len):
            wt_lp = lp[pos0, aa_to_alpha_idx[seq[pos0]]]
            for a_i, aa in enumerate(STANDARD_AA):
                mat[pos0, a_i] = lp[pos0, aa_to_alpha_idx[aa]] - wt_lp
        if1_saturation[sym] = mat

        subs = tcga_variants_by_gene.get(sym, [])
        if subs:
            all_seqs = [("WT", seq)] + [
                (pch, seq[: pos - 1] + mut + seq[pos:])
                for pch, wt, pos, mut in subs
                if 1 <= pos <= seq_len and seq[pos - 1] == wt
            ]
            seq_total_ll = {}
            with torch.no_grad():
                for lbl, s_str in all_seqs:
                    _, _, _, tok_t, _ = coord_converter([(coords, None, s_str)], device=device)
                    lg, _ = model.decoder(tok_t[:, :-1], encoder_out=enc_out)
                    if lg.shape[0] == len(alphabet) and lg.shape[1] == 1:
                        lg_btv = lg[:, :, :seq_len].permute(1, 2, 0).float()
                    elif lg.shape[1] == len(alphabet):
                        lg_btv = lg[:, :, :seq_len].permute(0, 2, 1).float()
                    else:
                        lg_btv = lg[:, :seq_len, :].float()
                    lprobs = torch.log_softmax(lg_btv, dim=-1)
                    gathered = lprobs.gather(-1, tok_t[:, 1: seq_len + 1].unsqueeze(-1)).squeeze(-1)
                    seq_total_ll[lbl] = float((gathered * valid_mask_1d).sum().item())

            wt_ll = seq_total_ll["WT"]
            for pch, _, _, _ in subs:
                if pch in seq_total_ll:
                    if1_tcga_full[(sym, pch)] = round(seq_total_ll[pch] - wt_ll, 5)

    cache_dir.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(ckpt_npz, **if1_saturation)
    ckpt_json.write_text(json.dumps({f"{k[0]}|{k[1]}": v for k, v in if1_tcga_full.items()}), encoding="utf-8")
    del model
    return if1_saturation, if1_tcga_full


def run_sequence_plm(seqs: dict, model_name: str, loader_fn, device: torch.device, batch_size: int, cache_dir: pathlib.Path):
    ckpt_npz = cache_dir / f"ckpt_{model_name}.npz"
    sat_matrices = {}
    if ckpt_npz.exists():
        z = np.load(ckpt_npz)
        if all(sym in z for sym in GENE_ORDER):
            return {sym: z[sym] for sym in GENE_ORDER}
        for sym in z.files:
            sat_matrices[sym] = z[sym]

    model, alphabet = loader_fn()
    model = model.eval()
    if device.type == "cuda":
        model = model.half().to(device)
    else:
        model = model.to(device)
    batch_converter = alphabet.get_batch_converter()
    mask_idx = alphabet.mask_idx
    aa_to_tok = [alphabet.get_idx(aa) for aa in STANDARD_AA]

    cache_dir.mkdir(parents=True, exist_ok=True)
    for sym in GENE_ORDER:
        if sym in sat_matrices:
            continue
        seq = seqs[sym]
        seq_len = len(seq)
        mat = np.zeros((seq_len, 20), dtype=np.float32)
        positions = list(range(1, seq_len + 1))

        for b0 in range(0, seq_len, batch_size):
            b_pos = positions[b0: b0 + batch_size]
            b_data, w_idx0 = [], []
            for p in b_pos:
                wseq, wi = extract_context_window(seq, p, max_residues=1022)
                b_data.append((f"{sym}_{p}", wseq))
                w_idx0.append(wi)
            _, _, tokens = batch_converter(b_data)
            for r_i, wi in enumerate(w_idx0):
                tokens[r_i, wi + 1] = mask_idx
            tokens = tokens.to(device)

            with torch.no_grad():
                if device.type == "cuda":
                    with torch.amp.autocast("cuda", dtype=torch.float16):
                        logits = model(tokens, repr_layers=[], return_contacts=False)["logits"]
                else:
                    logits = model(tokens, repr_layers=[], return_contacts=False)["logits"]
            lprobs = torch.log_softmax(logits.float(), dim=-1).cpu().numpy()

            for r_i, p in enumerate(b_pos):
                lp_vec = lprobs[r_i, w_idx0[r_i] + 1]
                wt_lp = lp_vec[alphabet.get_idx(seq[p - 1])]
                mat[p - 1, :] = lp_vec[aa_to_tok] - wt_lp

        sat_matrices[sym] = mat
        np.savez_compressed(ckpt_npz, **sat_matrices)

    del model
    return sat_matrices


def download_if_needed(urls, dest_gz: pathlib.Path, label: str):
    if dest_gz.exists() and dest_gz.stat().st_size > 5_000_000:
        return
    for url in urls:
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=600, context=CTX) as resp:
                head = resp.read(2)
                if head == b"\x1f\x8b":
                    with open(dest_gz, "wb") as out:
                        out.write(head)
                        while True:
                            chunk = resp.read(1024 * 1024)
                            if not chunk:
                                break
                            out.write(chunk)
                else:
                    with gzip.open(dest_gz, "wb", compresslevel=4) as out:
                        out.write(head)
                        while True:
                            chunk = resp.read(1024 * 1024)
                            if not chunk:
                                break
                            out.write(chunk)
            if dest_gz.stat().st_size > 5_000_000:
                return
        except Exception:
            continue
    raise RuntimeError(f"Could not download {label}")


def run_mirna_mrna_couplings(cohort_samples: set, device: torch.device, cache_dir: pathlib.Path):
    cache_dir.mkdir(parents=True, exist_ok=True)
    eb_urls = [
        "https://api.gdc.cancer.gov/data/3586c0da-64d0-4b74-a449-5ff4d9136611",
        "https://pancanatlas.xenahubs.net/download/EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz",
    ]
    mir_urls = [
        "https://api.gdc.cancer.gov/data/1c6174d9-8ffb-466e-b5ee-07b204c15cf8",
        "https://pancanatlas.xenahubs.net/download/pancanMiRs_EBadjOnProtocolPlatformWithoutRepsWithUnCorrectMiRs_08_04_16.xena.gz",
    ]

    mrna_gz = cache_dir / "EBPlusPlusAdjustPANCAN_RNASeqV2.tsv.gz"
    mir_gz = cache_dir / "pancanMiRs_08_04_16.csv.gz"
    download_if_needed(mir_urls, mir_gz, "TCGA PanCanAtlas Mature miRNA Matrix")
    download_if_needed(eb_urls, mrna_gz, "TCGA PanCanAtlas EB++ mRNA Matrix")

    with gzip.open(mir_gz, "rt", encoding="utf-8", errors="replace") as f:
        first = f.readline().rstrip("\n\r")
        sep_m = "," if first.count(",") > first.count("\t") else "\t"
        hdr_m = [x.strip().strip('"') for x in first.split(sep_m)]
        has_corr = hdr_m[1].lower() == "correction"
        start_c = 2 if has_corr else 1
        mir_col_sid, seen_mir_sids = {}, set()
        for ci in range(start_c, len(hdr_m)):
            s15 = hdr_m[ci][:15]
            if s15 in cohort_samples and s15 not in seen_mir_sids:
                mir_col_sid[ci] = s15
                seen_mir_sids.add(s15)

        mir_sids_avail = set(mir_col_sid.values())
        col_indices_m = sorted(mir_col_sid.keys())
        mir_sid_order = [mir_col_sid[ci] for ci in col_indices_m]
        mir_names, mir_rows_raw = [], []

        for ln in f:
            p = ln.rstrip("\n\r").split(sep_m)
            if len(p) < start_c + 10 or (has_corr and p[1].strip().strip('"') != "Corrected"):
                continue
            vals = np.zeros(len(col_indices_m), dtype=np.float32)
            for k_i, ci in enumerate(col_indices_m):
                v_str = p[ci].strip().strip('"') if ci < len(p) else ""
                if v_str not in ("", "NA", "NaN"):
                    try:
                        v = float(v_str)
                        if not math.isnan(v) and v > 0:
                            vals[k_i] = math.log2(v + 1.0)
                    except Exception:
                        pass
            mir_names.append(p[0].strip().strip('"'))
            mir_rows_raw.append(vals)

    m_mir_full = np.stack(mir_rows_raw, axis=1)
    mir_sid_to_row = {sid: idx for idx, sid in enumerate(mir_sid_order)}

    with gzip.open(mrna_gz, "rt", encoding="utf-8", errors="replace") as f:
        hdr_g = [x.strip().strip('"') for x in f.readline().rstrip("\n\r").split("\t")]
        mrna_col_sid, seen_mrna_sids = {}, set()
        for ci in range(1, len(hdr_g)):
            s15 = hdr_g[ci][:15]
            if s15 in mir_sids_avail and s15 not in seen_mrna_sids:
                mrna_col_sid[ci] = s15
                seen_mrna_sids.add(s15)

        col_indices_g = sorted(mrna_col_sid.keys())
        dual_sids = [mrna_col_sid[ci] for ci in col_indices_g]
        n_dual = len(dual_sids)
        gene_names, gene_rows_raw = [], []
        scale_checked = False
        needs_log2 = True
        min_valid = int(0.80 * n_dual)
        bad_tokens = {"", "NA", "NaN", "nan", "null", '"NA"', '"NaN"'}

        for ln in f:
            p = ln.rstrip("\n\r").split("\t")
            if len(p) < 50:
                continue
            gsym = p[0].strip().strip('"').split("|")[0].strip()
            if not gsym or gsym == "?":
                continue
            vals = np.array(
                [float(p[ci].strip('"')) if (ci < len(p) and p[ci] not in bad_tokens) else np.nan for ci in col_indices_g],
                dtype=np.float32,
            )
            valid_mask = ~np.isnan(vals)
            valid_c = int(valid_mask.sum())
            if valid_c >= min_valid:
                if not scale_checked and float(np.nanmax(vals)) > 50.0:
                    needs_log2, scale_checked = True, True
                elif not scale_checked and len(gene_rows_raw) >= 25:
                    needs_log2, scale_checked = False, True
                vals = np.where(valid_mask, np.log2(np.maximum(0.0, vals) + 1.0) if needs_log2 else np.maximum(0.0, vals), np.nan)
                if valid_c < n_dual:
                    vals = np.where(valid_mask, vals, float(np.nanmedian(vals)))
                if float(np.std(vals)) > 0.05:
                    gene_names.append(gsym)
                    gene_rows_raw.append(vals.astype(np.float32, copy=False))

    m_mrna = np.stack(gene_rows_raw, axis=1)
    m_mir = np.stack([m_mir_full[mir_sid_to_row[sid]] for sid in dual_sids], axis=0)

    x_mir = torch.tensor(m_mir, device=device, dtype=torch.float32)
    x_mrna = torch.tensor(m_mrna, device=device, dtype=torch.float32)
    z_mir = (x_mir - x_mir.mean(dim=0, keepdim=True)) / x_mir.std(dim=0, keepdim=True).clamp(min=1e-6)
    z_mrna = (x_mrna - x_mrna.mean(dim=0, keepdim=True)) / x_mrna.std(dim=0, keepdim=True).clamp(min=1e-6)
    r_mat = (z_mir.T @ z_mrna) / float(n_dual - 1)
    n_pairs_total = int(r_mat.shape[0] * r_mat.shape[1])

    r_np = r_mat.cpu().numpy()
    z_mrna_np = z_mrna.cpu().numpy()

    flat_idx = np.argpartition(r_np.ravel(), 250)[:250]
    flat_idx = flat_idx[np.argsort(r_np.ravel()[flat_idx])]
    top_neg_pairs = []
    for rank_i, f_i in enumerate(flat_idx, 1):
        mi_i, ge_i = divmod(int(f_i), r_np.shape[1])
        r_val = float(r_np[mi_i, ge_i])
        t_stat = r_val * math.sqrt((n_dual - 2) / max(1e-12, 1.0 - r_val * r_val))
        top_neg_pairs.append({
            "rank": rank_i,
            "matureMiRNA": mir_names[mi_i],
            "targetGeneSymbol": gene_names[ge_i],
            "n_dual_tumours": n_dual,
            "pearson_r_log2": round(r_val, 5),
            "t_statistic": round(t_stat, 2),
        })

    gene_to_col = {}
    for idx_g, gsym in enumerate(gene_names):
        if gsym not in gene_to_col:
            gene_to_col[gsym] = idx_g

    let7_targets = [g for g in ("HMGA2", "IGF2BP1", "IGF2BP2", "IGF2BP3", "LIN28B", "MYCN", "NRAS", "CCND2", "CDC25A", "CDK6") if g in gene_to_col]
    mir200_targets = [g for g in ("ZEB1", "ZEB2", "FN1", "VIM", "SNAI2", "TWIST1") if g in gene_to_col]

    idx_5p = [i for i, m in enumerate(mir_names) if m.endswith("-5p")]
    idx_3p = [i for i, m in enumerate(mir_names) if m.endswith("-3p")]
    mean_r_5p = r_np[idx_5p, :].mean(axis=0)
    pref_5p_score = mean_r_5p - r_np[idx_3p, :].mean(axis=0)
    cand_mask = np.where(mean_r_5p < -0.08)[0]
    top100_5p_idx = cand_mask[np.argsort(pref_5p_score[cand_mask])[:100]] if len(cand_mask) >= 100 else np.argsort(pref_5p_score)[:100]
    top200_risc_idx = np.argsort(r_np.mean(axis=0))[:200]

    let7_tgt_z = z_mrna_np[:, [gene_to_col[g] for g in let7_targets]].mean(axis=1)
    m200_tgt_z = z_mrna_np[:, [gene_to_col[g] for g in mir200_targets]].mean(axis=1)
    emp_5p_derep_z = z_mrna_np[:, top100_5p_idx].mean(axis=1) - z_mrna_np[:, np.argsort(-pref_5p_score)[:100]].mean(axis=1)
    global_risc_derep_z = z_mrna_np[:, top200_risc_idx].mean(axis=1)

    u_mat, s_vec, _ = torch.linalg.svd(z_mir, full_matrices=False)
    mir_pc1 = (u_mat[:, 0] * s_vec[0]).cpu().numpy()

    derep_rows = []
    for i, sid in enumerate(dual_sids):
        derep_rows.append({
            "sampleId": sid,
            "let7_target_oncofetal_mrna_z": round(float(let7_tgt_z[i]), 5),
            "mir200_target_emt_mrna_z": round(float(m200_tgt_z[i]), 5),
            "empirical_5p_target_derepression_z": round(float(emp_5p_derep_z[i]), 5),
            "global_risc_target_derepression_z": round(float(global_risc_derep_z[i]), 5),
            "full662_mirnaome_pc1": round(float(mir_pc1[i]), 5),
        })

    return {
        "n_dual_tumours": n_dual,
        "n_mirnas": len(mir_names),
        "n_mrnas": len(gene_names),
        "n_pairs_evaluated": n_pairs_total,
        "let7_targets_used": let7_targets,
        "mir200_targets_used": mir200_targets,
        "top_neg_pairs": top_neg_pairs,
        "derep_rows": derep_rows,
    }


def main():
    ap = argparse.ArgumentParser(description="Recompute PLM/ESM-IF1 saturation scores and miRNA-mRNA target couplings.")
    ap.add_argument("--device", default="auto", help="Torch device: auto, cuda, mps, or cpu")
    ap.add_argument("--batch-size", type=int, default=16, help="Batch size for masked-marginal inference")
    args = ap.parse_args()

    device = select_device(args.device)
    print(f"Running PLM/ESM-IF1 and miRNA-mRNA coupling inference on device={device} (PyTorch {torch.__version__})")

    with gzip.open(EXT_DIR / "structural_3d_and_uniprot_15genes.json.gz", "rt", encoding="utf-8") as f:
        struct_snap = json.load(f)
    with gzip.open(EXT_DIR / "alphamissense_15genes_lookup.json.gz", "rt", encoding="utf-8") as f:
        am_snap = json.load(f)
    seqs = struct_snap["uniprot_sequences"]

    tcga_variants_by_gene = defaultdict(list)
    tcga_obs_counts = Counter()
    tcga_codon_counts = Counter()
    seen_vars = set()
    with open(DATA_DIR / "tcga_pancan_miRNA_biogenesis_somatic_mutations_annotated.csv", "r", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            sym = r["geneSymbol"]
            pch = r["proteinChange"].strip()
            pos_s = str(r.get("proteinPosStart") or "").strip()
            if pos_s.isdigit():
                tcga_codon_counts[(sym, int(pos_s))] += 1
            m = MISSENSE_RE.match(pch)
            if m:
                wt, pos, mut = m.group(1), int(m.group(2)), m.group(3)
                tcga_obs_counts[(sym, f"{wt}{pos}{mut}")] += 1
                if wt in STANDARD_AA_SET and mut in STANDARD_AA_SET and (sym, pch) not in seen_vars:
                    seen_vars.add((sym, pch))
                    tcga_variants_by_gene[sym].append((pch, wt, pos, mut))

    with open(DATA_DIR / "tcga_pancan_miRNA_biogenesis_patient_level_burden.csv", "r", encoding="utf-8") as f:
        cohort_samples = {r["sampleId"] for r in csv.DictReader(f)}

    cache_dir = WORK_ROOT / "gdc_cache"
    out_derep_csv = DATA_DIR / "tcga_pancan_sample_risc_target_derepression.csv"
    out_top_pairs_csv = DATA_DIR / "tcga_pancan_mirna_mrna_top250_repressive_axes.csv"

    m3_out = run_mirna_mrna_couplings(cohort_samples, device, cache_dir=cache_dir)
    with open(out_derep_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(m3_out["derep_rows"][0].keys()))
        w.writeheader()
        w.writerows(m3_out["derep_rows"])
    with open(out_top_pairs_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(m3_out["top_neg_pairs"][0].keys()))
        w.writeheader()
        w.writerows(m3_out["top_neg_pairs"])

    if1_sat, if1_tcga_full = run_esm_if1(seqs, tcga_variants_by_gene, device, cache_dir=cache_dir)

    import esm
    esm2_sat = run_sequence_plm(seqs, "esm2_t36_3B_UR50D", esm.pretrained.esm2_t36_3B_UR50D, device, batch_size=args.batch_size, cache_dir=cache_dir)
    esm1v_sat = run_sequence_plm(seqs, "esm1v_t33_650M_UR90S_1", esm.pretrained.esm1v_t33_650M_UR90S_1, device, batch_size=args.batch_size * 2, cache_dir=cache_dir)

    aa_idx_map = {aa: i for i, aa in enumerate(STANDARD_AA)}
    out_tcga_var_csv = DATA_DIR / "tcga_somatic_variant_plm_if1_scores.csv"
    tcga_var_rows = []
    for sym, subs in sorted(tcga_variants_by_gene.items()):
        seq = seqs[sym]
        for pch, wt, pos, mut in sorted(subs, key=lambda t: (t[2], t[0])):
            if not (1 <= pos <= len(seq)) or seq[pos - 1] != wt:
                continue
            mi = aa_idx_map[mut]
            e2 = round(float(esm2_sat[sym][pos - 1, mi]), 5)
            e1 = round(float(esm1v_sat[sym][pos - 1, mi]), 5)
            if1_m = round(float(if1_sat[sym][pos - 1, mi]), 5)
            if1_f = if1_tcga_full.get((sym, pch), if1_m)
            esm_tier = "Deleterious_ESM2_LLR_le_m7" if e2 <= -7.0 else ("Intermediate_ESM2" if e2 <= -4.0 else "Neutral_ESM2_LLR_gt_m4")
            if e2 <= -5.0 or e1 <= -5.0:
                mech = "Thermodynamic_Fold_Destabilizing" if float(if1_f) <= -1.5 else "Fold_Preserving_Catalytic_or_Interface"
            else:
                mech = "Thermodynamic_Fold_Destabilizing" if float(if1_f) <= -2.0 else "Tolerated_Neutral"
            tcga_var_rows.append({
                "geneSymbol": sym,
                "proteinChange": pch,
                "esm2_3b_llr": e2,
                "esm1v_650m_llr": e1,
                "esm2_3b_effect_tier": esm_tier,
                "esm_if1_3d_dlogp": if1_f,
                "esm_if1_3d_marginal_llr": if1_m,
                "biophysicalMechanismClass": mech,
            })

    with open(out_tcga_var_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(tcga_var_rows[0].keys()))
        w.writeheader()
        w.writerows(tcga_var_rows)

    sat_gz_path = DATA_DIR / "miRNA_biogenesis_saturation_mutagenesis_213940_atlas.csv.gz"
    codon_gz_path = DATA_DIR / "miRNA_biogenesis_residue_biophysical_11260_profile.csv.gz"
    sat_fields = [
        "geneSymbol", "uniprotAccession", "codonPos", "wtAA", "mutAA", "variant",
        "pdbComplexId", "distSubstrateA", "distMetalA", "distPartnerA", "afdbPlddt",
        "alphamissenseScore", "alphamissenseClass", "esm2_3b_llr", "esm1v_650m_llr",
        "esm_if1_3d_llr", "observedTcgaSomaticCount",
    ]
    codon_fields = [
        "geneSymbol", "uniprotAccession", "codonPos", "wtAA", "afdbPlddt",
        "pdbComplexId", "distSubstrateA", "distMetalA", "distPartnerA",
        "mean_alphamissense", "mean_esm2_3b_llr", "mean_esm1v_650m_llr", "mean_esm_if1_3d_llr",
        "observedTcgaCodonMutCount",
    ]

    total_sat_rows, total_codon_rows = 0, 0
    with gzip.open(sat_gz_path, "wt", newline="", encoding="utf-8", compresslevel=6) as f_sat, \
         gzip.open(codon_gz_path, "wt", newline="", encoding="utf-8", compresslevel=6) as f_cod:
        w_sat, w_cod = csv.writer(f_sat), csv.writer(f_cod)
        w_sat.writerow(sat_fields)
        w_cod.writerow(codon_fields)

        for sym in GENE_ORDER:
            acc = UNIPROT_MAP[sym]
            seq = seqs[sym]
            per_res = struct_snap["per_residue"][sym]
            am_gene = am_snap["genes"][sym]["subs"]
            am_rmean = am_snap["genes"][sym]["res_mean"]
            m_e2, m_e1, m_if = esm2_sat[sym], esm1v_sat[sym], if1_sat[sym]

            for pos in range(1, len(seq) + 1):
                wt = seq[pos - 1]
                rinfo = per_res.get(str(pos), [0.0, "", -1.0, -1.0, -1.0])
                plddt, pdb_id = rinfo[0], rinfo[1]
                d_na = rinfo[2] if rinfo[2] >= 0 else ""
                d_met = rinfo[3] if rinfo[3] >= 0 else ""
                d_prt = rinfo[4] if rinfo[4] >= 0 else ""

                e2_list, e1_list, if_list = [], [], []
                for a_i, mut in enumerate(STANDARD_AA):
                    if mut == wt:
                        continue
                    vkey = f"{wt}{pos}{mut}"
                    am_hit = am_gene.get(vkey, ["", ""])
                    e2_v = round(float(m_e2[pos - 1, a_i]), 4)
                    e1_v = round(float(m_e1[pos - 1, a_i]), 4)
                    if_v = round(float(m_if[pos - 1, a_i]), 4)
                    e2_list.append(e2_v)
                    e1_list.append(e1_v)
                    if_list.append(if_v)
                    w_sat.writerow([
                        sym, acc, pos, wt, mut, vkey,
                        pdb_id, d_na, d_met, d_prt, plddt,
                        am_hit[0], am_hit[1], e2_v, e1_v, if_v,
                        tcga_obs_counts.get((sym, vkey), 0),
                    ])
                    total_sat_rows += 1

                w_cod.writerow([
                    sym, acc, pos, wt, plddt,
                    pdb_id, d_na, d_met, d_prt,
                    am_rmean.get(str(pos), ""),
                    round(sum(e2_list) / len(e2_list), 4),
                    round(sum(e1_list) / len(e1_list), 4),
                    round(sum(if_list) / len(if_list), 4),
                    tcga_codon_counts.get((sym, pos), 0),
                ])
                total_codon_rows += 1

    prov = {
        "completed_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "software": {"pytorch_version": torch.__version__.split("+")[0]},
        "esm_if1_inverse_folding": {
            "model": "esm_if1_gvp4_t16_142M_UR50",
            "saturation_variants_scored": total_sat_rows,
            "tcga_full_autoregressive_variants_scored": len(if1_tcga_full),
        },
        "sequence_saturation_mutagenesis": {
            "models": ["esm2_t36_3B_UR50D", "esm1v_t33_650M_UR90S_1"],
            "total_codons_scored": total_codon_rows,
            "total_single_aa_variants_scored": total_sat_rows,
        },
        "mirna_mrna_target_coupling": {
            "n_dual_profiled_tcga_tumours": m3_out["n_dual_tumours"],
            "n_mature_mirnas": m3_out["n_mirnas"],
            "n_expressed_mrnas": m3_out["n_mrnas"],
            "total_mirna_mrna_pairs_evaluated": m3_out["n_pairs_evaluated"],
            "let7_targets_used": m3_out["let7_targets_used"],
            "mir200_targets_used": m3_out["mir200_targets_used"],
        },
        "outputs": {
            out_tcga_var_csv.name: sha256_file(out_tcga_var_csv),
            sat_gz_path.name: sha256_file(sat_gz_path),
            codon_gz_path.name: sha256_file(codon_gz_path),
            out_derep_csv.name: sha256_file(out_derep_csv),
            out_top_pairs_csv.name: sha256_file(out_top_pairs_csv),
        },
    }
    (DATA_DIR / "plm_inference_provenance.json").write_text(json.dumps(prov, indent=2) + "\n", encoding="utf-8")
    print("Completed PLM/ESM-IF1 saturation and miRNA-mRNA coupling outputs.")


if __name__ == "__main__":
    main()
