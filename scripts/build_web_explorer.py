#!/usr/bin/env python3
"""Compile CSV/CSV.GZ tables and 3D atomic coordinates into docs/index.html, docs/sat/*.js, and docs/struct3d/*.js."""

import csv
import gzip
import json
import math
import pathlib
import sys
from collections import Counter, defaultdict

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from explainer_content import (
    DATASET_DESCRIPTIONS,
    FORMULA_CARDS,
    SOURCE_CONSORTIA,
    TERM_GLOSSARY,
)

REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent
DATA_DIR = REPO_ROOT / "data"
DOCS_DIR = REPO_ROOT / "docs"
SAT_DIR = DOCS_DIR / "sat"
STRUCT3D_DIR = DOCS_DIR / "struct3d"


def read_csv_rows(path: pathlib.Path):
    if not path.exists():
        return []
    opener = gzip.open if path.name.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def _sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _cross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def _dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _norm(a):
    length = math.sqrt(_dot(a, a)) or 1.0
    return (a[0] / length, a[1] / length, a[2] / length)


def _rotamer_sq_dist(sc_a, sc_b):
    map_b = {x[0]: x[1:] for x in sc_b}
    tot = 0.0
    n = 0
    for name, u, v, w in sc_a:
        if name in map_b:
            ub, vb, wb = map_b[name]
            tot += (u - ub) ** 2 + (v - vb) ** 2 + (w - wb) ** 2
            n += 1
    return tot / max(1, n)


def build_struct3d_chunks_and_rotamers():
    """Export lazy-loaded 3D heavy-atom JS files per gene and a compact 20-AA rotamer template library."""
    snap_gz = DATA_DIR / "raw_snapshot" / "external" / "afdb_v6_15genes_heavy_atoms.json.gz"
    if not snap_gz.exists():
        return {}

    STRUCT3D_DIR.mkdir(parents=True, exist_ok=True)
    with gzip.open(snap_gz, "rt", encoding="utf-8") as f:
        snap = json.load(f)

    candidates_by_aa = defaultdict(list)
    total_kb = 0.0

    for gene, gd in sorted(snap.get("genes", {}).items()):
        out_obj = {
            "u": gd["uniprot"],
            "s": gd["sequence"],
            "a": gd["atoms"],
        }
        js_str = (
            "window.STRUCT3D_CACHE=window.STRUCT3D_CACHE||{};window.STRUCT3D_CACHE["
            + json.dumps(gene)
            + "]="
            + json.dumps(out_obj, separators=(",", ":"))
            + ";\n"
        )
        out_f = STRUCT3D_DIR / f"{gene}.js"
        out_f.write_text(js_str, encoding="utf-8")
        total_kb += out_f.stat().st_size / 1024.0

        seq = gd["sequence"]
        res_atoms = defaultdict(dict)
        for seq_id, aname, x, y, z, plddt in gd["atoms"]:
            res_atoms[seq_id][aname] = (x, y, z, plddt)

        for pos, amap in res_atoms.items():
            aa = seq[pos - 1]
            if aa == "G" or not all(k in amap for k in ("N", "CA", "C", "CB")):
                continue
            if amap["CA"][3] < 88.0:
                continue
            r_n, r_ca, r_c = amap["N"][:3], amap["CA"][:3], amap["C"][:3]
            u = _norm(_sub(r_n, r_ca))
            w = _norm(_cross(u, _sub(r_c, r_ca)))
            v = _cross(w, u)
            sc = []
            for aname in sorted(amap.keys()):
                if aname in ("N", "CA", "C", "O", "OXT"):
                    continue
                d = _sub(amap[aname][:3], r_ca)
                sc.append([aname, round(_dot(d, u), 3), round(_dot(d, v), 3), round(_dot(d, w), 3)])
            candidates_by_aa[aa].append(sc)

    rotamers = {"G": [[]]}
    for aa, pool in sorted(candidates_by_aa.items()):
        k_keep = 1 if aa == "A" else min(8, len(pool))
        chosen = [pool[0]]
        while len(chosen) < k_keep:
            best_cand = max(pool, key=lambda c: min(_rotamer_sq_dist(c, prev) for prev in chosen))
            chosen.append(best_cand)
        rotamers[aa] = chosen

    print(f"   Exported 15 lazy-loaded 3D structure chunks in docs/struct3d/*.js ({total_kb:.1f} KB total)")
    return rotamers


def build_saturation_js_chunks(clinvar_rows):
    """Export 15 compact docs/sat/<GENE>.js files covering all 213,940 variants and codon tracks."""
    sat_gz = DATA_DIR / "miRNA_biogenesis_saturation_mutagenesis_213940_atlas.csv.gz"
    res_gz = DATA_DIR / "miRNA_biogenesis_residue_biophysical_11260_profile.csv.gz"
    if not (sat_gz.exists() and res_gz.exists()):
        return {}

    SAT_DIR.mkdir(parents=True, exist_ok=True)

    cv_pos_path = Counter()
    cv_var_sig = {}
    for r in clinvar_rows:
        g = r["geneSymbol"]
        pos = int(r["codonPos"]) if str(r.get("codonPos", "")).isdigit() else 0
        if r.get("significanceCategory") == "Pathogenic_or_Likely_Pathogenic" and pos > 0:
            cv_pos_path[(g, pos)] += 1
        cv_var_sig[(g, r["proteinChange"])] = r["clinicalSignificance"]

    tracks = defaultdict(list)
    with gzip.open(res_gz, "rt", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            g = r["geneSymbol"]
            pos = int(r["codonPos"])
            is_pkt = int(any(r.get(k, "") != "" and float(r[k]) <= 6.0 for k in ("distSubstrateA", "distMetalA", "distPartnerA")))
            tracks[g].append([
                pos,
                r["wtAA"],
                round(float(r["afdbPlddt"]), 1) if r["afdbPlddt"] != "" else None,
                is_pkt,
                round(float(r["mean_alphamissense"]), 3) if r["mean_alphamissense"] != "" else None,
                round(float(r["mean_esm2_3b_llr"]), 2) if r["mean_esm2_3b_llr"] != "" else None,
                round(float(r["mean_esm_if1_3d_llr"]), 2) if r["mean_esm_if1_3d_llr"] != "" else None,
                int(r["observedTcgaCodonMutCount"]),
                cv_pos_path.get((g, pos), 0),
                r["pdbComplexId"] or "",
            ])

    per_gene_codons = defaultdict(lambda: defaultdict(list))
    per_gene_meta = {}
    with gzip.open(sat_gz, "rt", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            g = r["geneSymbol"]
            pos = int(r["codonPos"])
            if (g, pos) not in per_gene_meta:
                per_gene_meta[(g, pos)] = [
                    r["wtAA"],
                    r["pdbComplexId"] or "",
                    round(float(r["distSubstrateA"]), 2) if r["distSubstrateA"] != "" else None,
                    round(float(r["distMetalA"]), 2) if r["distMetalA"] != "" else None,
                    round(float(r["distPartnerA"]), 2) if r["distPartnerA"] != "" else None,
                    round(float(r["afdbPlddt"]), 1) if r["afdbPlddt"] != "" else None,
                ]
            vname = r["variant"]
            per_gene_codons[g][pos].append([
                r["mutAA"],
                round(float(r["alphamissenseScore"]), 3) if r["alphamissenseScore"] != "" else None,
                round(float(r["esm2_3b_llr"]), 2) if r["esm2_3b_llr"] != "" else None,
                round(float(r["esm1v_650m_llr"]), 2) if r["esm1v_650m_llr"] != "" else None,
                round(float(r["esm_if1_3d_llr"]), 2) if r["esm_if1_3d_llr"] != "" else None,
                int(r["observedTcgaSomaticCount"]),
                cv_var_sig.get((g, vname), ""),
            ])

    total_kb = 0.0
    for g, cmap in sorted(per_gene_codons.items()):
        out_obj = {}
        for pos in sorted(cmap.keys()):
            out_obj[str(pos)] = {
                "m": per_gene_meta[(g, pos)],
                "s": cmap[pos],
            }
        js_payload = "window.SAT_CACHE=window.SAT_CACHE||{};window.SAT_CACHE[" + json.dumps(g) + "]=" + json.dumps(out_obj, separators=(",", ":")) + ";\n"
        out_f = SAT_DIR / f"{g}.js"
        out_f.write_text(js_payload, encoding="utf-8")
        total_kb += out_f.stat().st_size / 1024.0
    print(f"   Exported 15 lazy-loaded saturation mutagenesis chunks in docs/sat/*.js ({total_kb:.1f} KB total across 213,940 variants)")
    return dict(tracks)


def build_web_explorer():
    DOCS_DIR.mkdir(parents=True, exist_ok=True)

    panel = read_csv_rows(DATA_DIR / "miRNA_biogenesis_gene_panel_metadata.csv")
    long_df = read_csv_rows(DATA_DIR / "tcga_pancan_miRNA_biogenesis_alteration_freqs_long_FINAL.csv")
    muts = read_csv_rows(DATA_DIR / "tcga_pancan_miRNA_biogenesis_somatic_mutations_annotated.csv")
    mir_bench = read_csv_rows(DATA_DIR / "validation_mirnaseq_enzymatic_benchmarks.csv")
    mir_rows = read_csv_rows(DATA_DIR / "tcga_pancan_miRNA_biogenesis_mirnaseq_functional_readouts.csv")
    cptac = read_csv_rows(DATA_DIR / "cptac_pancan_miRNA_biogenesis_proteogenomic_validation.csv")
    depmap = read_csv_rows(DATA_DIR / "depmap_crispr_miRNA_biogenesis_essentiality.csv")
    collat = read_csv_rows(DATA_DIR / "validation_depmap_collateral_vulnerability.csv")
    top_axes = read_csv_rows(DATA_DIR / "tcga_pancan_mirna_mrna_top250_repressive_axes.csv")[:60]
    gtex_rows = read_csv_rows(DATA_DIR / "gtex_v10_normal_54tissues_miRNA_biogenesis_tpm.csv")
    clinvar_rows = read_csv_rows(DATA_DIR / "clinvar_clinical_variants_15genes_annotated.csv")
    clinvar_roc = read_csv_rows(DATA_DIR / "validation_clinvar_roc_auc_benchmarks.csv")
    msk_rep = read_csv_rows(DATA_DIR / "msk_impact_10945_clinical_cohort_replication.csv")
    msk_dicer = read_csv_rows(DATA_DIR / "msk_impact_dicer1_cancer_type_enrichment.csv")[:20]
    surv_bench = read_csv_rows(DATA_DIR / "validation_clinical_survival_benchmarks.csv")
    imm_bench = read_csv_rows(DATA_DIR / "validation_panimmune_and_hrd_benchmarks.csv")
    study_summary = read_csv_rows(DATA_DIR / "tcga_pancan_atlas_sample_set_summary.csv")
    excluded_cands = read_csv_rows(DATA_DIR / "miRNA_biogenesis_excluded_candidates_rationale.csv")
    top_tumours = read_csv_rows(DATA_DIR / "miRNA_biogenesis_top_tumours_excluding_TP53.csv")

    derived_catalogue = []
    for fpath in sorted(DATA_DIR.iterdir(), key=lambda p: p.name):
        if not fpath.is_file() or not (fpath.name.endswith(".csv") or fpath.name.endswith(".csv.gz")):
            continue
        opener = gzip.open if fpath.name.endswith(".gz") else open
        with opener(fpath, "rt", encoding="utf-8") as fh:
            rdr = csv.reader(fh)
            hdr = next(rdr, [])
            n_rows = sum(1 for _ in rdr)
        sz_bytes = fpath.stat().st_size
        sz_str = f"{sz_bytes / (1024 * 1024):.2f} MB" if sz_bytes >= 512 * 1024 else f"{max(1.0, sz_bytes / 1024.0):.1f} KB"
        pk, desc = DATASET_DESCRIPTIONS.get(fpath.name, ("Key columns", "Derived analytical table in data/."))
        derived_catalogue.append({
            "file": fpath.name,
            "rows": n_rows,
            "cols": len(hdr),
            "size": sz_str,
            "pk": pk,
            "desc": desc,
        })

    residue_tracks = build_saturation_js_chunks(clinvar_rows)
    rotamer_lib = build_struct3d_chunks_and_rotamers()

    compact_long = [
        {
            "t": r["tumour_acronym"],
            "g": r["geneSymbol"],
            "N": int(r["N"]),
            "any": round(float(r["any_alter_freq"]), 4),
            "hc": round(float(r["high_conf_alter_freq"]), 4),
            "fcna": round(float(r["focal_cna_freq"]), 4),
            "acna": round(float(r["any_cna_freq"]), 4),
            "mut": round(float(r["mut_freq"]), 4),
            "nhm": round(float(r["nonhypermut_mut_freq"]), 4),
            "amp_n": int(r["amp_n"]),
            "dd_n": int(r["deep_del_n"]),
            "mut_n": int(r["mut_n"]),
        }
        for r in long_df
    ]

    var_agg = {}
    var_cancers = defaultdict(Counter)
    for m in muts:
        key = (m["geneSymbol"], m["proteinChange"] or m["mutationType"])
        var_cancers[key][m["tumour_acronym"]] += 1
        if key not in var_agg:
            var_agg[key] = {
                "g": m["geneSymbol"],
                "p": m["proteinChange"] or m["mutationType"],
                "pos": int(m["proteinPosStart"]) if str(m["proteinPosStart"]).isdigit() else 0,
                "mtype": m["mutationType"],
                "tier": m["functionalVariantTierV2"],
                "pdb": m["pdbComplexId"],
                "dna": float(m["distToSubstrateA"]) if m["distToSubstrateA"] != "" else None,
                "dmet": float(m["distToCatalyticMetalA"]) if m["distToCatalyticMetalA"] != "" else None,
                "plddt": float(m["afdbResiduePlddt"]) if m["afdbResiduePlddt"] != "" else None,
                "am": float(m["alphamissenseScore"]) if m["alphamissenseScore"] != "" else None,
                "amc": m["alphamissenseClass"],
                "esm": float(m["esm2_3b_llr"]) if m.get("esm2_3b_llr", "") not in ("", None) else None,
                "esm1v": float(m["esm1v_650m_llr"]) if m.get("esm1v_650m_llr", "") not in ("", None) else None,
                "if1": float(m["esm_if1_3d_dlogp"]) if m.get("esm_if1_3d_dlogp", "") not in ("", None) else None,
                "mech": m.get("biophysicalMechanismClass", "") or "",
                "cv": m.get("clinvarClinicalSignificance", "") or "",
                "n": 0,
            }
        var_agg[key]["n"] += 1

    compact_vars = []
    for key, rec in var_agg.items():
        rec["cohorts"] = ", ".join(f"{ac}({cnt})" for ac, cnt in var_cancers[key].most_common(4))
        compact_vars.append(rec)
    compact_vars.sort(key=lambda x: (-x["n"], x["g"], x["pos"]))

    dicer_pts = []
    for r in mir_rows:
        if r["mirna_duplex_5p_3p_log2_ratio"] != "" and r["dicer1FunctionalClass"] != "DICER1_WildType":
            dicer_pts.append({
                "s": r["sampleId"],
                "t": r["tumour_acronym"],
                "c": r["dicer1FunctionalClass"],
                "d53": round(float(r["mirna_duplex_5p_3p_log2_ratio"]), 4),
                "f5": round(float(r["mirna_5p_arm_fraction"]), 4) if r["mirna_5p_arm_fraction"] != "" else None,
                "l7z": round(float(r["let7_target_oncofetal_mrna_z"]), 4) if r.get("let7_target_oncofetal_mrna_z", "") not in ("", None) else None,
            })

    gtex_tissues = sorted({r["tissueSiteDetailId"] for r in gtex_rows})
    t_idx = {t: i for i, t in enumerate(gtex_tissues)}
    gtex_by_gene = {p["geneSymbol"]: [0.0] * len(gtex_tissues) for p in panel}
    for r in gtex_rows:
        if r["geneSymbol"] in gtex_by_gene:
            gtex_by_gene[r["geneSymbol"]][t_idx[r["tissueSiteDetailId"]]] = round(float(r["medianTPM"]), 2)

    bundle = {
        "meta": {
            "version": "2.0.2",
            "doi": "10.5281/zenodo.24455067",
            "tcga_n": 10195,
            "msk_n": 10945,
            "combined_tumours_n": 21140,
            "mirna_n": len(mir_rows),
            "dual_mrna_mirna_n": 9353,
            "sat_mut_n": 213940,
            "clinvar_n": 16569,
            "cptac_n": 1002,
            "depmap_n": 1208,
        },
        "panel": panel,
        "long": compact_long,
        "variants": compact_vars,
        "tracks": residue_tracks,
        "rotamers": rotamer_lib,
        "mir_bench": mir_bench,
        "dicer_pts": dicer_pts,
        "top_axes": top_axes,
        "cptac": cptac,
        "depmap": depmap,
        "collateral": collat,
        "clinvar_roc": clinvar_roc,
        "msk_rep": msk_rep,
        "msk_dicer": msk_dicer,
        "gtex_tissues": gtex_tissues,
        "gtex_tpm": gtex_by_gene,
        "surv_bench": surv_bench,
        "imm_bench": imm_bench,
        "study_summary": study_summary,
        "excluded_cands": excluded_cands,
        "top_tumours": top_tumours,
        "datasets": derived_catalogue,
        "consortia": SOURCE_CONSORTIA,
        "formulas": FORMULA_CARDS,
        "glossary": TERM_GLOSSARY,
    }

    data_json = json.dumps(bundle, separators=(",", ":"))

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>miRBioGenAtlas — Human MicroRNA Biogenesis Somatic &amp; Structural Reference</title>
<script src="vendor/3Dmol-min.js"></script>
<style>
  :root {{
    --canvas: #f5f5f7;
    --surface: #ffffff;
    --surface-subtle: #fafafa;
    --ink: #1d1d1f;
    --ink-secondary: #515154;
    --ink-tertiary: #6e6e73;
    --line: #d2d2d7;
    --line-light: #e8e8ed;
    --embl-petrol: #00709b;
    --embl-petrol-dark: #004f6e;
    --embl-teal: #007c82;
    --crimson: #b91c1c;
    --amber: #b45309;
    --emerald: #047857;
    --mono: "SF Mono", "JetBrains Mono", ui-monospace, Menlo, Consolas, monospace;
    --sans: -apple-system, BlinkMacSystemFont, "SF Pro Text", "IBM Plex Sans", "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
  }}
  * {{ box-sizing: border-box; }}
  html {{ -webkit-font-smoothing: antialiased; -moz-osx-font-smoothing: grayscale; }}
  body {{
    margin: 0;
    font-family: var(--sans);
    background: var(--canvas);
    color: var(--ink);
    font-size: 14px;
    line-height: 1.47;
  }}
  .container {{ max-width: 1440px; margin: 0 auto; padding: 0 24px; }}

  /* EMBL-EBI style top institutional utility bar */
  .global-bar {{
    background: #1d1d1f;
    color: #d2d2d7;
    font-size: 12px;
    height: 36px;
    display: flex;
    align-items: center;
    border-bottom: 1px solid #2c2c2e;
  }}
  .global-bar .container {{
    width: 100%;
    display: flex;
    justify-content: space-between;
    align-items: center;
    gap: 16px;
  }}
  .global-brand {{
    display: flex;
    align-items: center;
    gap: 10px;
    color: #f5f5f7;
    font-weight: 500;
    letter-spacing: -0.01em;
  }}
  .embl-mark {{
    width: 8px;
    height: 8px;
    border-radius: 50%;
    background: #009f4d;
    box-shadow: 12px 0 0 #00709b;
    margin-right: 12px;
    display: inline-block;
  }}
  .global-links {{
    display: flex;
    align-items: center;
    gap: 18px;
    font-family: var(--mono);
    font-size: 11.5px;
  }}
  .global-links a {{
    color: #d2d2d7;
    text-decoration: none;
    transition: color 0.15s ease;
  }}
  .global-links a:hover {{ color: #ffffff; }}

  /* Architectural masthead */
  .masthead {{
    background: var(--surface);
    border-bottom: 1px solid var(--line-light);
    padding: 22px 0 0;
  }}
  .masthead-top {{
    display: flex;
    justify-content: space-between;
    align-items: flex-end;
    flex-wrap: wrap;
    gap: 20px;
    padding-bottom: 18px;
  }}
  .title-block h1 {{
    margin: 0;
    font-size: 26px;
    font-weight: 650;
    letter-spacing: -0.028em;
    color: var(--ink);
    display: flex;
    align-items: baseline;
    gap: 10px;
  }}
  .version-tag {{
    font-family: var(--mono);
    font-size: 11.5px;
    font-weight: 500;
    color: var(--embl-petrol);
    border: 1px solid rgba(0, 112, 155, 0.3);
    padding: 1px 7px;
    border-radius: 4px;
    letter-spacing: 0;
  }}
  .title-block p {{
    margin: 4px 0 0;
    color: var(--ink-secondary);
    font-size: 13.5px;
    max-width: 780px;
  }}
  .omni-search {{
    position: relative;
    width: 350px;
    max-width: 100%;
  }}
  .omni-search input {{
    width: 100%;
    height: 36px;
    padding: 0 34px 0 12px;
    background: var(--canvas);
    border: 1px solid var(--line);
    border-radius: 7px;
    font-size: 13px;
    color: var(--ink);
    font-family: var(--sans);
    transition: border-color 0.15s ease, background 0.15s ease;
  }}
  .omni-search input:focus {{
    outline: none;
    background: var(--surface);
    border-color: var(--embl-petrol);
    box-shadow: 0 0 0 3px rgba(0, 112, 155, 0.12);
  }}
  .omni-kbd {{
    position: absolute;
    right: 9px;
    top: 8px;
    font-family: var(--mono);
    font-size: 11px;
    color: var(--ink-tertiary);
    border: 1px solid var(--line);
    border-radius: 4px;
    padding: 0 5px;
    background: var(--surface);
    pointer-events: none;
  }}

  /* Restrained specification bar (replaces KPI cards) */
  .spec-strip {{
    display: grid;
    grid-template-columns: repeat(6, 1fr);
    border-top: 1px solid var(--line-light);
    padding: 12px 0;
    gap: 16px;
  }}
  @media (max-width: 1024px) {{
    .spec-strip {{ grid-template-columns: repeat(3, 1fr); }}
  }}
  @media (max-width: 640px) {{
    .spec-strip {{ grid-template-columns: repeat(2, 1fr); }}
  }}
  .spec-item {{
    border-left: 1px solid var(--line-light);
    padding-left: 14px;
  }}
  .spec-item:first-child {{ border-left: none; padding-left: 0; }}
  .spec-val {{
    font-size: 16px;
    font-weight: 600;
    letter-spacing: -0.02em;
    color: var(--ink);
    font-variant-numeric: tabular-nums;
  }}
  .spec-lbl {{
    font-size: 11.5px;
    color: var(--ink-tertiary);
    margin-top: 1px;
  }}

  /* Sticky frosted sub-navigation */
  .subnav-wrap {{
    position: sticky;
    top: 0;
    z-index: 30;
    background: rgba(255, 255, 255, 0.88);
    backdrop-filter: saturate(180%) blur(18px);
    -webkit-backdrop-filter: saturate(180%) blur(18px);
    border-bottom: 1px solid var(--line);
  }}
  .tabs {{
    display: flex;
    gap: 24px;
    overflow-x: auto;
    scrollbar-width: none;
  }}
  .tabs::-webkit-scrollbar {{ display: none; }}
  .tab-btn {{
    border: none;
    background: transparent;
    padding: 12px 2px 11px;
    font-family: var(--sans);
    font-size: 13px;
    font-weight: 500;
    color: var(--ink-secondary);
    cursor: pointer;
    border-bottom: 2px solid transparent;
    white-space: nowrap;
    transition: color 0.15s ease, border-color 0.15s ease;
  }}
  .tab-btn:hover {{ color: var(--ink); }}
  .tab-btn.active {{
    color: var(--embl-petrol);
    font-weight: 600;
    border-bottom-color: var(--embl-petrol);
  }}

  /* Content layout & surfaces */
  main.workspace {{ padding-top: 22px; padding-bottom: 36px; }}
  .panel {{ display: none; }}
  .panel.active {{ display: block; }}
  .surface {{
    background: var(--surface);
    border: 1px solid var(--line-light);
    border-radius: 10px;
    padding: 20px 22px;
    margin-bottom: 18px;
  }}
  .section-head {{
    display: flex;
    justify-content: space-between;
    align-items: baseline;
    flex-wrap: wrap;
    gap: 10px;
    margin-bottom: 14px;
    padding-bottom: 10px;
    border-bottom: 1px solid var(--line-light);
  }}
  .section-head h2 {{
    margin: 0;
    font-size: 16px;
    font-weight: 600;
    letter-spacing: -0.018em;
    color: var(--ink);
  }}
  .section-meta {{
    font-size: 12px;
    color: var(--ink-tertiary);
    font-variant-numeric: tabular-nums;
  }}

  /* Clean form controls & segmented groups */
  .toolbar {{
    display: flex;
    flex-wrap: wrap;
    gap: 12px;
    align-items: center;
    margin-bottom: 14px;
  }}
  .field-group {{
    display: inline-flex;
    align-items: center;
    gap: 7px;
    font-size: 12.5px;
    color: var(--ink-secondary);
  }}
  .field-label {{
    font-weight: 500;
    color: var(--ink-tertiary);
    text-transform: uppercase;
    font-size: 11px;
    letter-spacing: 0.04em;
  }}
  select, input[type="text"] {{
    height: 32px;
    padding: 0 10px;
    border: 1px solid var(--line);
    border-radius: 6px;
    font-size: 12.5px;
    font-family: var(--sans);
    color: var(--ink);
    background: var(--surface);
  }}
  select:focus, input[type="text"]:focus {{
    outline: none;
    border-color: var(--embl-petrol);
  }}
  .btn-primary {{
    height: 32px;
    background: var(--embl-petrol);
    color: #ffffff;
    border: none;
    padding: 0 13px;
    border-radius: 6px;
    font-size: 12.5px;
    font-weight: 500;
    cursor: pointer;
    transition: background 0.15s ease;
  }}
  .btn-primary:hover {{ background: var(--embl-petrol-dark); }}
  .btn-quiet {{
    height: 28px;
    background: var(--surface);
    color: var(--ink);
    border: 1px solid var(--line);
    padding: 0 10px;
    border-radius: 6px;
    font-size: 12px;
    font-weight: 500;
    cursor: pointer;
    transition: background 0.12s ease;
  }}
  .btn-quiet:hover {{ background: var(--canvas); }}

  /* Reference loci inline strip */
  .ref-strip {{
    display: flex;
    flex-wrap: wrap;
    align-items: center;
    gap: 6px;
    font-size: 12px;
    color: var(--ink-tertiary);
    margin-bottom: 12px;
  }}
  .ref-link {{
    background: transparent;
    border: 1px solid var(--line-light);
    color: var(--ink-secondary);
    font-family: var(--mono);
    font-size: 11.5px;
    padding: 2px 8px;
    border-radius: 5px;
    cursor: pointer;
    transition: all 0.12s ease;
  }}
  .ref-link:hover {{
    border-color: var(--embl-petrol);
    color: var(--embl-petrol);
    background: rgba(0, 112, 155, 0.04);
  }}

  /* Linear sequence & pocket SVG track */
  .svg-box {{
    width: 100%;
    overflow-x: auto;
    border: 1px solid var(--line-light);
    border-radius: 8px;
    background: var(--surface-subtle);
    padding: 8px 10px 4px;
  }}

  /* 3D Structure & Residue Specification Workbench */
  .workbench {{
    display: grid;
    grid-template-columns: 1.18fr 1fr;
    gap: 18px;
    margin-top: 16px;
  }}
  @media (max-width: 1060px) {{
    .workbench {{ grid-template-columns: 1fr; }}
  }}
  .viewport-shell {{
    border: 1px solid var(--line);
    border-radius: 8px;
    overflow: hidden;
    background: #111315;
    display: flex;
    flex-direction: column;
  }}
  .viewport-bar {{
    background: #1c1e21;
    border-bottom: 1px solid #2d3035;
    padding: 8px 10px;
    display: flex;
    flex-wrap: wrap;
    gap: 6px;
    align-items: center;
    justify-content: space-between;
  }}
  .viewport-bar select {{
    height: 27px;
    background: #282b30;
    color: #f5f5f7;
    border: 1px solid #3c4047;
    font-size: 11.5px;
    border-radius: 5px;
    padding: 0 8px;
  }}
  .vp-btn {{
    height: 27px;
    background: #282b30;
    color: #e5e5ea;
    border: 1px solid #3c4047;
    border-radius: 5px;
    padding: 0 9px;
    font-size: 11.5px;
    cursor: pointer;
  }}
  .vp-btn:hover {{ background: #34383e; color: #ffffff; }}
  .mol3d-stage {{
    position: relative;
    width: 100%;
    height: 470px;
    background: #111315;
  }}
  .mol3d-hud {{
    position: absolute;
    top: 10px;
    left: 12px;
    color: rgba(245, 245, 247, 0.85);
    font-family: var(--mono);
    font-size: 11px;
    background: rgba(28, 30, 33, 0.78);
    border: 1px solid rgba(255, 255, 255, 0.1);
    padding: 4px 9px;
    border-radius: 5px;
    pointer-events: none;
    z-index: 5;
  }}
  .mol3d-legend {{
    position: absolute;
    bottom: 10px;
    left: 12px;
    right: 12px;
    background: rgba(22, 24, 27, 0.84);
    border: 1px solid rgba(255, 255, 255, 0.09);
    color: #d2d2d7;
    padding: 6px 10px;
    border-radius: 5px;
    font-size: 11px;
    display: flex;
    flex-wrap: wrap;
    gap: 14px;
    pointer-events: none;
    z-index: 5;
  }}
  .mol3d-swatch {{
    display: inline-block;
    width: 8px;
    height: 8px;
    border-radius: 50%;
    margin-right: 5px;
    vertical-align: middle;
  }}
  .vp-footnote {{
    margin: 7px 0 0;
    font-size: 11.5px;
    color: var(--ink-tertiary);
  }}

  /* Residue specification sheet (Right column) */
  .spec-sheet {{
    border: 1px solid var(--line-light);
    border-radius: 8px;
    padding: 16px 18px;
    background: var(--surface);
  }}
  .aa-strip {{
    display: grid;
    grid-template-columns: repeat(19, minmax(24px, 1fr));
    gap: 3px;
    margin: 8px 0 14px;
  }}
  .aa-cell {{
    background: var(--surface-subtle);
    border: 1px solid var(--line);
    border-top-width: 2.5px;
    border-radius: 4px;
    height: 28px;
    font-family: var(--mono);
    font-size: 11.5px;
    font-weight: 600;
    color: var(--ink);
    cursor: pointer;
    padding: 0;
    transition: all 0.1s ease;
  }}
  .aa-cell:hover {{ background: #ededf2; }}
  .aa-cell.active {{
    background: var(--ink);
    color: #ffffff;
    border-color: var(--ink);
  }}
  .interp-box {{
    margin-top: 12px;
    padding: 10px 12px;
    background: var(--surface-subtle);
    border: 1px solid var(--line-light);
    border-radius: 6px;
    font-size: 12.5px;
    color: var(--ink-secondary);
    line-height: 1.5;
  }}

  /* Clean institutional data tables & subtle status dots */
  .table-wrap {{
    overflow-x: auto;
    max-height: 500px;
    overflow-y: auto;
    border: 1px solid var(--line-light);
    border-radius: 8px;
  }}
  table {{
    width: 100%;
    border-collapse: collapse;
    font-size: 12.5px;
    font-variant-numeric: tabular-nums;
  }}
  th, td {{
    padding: 8px 11px;
    border-bottom: 1px solid var(--line-light);
    text-align: left;
    white-space: nowrap;
  }}
  th {{
    background: var(--surface-subtle);
    position: sticky;
    top: 0;
    z-index: 2;
    font-size: 11px;
    font-weight: 600;
    text-transform: uppercase;
    letter-spacing: 0.035em;
    color: var(--ink-tertiary);
    border-bottom: 1px solid var(--line);
  }}
  tbody tr:hover td {{ background: #f8f9fb; }}
  code {{
    font-family: var(--mono);
    font-size: 12px;
    color: var(--ink);
  }}

  /* Restrained EMBL-style status tags (replaces loud pastel pills) */
  .tag {{
    display: inline-flex;
    align-items: center;
    gap: 5px;
    font-size: 11.5px;
    font-weight: 500;
    color: var(--ink);
    padding: 1px 7px;
    border-radius: 4px;
    border: 1px solid var(--line-light);
    background: var(--surface-subtle);
  }}
  .dot {{
    width: 6px;
    height: 6px;
    border-radius: 50%;
    display: inline-block;
    flex-shrink: 0;
  }}
  .dot-crimson {{ background: var(--crimson); }}
  .dot-amber {{ background: var(--amber); }}
  .dot-petrol {{ background: var(--embl-petrol); }}
  .dot-emerald {{ background: var(--emerald); }}
  .dot-slate {{ background: #6e6e73; }}

  .heatmap-grid {{ overflow-x: auto; border: 1px solid var(--line-light); border-radius: 8px; }}
  .hm-cell {{
    text-align: center;
    font-family: var(--mono);
    font-size: 11px;
    font-variant-numeric: tabular-nums;
    padding: 6px 4px;
    cursor: default;
  }}

  footer {{
    border-top: 1px solid var(--line-light);
    background: var(--surface);
    padding: 24px 0 32px;
    color: var(--ink-tertiary);
    font-size: 12px;
  }}
  .footer-inner {{
    display: flex;
    justify-content: space-between;
    flex-wrap: wrap;
    gap: 12px;
  }}
  .footer-inner a {{ color: var(--embl-petrol); text-decoration: none; }}
  .footer-inner a:hover {{ text-decoration: underline; }}

  /* Interactive Term Tooltips, Section Explainers & Reference Hub */
  [data-term] {{
    border-bottom: 1px dotted rgba(0, 112, 155, 0.58);
    cursor: help;
    transition: color 0.12s ease, border-color 0.12s ease;
  }}
  [data-term]:hover {{
    color: var(--embl-petrol);
    border-bottom-color: var(--embl-petrol);
  }}
  th [data-term], th[data-term] {{
    border-bottom: 1px dotted rgba(0, 112, 155, 0.65);
  }}
  .spec-item[data-term] {{
    border-bottom: none;
    border-radius: 6px;
    padding: 4px 10px 4px 14px;
    transition: background 0.12s ease;
  }}
  .spec-item[data-term]:first-child {{ padding-left: 6px; }}
  .spec-item[data-term]:hover {{ background: var(--surface-subtle); }}
  .info-badge {{
    display: inline-flex;
    align-items: center;
    justify-content: center;
    width: 14px;
    height: 14px;
    border-radius: 50%;
    font-family: var(--sans);
    font-size: 10px;
    font-weight: 650;
    color: var(--embl-petrol);
    background: rgba(0, 112, 155, 0.10);
    margin-left: 4px;
    vertical-align: middle;
    text-transform: none;
  }}
  .explainer-banner {{
    background: var(--surface);
    border: 1px solid var(--line-light);
    border-left: 3px solid var(--embl-petrol);
    border-radius: 8px;
    margin-bottom: 16px;
    overflow: hidden;
  }}
  .explainer-banner summary {{
    list-style: none;
    cursor: pointer;
    padding: 11px 16px;
    display: flex;
    justify-content: space-between;
    align-items: center;
    gap: 12px;
    font-size: 13px;
    color: var(--ink);
    user-select: none;
  }}
  .explainer-banner summary::-webkit-details-marker {{ display: none; }}
  .explainer-banner summary:hover {{ background: var(--surface-subtle); }}
  .exp-sum-left {{
    display: flex;
    align-items: center;
    gap: 9px;
    flex-wrap: wrap;
  }}
  .exp-kicker {{
    font-family: var(--mono);
    font-size: 10.5px;
    font-weight: 600;
    text-transform: uppercase;
    letter-spacing: 0.04em;
    color: var(--embl-petrol);
    background: rgba(0, 112, 155, 0.08);
    padding: 2px 7px;
    border-radius: 4px;
  }}
  .exp-toggle-lbl {{
    font-size: 12px;
    font-weight: 500;
    color: var(--embl-petrol);
    white-space: nowrap;
  }}
  .explainer-banner[open] .exp-toggle-lbl::after {{ content: " ▴"; }}
  .explainer-banner:not([open]) .exp-toggle-lbl::after {{ content: " ▾"; }}
  .exp-body {{
    padding: 14px 18px 16px;
    border-top: 1px solid var(--line-light);
    background: var(--surface-subtle);
    display: grid;
    grid-template-columns: repeat(3, 1fr);
    gap: 16px;
    font-size: 12.5px;
    color: var(--ink-secondary);
    line-height: 1.52;
  }}
  @media (max-width: 960px) {{
    .exp-body {{ grid-template-columns: 1fr; }}
  }}
  .exp-col h4 {{
    margin: 0 0 5px;
    font-size: 11.5px;
    font-weight: 600;
    text-transform: uppercase;
    letter-spacing: 0.035em;
    color: var(--ink);
  }}
  .exp-col p {{ margin: 0; }}
  .exp-chips {{
    grid-column: 1 / -1;
    display: flex;
    flex-wrap: wrap;
    align-items: center;
    gap: 6px;
    padding-top: 8px;
    border-top: 1px solid var(--line-light);
    font-size: 11.5px;
    color: var(--ink-tertiary);
  }}
  .term-chip {{
    background: var(--surface);
    border: 1px solid var(--line);
    border-radius: 5px;
    padding: 2px 8px;
    font-family: var(--mono);
    font-size: 11px;
    color: var(--ink);
    cursor: help;
  }}
  .term-chip:hover {{
    border-color: var(--embl-petrol);
    color: var(--embl-petrol);
  }}

  /* Floating Term Specification Popover */
  .term-popover {{
    position: fixed;
    z-index: 9999;
    width: 390px;
    max-width: calc(100vw - 24px);
    background: #1d1d1f;
    color: #f5f5f7;
    border: 1px solid #3a3a3c;
    border-radius: 10px;
    padding: 13px 15px;
    box-shadow: 0 14px 36px rgba(0, 0, 0, 0.28);
    font-size: 12.5px;
    line-height: 1.48;
    display: none;
  }}
  .term-popover.visible {{ display: block; }}
  .tp-head {{
    display: flex;
    justify-content: space-between;
    align-items: flex-start;
    gap: 8px;
    margin-bottom: 6px;
  }}
  .tp-cat {{
    font-family: var(--mono);
    font-size: 10px;
    text-transform: uppercase;
    letter-spacing: 0.05em;
    color: #5ac8fa;
    display: block;
    margin-bottom: 2px;
  }}
  .tp-title {{
    font-size: 13.5px;
    font-weight: 600;
    color: #ffffff;
  }}
  .tp-close {{
    background: transparent;
    border: none;
    color: #98989d;
    font-size: 15px;
    cursor: pointer;
    padding: 0 2px;
    line-height: 1;
  }}
  .tp-close:hover {{ color: #ffffff; }}
  .tp-meaning {{ color: #e5e5ea; margin-bottom: 8px; }}
  .tp-formula {{
    font-family: var(--mono);
    font-size: 11.2px;
    background: rgba(255, 255, 255, 0.07);
    border: 1px solid rgba(255, 255, 255, 0.11);
    border-radius: 6px;
    padding: 6px 9px;
    color: #f5f5f7;
    margin-bottom: 8px;
    overflow-x: auto;
    white-space: pre-wrap;
    word-break: break-word;
  }}
  .tp-thresh {{
    font-size: 11.5px;
    color: #AEAEB2;
    margin-bottom: 8px;
  }}
  .tp-link {{
    background: transparent;
    border: none;
    padding: 0;
    color: #5ac8fa;
    font-size: 11.5px;
    font-weight: 500;
    cursor: pointer;
    font-family: var(--sans);
  }}
  .tp-link:hover {{ text-decoration: underline; }}

  /* Reference Hub (Tab 7: Pathway, Formulas, Datasets & Glossary) */
  .doc-jumpbar {{
    display: flex;
    flex-wrap: wrap;
    gap: 8px;
    margin-bottom: 16px;
  }}
  .pathway-flow {{
    display: grid;
    grid-template-columns: repeat(5, 1fr);
    gap: 10px;
    margin-bottom: 16px;
  }}
  @media (max-width: 1060px) {{
    .pathway-flow {{ grid-template-columns: repeat(2, 1fr); }}
  }}
  @media (max-width: 600px) {{
    .pathway-flow {{ grid-template-columns: 1fr; }}
  }}
  .pf-step {{
    border: 1px solid var(--line-light);
    border-top: 3px solid var(--embl-petrol);
    border-radius: 8px;
    padding: 12px 13px;
    background: var(--surface-subtle);
  }}
  .pf-num {{
    font-family: var(--mono);
    font-size: 10.5px;
    color: var(--embl-petrol);
    font-weight: 600;
    text-transform: uppercase;
  }}
  .pf-title {{
    font-size: 13px;
    font-weight: 600;
    margin: 3px 0 5px;
    color: var(--ink);
  }}
  .pf-desc {{
    font-size: 12px;
    color: var(--ink-secondary);
    line-height: 1.45;
  }}
  .formula-grid {{
    display: grid;
    grid-template-columns: repeat(2, 1fr);
    gap: 14px;
  }}
  @media (max-width: 1024px) {{
    .formula-grid {{ grid-template-columns: 1fr; }}
  }}
  .formula-card {{
    border: 1px solid var(--line-light);
    border-radius: 8px;
    padding: 14px 16px;
    background: var(--surface);
    display: flex;
    flex-direction: column;
    gap: 7px;
  }}
  .formula-card:target {{
    border-color: var(--embl-petrol);
    box-shadow: 0 0 0 3px rgba(0, 112, 155, 0.14);
  }}
  .fc-top {{
    display: flex;
    justify-content: space-between;
    align-items: center;
    gap: 8px;
  }}
  .fc-id {{
    font-family: var(--mono);
    font-size: 11px;
    font-weight: 600;
    color: var(--embl-petrol);
    background: rgba(0, 112, 155, 0.08);
    padding: 2px 7px;
    border-radius: 4px;
  }}
  .fc-cat {{
    font-size: 11px;
    color: var(--ink-tertiary);
    font-family: var(--mono);
  }}
  .fc-title {{
    font-size: 13.5px;
    font-weight: 600;
    color: var(--ink);
  }}
  .fc-eq {{
    font-family: var(--mono);
    font-size: 11.8px;
    background: var(--surface-subtle);
    border: 1px solid var(--line-light);
    border-radius: 6px;
    padding: 8px 10px;
    color: var(--ink);
    overflow-x: auto;
    white-space: pre-wrap;
    word-break: break-word;
  }}
  .fc-vars {{
    font-size: 11.8px;
    color: var(--ink-tertiary);
  }}
  .fc-meaning {{
    font-size: 12.3px;
    color: var(--ink-secondary);
    line-height: 1.48;
  }}
  .fc-script {{
    font-family: var(--mono);
    font-size: 11px;
    color: var(--embl-teal);
  }}
  .glossary-grid {{
    display: grid;
    grid-template-columns: repeat(2, 1fr);
    gap: 12px;
  }}
  @media (max-width: 960px) {{
    .glossary-grid {{ grid-template-columns: 1fr; }}
  }}
  .gloss-card {{
    border: 1px solid var(--line-light);
    border-radius: 8px;
    padding: 12px 14px;
    background: var(--surface-subtle);
  }}
</style>
</head>
<body>

<div class="global-bar">
  <div class="container">
    <div class="global-brand">
      <span class="embl-mark"></span>
      <span>miRBioGenAtlas &middot; Multi-Consortium Cancer Genomics &amp; Structural Reference</span>
    </div>
    <div class="global-links">
      <a href="#tab-docs" onclick="event.preventDefault();openDocsSection('sec-pathway')" style="color:#5ac8fa;font-weight:600;">Guide, Formulas &amp; Glossary</a>
      <a href="https://doi.org/10.5281/zenodo.24455067" target="_blank" rel="noopener">DOI:10.5281/zenodo.24455067</a>
      <a href="https://github.com/siddugayam/miRBioGenAtlas" target="_blank" rel="noopener">GitHub</a>
      <span>CC BY 4.0 / MIT</span>
    </div>
  </div>
</div>

<header class="masthead">
  <div class="container">
    <div class="masthead-top">
      <div class="title-block">
        <h1>miRBioGenAtlas <span class="version-tag">v2.0.2</span></h1>
        <p>Sample-, variant-, and atomic-level reference of the human microRNA biogenesis machinery across TCGA, MSK-IMPACT, GDC, PDB, AlphaFold, ClinVar, gnomAD, GTEx, CPTAC, and DepMap. <span style="color:var(--embl-petrol);cursor:pointer;font-weight:500;" onclick="openDocsSection('sec-pathway')">Hover or click any underlined metric for definitions &amp; formulas &rarr;</span></p>
      </div>
      <div class="omni-search">
        <input type="text" id="global-search" placeholder="Search gene, variant, term or formula (e.g. DICER1 E1705K, LOEUF, 5p:3p)..." onkeydown="if(event.key==='Enter') runGlobalSearch(this.value)">
        <span class="omni-kbd">/</span>
      </div>
    </div>

    <div class="spec-strip">
      <div class="spec-item" data-term="kpi_tumours">
        <div class="spec-val">21,140 <span class="info-badge">?</span></div>
        <div class="spec-lbl">Sequenced Tumours (TCGA &amp; MSK)</div>
      </div>
      <div class="spec-item" data-term="kpi_proteins">
        <div class="spec-val">15 Proteins <span class="info-badge">?</span></div>
        <div class="spec-lbl">11,260 Codons &middot; 89,132 Heavy Atoms</div>
      </div>
      <div class="spec-item" data-term="kpi_sat">
        <div class="spec-val">213,940 <span class="info-badge">?</span></div>
        <div class="spec-lbl">Missense Substitutions &amp; 3D Rotamers</div>
      </div>
      <div class="spec-item" data-term="kpi_pdb">
        <div class="spec-val">16 Complexes <span class="info-badge">?</span></div>
        <div class="spec-lbl">Experimental Cryo-EM &amp; X-ray PDBs</div>
      </div>
      <div class="spec-item" data-term="kpi_mirna">
        <div class="spec-val">9,578 <span class="info-badge">?</span></div>
        <div class="spec-lbl">Mature miRNA-seq &amp; 13.37M Pairs</div>
      </div>
      <div class="spec-item" data-term="kpi_auc">
        <div class="spec-val">0.978 AUC <span class="info-badge">?</span></div>
        <div class="spec-lbl">ClinVar Benchmark (16,569 Records)</div>
      </div>
    </div>
  </div>
</header>

<nav class="subnav-wrap" aria-label="Primary sections">
  <div class="container">
    <div class="tabs">
      <button class="tab-btn active" id="btn-tab-variants" onclick="switchTab('tab-variants', this)">Structure &amp; Saturation</button>
      <button class="tab-btn" id="btn-tab-heatmap" onclick="switchTab('tab-heatmap', this)">Pan-Cancer Matrix</button>
      <button class="tab-btn" id="btn-tab-mirna" onclick="switchTab('tab-mirna', this)">Small-RNA &amp; Targets</button>
      <button class="tab-btn" id="btn-tab-clinvar" onclick="switchTab('tab-clinvar', this)">Clinical &amp; Germline</button>
      <button class="tab-btn" id="btn-tab-gtex" onclick="switchTab('tab-gtex', this)">Tissue &amp; Tractability</button>
      <button class="tab-btn" id="btn-tab-consortia" onclick="switchTab('tab-consortia', this)">Proteomics, CRISPR &amp; Survival</button>
      <button class="tab-btn" id="btn-tab-docs" onclick="switchTab('tab-docs', this)">Methods, Formulas &amp; Datasets</button>
    </div>
  </div>
</nav>

<main class="container workspace">

  <!-- SECTION 1: STRUCTURE, LOLLIPOP & 213,940 SATURATION WORKBENCH -->
  <section id="tab-variants" class="panel active">
    <details class="explainer-banner">
      <summary>
        <div class="exp-sum-left">
          <span class="exp-kicker">Section Guide</span>
          <span><strong>How to interpret 3D mutant side-chain rotamers, active-pocket distances &amp; saturation scores</strong> across all 213,940 substitutions</span>
        </div>
        <span class="exp-toggle-lbl">Expand explainer</span>
      </summary>
      <div class="exp-body">
        <div class="exp-col">
          <h4>1. Biological &amp; Structural Question</h4>
          <p>Large biogenesis ribonucleases (<em>DICER1</em> 1,922 aa; <em>DROSHA</em> 1,374 aa) accumulate both true driver mutations and neutral hypermutator passengers. This workbench maps every codon onto 16 human Cryo-EM/X-ray complexes and AlphaFold v6 coordinates to resolve atomic mechanism.</p>
        </div>
        <div class="exp-col">
          <h4>2. Three Biophysical Regimes</h4>
          <p>By crossing sequence constraint (<strong>AlphaMissense &ge; 0.564</strong> or <strong>ESM-2 &le; &minus;5.0</strong>) with 3D inverse folding (<strong>ESM-IF1 &Delta;logP</strong>), variants separate into <strong>Fold-preserving active/interface (&Delta;logP &gt; &minus;1.5)</strong>, <strong>Fold-destabilising (&le; &minus;1.5)</strong>, and <strong>Tolerated neutral</strong>.</p>
        </div>
        <div class="exp-col">
          <h4>3. Reading the Lollipop &amp; 3D Stage</h4>
          <p>Click any lollipop stem, 19-amino-acid button, or 3D residue. Cyan sticks show the wild-type side chain; crimson sticks show the lowest-clash mutant rotamer; amber dashed lines mark contacts (&le;3.8 &Aring;) and red cylinders flag van der Waals clashes (&lt;2.4 &Aring;).</p>
        </div>
        <div class="exp-chips">
          <span>Click or hover any term:</span>
          <span class="term-chip" data-term="alphamissense">AlphaMissense</span>
          <span class="term-chip" data-term="esm2_llr">ESM-2 3B &Delta;&ell;</span>
          <span class="term-chip" data-term="esm1v_llr">ESM-1v &Delta;&ell;</span>
          <span class="term-chip" data-term="esm_if1">ESM-IF1 3D &Delta;logP</span>
          <span class="term-chip" data-term="mech_regime">3D Mechanism Regimes</span>
          <span class="term-chip" data-term="rotamer_clash">Rotamer &amp; Steric Clash (&lt;2.4&Aring;)</span>
          <span class="term-chip" data-term="grantham">Grantham Distance</span>
          <span class="term-chip" data-term="plddt">AlphaFold pLDDT</span>
        </div>
      </div>
    </details>

    <div class="surface" id="struct3d-card">
      <div class="section-head">
        <h2>3D Structural, Side-Chain Rotamer &amp; Saturation Mutagenesis Workbench</h2>
        <span class="section-meta" id="sat-status">Initialising structural coordinates...</span>
      </div>

      <div class="toolbar">
        <div class="field-group">
          <span class="field-label">Protein</span>
          <select id="lol-gene" onchange="renderLollipopAndLoadSat()"></select>
        </div>
        <div class="field-group">
          <span class="field-label">Residue / Variant</span>
          <input type="text" id="sat-query" placeholder="e.g. E1705K, 1705" oninput="inspectSaturationCodon()" style="width:148px;font-family:var(--mono);font-weight:500;">
        </div>
      </div>

      <div class="ref-strip">
        <span>Reference loci:</span>
        <button class="ref-link" onclick="loadPresetMutation('DICER1','E1705K')">DICER1 E1705K</button>
        <button class="ref-link" onclick="loadPresetMutation('DICER1','D1810Y')">DICER1 D1810Y</button>
        <button class="ref-link" onclick="loadPresetMutation('DROSHA','E1147K')">DROSHA E1147K</button>
        <button class="ref-link" onclick="loadPresetMutation('TP53','R175H')">TP53 R175H</button>
        <button class="ref-link" onclick="loadPresetMutation('TP53','R248Q')">TP53 R248Q</button>
        <button class="ref-link" onclick="loadPresetMutation('SMAD4','R361H')">SMAD4 R361H</button>
        <button class="ref-link" onclick="loadPresetMutation('AGO2','M364T')">AGO2 M364T</button>
        <button class="ref-link" onclick="loadPresetMutation('DICER1','R944Q')">DICER1 R944Q (neutral)</button>
      </div>

      <div id="lollipop-svg-wrap" class="svg-box"></div>

      <div class="workbench">
        <div>
          <div class="viewport-shell">
            <div class="viewport-bar">
              <div style="display:flex;flex-wrap:wrap;gap:6px;">
                <select id="mol3d-mode" onchange="refresh3DScene()" aria-label="Structure comparison mode">
                  <option value="overlay">Overlay: Wild-Type + Mutant</option>
                  <option value="mutant">Mutant Side-Chain Only</option>
                  <option value="wt">Wild-Type Only</option>
                  <option value="rcsb">Experimental PDB Complex</option>
                </select>
                <select id="mol3d-color" onchange="refresh3DScene()" aria-label="Backbone coloring">
                  <option value="plddt">Color: AlphaFold pLDDT</option>
                  <option value="am">Color: AlphaMissense</option>
                  <option value="pocket">Color: 3D Active Pocket</option>
                </select>
              </div>
              <div style="display:flex;flex-wrap:wrap;gap:5px;">
                <button class="vp-btn" onclick="zoomToMutatedResidue()">Focus 8&Aring;</button>
                <button class="vp-btn" onclick="resetFullProtein3D()">Full Chain</button>
                <button class="vp-btn" id="btn-nbrs" onclick="toggleNeighbors3D()">Shell: On</button>
                <button class="vp-btn" id="btn-spin" onclick="toggleSpin3D()">Spin</button>
                <button class="vp-btn" onclick="downloadMutatedPDB()">PDB</button>
              </div>
            </div>
            <div class="mol3d-stage">
              <div class="mol3d-hud" id="mol3d-hud">DICER1 &middot; Q9UPY3</div>
              <div id="mol3d-canvas" style="width:100%;height:100%;position:relative;"></div>
              <div class="mol3d-legend">
                <span><i class="mol3d-swatch" style="background:#0095c8;"></i>Wild-type side chain</span>
                <span data-term="rotamer_clash"><i class="mol3d-swatch" style="background:#e11d48;"></i>Mutant 3D rotamer</span>
                <span data-term="burial_8a"><i class="mol3d-swatch" style="background:#94a3b8;"></i>8&Aring; shell</span>
                <span data-term="dist_na"><i class="mol3d-swatch" style="background:#f59e0b;"></i>Contact (&le;3.8&Aring;)</span>
                <span data-term="rotamer_clash"><i class="mol3d-swatch" style="background:#ff3b30;"></i>Steric clash (&lt;2.4&Aring;)</span>
              </div>
            </div>
          </div>
          <p class="vp-footnote">Coordinates: AlphaFold DB v6 canonical heavy-atom models and experimental RCSB PDB macromolecular complexes. Mutant side chains represent the lowest-clash AMBER-relaxed rotamer projected onto the local N&ndash;C&alpha;&ndash;C frame.</p>
        </div>

        <div class="spec-sheet" id="prop-delta-panel">
          <div style="color:var(--ink-tertiary);font-size:13px;">Loading residue specification...</div>
        </div>
      </div>

      <div id="sat-codon-detail" style="margin-top:18px;"></div>
    </div>

    <div class="surface">
      <div class="section-head">
        <h2>Observed TCGA PanCanAtlas Somatic Coding Variants</h2>
        <span class="section-meta" id="var-count">6,292 coding calls across 10,195 tumours</span>
      </div>
      <div class="toolbar">
        <div class="field-group">
          <span class="field-label">Gene</span>
          <select id="var-gene" onchange="renderVariants()">
            <option value="ALL">All 15 genes</option>
          </select>
        </div>
        <div class="field-group">
          <span class="field-label">Tier</span>
          <select id="var-tier" onchange="renderVariants()">
            <option value="ALL">All functional tiers</option>
            <option value="3D_Catalytic_or_Substrate_Pocket_Driver">3D pocket (&le;6 &Aring;)</option>
            <option value="Pathogenic_Fold_or_Interface_Missense">Pathogenic fold / interface</option>
            <option value="Truncating_LoF">Truncating LoF</option>
            <option value="Ambiguous_Missense_VUS">Ambiguous VUS</option>
            <option value="Likely_Benign_Passenger_Missense">Likely benign VUS</option>
          </select>
        </div>
        <div class="field-group">
          <span class="field-label">Mechanism</span>
          <select id="var-mech" onchange="renderVariants()">
            <option value="ALL">All mechanisms</option>
            <option value="Fold_Preserving_Catalytic_or_Interface">Fold-preserving active / interface</option>
            <option value="Thermodynamic_Fold_Destabilizing">Fold-destabilising</option>
            <option value="Tolerated_Neutral">Tolerated neutral</option>
          </select>
        </div>
        <div class="field-group">
          <span class="field-label">Filter</span>
          <input type="text" id="var-q" placeholder="HGVS, study, ClinVar..." oninput="renderVariants()">
        </div>
        <button class="btn-primary" onclick="exportFilteredVariantsCSV()">Export CSV</button>
      </div>
      <div class="table-wrap">
        <table id="var-table">
          <thead>
            <tr>
              <th>3D</th>
              <th>Gene</th>
              <th>Variant</th>
              <th data-term="kpi_tumours">TCGA</th>
              <th data-term="kpi_auc">ClinVar</th>
              <th data-term="struct_tier">Structural Tier</th>
              <th data-term="mech_regime">3D Mechanism</th>
              <th data-term="kpi_pdb">PDB</th>
              <th data-term="dist_na">RNA/DNA (&Aring;)</th>
              <th data-term="dist_metal">Metal (&Aring;)</th>
              <th data-term="alphamissense">AlphaMissense</th>
              <th data-term="esm2_llr">ESM-2 3B</th>
              <th data-term="esm1v_llr">ESM-1v</th>
              <th data-term="esm_if1">ESM-IF1 &Delta;logP</th>
              <th>Studies</th>
            </tr>
          </thead>
          <tbody></tbody>
        </table>
      </div>
    </div>
  </section>

  <!-- SECTION 2: PAN-CANCER MATRIX -->
  <section id="tab-heatmap" class="panel">
    <details class="explainer-banner">
      <summary>
        <div class="exp-sum-left">
          <span class="exp-kicker">Section Guide</span>
          <span><strong>Why strict vs broad alteration rates differ dramatically</strong>: whole-arm aneuploidy (&plusmn;1), hypermutators &amp; syntenic gene clusters</span>
        </div>
        <span class="exp-toggle-lbl">Expand explainer</span>
      </summary>
      <div class="exp-body">
        <div class="exp-col">
          <h4>1. Removing Whole-Arm Aneuploidy Noise</h4>
          <p>Standard portal queries include shallow <strong>GISTIC2 &plusmn;1</strong> gains/losses, which account for <strong>81%</strong> of raw calls due to broad chromosome <code>8q</code>, <code>1p</code>, <code>18q</code>, and <code>17p</code> aneuploidy. Our primary <strong>High-Confidence (Strict)</strong> metric restricts CNAs to focal <strong>GISTIC2 &plusmn;2</strong> plus coding mutations.</p>
        </div>
        <div class="exp-col">
          <h4>2. Exact Denominators &amp; Union Identity</h4>
          <p>Every study uses the exact matched intersection <code>S_seq &cap; S_CNA</code> (<em>N</em> = 10,195) and satisfies <code>N_any = N_mut-only + N_CNA-only + N_both</code> (verified with 0.0000 discrepancy across 20 cBioPortal API spot-checks).</p>
        </div>
        <div class="exp-col">
          <h4>3. Synteny-Aware Patient Pathway Tiers</h4>
          <p>Because <em>AGO1/3/4</em> form a tandem cluster at <code>1p34.3</code> and <em>SMAD2/4</em> colocalise at <code>18q21.1</code>, patient-level burden collapses co-deleted clusters and stratifies genes into <strong>Core-10</strong>, <strong>Extended-12</strong>, and <strong>Panel-14</strong> (excluding <em>TP53</em>).</p>
        </div>
        <div class="exp-chips">
          <span>Key definitions:</span>
          <span class="term-chip" data-term="strict_hc">High-Confidence (Strict)</span>
          <span class="term-chip" data-term="focal_cna">Focal CNA (&plusmn;2)</span>
          <span class="term-chip" data-term="gistic2">GISTIC2 (&minus;2..+2)</span>
          <span class="term-chip" data-term="hypermutator">Non-Hypermutated</span>
          <span class="term-chip" data-term="broad_any">Broad Alteration</span>
          <span class="term-chip" data-term="synteny">1p34.3 / 18q21.1 Synteny</span>
          <span class="term-chip" data-term="tier_core10">Core-10</span>
          <span class="term-chip" data-term="tier_ext12">Extended-12</span>
          <span class="term-chip" data-term="tier_panel14">Panel-14</span>
        </div>
      </div>
    </details>

    <div class="surface">
      <div class="section-head">
        <h2>Pan-Cancer Alteration Frequency Matrix (15 Genes &times; 32 TCGA Studies)</h2>
        <span class="section-meta">Matched sequencing &amp; GISTIC2 cohort (N = 10,195)</span>
      </div>
      <div class="toolbar">
        <div class="field-group">
          <span class="field-label">Metric</span>
          <select id="hm-metric" onchange="renderHeatmap()">
            <option value="hc">High-confidence alteration (nonsynonymous mutation or GISTIC2 &plusmn;2)</option>
            <option value="fcna">Focal high-level CNA only (GISTIC2 &plusmn;2)</option>
            <option value="mut">Somatic nonsynonymous mutation (all samples)</option>
            <option value="nhm">Somatic nonsynonymous mutation (non-hypermutated, n = 8,810)</option>
            <option value="any">Broad alteration (mutation or any GISTIC2 &ne; 0)</option>
          </select>
        </div>
        <button class="btn-primary" onclick="exportLongCSV()">Download 480-Row CSV</button>
        <button class="btn-quiet" onclick="openDocsSection('sec-tcga-codes')">Decode 32 TCGA Study Acronyms &rarr;</button>
      </div>
      <div id="hm-container" class="heatmap-grid"></div>
      <p class="vp-footnote">Cells report exact study proportions in [0.00, 1.00]. Hover any cell for denominator N, somatic mutation count, +2 amplification count, and &minus;2 deep-deletion count. Studies marked with * have N &lt; 100.</p>
    </div>

    <div class="surface">
      <div class="section-head">
        <h2>Synteny-Aware Patient Pathway Burden Ranking Across 32 TCGA Studies (Excluding TP53)</h2>
        <span class="section-meta">Boolean patient union across nested gene-set tiers</span>
      </div>
      <div class="table-wrap" style="max-height:360px;">
        <table id="cohort-burden-table">
          <thead>
            <tr>
              <th>Rank</th>
              <th>Study</th>
              <th>Matched N</th>
              <th data-term="tier_core10">Core-10 Strict % (n)</th>
              <th data-term="tier_ext12">Extended-12 Strict %</th>
              <th data-term="tier_panel14">Panel-14 Strict %</th>
              <th data-term="broad_any">Core-10 Broad % (incl. &plusmn;1)</th>
              <th data-term="synteny">1p34.3 AGO1/3/4 Co-CNA %</th>
            </tr>
          </thead>
          <tbody></tbody>
        </table>
      </div>
    </div>
  </section>

  <!-- SECTION 3: MIRNA-SEQ & TARGET COUPLINGS -->
  <section id="tab-mirna" class="panel">
    <details class="explainer-banner">
      <summary>
        <div class="exp-sum-left">
          <span class="exp-kicker">Section Guide</span>
          <span><strong>How endogenous mature miRNA-seq proves enzymatic cleavage asymmetry</strong>: 5p:3p duplex ratios, AGO2-Slicer control &amp; target mRNA derepression</span>
        </div>
        <span class="exp-toggle-lbl">Expand explainer</span>
      </summary>
      <div class="exp-body">
        <div class="exp-col">
          <h4>1. Paired 5p:3p Hairpin Enzymology</h4>
          <p>Inside cytoplasmic <em>DICER1</em>, the <strong>RNase IIIb</strong> domain cleaves the <code>-5p</code> hairpin strand while <strong>RNase IIIa</strong> cleaves the <code>-3p</code> strand. Pairing <code>log2((RPM_5p+1)/(RPM_3p+1))</code> within 192 co-transcribed duplexes cancels transcriptional noise and isolates catalytic RNase IIIb failure (<em>p</em> = 1.46 &times; 10<sup>&minus;15</sup>).</p>
        </div>
        <div class="exp-col">
          <h4>2. Endogenous AGO2-Slicer Control</h4>
          <p>Bicistronic <code>hsa-miR-451a</code> has a short 17-bp stem that completely bypasses <em>DICER1</em> and is cleaved directly by <em>AGO2</em> PIWI (&lsquo;Slicer&rsquo;), whereas co-transcribed <code>hsa-miR-144-5p</code> requires <em>DICER1</em>—providing an internal orthogonal control (<em>p</em> = 0.967).</p>
        </div>
        <div class="exp-col">
          <h4>3. Downstream 13.37M mRNA Couplings</h4>
          <p>Across 9,353 dual-profiled tumours (662 miRNAs &times; 20,196 mRNAs), selective collapse of <code>5p</code> <em>let-7</em> strands in <em>DICER1</em> pocket mutants derepresses oncofetal <em>let-7</em> target mRNAs (<em>p</em> = 6.53 &times; 10<sup>&minus;5</sup>) while sparing <code>3p</code>-dominant <em>miR-200</em> EMT targets.</p>
        </div>
        <div class="exp-chips">
          <span>Key definitions:</span>
          <span class="term-chip" data-term="duplex_5p3p">Paired 5p:3p Duplex Ratio</span>
          <span class="term-chip" data-term="frac_5p">Linear 5p Fraction</span>
          <span class="term-chip" data-term="mir451_slicer">miR-451a / miR-144 Slicer Control</span>
          <span class="term-chip" data-term="let7_target_z">let-7 Target mRNA Z</span>
          <span class="term-chip" data-term="pearson_couplings">13.37M Pearson Couplings</span>
          <span class="term-chip" data-term="rank_shift_p">Mann&ndash;Whitney U Shift</span>
        </div>
      </div>
    </details>

    <div class="surface">
      <div class="section-head">
        <h2>Mature miRNA-seq Enzymatic &amp; Target mRNA Benchmarks</h2>
        <span class="section-meta">9,578 miRNA-seq tumours &middot; 192 canonical 5p/3p duplexes &middot; 9,353 dual-profiled tumours</span>
      </div>
      <div class="table-wrap">
        <table id="mir-bench-table">
          <thead>
            <tr>
              <th>Benchmark</th>
              <th data-term="duplex_5p3p">Readout Metric</th>
              <th>Group A (n)</th>
              <th>Median A</th>
              <th>Group B (n)</th>
              <th>Median B</th>
              <th data-term="rank_shift_p">Median Shift</th>
              <th data-term="rank_shift_p">Two-Sided p</th>
            </tr>
          </thead>
          <tbody></tbody>
        </table>
      </div>
    </div>
    <div class="surface">
      <div class="section-head">
        <h2>DICER1-Altered TCGA Tumours: 5p:3p Strand Balance &amp; let-7 Target Derepression</h2>
        <button class="btn-primary" onclick="exportDicerMirnaCSV()">Export CSV</button>
      </div>
      <div class="table-wrap" style="max-height:340px;">
        <table id="dicer-mir-table">
          <thead>
            <tr>
              <th>Sample ID</th>
              <th>Study</th>
              <th data-term="struct_tier">DICER1 Class</th>
              <th data-term="duplex_5p3p">Paired Duplex log2((5p+1)/(3p+1))</th>
              <th data-term="frac_5p">Linear 5p Fraction</th>
              <th data-term="let7_target_z">let-7 Target mRNA Z</th>
            </tr>
          </thead>
          <tbody></tbody>
        </table>
      </div>
    </div>
    <div class="surface">
      <div class="section-head">
        <h2>Top Pan-Cancer Repressive miRNA &times; mRNA Couplings</h2>
        <span class="section-meta">13,369,752 transcriptome-wide pairs evaluated across N = 9,353 dual-profiled tumours</span>
      </div>
      <div class="table-wrap" style="max-height:320px;">
        <table id="top-axes-table">
          <thead>
            <tr>
              <th>Rank</th>
              <th>Mature miRNA</th>
              <th>Target mRNA</th>
              <th>Dual Tumours (n)</th>
              <th data-term="pearson_couplings">Pearson r (log2)</th>
              <th data-term="pearson_couplings">t-Statistic</th>
            </tr>
          </thead>
          <tbody></tbody>
        </table>
      </div>
    </div>
  </section>

  <!-- SECTION 4: CLINVAR, GNOMAD & MSK-IMPACT -->
  <section id="tab-clinvar" class="panel">
    <details class="explainer-banner">
      <summary>
        <div class="exp-sum-left">
          <span class="exp-kicker">Section Guide</span>
          <span><strong>Non-circular clinical calibration &amp; independent clinical cohort replication</strong>: NCBI ClinVar (16,569 records) &amp; MSK-IMPACT (N = 10,945)</span>
        </div>
        <span class="exp-toggle-lbl">Expand explainer</span>
      </summary>
      <div class="exp-body">
        <div class="exp-col">
          <h4>1. Non-Circular ClinVar Benchmarking</h4>
          <p>None of our four sequence/3D predictors (<code>AlphaMissense</code>, <code>ESM-2 3B</code>, <code>ESM-1v</code>, <code>ESM-IF1</code>) were trained on ClinVar clinical labels. Evaluating them against 612 expert-curated Pathogenic/LP vs Benign/LB missense variants yields an unbiased clinical calibration (<code>ROC-AUC = 0.9784</code>).</p>
        </div>
        <div class="exp-col">
          <h4>2. Reading ROC-AUC vs PR-AUC</h4>
          <p><strong>ROC-AUC</strong> (0.5 random &rarr; 1.0 perfect; Hanley&ndash;McNeil 95% CI) measures rank separation; <strong>PR-AUC</strong> measures precision-recall performance and must be compared against the positive-class prevalence baseline shown in parentheses.</p>
        </div>
        <div class="exp-col">
          <h4>3. MSK-IMPACT Replication (N = 10,945)</h4>
          <p>Replicates somatic mutation, non-hypermutated, truncating LoF, and 3D RNase IIIb hotspot frequencies in an independent prospective clinical sequencing cohort spanning primary and metastatic tumours.</p>
        </div>
        <div class="exp-chips">
          <span>Key definitions:</span>
          <span class="term-chip" data-term="roc_pr_auc">ROC-AUC &amp; PR-AUC</span>
          <span class="term-chip" data-term="alphamissense">AlphaMissense</span>
          <span class="term-chip" data-term="esm2_llr">ESM-2 3B</span>
          <span class="term-chip" data-term="esm_if1">ESM-IF1 3D</span>
          <span class="term-chip" data-term="hypermutator">Non-Hypermutated %</span>
        </div>
      </div>
    </details>

    <div class="surface">
      <div class="section-head">
        <h2>NCBI ClinVar Pathogenic vs Benign Missense Discrimination</h2>
        <span class="section-meta">Non-circular human clinical assertions joined to the 213,940-variant saturation atlas</span>
      </div>
      <div class="table-wrap">
        <table id="clinvar-roc-table">
          <thead>
            <tr>
              <th>Benchmark Subset</th>
              <th>Predictor</th>
              <th>Pathogenic (n)</th>
              <th>Benign (n)</th>
              <th data-term="roc_pr_auc">ROC-AUC (95% CI)</th>
              <th data-term="roc_pr_auc">PR-AUC (Baseline)</th>
              <th data-term="rank_shift_p">Mann&ndash;Whitney p</th>
            </tr>
          </thead>
          <tbody></tbody>
        </table>
      </div>
    </div>
    <div class="surface">
      <div class="section-head">
        <h2>Prospective Clinical Sequencing Replication: MSK-IMPACT (N = 10,945) vs TCGA (N = 10,195)</h2>
        <span class="section-meta">21,140 combined tumours</span>
      </div>
      <div class="table-wrap">
        <table id="msk-rep-table">
          <thead>
            <tr>
              <th>Gene</th>
              <th>MSK Profiled N</th>
              <th>MSK Mut % (n)</th>
              <th>TCGA Mut %</th>
              <th data-term="hypermutator">MSK Non-Hypermut %</th>
              <th>MSK Truncating % (n)</th>
              <th data-term="struct_tier">MSK 3D Pocket / Hotspot % (n)</th>
              <th>Primary vs Metastasis %</th>
            </tr>
          </thead>
          <tbody></tbody>
        </table>
      </div>
    </div>
    <div class="surface">
      <div class="section-head">
        <h2>MSK-IMPACT Cancer-Type Distribution of Somatic DICER1 Variants</h2>
        <span class="section-meta">RNase IIIb catalytic hotspots vs truncating variants</span>
      </div>
      <div class="table-wrap" style="max-height:300px;">
        <table id="msk-dicer-table">
          <thead>
            <tr>
              <th>Clinical Cancer Type</th>
              <th>Samples (N)</th>
              <th>DICER1 Mutated (n, %)</th>
              <th data-term="struct_tier">RNase IIIb Hotspot (n, %)</th>
              <th>Truncating LoF (n, %)</th>
            </tr>
          </thead>
          <tbody></tbody>
        </table>
      </div>
    </div>
  </section>

  <!-- SECTION 5: GTEX, HPA & TRACTABILITY -->
  <section id="tab-gtex" class="panel">
    <details class="explainer-banner">
      <summary>
        <div class="exp-sum-left">
          <span class="exp-kicker">Section Guide</span>
          <span><strong>Healthy human germline constraint, adult tissue specificity (&tau;) &amp; therapeutic tractability</strong> across gnomAD v4, GTEx v10, HPA &amp; Open Targets</span>
        </div>
        <span class="exp-toggle-lbl">Expand explainer</span>
      </summary>
      <div class="exp-body">
        <div class="exp-col">
          <h4>1. Human Germline Purifying Selection</h4>
          <p><strong>gnomAD v4</strong> (&gt;800,000 individuals) shows extreme intolerance to germline truncation (<code>pLI &approx; 1.0</code>, <code>LOEUF &lt; 0.35</code>) and missense variation across core biogenesis enzymes—peaking in <em>AGO2</em> (<code>Missense Z = +7.97</code>, <code>LOEUF = 0.103</code>).</p>
        </div>
        <div class="exp-col">
          <h4>2. Housekeeping vs Oncofetal Specificity (&tau;)</h4>
          <p>Yanai&rsquo;s specificity index <strong>&tau; &isin; [0, 1]</strong> across 54 adult GTEx tissues separates ubiquitous housekeeping Core-10 enzymes (<code>&tau; = 0.22&ndash;0.38</code>) from testis-restricted <em>LIN28A</em> (<code>&tau; = 0.997</code>) and <em>LIN28B</em> (<code>&tau; = 0.978</code>), whose expression in somatic tumours represents true oncofetal reactivation.</p>
        </div>
        <div class="exp-col">
          <h4>3. Drug Discovery Whitespace</h4>
          <p><strong>Open Targets</strong> reveals that <strong>0 of 14</strong> biogenesis genes (excluding <em>TP53</em>) have clinical-stage drugs, despite structural small-molecule ligandability in <em>AGO2</em>, <em>XPO5</em>, and <em>SMAD2/4</em>.</p>
        </div>
        <div class="exp-chips">
          <span>Key definitions:</span>
          <span class="term-chip" data-term="gnomad_pli">gnomAD pLI</span>
          <span class="term-chip" data-term="gnomad_loeuf">LOEUF</span>
          <span class="term-chip" data-term="gnomad_misz">Missense Z</span>
          <span class="term-chip" data-term="gtex_tpm">GTEx Median TPM</span>
          <span class="term-chip" data-term="tau_spec">Specificity (&tau;)</span>
          <span class="term-chip" data-term="tractability">Open Targets Tractability</span>
        </div>
      </div>
    </details>

    <div class="surface">
      <div class="section-head">
        <h2>Germline Constraint (gnomAD v4), Subcellular Localisation (HPA) &amp; Tractability (Open Targets)</h2>
        <span class="section-meta">15-gene reference summary</span>
      </div>
      <div class="table-wrap">
        <table id="gene-ref-table">
          <thead>
            <tr>
              <th>Gene</th>
              <th data-term="gnomad_pli">gnomAD pLI</th>
              <th data-term="gnomad_loeuf">LOEUF</th>
              <th data-term="gnomad_misz">Missense Z</th>
              <th>HPA Subcellular Compartment</th>
              <th data-term="gtex_tpm">GTEx Pan-54 Median TPM</th>
              <th data-term="gtex_tpm">GTEx Peak Tissue (TPM)</th>
              <th data-term="tau_spec">Specificity (&tau;)</th>
              <th data-term="tractability">Tractability</th>
              <th data-term="kpi_auc">ClinVar Pathogenic / Total</th>
            </tr>
          </thead>
          <tbody></tbody>
        </table>
      </div>
    </div>
    <div class="surface">
      <div class="section-head">
        <h2>GTEx v10 Healthy Adult Expression Baseline Across 54 Tissues</h2>
        <div class="field-group">
          <span class="field-label">Gene</span>
          <select id="gtex-gene" onchange="renderGtexBars()"></select>
        </div>
      </div>
      <div id="gtex-bars" class="table-wrap" style="max-height:420px;"></div>
    </div>
  </section>

  <!-- SECTION 6: CPTAC, DEPMAP, SURVIVAL & PANIMMUNE -->
  <section id="tab-consortia" class="panel">
    <details class="explainer-banner">
      <summary>
        <div class="exp-sum-left">
          <span class="exp-kicker">Section Guide</span>
          <span><strong>Proteogenomic buffering, CRISPR collateral synthetic lethality, multivariable survival &amp; HRD</strong> across CPTAC, DepMap &amp; TCGA-CDR</span>
        </div>
        <span class="exp-toggle-lbl">Expand explainer</span>
      </summary>
      <div class="exp-body">
        <div class="exp-col">
          <h4>1. CPTAC Protein Buffering (&Delta;&rho;)</h4>
          <p>Comparing <code>&rho;(CNA, mRNA)</code> vs mass-spectrometry <code>&rho;(CNA, Protein)</code> across 1,002 tumours shows strong <em>AGO2</em> amplification-to-protein transmission alongside stoichiometry buffering (<code>&Delta;&rho; &gt; 0</code>) in heterodimers (<em>DROSHA&ndash;DGCR8</em>, <em>DICER1&ndash;TARBP2</em>).</p>
        </div>
        <div class="exp-col">
          <h4>2. DepMap Chronos &amp; 1p34.3 Synthetic Lethality</h4>
          <p><strong>Chronos &le; &minus;0.5</strong> marks cell essentiality (<em>DROSHA</em>, <em>DGCR8</em>, <em>XPO5</em>, <em>DICER1</em>). Argonaute paralogues buffer single knockouts (~0.0) until coordinate <code>1p34.3</code> <em>AGO1/3/4</em> depletion exposes collateral dependency on <em>AGO2</em> (<em>p</em> = 5.40 &times; 10<sup>&minus;4</sup>).</p>
        </div>
        <div class="exp-col">
          <h4>3. Stratified Survival &amp; Extreme HRD</h4>
          <p>Survival models stratify by cancer type (Mantel&ndash;Haenszel) and adjust for age, sex, AJCC stage, and hypermutation (Breslow multivariable Cox). Focal <code>8q24.3</code> <em>AGO2</em> <code>+2</code> marks severe homologous recombination deficiency (median HRD 33 vs 7, <em>p</em> = 1.43 &times; 10<sup>&minus;183</sup>).</p>
        </div>
        <div class="exp-chips">
          <span>Key definitions:</span>
          <span class="term-chip" data-term="cptac_rho">CPTAC CNA&rarr;Protein &rho;</span>
          <span class="term-chip" data-term="cptac_atten">Attenuation (&Delta;&rho;)</span>
          <span class="term-chip" data-term="depmap_chronos">DepMap Chronos</span>
          <span class="term-chip" data-term="strat_hr">Stratified Log-Rank HR</span>
          <span class="term-chip" data-term="multivar_cox">Multivariable Breslow Cox HR</span>
          <span class="term-chip" data-term="hrd_score">HRD Genomic-Scar Score</span>
        </div>
      </div>
    </details>

    <div class="surface">
      <div class="section-head">
        <h2>CPTAC Mass-Spectrometry Proteogenomics (N = 1,002) &amp; Broad DepMap CRISPR Essentiality (N = 1,208)</h2>
        <span class="section-meta">8 CPTAC cohorts &middot; 29 DepMap lineages</span>
      </div>
      <div class="table-wrap">
        <table id="consortia-table">
          <thead>
            <tr>
              <th>Gene</th>
              <th>Role</th>
              <th data-term="kpi_pdb">PDB Complexes</th>
              <th data-term="cptac_rho">CPTAC CNA&rarr;mRNA &rho;</th>
              <th data-term="cptac_rho">CPTAC CNA&rarr;Protein &rho;</th>
              <th data-term="cptac_atten">Attenuation (&Delta;&rho;)</th>
              <th data-term="depmap_chronos">DepMap Median Chronos</th>
              <th data-term="depmap_chronos">% Dependent (&le;&minus;0.5)</th>
              <th data-term="depmap_chronos">Dependency Class</th>
            </tr>
          </thead>
          <tbody></tbody>
        </table>
      </div>
    </div>

    <div class="surface">
      <div class="section-head">
        <h2>DepMap CRISPR-Cas9 Collateral Synthetic Lethality, CYCLOPS &amp; Oncofetal Addiction Benchmarks</h2>
        <span class="section-meta">1,208 CRISPR-screened cancer cell lines</span>
      </div>
      <div class="table-wrap">
        <table id="collat-table">
          <thead>
            <tr>
              <th>Test ID</th>
              <th>CRISPR Target</th>
              <th>Biological Hypothesis</th>
              <th>Altered Lines (n, Median Chronos)</th>
              <th>Control Lines (n, Median Chronos)</th>
              <th data-term="depmap_chronos">&Delta;Chronos</th>
              <th data-term="rank_shift_p">Two-Sided p</th>
            </tr>
          </thead>
          <tbody></tbody>
        </table>
      </div>
    </div>

    <div class="surface">
      <div class="section-head">
        <h2>Cancer-Type-Stratified Log-Rank &amp; Multivariable Breslow Cox Survival Benchmarks</h2>
        <span class="section-meta">Adjusted for age, sex, AJCC stage, hypermutator status, and tumour strata</span>
      </div>
      <div class="table-wrap">
        <table id="surv-table">
          <thead>
            <tr>
              <th>Benchmark ID</th>
              <th>Cohort</th>
              <th data-term="surv_endpoints">Endpoint</th>
              <th>Altered (Events / n)</th>
              <th>Control (Events / n)</th>
              <th data-term="strat_hr">Stratified HR (95% CI)</th>
              <th data-term="strat_hr">Log-Rank p</th>
              <th data-term="multivar_cox">Multivariable Cox HR (95% CI)</th>
              <th data-term="multivar_cox">Multivariable p</th>
            </tr>
          </thead>
          <tbody></tbody>
        </table>
      </div>
    </div>
    <div class="surface">
      <div class="section-head">
        <h2>Thorsson Pan-Cancer Immune &amp; Homologous Recombination Deficiency (HRD) Benchmarks</h2>
        <span class="section-meta">N = 9,899 immune-profiled TCGA tumours</span>
      </div>
      <div class="table-wrap">
        <table id="imm-table">
          <thead>
            <tr>
              <th>Benchmark ID</th>
              <th data-term="hrd_score">Metric</th>
              <th>Group A (n)</th>
              <th>Median A</th>
              <th>Group B (n)</th>
              <th>Median B</th>
              <th data-term="rank_shift_p">Median Shift</th>
              <th data-term="rank_shift_p">Two-Sided p</th>
            </tr>
          </thead>
          <tbody></tbody>
        </table>
      </div>
    </div>
  </section>

  <!-- SECTION 7: METHODS, MATHEMATICAL FORMULAS, DATASETS & A-Z GLOSSARY -->
  <section id="tab-docs" class="panel">
    <div class="doc-jumpbar">
      <button class="btn-quiet" onclick="document.getElementById('sec-pathway').scrollIntoView({{behavior:'smooth'}})">1. Pathway &amp; Gene Tiers</button>
      <button class="btn-quiet" onclick="document.getElementById('sec-formulas').scrollIntoView({{behavior:'smooth'}})">2. Mathematical Formulas (F1&ndash;F16)</button>
      <button class="btn-quiet" onclick="document.getElementById('sec-datasets').scrollIntoView({{behavior:'smooth'}})">3. Source Consortia &amp; 41 Datasets</button>
      <button class="btn-quiet" onclick="document.getElementById('sec-tcga-codes').scrollIntoView({{behavior:'smooth'}})">4. 32 TCGA Study Codes &amp; N</button>
      <button class="btn-quiet" onclick="document.getElementById('sec-glossary').scrollIntoView({{behavior:'smooth'}})">5. Searchable A&ndash;Z Glossary</button>
    </div>

    <!-- 7.1 Pathway Architecture & Nested Tiers -->
    <div class="surface" id="sec-pathway">
      <div class="section-head">
        <h2>1. Canonical Human miRNA Biogenesis Cascade, Nested Gene-Set Tiers &amp; Exclusion Rationale</h2>
        <span class="section-meta">15 included proteins &middot; 4 nested analysis tiers &middot; 8 excluded pleiotropic factors</span>
      </div>

      <div class="pathway-flow">
        <div class="pf-step">
          <div class="pf-num">Step 1 &middot; Nucleus</div>
          <div class="pf-title">Microprocessor Cropping</div>
          <div class="pf-desc"><strong>DROSHA</strong> (RNase III, 5p13.3) + <strong>DGCR8</strong> (dsRBD dimer, 22q11.21) cleave capped/polyadenylated pri-miRNAs into ~65&ndash;70 nt hairpin pre-miRNAs (PDB <code>9ASM</code>, <code>6V5B</code>).</div>
        </div>
        <div class="pf-step">
          <div class="pf-num">Step 2 &middot; Nuclear Pore</div>
          <div class="pf-title">Pre-miRNA Nuclear Export</div>
          <div class="pf-desc"><strong>XPO5</strong> (Exportin-5, 6p21.1) recognises the 2-nt 3&prime; overhang and dsRNA stem of pre-miRNAs in a RanGTP-dependent tunnel (PDB <code>3A6P</code>).</div>
        </div>
        <div class="pf-step">
          <div class="pf-num">Step 3 &middot; Cytoplasm</div>
          <div class="pf-title">Dicing into 5p / 3p Duplex</div>
          <div class="pf-desc"><strong>DICER1</strong> (14q32.13; RNase IIIb cleaves <code>5p</code>, RNase IIIa cleaves <code>3p</code>) with dsRNA cofactors <strong>TARBP2</strong> (12q13.13) &amp; <strong>PRKRA</strong> (2q31.2) produces ~22-nt duplexes (PDB <code>7XW2</code>, <code>5ZAK</code>).</div>
        </div>
        <div class="pf-step">
          <div class="pf-num">Step 4 &middot; RISC Complex</div>
          <div class="pf-title">Guide Loading &amp; Slicing</div>
          <div class="pf-desc"><strong>AGO2</strong> (8q24.3, catalytic PIWI Slicer; cleaves <em>miR-451a</em> &amp; aids RAD51 HR repair; PDB <code>4W5N</code>) + syntenic <code>1p34.3</code> paralogues <strong>AGO1, AGO3, AGO4</strong> silence target mRNAs.</div>
        </div>
        <div class="pf-step">
          <div class="pf-num">Step 5 &middot; Regulation</div>
          <div class="pf-title">Let-7 Gate &amp; Accessories</div>
          <div class="pf-desc">Oncofetal <strong>LIN28A/B</strong> block pri/pre-<em>let-7</em> maturation (PDB <code>5UDZ</code>, <code>4A4I</code>); <strong>SMAD2/4</strong> (18q21.1) &amp; <strong>TP53</strong> (17p13.1) modulate pri-miRNA cropping &amp; <em>miR-34</em> transcription.</div>
        </div>
      </div>

      <div class="table-wrap" style="max-height:none;margin-bottom:14px;">
        <table>
          <thead>
            <tr>
              <th>Nested Tier</th>
              <th>Gene Count</th>
              <th>Included Symbols</th>
              <th>Syntenic Clusters Collapsed in Patient Union</th>
              <th>When to Use This Tier</th>
            </tr>
          </thead>
          <tbody>
            <tr>
              <td><strong>Core-10</strong></td>
              <td>10</td>
              <td><code>DROSHA, DGCR8, XPO5, DICER1, TARBP2, PRKRA, AGO1, AGO2, AGO3, AGO4</code></td>
              <td><code>1p34.3</code> (<em>AGO1/3/4</em>)</td>
              <td>Primary specific burden of dedicated canonical miRNA maturation &amp; RISC machinery.</td>
            </tr>
            <tr>
              <td><strong>Extended-12</strong></td>
              <td>12</td>
              <td>Core-10 + <code>LIN28A, LIN28B</code></td>
              <td><code>1p34.3</code> (<em>AGO1/3/4</em>)</td>
              <td>Captures direct pri/pre-<em>let-7</em> oncofetal blockade alongside core enzymatic defects.</td>
            </tr>
            <tr>
              <td><strong>Panel-14</strong></td>
              <td>14</td>
              <td>Extended-12 + <code>SMAD2, SMAD4</code></td>
              <td><code>1p34.3</code> (<em>AGO1/3/4</em>) &amp; <code>18q21.1</code> (<em>SMAD2/4</em>)</td>
              <td>Includes TGF-&beta;/SMAD Microprocessor accessory cofactors while still excluding <em>TP53</em>.</td>
            </tr>
            <tr>
              <td><strong>All-15 (Reference)</strong></td>
              <td>15</td>
              <td>Panel-14 + <code>TP53</code></td>
              <td>All above + <code>17p13.1</code> (<em>TP53</em>)</td>
              <td>Gene-level &amp; structural comparison only; <em>TP53</em> is excluded from pathway union summaries to prevent pan-cancer tumour-suppressor dominance.</td>
            </tr>
          </tbody>
        </table>
      </div>

      <details style="margin-top:8px;">
        <summary style="cursor:pointer;font-size:12.5px;font-weight:600;color:var(--embl-petrol);">View 8 Pleiotropic / Downstream RNA Factors Explicitly Excluded from the Core Panel &darr;</summary>
        <div class="table-wrap" style="margin-top:10px;max-height:none;">
          <table id="excluded-cands-table">
            <thead>
              <tr><th>Candidate Gene</th><th>Cytoband</th><th>Molecular Function</th><th>Exclusion Rationale</th></tr>
            </thead>
            <tbody></tbody>
          </table>
        </div>
      </details>
    </div>

    <!-- 7.2 Complete Mathematical Formulas (F1-F16) -->
    <div class="surface" id="sec-formulas">
      <div class="section-head">
        <h2>2. Mathematical, Biophysical &amp; Biostatistical Formulas Used Across miRBioGenAtlas</h2>
        <span class="section-meta">16 exact equations verified against pipeline source code</span>
      </div>
      <div class="toolbar">
        <div class="field-group">
          <span class="field-label">Category</span>
          <select id="formula-cat" onchange="renderFormulas()">
            <option value="ALL">All 16 formulas</option>
            <option value="Genomics & Cohorts">Genomics &amp; Cohorts (F1&ndash;F2)</option>
            <option value="3D Structural Biophysics">3D Structural Biophysics (F3, F5&ndash;F7)</option>
            <option value="Deep Mutational Scoring">Deep Mutational Scoring (F4)</option>
            <option value="Small-RNA Enzymology">Small-RNA Enzymology &amp; Couplings (F8&ndash;F11)</option>
            <option value="Tissue & Proteogenomics">Tissue Specificity &amp; Proteogenomics (F12&ndash;F13)</option>
            <option value="Biostatistics & Survival">Biostatistics, ROC-AUC &amp; Survival (F14&ndash;F16)</option>
          </select>
        </div>
        <div class="field-group">
          <span class="field-label">Filter</span>
          <input type="text" id="formula-q" placeholder="e.g. Cox, 5p:3p, Grantham, ESM-IF1..." oninput="renderFormulas()" style="width:240px;">
        </div>
      </div>
      <div id="formula-cards-wrap" class="formula-grid"></div>
    </div>

    <!-- 7.3 Source Consortia & 41 Derived Datasets -->
    <div class="surface" id="sec-datasets">
      <div class="section-head">
        <h2>3. Integrated Public Consortia (10 Sources) &amp; Derived Dataset Files Catalogue</h2>
        <span class="section-meta">46 SHA-256-verified raw snapshots &rarr; 41 derived CSV / CSV.GZ tables</span>
      </div>
      <h3 style="font-size:13.5px;margin:0 0 8px;color:var(--ink);">3.1 Upstream Multi-Consortium Data Sources &amp; Harmonisation Keys</h3>
      <div class="table-wrap" style="max-height:none;margin-bottom:18px;">
        <table id="consortia-sources-table">
          <thead>
            <tr><th>Consortium / Resource</th><th>Cohort / Variant Scale</th><th>Data Layers Extracted</th><th>Harmonisation Join Key</th></tr>
          </thead>
          <tbody></tbody>
        </table>
      </div>

      <div style="display:flex;justify-content:space-between;align-items:center;flex-wrap:wrap;gap:10px;margin-bottom:10px;">
        <h3 style="font-size:13.5px;margin:0;color:var(--ink);">3.2 Derived CSV / CSV.GZ Tables in <code>data/</code> (Zenodo DOI: 10.5281/zenodo.24455067)</h3>
        <div class="field-group">
          <span class="field-label">Search files</span>
          <input type="text" id="ds-q" placeholder="Filter filename, key, or metric..." oninput="renderDatasetsCatalogue()" style="width:240px;">
        </div>
      </div>
      <div class="table-wrap" style="max-height:420px;">
        <table id="datasets-cat-table">
          <thead>
            <tr><th>Dataset Filename (<code>data/</code>)</th><th>Rows &times; Cols</th><th>Size</th><th>Primary Key / Grain</th><th>Contents &amp; Analytical Purpose</th></tr>
          </thead>
          <tbody></tbody>
        </table>
      </div>
    </div>

    <!-- 7.4 32 TCGA Study Codes Decoder -->
    <div class="surface" id="sec-tcga-codes">
      <div class="section-head">
        <h2>4. TCGA PanCancer Atlas 32-Study Acronym Decoder &amp; Matched Sample Denominators</h2>
        <span class="section-meta">Total matched N = 10,195 (33 cancer types; COAD + READ combined in COADREAD)</span>
      </div>
      <div class="table-wrap" style="max-height:380px;">
        <table id="study-decoder-table">
          <thead>
            <tr>
              <th>Code</th>
              <th>Full Tumour Study Name</th>
              <th>Matched N (Seq &cap; CNA)</th>
              <th>Exome Seq</th>
              <th>CNA SNP6</th>
              <th>mRNA-seq</th>
              <th>miRNA-seq</th>
              <th data-term="hypermutator">Hypermutated (n)</th>
              <th data-term="hypermutator">Non-Hypermut (n)</th>
            </tr>
          </thead>
          <tbody></tbody>
        </table>
      </div>
    </div>

    <!-- 7.5 Searchable A-Z Glossary -->
    <div class="surface" id="sec-glossary">
      <div class="section-head">
        <h2>5. Interactive A&ndash;Z Glossary of Terms, Biophysical Scores &amp; Statistical Metrics</h2>
        <span class="section-meta" id="gloss-count">All terms</span>
      </div>
      <div class="toolbar">
        <div class="field-group">
          <span class="field-label">Search Glossary</span>
          <input type="text" id="gloss-q" placeholder="Type any term (e.g. LOEUF, GISTIC2, RNase IIIb, Chronos, HRD, pLDDT)..." oninput="renderGlossaryCards()" style="width:340px;">
        </div>
      </div>
      <div id="glossary-cards-wrap" class="glossary-grid"></div>
    </div>
  </section>
</main>

<div id="term-popover" class="term-popover" role="dialog" aria-live="polite">
  <div class="tp-head">
    <div>
      <span class="tp-cat" id="tp-cat"></span>
      <div class="tp-title" id="tp-title"></div>
    </div>
    <button class="tp-close" onclick="hideTermPopover(true)" aria-label="Close definition">&times;</button>
  </div>
  <div class="tp-meaning" id="tp-meaning"></div>
  <div class="tp-formula" id="tp-formula"></div>
  <div class="tp-thresh" id="tp-thresh"></div>
  <button class="tp-link" id="tp-link">Open full derivation in Methods, Formulas &amp; Datasets &rarr;</button>
</div>

<footer>
  <div class="container footer-inner">
    <div>miRBioGenAtlas v2.0.2 &middot; Manipal Academy of Higher Education (MAHE) &middot; Built deterministically from 46 SHA-256-verified multi-consortium archives.</div>
    <div>
      <a href="https://doi.org/10.5281/zenodo.24455067" target="_blank" rel="noopener">Zenodo Deposit (CC BY 4.0)</a> &middot;
      <a href="https://github.com/siddugayam/miRBioGenAtlas" target="_blank" rel="noopener">Source &amp; Documentation (MIT)</a>
    </div>
  </div>
</footer>

<script>
window.SAT_CACHE = window.SAT_CACHE || {{}};
window.STRUCT3D_CACHE = window.STRUCT3D_CACHE || {{}};
window.RCSB_CACHE = window.RCSB_CACHE || {{}};
const ATLAS = {data_json};
const LOW_N = new Set(["ACC","CHOL","DLBC","KICH","MESO","UCS","UVM"]);
const GENE_SET = new Set(ATLAS.panel.map(p => p.geneSymbol));

const AA_PROPS = {{
  A: {{three:'ALA', name:'Alanine', cls:'Aliphatic nonpolar', q:0.0, vol:88.6, mw:89.09, kd:1.8, hd:0, ha:0, arom:0, metal:0}},
  R: {{three:'ARG', name:'Arginine', cls:'Basic guanidinium (+)', q:1.0, vol:173.4, mw:174.20, kd:-4.5, hd:5, ha:0, arom:0, metal:0}},
  N: {{three:'ASN', name:'Asparagine', cls:'Polar carboxamide', q:0.0, vol:114.1, mw:132.12, kd:-3.5, hd:2, ha:1, arom:0, metal:1}},
  D: {{three:'ASP', name:'Aspartate', cls:'Acidic carboxylate (-)', q:-1.0, vol:111.1, mw:133.10, kd:-3.5, hd:0, ha:2, arom:0, metal:1}},
  C: {{three:'CYS', name:'Cysteine', cls:'Thiol / sulfhydryl', q:0.0, vol:108.5, mw:121.16, kd:2.5, hd:1, ha:1, arom:0, metal:1}},
  Q: {{three:'GLN', name:'Glutamine', cls:'Polar carboxamide', q:0.0, vol:143.8, mw:146.15, kd:-3.5, hd:2, ha:1, arom:0, metal:1}},
  E: {{three:'GLU', name:'Glutamate', cls:'Acidic carboxylate (-)', q:-1.0, vol:138.4, mw:147.13, kd:-3.5, hd:0, ha:2, arom:0, metal:1}},
  G: {{three:'GLY', name:'Glycine', cls:'Achiral H-only', q:0.0, vol:60.1, mw:75.07, kd:-0.4, hd:0, ha:0, arom:0, metal:0}},
  H: {{three:'HIS', name:'Histidine', cls:'Imidazole / Zn-ligand', q:0.1, vol:153.2, mw:155.16, kd:-3.2, hd:1, ha:1, arom:1, metal:1}},
  I: {{three:'ILE', name:'Isoleucine', cls:'Branched hydrophobic', q:0.0, vol:166.7, mw:131.17, kd:4.5, hd:0, ha:0, arom:0, metal:0}},
  L: {{three:'LEU', name:'Leucine', cls:'Branched hydrophobic', q:0.0, vol:166.7, mw:131.17, kd:3.8, hd:0, ha:0, arom:0, metal:0}},
  K: {{three:'LYS', name:'Lysine', cls:'Basic primary amine (+)', q:1.0, vol:168.6, mw:146.19, kd:-3.9, hd:3, ha:0, arom:0, metal:0}},
  M: {{three:'MET', name:'Methionine', cls:'Thioether nonpolar', q:0.0, vol:162.9, mw:149.21, kd:1.9, hd:0, ha:1, arom:0, metal:0}},
  F: {{three:'PHE', name:'Phenylalanine', cls:'Aromatic nonpolar', q:0.0, vol:189.9, mw:165.19, kd:2.8, hd:0, ha:0, arom:1, metal:0}},
  P: {{three:'PRO', name:'Proline', cls:'Cyclic imino (rigid)', q:0.0, vol:112.7, mw:115.13, kd:-1.6, hd:0, ha:0, arom:0, metal:0}},
  S: {{three:'SER', name:'Serine', cls:'Polar hydroxyl', q:0.0, vol:89.0, mw:105.09, kd:-0.8, hd:1, ha:1, arom:0, metal:0}},
  T: {{three:'THR', name:'Threonine', cls:'Beta-branched hydroxyl', q:0.0, vol:116.1, mw:119.12, kd:-0.7, hd:1, ha:1, arom:0, metal:0}},
  W: {{three:'TRP', name:'Tryptophan', cls:'Indole aromatic', q:0.0, vol:227.8, mw:204.23, kd:-0.9, hd:1, ha:0, arom:1, metal:0}},
  Y: {{three:'TYR', name:'Tyrosine', cls:'Phenolic aromatic', q:0.0, vol:193.6, mw:181.19, kd:-1.3, hd:1, ha:1, arom:1, metal:0}},
  V: {{three:'VAL', name:'Valine', cls:'Beta-branched hydrophobic', q:0.0, vol:140.0, mw:117.15, kd:4.2, hd:0, ha:0, arom:0, metal:0}}
}};

const GRANTHAM_CPV = {{
  S:[1.42,9.2,32], R:[0.65,10.5,124], L:[0.00,4.9,111], P:[0.39,8.0,32.5], T:[0.71,8.6,61],
  A:[0.00,8.1,31], V:[0.00,5.9,84], G:[0.74,9.0,3], I:[0.00,5.2,111], F:[0.00,5.2,132],
  Y:[0.20,6.2,136], C:[2.75,5.5,55], H:[0.58,10.4,96], Q:[0.89,10.5,85], N:[1.33,11.6,56],
  K:[0.33,11.3,119], D:[1.38,13.0,54], E:[0.92,12.3,83], M:[0.00,5.7,105], W:[0.13,5.4,170]
}};
function granthamDist(a, b) {{
  if (a === b || !GRANTHAM_CPV[a] || !GRANTHAM_CPV[b]) return 0;
  const [ca, pa, va] = GRANTHAM_CPV[a], [cb, pb, vb] = GRANTHAM_CPV[b];
  return Math.round(50.723 * Math.sqrt(1.833*(ca-cb)*(ca-cb) + 0.1018*(pa-pb)*(pa-pb) + 0.000399*(va-vb)*(va-vb)));
}}

function tagDot(cls, label) {{
  return `<span class="tag"><i class="dot ${{cls}}"></i>${{label}}</span>`;
}}

let viewer3d = null;
let spinActive = false;
let showNeighbors = true;
let active3DState = {{ gene: 'DICER1', pos: 1705, wt: 'E', mut: 'K', lastPdbText: '' }};

function switchTab(id, btn) {{
  document.querySelectorAll('.panel').forEach(p => p.classList.remove('active'));
  document.querySelectorAll('.tab-btn').forEach(b => b.classList.remove('active'));
  document.getElementById(id).classList.add('active');
  if (btn) btn.classList.add('active');
  if (id === 'tab-variants') {{
    setTimeout(() => {{
      if (viewer3d) viewer3d.resize();
      refresh3DScene(true);
    }}, 40);
  }}
}}

function openDocsSection(anchorId) {{
  switchTab('tab-docs', document.getElementById('btn-tab-docs'));
  if (anchorId) {{
    setTimeout(() => {{
      const el = document.getElementById(anchorId);
      if (el) el.scrollIntoView({{behavior: 'smooth', block: 'start'}});
    }}, 40);
  }}
}}

let tpPinned = false;
let tpHideTimer = null;
function showTermPopover(targetEl, pin) {{
  const key = targetEl.getAttribute('data-term');
  const item = ATLAS.glossary && ATLAS.glossary[key];
  if (!item) return;
  clearTimeout(tpHideTimer);
  if (pin) tpPinned = true;

  const pop = document.getElementById('term-popover');
  document.getElementById('tp-cat').textContent = item.cat || 'Definition';
  document.getElementById('tp-title').textContent = item.title || key;
  document.getElementById('tp-meaning').textContent = item.meaning || '';
  const fEl = document.getElementById('tp-formula');
  if (item.formula) {{ fEl.style.display = 'block'; fEl.textContent = item.formula; }}
  else {{ fEl.style.display = 'none'; }}
  document.getElementById('tp-thresh').textContent = item.thresholds || '';
  const lnk = document.getElementById('tp-link');
  lnk.onclick = () => {{
    hideTermPopover(true);
    openDocsSection(item.anchor || 'sec-formulas');
  }};

  pop.classList.add('visible');
  const rect = targetEl.getBoundingClientRect();
  const pw = pop.offsetWidth || 380;
  const ph = pop.offsetHeight || 210;
  let left = rect.left;
  if (left + pw > window.innerWidth - 14) left = Math.max(12, window.innerWidth - pw - 14);
  let top = rect.bottom + 8;
  if (top + ph > window.innerHeight - 12 && rect.top - ph - 8 > 8) {{
    top = rect.top - ph - 8;
  }}
  pop.style.left = left + 'px';
  pop.style.top = top + 'px';
}}

function hideTermPopover(force) {{
  if (tpPinned && !force) return;
  tpPinned = false;
  const pop = document.getElementById('term-popover');
  if (pop) pop.classList.remove('visible');
}}

document.addEventListener('mouseover', e => {{
  const t = e.target.closest('[data-term]');
  if (t && !tpPinned) showTermPopover(t, false);
}});
document.addEventListener('mouseout', e => {{
  const t = e.target.closest('[data-term]');
  if (t && !tpPinned) {{
    tpHideTimer = setTimeout(() => hideTermPopover(false), 220);
  }}
}});
document.addEventListener('click', e => {{
  const pop = document.getElementById('term-popover');
  if (pop && pop.contains(e.target)) return;
  const t = e.target.closest('[data-term]');
  if (t) {{
    e.stopPropagation();
    showTermPopover(t, true);
    return;
  }}
  if (tpPinned) hideTermPopover(true);
}});
document.getElementById('term-popover').addEventListener('mouseenter', () => clearTimeout(tpHideTimer));
document.getElementById('term-popover').addEventListener('mouseleave', () => {{
  if (!tpPinned) tpHideTimer = setTimeout(() => hideTermPopover(false), 180);
}});

function runGlobalSearch(raw) {{
  const q = (raw || '').trim().toUpperCase().replace(/^P\\./, '');
  if (!q) return;
  const parts = q.split(/[\\s:]+/);
  if (parts.length >= 2 && GENE_SET.has(parts[0])) {{
    loadPresetMutation(parts[0], parts[1]);
    return;
  }}
  if (GENE_SET.has(q)) {{
    document.getElementById('lol-gene').value = q;
    document.getElementById('sat-query').value = '';
    switchTab('tab-variants', document.getElementById('btn-tab-variants'));
    renderLollipopAndLoadSat();
    return;
  }}
  if (/^[A-Z]\\d+[A-Z]$/.test(q)) {{
    document.getElementById('sat-query').value = q;
    switchTab('tab-variants', document.getElementById('btn-tab-variants'));
    inspectSaturationCodon();
    return;
  }}
  const studyHit = (ATLAS.study_summary || []).find(s => s.tumour_acronym === q);
  if (studyHit) {{
    openDocsSection('sec-tcga-codes');
    return;
  }}
  const glossEntries = Object.values(ATLAS.glossary || {{}});
  const qLow = raw.trim().toLowerCase();
  const glossHit = glossEntries.find(g => (g.title + ' ' + g.cat + ' ' + g.formula).toLowerCase().includes(qLow));
  if (glossHit) {{
    document.getElementById('gloss-q').value = raw.trim();
    renderGlossaryCards();
    openDocsSection('sec-glossary');
    return;
  }}
  document.getElementById('var-q').value = raw.trim();
  switchTab('tab-variants', document.getElementById('btn-tab-variants'));
  renderVariants();
  document.getElementById('var-table').scrollIntoView({{behavior:'smooth'}});
}}

document.addEventListener('keydown', e => {{
  if (e.key === 'Escape') hideTermPopover(true);
  if (e.key === '/' && document.activeElement.tagName !== 'INPUT' && document.activeElement.tagName !== 'SELECT') {{
    e.preventDefault();
    document.getElementById('global-search').focus();
  }}
}});

function loadPresetMutation(gene, mutStr) {{
  document.getElementById('lol-gene').value = gene;
  document.getElementById('sat-query').value = mutStr;
  if (!document.getElementById('tab-variants').classList.contains('active')) {{
    switchTab('tab-variants', document.getElementById('btn-tab-variants'));
  }}
  renderLollipopAndLoadSat();
}}

function openVariantIn3D(gene, mutStr) {{
  loadPresetMutation(gene, mutStr);
  document.getElementById('struct3d-card').scrollIntoView({{behavior: 'smooth', block: 'start'}});
}}

// EMBL-petrol monochrome ramp
function heatColor(val, maxScale) {{
  const t = Math.min(1, Math.max(0, val / maxScale));
  const r = Math.round(255 - t * (255 - 0));
  const g = Math.round(255 - t * (255 - 88));
  const b = Math.round(255 - t * (255 - 124));
  const txt = t > 0.52 ? '#ffffff' : '#1d1d1f';
  return `background:rgb(${{r}},${{g}},${{b}});color:${{txt}};`;
}}

function renderHeatmap() {{
  const metric = document.getElementById('hm-metric').value;
  const maxScale = (metric === 'mut' || metric === 'nhm') ? 0.12 : (metric === 'any' ? 0.85 : 0.40);
  const studies = [...new Set(ATLAS.long.map(r => r.t))].sort();
  const genes = ATLAS.panel.map(p => p.geneSymbol);
  const map = {{}};
  ATLAS.long.forEach(r => {{ map[r.g + '|' + r.t] = r; }});

  let html = '<table><thead><tr><th>Gene</th>';
  studies.forEach(s => {{ html += `<th style="text-align:center">${{s}}${{LOW_N.has(s)?'*':''}}</th>`; }});
  html += '</tr></thead><tbody>';
  genes.forEach(g => {{
    html += `<tr><td><strong>${{g}}</strong></td>`;
    studies.forEach(s => {{
      const r = map[g + '|' + s];
      const v = r ? r[metric] : 0;
      const tip = r ? `${{g}} in ${{s}} (N=${{r.N}}): ${{v.toFixed(3)}} | mut=${{r.mut_n}}, amp(+2)=${{r.amp_n}}, deep_del(-2)=${{r.dd_n}}` : '';
      html += `<td class="hm-cell" style="${{heatColor(v, maxScale)}}" title="${{tip}}">${{v.toFixed(2)}}</td>`;
    }});
    html += '</tr>';
  }});
  html += '</tbody></table>';
  document.getElementById('hm-container').innerHTML = html;
}}

function ensureSatLoaded(gene, cb) {{
  if (window.SAT_CACHE && window.SAT_CACHE[gene]) {{
    const nCod = Object.keys(window.SAT_CACHE[gene]).length;
    document.getElementById('sat-status').textContent = `${{gene}} · ${{nCod.toLocaleString()}} codons · ${{(nCod * 19).toLocaleString()}} missense substitutions`;
    if (cb) cb();
    return;
  }}
  document.getElementById('sat-status').textContent = `Loading ${{gene}} saturation & structural profile...`;
  const sc = document.createElement('script');
  sc.src = `sat/${{gene}}.js`;
  sc.onload = () => {{
    const nCod = window.SAT_CACHE[gene] ? Object.keys(window.SAT_CACHE[gene]).length : 0;
    document.getElementById('sat-status').textContent = `${{gene}} · ${{nCod.toLocaleString()}} codons · ${{(nCod * 19).toLocaleString()}} missense substitutions`;
    if (cb) cb();
  }};
  sc.onerror = () => {{
    document.getElementById('sat-status').textContent = `${{gene}} codon track`;
    if (cb) cb();
  }};
  document.head.appendChild(sc);
}}

function ensureStruct3DLoaded(gene, cb) {{
  if (window.STRUCT3D_CACHE && window.STRUCT3D_CACHE[gene]) {{
    if (!window.STRUCT3D_CACHE[gene]._indexed) indexStruct3D(window.STRUCT3D_CACHE[gene]);
    if (cb) cb(window.STRUCT3D_CACHE[gene]);
    return;
  }}
  const sc = document.createElement('script');
  sc.src = `struct3d/${{gene}}.js`;
  sc.onload = () => {{
    const st = window.STRUCT3D_CACHE[gene];
    if (st && !st._indexed) indexStruct3D(st);
    if (cb) cb(st);
  }};
  sc.onerror = () => {{ if (cb) cb(null); }};
  document.head.appendChild(sc);
}}

function indexStruct3D(st) {{
  const byRes = {{}};
  const caList = [];
  st.a.forEach(row => {{
    const [pos, aname, x, y, z, plddt] = row;
    if (!byRes[pos]) byRes[pos] = [];
    byRes[pos].push(row);
    if (aname === 'CA') caList.push([pos, x, y, z, plddt]);
  }});
  st._byRes = byRes;
  st._caList = caList;
  st._indexed = true;
}}

function selectCodonFromTrack(pos) {{
  const gene = document.getElementById('lol-gene').value;
  const gCache = window.SAT_CACHE ? window.SAT_CACHE[gene] : null;
  let q = String(pos);
  if (gCache && gCache[q]) {{
    const subs = gCache[q].s;
    const obs = subs.slice().sort((a,b) => (b[5] - a[5]) || ((b[1]||0) - (a[1]||0)))[0];
    if (obs) q = gCache[q].m[0] + pos + obs[0];
  }}
  document.getElementById('sat-query').value = q;
  inspectSaturationCodon();
}}

function selectMutAA(pos, wt, mutAA) {{
  document.getElementById('sat-query').value = wt + String(pos) + mutAA;
  inspectSaturationCodon();
}}

function renderLollipopAndLoadSat() {{
  const gene = document.getElementById('lol-gene').value;
  const tr = ATLAS.tracks[gene] || [];
  if (!tr.length) return;
  const L = tr.length;
  const W = 1320, H = 178, padL = 96, padR = 24;
  const plotW = W - padL - padR;
  const xOf = pos => padL + ((pos - 1) / Math.max(1, L - 1)) * plotW;

  let maxCnt = 1;
  tr.forEach(c => {{
    const tot = c[7] + c[8];
    if (tot > maxCnt) maxCnt = tot;
  }});

  let pockets = '';
  let amPath = '';
  let stems = '';
  tr.forEach(c => {{
    const [pos, wt, plddt, pkt, am, e2, if1, tcgaN, cvN, pdb] = c;
    const x = xOf(pos);
    if (pkt === 1) {{
      pockets += `<rect x="${{x-1.2}}" y="108" width="2.4" height="14" fill="#b91c1c" opacity="0.78" style="cursor:pointer" onclick="selectCodonFromTrack(${{pos}})"><title>${{gene}} ${{wt}}${{pos}} · 3D active site (&le;6&Aring; in ${{pdb}})</title></rect>`;
    }} else if (pdb) {{
      pockets += `<rect x="${{x-0.8}}" y="111" width="1.6" height="8" fill="#00709b" opacity="0.30"></rect>`;
    }}
    if (am !== null) {{
      const yAm = 160 - am * 26;
      amPath += (amPath ? ' L' : 'M') + x.toFixed(1) + ',' + yAm.toFixed(1);
    }}
    if (tcgaN > 0 || cvN > 0) {{
      const h = 14 + (Math.log2(tcgaN + cvN + 1) / Math.log2(maxCnt + 1)) * 74;
      const yTop = 104 - h;
      const col = pkt === 1 ? '#b91c1c' : (cvN > 0 ? '#007c82' : '#515154');
      const r = Math.min(5.8, 2.5 + Math.log2(tcgaN + cvN + 1));
      stems += `<g style="cursor:pointer" onclick="selectCodonFromTrack(${{pos}})">
        <line x1="${{x.toFixed(1)}}" y1="104" x2="${{x.toFixed(1)}}" y2="${{yTop.toFixed(1)}}" stroke="${{col}}" stroke-width="1.15" opacity="0.72"/>
        <circle cx="${{x.toFixed(1)}}" cy="${{yTop.toFixed(1)}}" r="${{r.toFixed(1)}}" fill="${{col}}" stroke="#ffffff" stroke-width="0.7">
          <title>${{gene}} ${{wt}}${{pos}} · TCGA: ${{tcgaN}} · ClinVar P/LP: ${{cvN}} · AlphaMissense: ${{am ?? 'NA'}}</title>
        </circle>
      </g>`;
    }}
  }});

  let ticks = '';
  const step = L > 1200 ? 200 : (L > 500 ? 100 : 50);
  for (let p = 1; p <= L; p += step) {{
    const tx = xOf(p);
    ticks += `<text x="${{tx.toFixed(1)}}" y="174" font-size="10" font-family="-apple-system,sans-serif" fill="#6e6e73" text-anchor="middle">${{p}}</text>`;
  }}
  ticks += `<text x="${{xOf(L).toFixed(1)}}" y="174" font-size="10" font-family="-apple-system,sans-serif" fill="#6e6e73" text-anchor="end">${{L}} aa</text>`;

  const svg = `<svg viewBox="0 0 ${{W}} ${{H}}" width="100%" height="178">
    <text x="6" y="22" font-size="10" font-weight="600" fill="#6e6e73" letter-spacing="0.03em">VARIANTS</text>
    <text x="6" y="119" font-size="10" font-weight="600" fill="#6e6e73" letter-spacing="0.03em">3D COVERAGE</text>
    <text x="6" y="152" font-size="10" font-weight="600" fill="#6e6e73" letter-spacing="0.03em">CONSTRAINT</text>
    <rect x="${{padL}}" y="108" width="${{plotW}}" height="14" rx="3" fill="#e8e8ed"/>
    ${{pockets}}
    ${{stems}}
    <path d="${{amPath}}" fill="none" stroke="#00709b" stroke-width="1.05" opacity="0.82"/>
    <line x1="${{padL}}" y1="${{160 - 0.564*26}}" x2="${{W-padR}}" y2="${{160 - 0.564*26}}" stroke="#b91c1c" stroke-dasharray="3,3" stroke-width="0.75" opacity="0.65"/>
    ${{ticks}}
    <g transform="translate(${{W - 445}}, 14)" font-size="10.5" fill="#515154">
      <circle cx="0" cy="-3" r="3.5" fill="#b91c1c"/><text x="7" y="0">Active pocket (&le;6&Aring;)</text>
      <circle cx="122" cy="-3" r="3.5" fill="#007c82"/><text x="129" y="0">ClinVar P/LP</text>
      <circle cx="212" cy="-3" r="3.5" fill="#515154"/><text x="219" y="0">TCGA somatic</text>
      <rect x="308" y="-7" width="7" height="7" fill="#00709b" opacity="0.4"/><text x="319" y="0">Experimental PDB</text>
    </g>
  </svg>`;
  document.getElementById('lollipop-svg-wrap').innerHTML = svg;

  ensureSatLoaded(gene, () => {{
    const rawQ = document.getElementById('sat-query').value.trim().toUpperCase();
    const gCache = window.SAT_CACHE[gene] || {{}};
    const m = rawQ.match(/^([A-Z])?(\\d+)([A-Z])?$/);
    if (!m || !gCache[m[2]] || (m[1] && gCache[m[2]].m[0] !== m[1])) {{
      const bestTrack = tr.slice().sort((a,b) => ((b[7]+b[8])*2 + b[3]*5) - ((a[7]+a[8])*2 + a[3]*5))[0];
      const defPos = bestTrack ? bestTrack[0] : 1;
      const entry = gCache[String(defPos)];
      if (entry) {{
        const topSub = entry.s.slice().sort((a,b) => (b[5]-a[5]) || ((b[1]||0)-(a[1]||0)))[0];
        document.getElementById('sat-query').value = entry.m[0] + defPos + (topSub ? topSub[0] : '');
      }}
    }}
    ensureStruct3DLoaded(gene, () => inspectSaturationCodon());
  }});
}}

function vecSub(a, b) {{ return [a[0]-b[0], a[1]-b[1], a[2]-b[2]]; }}
function vecCross(a, b) {{ return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]; }}
function vecDot(a, b) {{ return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]; }}
function vecNorm(a) {{ const l = Math.hypot(a[0], a[1], a[2]) || 1; return [a[0]/l, a[1]/l, a[2]/l]; }}

function rotateAroundAxis(pt, origin, axisUnit, rad) {{
  const p = [pt[0]-origin[0], pt[1]-origin[1], pt[2]-origin[2]];
  const c = Math.cos(rad), s = Math.sin(rad);
  const kxp = vecCross(axisUnit, p);
  const kdp = vecDot(axisUnit, p);
  return [
    origin[0] + p[0]*c + kxp[0]*s + axisUnit[0]*kdp*(1-c),
    origin[1] + p[1]*c + kxp[1]*s + axisUnit[1]*kdp*(1-c),
    origin[2] + p[2]*c + kxp[2]*s + axisUnit[2]*kdp*(1-c)
  ];
}}

function formatPdbAtom(serial, atomName, res3, chain, resSeq, x, y, z, bfac) {{
  const an = atomName.length < 4 ? (' ' + atomName).padEnd(4, ' ') : atomName.slice(0, 4);
  const el = atomName.replace(/[0-9]/g, '').slice(0, 1);
  return (
    'ATOM  ' +
    String(serial).padStart(5, ' ') + ' ' +
    an + ' ' +
    res3.padStart(3, ' ') + ' ' +
    chain +
    String(resSeq).padStart(4, ' ') + '    ' +
    x.toFixed(3).padStart(8, ' ') +
    y.toFixed(3).padStart(8, ' ') +
    z.toFixed(3).padStart(8, ' ') +
    '  1.00' +
    bfac.toFixed(2).padStart(6, ' ') +
    '          ' +
    el.padStart(2, ' ')
  );
}}

function computeMutant3DGeometry(st, pos, wtAA, mutAA) {{
  if (!st || !st._byRes[pos]) return null;
  const resRows = st._byRes[pos];
  const bb = {{}};
  resRows.forEach(r => {{ bb[r[1]] = [r[2], r[3], r[4], r[5]]; }});
  if (!bb.N || !bb.CA || !bb.C) return null;

  const ca = bb.CA;
  const envAtoms = [];
  let caNeighbors8A = 0;
  st._caList.forEach(c => {{
    if (c[0] === pos) return;
    const dCa = Math.hypot(c[1]-ca[0], c[2]-ca[1], c[3]-ca[2]);
    if (dCa <= 8.0) caNeighbors8A++;
    if (dCa <= 12.5) {{
      (st._byRes[c[0]] || []).forEach(ar => {{
        envAtoms.push([c[0], ar[1], ar[2], ar[3], ar[4]]);
      }});
    }}
  }});

  let wtMinDist = 99.9;
  resRows.forEach(r => {{
    if (['N','CA','C','O'].includes(r[1])) return;
    envAtoms.forEach(ea => {{
      const d = Math.hypot(r[2]-ea[2], r[3]-ea[3], r[4]-ea[4]);
      if (d < wtMinDist) wtMinDist = d;
    }});
  }});
  if (wtMinDist > 50) wtMinDist = null;

  const u = vecNorm(vecSub(bb.N, bb.CA));
  const w = vecNorm(vecCross(u, vecSub(bb.C, bb.CA)));
  const v = vecCross(w, u);

  const templates = (ATLAS.rotamers && ATLAS.rotamers[mutAA]) ? ATLAS.rotamers[mutAA] : [[]];
  const spinAngles = (mutAA === 'G' || mutAA === 'A' || mutAA === 'P') ? [0] : [0, -20, 20, -40, 40, 60, -60, 120, -120, 180];

  let bestAtoms = [];
  let bestScore = Infinity;
  let bestMinDist = 99.9;
  let bestClashes = [];
  let bestContacts = [];

  templates.forEach(tpl => {{
    const basePts = tpl.map(([aname, cu, cv, cw]) => [
      aname,
      ca[0] + cu*u[0] + cv*v[0] + cw*w[0],
      ca[1] + cu*u[1] + cv*v[1] + cw*w[1],
      ca[2] + cu*u[2] + cv*v[2] + cw*w[2]
    ]);
    const cbRow = basePts.find(p => p[0] === 'CB');
    const axis = cbRow ? vecNorm(vecSub([cbRow[1], cbRow[2], cbRow[3]], ca)) : u;

    spinAngles.forEach(deg => {{
      const rad = deg * Math.PI / 180.0;
      const cand = basePts.map(([aname, x, y, z]) => {{
        if (aname === 'CB' || deg === 0) return [aname, x, y, z];
        const rp = rotateAroundAxis([x, y, z], ca, axis, rad);
        return [aname, rp[0], rp[1], rp[2]];
      }});

      let clashScore = 0.0;
      let minD = 99.9;
      const clashes = [];
      const contacts = [];
      cand.forEach(([aname, x, y, z]) => {{
        envAtoms.forEach(([epos, eaname, ex, ey, ez]) => {{
          if (Math.abs(epos - pos) === 1 && aname === 'CB' && (eaname === 'N' || eaname === 'C')) return;
          const d = Math.hypot(x - ex, y - ey, z - ez);
          if (d < minD) minD = d;
          if (d < 2.40) {{
            const ov = 2.40 - d;
            clashScore += ov * ov * 25.0;
            clashes.push({{from:[x,y,z], to:[ex,ey,ez], d:d, atom:aname, nbr:`${{epos}}:${{eaname}}`}});
          }} else if (d < 3.0) {{
            clashScore += (3.0 - d) * 0.4;
          }}
          if (d >= 2.40 && d <= 3.75 && aname !== 'CB') {{
            contacts.push({{from:[x,y,z], to:[ex,ey,ez], d:d, atom:aname, nbr:`${{epos}}:${{eaname}}`}});
          }}
        }});
      }});
      clashScore += Math.abs(deg) * 0.0005;
      if (clashScore < bestScore) {{
        bestScore = clashScore;
        bestAtoms = cand;
        bestMinDist = minD > 50 ? null : minD;
        bestClashes = clashes.sort((a,b)=>a.d-b.d).slice(0, 6);
        bestContacts = contacts.sort((a,b)=>a.d-b.d).slice(0, 5);
      }}
    }});
  }});

  return {{
    bb: bb,
    mutSidechain: bestAtoms,
    caNeighbors8A: caNeighbors8A,
    wtMinDist: wtMinDist,
    mutMinDist: bestMinDist,
    clashes: bestClashes,
    contacts: bestContacts
  }};
}}

function buildPdbModels(gene, st, pos, wtAA, mutAA, geom, colorMode) {{
  const seq = st.s;
  const trMap = {{}};
  (ATLAS.tracks[gene] || []).forEach(c => {{ trMap[c[0]] = c; }});
  const mode = document.getElementById('mol3d-mode').value;

  const linesWT = [];
  let serial = 1;
  st.a.forEach(([rpos, aname, x, y, z, plddt]) => {{
    if (mode === 'mutant' && rpos === pos && !['N','CA','C','O'].includes(aname)) return;
    const aa1 = seq[rpos - 1] || 'A';
    const res3 = (AA_PROPS[aa1] ? AA_PROPS[aa1].three : 'ALA');
    let bval = plddt;
    if (colorMode === 'am') {{
      const tr = trMap[rpos];
      bval = (tr && tr[4] !== null) ? tr[4] * 100.0 : 20.0;
    }} else if (colorMode === 'pocket') {{
      const tr = trMap[rpos];
      bval = (tr && tr[3] === 1) ? 95.0 : 25.0;
    }}
    linesWT.push(formatPdbAtom(serial++, aname, res3, 'A', rpos, x, y, z, bval));
  }});
  linesWT.push('TER');
  linesWT.push('END');

  const linesMut = [];
  if (geom) {{
    const mut3 = AA_PROPS[mutAA] ? AA_PROPS[mutAA].three : 'ALA';
    let mser = 1;
    ['N','CA','C','O'].forEach(bn => {{
      if (geom.bb[bn]) {{
        const [bx, by, bz, bp] = geom.bb[bn];
        linesMut.push(formatPdbAtom(mser++, bn, mut3, 'M', pos, bx, by, bz, bp));
      }}
    }});
    geom.mutSidechain.forEach(([aname, mx, my, mz]) => {{
      linesMut.push(formatPdbAtom(mser++, aname, mut3, 'M', pos, mx, my, mz, 99.0));
    }});
    linesMut.push('TER');
    linesMut.push('END');
  }}
  return {{ wtPdb: linesWT.join('\\n'), mutPdb: linesMut.join('\\n') }};
}}

function ensureViewerInitialized() {{
  if (viewer3d) return viewer3d;
  const el = document.getElementById('mol3d-canvas');
  if (!el || typeof $3Dmol === 'undefined') return null;
  viewer3d = $3Dmol.createViewer(el, {{ backgroundColor: '#111315', antialias: true }});
  return viewer3d;
}}

function refresh3DScene(keepCamera) {{
  const v = ensureViewerInitialized();
  if (!v) return;
  const {{ gene, pos, wt, mut }} = active3DState;
  const mode = document.getElementById('mol3d-mode').value;
  const colorMode = document.getElementById('mol3d-color').value;
  const pInfo = ATLAS.panel.find(p => p.geneSymbol === gene) || {{}};

  if (mode === 'rcsb') {{
    renderExperimentalRCSB(gene, pos, wt, mut);
    return;
  }}

  ensureStruct3DLoaded(gene, st => {{
    if (!st) return;
    document.getElementById('mol3d-hud').textContent = `${{gene}} · ${{st.u}} · AlphaFold v6 (${{st.s.length}} aa)`;
    const geom = computeMutant3DGeometry(st, pos, wt, mut);
    const pdbs = buildPdbModels(gene, st, pos, wt, mut, geom, colorMode);
    active3DState.lastPdbText = pdbs.wtPdb + (mode !== 'wt' ? ('\\n' + pdbs.mutPdb) : '');

    v.removeAllModels();
    v.removeAllShapes();
    v.removeAllLabels();

    v.addModel(pdbs.wtPdb, 'pdb');
    const colorFunc = atom => {{
      const b = atom.b;
      if (colorMode === 'plddt') {{
        if (b >= 90) return '#0053d6';
        if (b >= 70) return '#65cbf3';
        if (b >= 50) return '#ffdb13';
        return '#ff7d45';
      }} else if (colorMode === 'am') {{
        if (b >= 56.4) return '#b91c1c';
        if (b >= 34.0) return '#d97706';
        return '#00709b';
      }} else {{
        return b >= 80 ? '#e11d48' : '#64748b';
      }}
    }};
    v.setStyle({{model: 0}}, {{cartoon: {{colorfunc: colorFunc, opacity: 0.86}}}});

    if (showNeighbors && geom && geom.bb.CA) {{
      const [cx, cy, cz] = geom.bb.CA;
      const nbrRes = [];
      st._caList.forEach(c => {{
        if (c[0] !== pos && Math.hypot(c[1]-cx, c[2]-cy, c[3]-cz) <= 8.5) nbrRes.push(c[0]);
      }});
      if (nbrRes.length) {{
        v.addStyle({{model: 0, resi: nbrRes}}, {{stick: {{radius: 0.10, colorscheme: 'Jmol', opacity: 0.52}}}});
      }}
    }}

    if (mode === 'overlay' || mode === 'wt') {{
      v.addStyle(
        {{model: 0, resi: pos}},
        {{
          stick: {{radius: 0.21, color: '#0095c8', opacity: mode === 'overlay' ? 0.72 : 1.0}},
          sphere: {{scale: 0.25, color: '#0095c8', opacity: mode === 'overlay' ? 0.72 : 1.0}}
        }}
      );
    }}

    if ((mode === 'overlay' || mode === 'mutant') && pdbs.mutPdb && wt !== mut) {{
      v.addModel(pdbs.mutPdb, 'pdb');
      v.setStyle(
        {{model: 1}},
        {{
          stick: {{radius: 0.24, color: '#e11d48'}},
          sphere: {{scale: 0.28, color: '#e11d48'}}
        }}
      );
      if (geom) {{
        geom.contacts.forEach(ct => {{
          v.addCylinder({{
            start: {{x:ct.from[0], y:ct.from[1], z:ct.from[2]}},
            end: {{x:ct.to[0], y:ct.to[1], z:ct.to[2]}},
            radius: 0.045, dashed: true, dashLength: 0.22, gapLength: 0.16, color: '#f59e0b'
          }});
        }});
        geom.clashes.forEach(cl => {{
          v.addCylinder({{
            start: {{x:cl.from[0], y:cl.from[1], z:cl.from[2]}},
            end: {{x:cl.to[0], y:cl.to[1], z:cl.to[2]}},
            radius: 0.075, dashed: true, dashLength: 0.18, gapLength: 0.12, color: '#ff3b30'
          }});
        }});
      }}
    }}

    if (geom && geom.bb.CA) {{
      const [lx, ly, lz] = geom.bb.CA;
      const lblTxt = wt === mut ? `${{gene}} ${{wt}}${{pos}}` : `${{gene}} ${{wt}}${{pos}}${{mut}}`;
      v.addLabel(lblTxt, {{
        position: {{x: lx, y: ly + 1.8, z: lz}},
        backgroundColor: '#1d1d1f',
        fontColor: '#f5f5f7',
        fontSize: 11.5,
        backgroundOpacity: 0.86,
        borderColor: '#3c4047',
        borderThickness: 1
      }});
    }}

    v.setClickable({{model: 0}}, true, atom => {{
      if (atom && atom.resi) selectCodonFromTrack(atom.resi);
    }});

    v.render();
    if (!keepCamera) zoomToMutatedResidue();
  }});
}}

function renderExperimentalRCSB(gene, pos, wt, mut) {{
  const v = ensureViewerInitialized();
  if (!v) return;
  const gCache = window.SAT_CACHE ? window.SAT_CACHE[gene] : null;
  const entry = gCache ? gCache[String(pos)] : null;
  const pInfo = ATLAS.panel.find(p => p.geneSymbol === gene);
  const pdbRaw = (entry && entry.m[1]) ? entry.m[1].split(':')[0] : ((pInfo && pInfo.pdbComplexIds) ? pInfo.pdbComplexIds.split(';')[0].trim() : '');
  if (!pdbRaw) {{
    document.getElementById('mol3d-mode').value = 'overlay';
    refresh3DScene();
    return;
  }}

  document.getElementById('mol3d-hud').textContent = `${{gene}} · Experimental RCSB PDB ${{pdbRaw}}`;
  const drawCif = cifTxt => {{
    v.removeAllModels();
    v.removeAllShapes();
    v.removeAllLabels();
    v.addModel(cifTxt, 'cif');
    v.setStyle({{}}, {{cartoon: {{color: '#64748b', opacity: 0.82}}}});
    v.addStyle({{resn: ['A','U','G','C','DA','DT','DG','DC']}}, {{cartoon: {{color: '#d97706'}}, stick: {{radius: 0.13, color: '#f59e0b'}}}});
    v.addStyle({{elem: ['MG','ZN','CA','MN']}}, {{sphere: {{radius: 1.3, color: '#10b981'}}}});
    v.addStyle({{resi: pos}}, {{stick: {{radius: 0.26, color: '#e11d48'}}, sphere: {{scale: 0.30, color: '#e11d48'}}}});
    v.addLabel(`PDB ${{pdbRaw}} · ${{gene}} ${{wt}}${{pos}}`, {{
      resi: pos, backgroundColor: '#1d1d1f', fontColor: '#f5f5f7', fontSize: 11.5, backgroundOpacity: 0.88
    }});
    v.zoomTo({{resi: pos}}, 420);
    v.render();
  }};

  if (window.RCSB_CACHE[pdbRaw]) {{
    drawCif(window.RCSB_CACHE[pdbRaw]);
    return;
  }}
  fetch(`https://files.rcsb.org/download/${{pdbRaw}}.cif`)
    .then(r => {{ if (!r.ok) throw new Error('RCSB fetch failed'); return r.text(); }})
    .then(txt => {{ window.RCSB_CACHE[pdbRaw] = txt; drawCif(txt); }})
    .catch(() => {{
      document.getElementById('mol3d-mode').value = 'overlay';
      refresh3DScene();
    }});
}}

function zoomToMutatedResidue() {{
  if (!viewer3d) return;
  const pos = active3DState.pos;
  viewer3d.zoomTo({{model: 0, resi: [Math.max(1, pos-4), pos, pos+4]}}, 380);
  viewer3d.zoom(0.82, 180);
}}

function resetFullProtein3D() {{
  if (!viewer3d) return;
  viewer3d.zoomTo({{model: 0}}, 420);
}}

function toggleNeighbors3D() {{
  showNeighbors = !showNeighbors;
  document.getElementById('btn-nbrs').textContent = `Shell: ${{showNeighbors ? 'On' : 'Off'}}`;
  refresh3DScene(true);
}}

function toggleSpin3D() {{
  if (!viewer3d) return;
  spinActive = !spinActive;
  viewer3d.spin(spinActive ? 'y' : false);
  document.getElementById('btn-spin').textContent = spinActive ? 'Stop' : 'Spin';
}}

function downloadMutatedPDB() {{
  if (!active3DState.lastPdbText) return;
  const {{ gene, wt, pos, mut }} = active3DState;
  const blob = new Blob([active3DState.lastPdbText], {{type: 'chemical/x-pdb'}});
  const a = document.createElement('a');
  a.href = URL.createObjectURL(blob);
  a.download = `${{gene}}_${{wt}}${{pos}}${{mut}}_3D_model.pdb`;
  a.click();
}}

function renderPropertyDeltaCard(gene, pos, wt, mut, entry, subRow) {{
  const box = document.getElementById('prop-delta-panel');
  const st = window.STRUCT3D_CACHE ? window.STRUCT3D_CACHE[gene] : null;
  const geom = st ? computeMutant3DGeometry(st, pos, wt, mut) : null;

  const wp = AA_PROPS[wt] || AA_PROPS.A;
  const mp = AA_PROPS[mut] || wp;
  const [_, pdb, dna, dmet, dprt, plddt] = entry.m;
  const [alt, am, e2, e1v, if1, tn, cv] = subRow || [mut, null, null, null, null, 0, ''];

  const dq = mp.q - wp.q;
  const dvol = mp.vol - wp.vol;
  const dkd = mp.kd - wp.kd;
  const gDist = granthamDist(wt, mut);
  const gLabel = gDist === 0 ? 'Reference' : (gDist <= 50 ? 'Conservative' : (gDist <= 100 ? 'Moderate' : (gDist <= 150 ? 'Moderately radical' : 'Radical')));
  const gDot = gDist <= 50 ? 'dot-petrol' : (gDist <= 100 ? 'dot-amber' : 'dot-crimson');

  let chargeNote = 'Preserved';
  let chargeDot = 'dot-slate';
  if (wp.q * mp.q < -0.2) {{ chargeNote = 'Charge inversion'; chargeDot = 'dot-crimson'; }}
  else if (Math.abs(dq) >= 0.8) {{ chargeNote = dq > 0 ? 'Gain of (+) charge' : 'Loss of (+) / gain of (-)'; chargeDot = 'dot-amber'; }}

  const alerts = [];
  if (mut === 'P' && wt !== 'P') alerts.push('Proline backbone rigidity');
  if (wt === 'G' && mut !== 'G') alerts.push('Loss of glycine flexibility');
  if ((wt === 'C') !== (mut === 'C')) alerts.push('Thiol / disulfide change');
  if (wp.metal && !mp.metal && dmet !== null && dmet <= 6.0) alerts.push('Metal-chelator loss');
  if (!wp.arom && mp.arom) alerts.push('Aromatic ring gain');

  const constrained = (am !== null && am >= 0.564) || (e2 !== null && e2 <= -5.0);
  const regimeTag = wt === mut
    ? tagDot('dot-emerald', 'Wild-type reference')
    : (constrained
        ? ((if1 !== null && if1 > -1.5) ? tagDot('dot-crimson', 'Fold-preserving active / interface') : tagDot('dot-amber', 'Predicted fold-destabilising'))
        : tagDot('dot-petrol', 'Tolerated / neutral'));

  const burialStr = geom
    ? (geom.caNeighbors8A >= 13 ? `Buried (${{geom.caNeighbors8A}} C&alpha; &le;8&Aring;)` : (geom.caNeighbors8A >= 8 ? `Intermediate / cleft (${{geom.caNeighbors8A}} C&alpha;)` : `Solvent-exposed (${{geom.caNeighbors8A}} C&alpha;)`))
    : '&mdash;';
  const clashCount = geom ? geom.clashes.length : 0;
  const clashTag = clashCount > 0
    ? tagDot('dot-crimson', `${{clashCount}} clash${{clashCount>1?'es':''}} (<2.4Å; min ${{geom.mutMinDist ? geom.mutMinDist.toFixed(2) : '?'}}Å)`)
    : tagDot('dot-emerald', `0 clashes (min ${{geom && geom.mutMinDist ? geom.mutMinDist.toFixed(2) : '>2.5'}}Å)`);

  const aaCells = entry.s.map(r => {{
    const a = r[0], rAm = r[1], rE2 = r[2], rIf1 = r[4];
    const cns = (rAm !== null && rAm >= 0.564) || (rE2 !== null && rE2 <= -5.0);
    const topCol = cns ? ((rIf1 !== null && rIf1 > -1.5) ? '#b91c1c' : '#d97706') : '#00709b';
    return `<button class="aa-cell ${{a===mut?'active':''}}" style="border-top-color:${{topCol}};" onclick="selectMutAA(${{pos}},'${{wt}}','${{a}}')" title="${{wt}}${{pos}}${{a}} · AM=${{rAm??'NA'}} · IF1=${{rIf1??'NA'}}">${{a}}</button>`;
  }}).join('');

  let mechTxt = '';
  if (wt === mut) {{
    mechTxt = `Wild-type ${{wp.name}} (${{wt}}${{pos}}) occupies a ${{burialStr.toLowerCase().replace(/<[^>]+>/g,'')}} site (AlphaFold pLDDT ${{plddt ?? 'NA'}}). Select any substitution above to compute the mutant 3D side-chain rotamer and physicochemical shifts.`;
  }} else {{
    const pocketParts = [];
    if (dmet !== null && dmet <= 6.0) pocketParts.push(`${{dmet}} Å from catalytic metal`);
    if (dna !== null && dna <= 6.0) pocketParts.push(`${{dna}} Å from bound RNA/DNA`);
    if (dprt !== null && dprt <= 6.0) pocketParts.push(`${{dprt}} Å from cofactor subunit`);
    const locPhrase = pocketParts.length ? `at the experimental active site (${{pocketParts.join(', ')}} in PDB ${{pdb}})` : `in a ${{burialStr.toLowerCase().replace(/<[^>]+>/g,'')}} environment`;
    const stabPhrase = (if1 !== null && if1 <= -1.5)
      ? `is predicted to destabilise local backbone fold thermodynamics (ESM-IF1 ΔlogP = ${{if1.toFixed(2)}}${{clashCount>0 ? ', ' + clashCount + ' steric clash(es) <2.4 Å' : ''}})`
      : (constrained
          ? `maintains backbone fold stability (ESM-IF1 ΔlogP = ${{if1 !== null ? if1.toFixed(2) : 'NA'}}, ${{clashCount}} steric clashes) while altering active-site or interface chemistry`
          : `is predicted tolerated across evolutionary conservation (AlphaMissense ${{am !== null ? am.toFixed(3) : 'NA'}}) and 3D inverse folding`);
    mechTxt = `<strong>${{gene}} p.${{wp.three}}${{pos}}${{mp.three}} (${{wt}}${{pos}}${{mut}})</strong> substitutes ${{wp.cls.toLowerCase()}} ${{wp.name}} with ${{mp.cls.toLowerCase()}} ${{mp.name}} ${{locPhrase}} (Δq = ${{dq>=0?'+':''}}${{dq.toFixed(1)}} e, ΔV = ${{dvol>=0?'+':''}}${{dvol.toFixed(1)}} Å³, Grantham D = ${{gDist}}), and ${{stabPhrase}}.`;
  }}

  box.innerHTML = `
    <div style="display:flex;justify-content:space-between;align-items:flex-start;flex-wrap:wrap;gap:8px;">
      <div>
        <div style="font-size:15.5px;font-weight:650;letter-spacing:-0.015em;">
          ${{gene}} <code>p.${{wp.three}}${{pos}}${{mp.three}}</code> <span style="color:var(--ink-tertiary);font-weight:400;">(${{wt}}${{pos}}${{mut}})</span>
        </div>
        <div style="font-size:11.5px;color:var(--ink-tertiary);margin-top:2px;">
          TCGA somatic: <strong>${{tn}}</strong> &middot; ClinVar: <strong>${{cv || 'Unreported'}}</strong> &middot; PDB: <strong>${{pdb || 'AlphaFold v6'}}</strong>
        </div>
      </div>
      <div>${{regimeTag}}</div>
    </div>

    <div style="margin-top:11px;font-size:11px;color:var(--ink-tertiary);display:flex;justify-content:space-between;flex-wrap:wrap;">
      <span>Substitute codon ${{wt}}${{pos}}:</span>
      <span><i class="dot dot-crimson"></i> Active/interface &nbsp;<i class="dot dot-amber"></i> Destabilising &nbsp;<i class="dot dot-petrol"></i> Neutral</span>
    </div>
    <div class="aa-strip">${{aaCells}}</div>

    <div class="table-wrap" style="max-height:none;">
      <table>
        <thead>
          <tr><th>Property</th><th>WT (${{wt}})</th><th>Mut (${{mut}})</th><th>Shift &amp; Assessment</th></tr>
        </thead>
        <tbody>
          <tr>
            <td><strong>Chemical group</strong></td>
            <td>${{wp.cls}}</td>
            <td>${{mp.cls}}</td>
            <td>${{wp.cls === mp.cls ? 'Conserved group' : 'Group transition'}}</td>
          </tr>
          <tr>
            <td><strong data-term="delta_q">Net charge (pH 7.4)</strong></td>
            <td>${{wp.q >= 0 ? '+' : ''}}${{wp.q.toFixed(1)}} e</td>
            <td>${{mp.q >= 0 ? '+' : ''}}${{mp.q.toFixed(1)}} e</td>
            <td>${{dq >= 0 ? '+' : ''}}${{dq.toFixed(1)}} e &middot; ${{tagDot(chargeDot, chargeNote)}}</td>
          </tr>
          <tr>
            <td><strong data-term="delta_vol">Side-chain volume</strong></td>
            <td>${{wp.vol.toFixed(1)}} &Aring;&sup3;</td>
            <td>${{mp.vol.toFixed(1)}} &Aring;&sup3;</td>
            <td>${{dvol >= 0 ? '+' : ''}}${{dvol.toFixed(1)}} &Aring;&sup3; (${{dvol >= 0 ? '+' : ''}}${{((dvol/wp.vol)*100).toFixed(0)}}%)</td>
          </tr>
          <tr>
            <td><strong data-term="hydropathy_kd">Hydropathy (KD)</strong></td>
            <td>${{wp.kd.toFixed(1)}}</td>
            <td>${{mp.kd.toFixed(1)}}</td>
            <td>${{dkd >= 0 ? '+' : ''}}${{dkd.toFixed(1)}} (${{dkd > 1.5 ? 'More hydrophobic' : (dkd < -1.5 ? 'More polar' : 'Comparable')}})</td>
          </tr>
          <tr>
            <td><strong data-term="dist_metal">H-bond / metal</strong></td>
            <td>${{wp.hd}}D / ${{wp.ha}}A${{wp.metal ? ' · Metal' : ''}}</td>
            <td>${{mp.hd}}D / ${{mp.ha}}A${{mp.metal ? ' · Metal' : ''}}</td>
            <td>&Delta;D=${{mp.hd-wp.hd>=0?'+':''}}${{mp.hd-wp.hd}}, &Delta;A=${{mp.ha-wp.ha>=0?'+':''}}${{mp.ha-wp.ha}}</td>
          </tr>
          <tr>
            <td><strong data-term="grantham">Grantham / alerts</strong></td>
            <td colspan="2">D = <strong>${{gDist}}</strong> &middot; ${{tagDot(gDot, gLabel)}}</td>
            <td>${{alerts.length ? alerts.join('; ') : 'None'}}</td>
          </tr>
          <tr>
            <td><strong data-term="burial_8a">3D environment</strong></td>
            <td>${{burialStr}}</td>
            <td>WT ${{geom && geom.wtMinDist ? geom.wtMinDist.toFixed(2)+'&Aring;' : '&mdash;'}}</td>
            <td data-term="rotamer_clash">${{clashTag}}</td>
          </tr>
          <tr>
            <td><strong data-term="esm_if1">Model scores</strong></td>
            <td colspan="2"><span data-term="alphamissense">AM ${{am !== null ? am.toFixed(3) : '&mdash;'}}</span> &middot; <span data-term="esm2_llr">ESM-2 ${{e2 !== null ? e2.toFixed(2) : '&mdash;'}}</span></td>
            <td data-term="esm_if1">ESM-IF1 &Delta;logP: <strong>${{if1 !== null ? if1.toFixed(2) : '&mdash;'}}</strong></td>
          </tr>
        </tbody>
      </table>
    </div>
    <div class="interp-box">${{mechTxt}}</div>
  `;
}}

function inspectSaturationCodon() {{
  const gene = document.getElementById('lol-gene').value;
  const rawQ = document.getElementById('sat-query').value.trim().toUpperCase();
  const box = document.getElementById('sat-codon-detail');
  if (!rawQ) {{
    box.innerHTML = '';
    return;
  }}
  const gCache = window.SAT_CACHE ? window.SAT_CACHE[gene] : null;
  if (!gCache) {{
    ensureSatLoaded(gene, inspectSaturationCodon);
    return;
  }}
  const mVar = rawQ.match(/^([A-Z])?(\\d+)([A-Z*])?$/);
  if (!mVar) {{
    box.innerHTML = `<p style="font-size:12.5px;color:var(--crimson);">Enter a residue position (e.g. <code>1705</code>) or missense change (e.g. <code>E1705K</code>).</p>`;
    return;
  }}
  const posStr = mVar[2];
  const pos = parseInt(posStr, 10);
  const entry = gCache[posStr];
  if (!entry) {{
    box.innerHTML = `<p style="font-size:12.5px;color:var(--crimson);">Position ${{pos}} exceeds ${{gene}} canonical length (1&ndash;${{Object.keys(gCache).length}} aa).</p>`;
    return;
  }}
  const [wt, pdb, dna, dmet, dprt, plddt] = entry.m;
  let targetAlt = (mVar[3] && mVar[3] !== '*') ? mVar[3] : '';

  let activeSubRow = null;
  if (targetAlt && targetAlt !== wt) {{
    activeSubRow = entry.s.find(r => r[0] === targetAlt) || null;
  }}
  if (!activeSubRow) {{
    activeSubRow = entry.s.slice().sort((a,b) => (b[5] - a[5]) || ((b[1]||0) - (a[1]||0)))[0];
  }}
  const mut3D = activeSubRow ? activeSubRow[0] : wt;

  active3DState = {{ gene: gene, pos: pos, wt: wt, mut: mut3D, lastPdbText: active3DState.lastPdbText }};
  ensureStruct3DLoaded(gene, () => {{
    renderPropertyDeltaCard(gene, pos, wt, mut3D, entry, activeSubRow);
    refresh3DScene(false);
  }});

  let rows = entry.s;
  if (targetAlt && targetAlt !== wt) rows = rows.filter(r => r[0] === targetAlt);

  let html = `<div>
    <div style="display:flex;justify-content:space-between;align-items:center;flex-wrap:wrap;gap:8px;margin-bottom:8px;font-size:12.5px;color:var(--ink-secondary);">
      <span>
        <strong>Codon ${{gene}} ${{wt}}${{pos}} saturation profile</strong> &middot;
        pLDDT ${{plddt ?? 'NA'}} &middot;
        PDB ${{pdb || 'AlphaFold v6'}} &middot;
        RNA/DNA ${{dna !== null ? dna + ' Å' : '—'}} &middot;
        Metal ${{dmet !== null ? dmet + ' Å' : '—'}} &middot;
        Partner ${{dprt !== null ? dprt + ' Å' : '—'}}
      </span>
      ${{targetAlt ? `<button class="btn-quiet" onclick="document.getElementById('sat-query').value='${{pos}}';inspectSaturationCodon();">Show all 19 substitutions</button>` : ''}}
    </div>
    <div class="table-wrap" style="max-height:240px;">
      <table>
        <thead><tr><th>3D</th><th>Variant</th><th data-term="grantham">Grantham</th><th data-term="delta_q">&Delta;q</th><th data-term="delta_vol">&Delta;Vol (&Aring;&sup3;)</th><th data-term="alphamissense">AlphaMissense</th><th data-term="esm2_llr">ESM-2 3B</th><th data-term="esm1v_llr">ESM-1v</th><th data-term="esm_if1">ESM-IF1 3D</th><th data-term="mech_regime">Regime</th><th>TCGA</th><th data-term="kpi_auc">ClinVar</th></tr></thead>
        <tbody>`;
  rows.forEach(r => {{
    const [alt, am, e2, e1v, if1, tn, cv] = r;
    const wp = AA_PROPS[wt] || AA_PROPS.A, mp = AA_PROPS[alt] || wp;
    const dq = mp.q - wp.q, dv = mp.vol - wp.vol, gd = granthamDist(wt, alt);
    const constrained = (am !== null && am >= 0.564) || (e2 !== null && e2 <= -5.0);
    const regime = constrained ? ((if1 !== null && if1 > -1.5) ? tagDot('dot-crimson','Active / interface') : tagDot('dot-amber','Fold-destabilising')) : tagDot('dot-petrol','Neutral');
    html += `<tr>
      <td><button class="btn-quiet" style="height:23px;padding:0 7px;font-size:11px;" onclick="selectMutAA(${{pos}},'${{wt}}','${{alt}}')">View</button></td>
      <td><code>${{wt}}${{pos}}${{alt}}</code></td>
      <td>${{gd}}</td>
      <td>${{dq>=0?'+':''}}${{dq.toFixed(1)}}</td>
      <td>${{dv>=0?'+':''}}${{dv.toFixed(0)}}</td>
      <td>${{am !== null ? am.toFixed(3) : '&mdash;'}}</td>
      <td>${{e2 !== null ? e2.toFixed(2) : '&mdash;'}}</td>
      <td>${{e1v !== null ? e1v.toFixed(2) : '&mdash;'}}</td>
      <td>${{if1 !== null ? if1.toFixed(2) : '&mdash;'}}</td>
      <td>${{regime}}</td>
      <td>${{tn > 0 ? '<strong>' + tn + '</strong>' : '0'}}</td>
      <td>${{cv || '&mdash;'}}</td>
    </tr>`;
  }});
  html += `</tbody></table></div></div>`;
  box.innerHTML = html;
}}

function tierPill(t) {{
  if (t === '3D_Catalytic_or_Substrate_Pocket_Driver') return tagDot('dot-crimson', '3D pocket (≤6Å)');
  if (t === 'Pathogenic_Fold_or_Interface_Missense') return tagDot('dot-amber', 'Pathogenic fold');
  if (t === 'Truncating_LoF') return tagDot('dot-slate', 'Truncating LoF');
  if (t === 'Ambiguous_Missense_VUS') return tagDot('dot-slate', 'Ambiguous VUS');
  if (t === 'Likely_Benign_Passenger_Missense') return tagDot('dot-petrol', 'Likely benign');
  return tagDot('dot-slate', t);
}}

function mechPill(m) {{
  if (m === 'Fold_Preserving_Catalytic_or_Interface') return tagDot('dot-crimson', 'Active / interface');
  if (m === 'Thermodynamic_Fold_Destabilizing') return tagDot('dot-amber', 'Fold-destabilising');
  if (m === 'Tolerated_Neutral') return tagDot('dot-petrol', 'Neutral');
  return '&mdash;';
}}

let currentFilteredVars = [];
function renderVariants() {{
  const g = document.getElementById('var-gene').value;
  const tier = document.getElementById('var-tier').value;
  const mech = document.getElementById('var-mech').value;
  const q = document.getElementById('var-q').value.trim().toUpperCase();
  currentFilteredVars = ATLAS.variants.filter(v => {{
    if (g !== 'ALL' && v.g !== g) return false;
    if (tier !== 'ALL' && v.tier !== tier) return false;
    if (mech !== 'ALL' && v.mech !== mech) return false;
    if (q && !(v.p.toUpperCase().includes(q) || v.g.toUpperCase().includes(q) || v.cohorts.toUpperCase().includes(q) || v.cv.toUpperCase().includes(q))) return false;
    return true;
  }});
  document.getElementById('var-count').textContent = `Showing ${{Math.min(250, currentFilteredVars.length)}} of ${{currentFilteredVars.length}} unique variants (${{currentFilteredVars.reduce((a,b)=>a+b.n,0)}} somatic calls)`;
  document.querySelector('#var-table tbody').innerHTML = currentFilteredVars.slice(0, 250).map(v => {{
    const cleanP = v.p.replace(/^p\\./i, '');
    const isMissense = /^[A-Z]\\d+[A-Z]$/.test(cleanP);
    const jumpArg = isMissense ? cleanP : (v.pos > 0 ? String(v.pos) : '');
    const btn3d = jumpArg ? `<button class="btn-quiet" style="height:23px;padding:0 7px;font-size:11px;" onclick="openVariantIn3D('${{v.g}}','${{jumpArg}}')">View</button>` : '&mdash;';
    return `
    <tr>
      <td>${{btn3d}}</td>
      <td><strong>${{v.g}}</strong></td>
      <td><code>${{v.p}}</code></td>
      <td>${{v.n}}</td>
      <td>${{v.cv || '&mdash;'}}</td>
      <td>${{tierPill(v.tier)}}</td>
      <td>${{mechPill(v.mech)}}</td>
      <td><code>${{v.pdb || '&mdash;'}}</code></td>
      <td>${{v.dna !== null ? v.dna.toFixed(2) : '&mdash;'}}</td>
      <td>${{v.dmet !== null ? v.dmet.toFixed(2) : '&mdash;'}}</td>
      <td>${{v.am !== null ? v.am.toFixed(3) : '&mdash;'}}</td>
      <td>${{v.esm !== null ? v.esm.toFixed(2) : '&mdash;'}}</td>
      <td>${{v.esm1v !== null ? v.esm1v.toFixed(2) : '&mdash;'}}</td>
      <td>${{v.if1 !== null ? v.if1.toFixed(2) : '&mdash;'}}</td>
      <td>${{v.cohorts}}</td>
    </tr>`;
  }}).join('');
}}

function renderGtexBars() {{
  const g = document.getElementById('gtex-gene').value;
  const vals = ATLAS.gtex_tpm[g] || [];
  const pairs = ATLAS.gtex_tissues.map((t, i) => [t, vals[i] || 0]).sort((a,b) => b[1] - a[1]);
  const mx = Math.max(1, ...pairs.map(x => x[1]));
  document.getElementById('gtex-bars').innerHTML = '<table><thead><tr><th>GTEx v10 Adult Tissue (54 Sites)</th><th data-term="gtex_tpm">Median TPM</th><th>Relative Expression</th></tr></thead><tbody>' +
    pairs.map(([t, v]) => `<tr>
      <td>${{t.replace(/_/g, ' ')}}</td>
      <td><strong>${{v.toFixed(2)}}</strong></td>
      <td><div style="width:${{Math.max(2, (v/mx)*280).toFixed(0)}}px;height:8px;background:var(--embl-petrol);border-radius:2px;"></div></td>
    </tr>`).join('') + '</tbody></table>';
}}

function renderFormulas() {{
  const cat = document.getElementById('formula-cat').value;
  const q = (document.getElementById('formula-q').value || '').trim().toLowerCase();
  const list = (ATLAS.formulas || []).filter(f => {{
    if (cat !== 'ALL' && f.cat !== cat) return false;
    if (q && !(`${{f.id}} ${{f.title}} ${{f.eq}} ${{f.vars}} ${{f.meaning}} ${{f.script}}`.toLowerCase().includes(q))) return false;
    return true;
  }});
  document.getElementById('formula-cards-wrap').innerHTML = list.map(f => `
    <div class="formula-card" id="${{f.id}}">
      <div class="fc-top">
        <span class="fc-id">${{f.id}}</span>
        <span class="fc-cat">${{f.cat}}</span>
      </div>
      <div class="fc-title">${{f.title}}</div>
      <div class="fc-eq">${{f.eq}}</div>
      <div class="fc-vars"><strong>Variables:</strong> ${{f.vars}}</div>
      <div class="fc-meaning">${{f.meaning}}</div>
      <div class="fc-script">Source: ${{f.script}}</div>
    </div>
  `).join('') || '<p style="font-size:12.5px;color:var(--ink-tertiary);">No formulas match the current filter.</p>';
}}

function renderDatasetsCatalogue() {{
  const q = (document.getElementById('ds-q').value || '').trim().toLowerCase();
  const list = (ATLAS.datasets || []).filter(d => !q || (`${{d.file}} ${{d.pk}} ${{d.desc}}`.toLowerCase().includes(q)));
  document.querySelector('#datasets-cat-table tbody').innerHTML = list.map(d => `
    <tr>
      <td><code>${{d.file}}</code></td>
      <td><strong>${{Number(d.rows).toLocaleString()}}</strong> &times; ${{d.cols}}</td>
      <td>${{d.size}}</td>
      <td><code>${{d.pk}}</code></td>
      <td style="white-space:normal;min-width:340px;">${{d.desc}}</td>
    </tr>
  `).join('');
}}

function renderGlossaryCards() {{
  const q = (document.getElementById('gloss-q').value || '').trim().toLowerCase();
  const entries = Object.entries(ATLAS.glossary || {{}})
    .map(([k, v]) => ({{key: k, ...v}}))
    .sort((a, b) => a.title.localeCompare(b.title))
    .filter(g => !q || (`${{g.title}} ${{g.cat}} ${{g.meaning}} ${{g.formula}} ${{g.thresholds}}`.toLowerCase().includes(q)));
  document.getElementById('gloss-count').textContent = `Showing ${{entries.length}} of ${{Object.keys(ATLAS.glossary || {{}}).length}} defined terms`;
  document.getElementById('glossary-cards-wrap').innerHTML = entries.map(g => `
    <div class="gloss-card">
      <div style="display:flex;justify-content:space-between;align-items:baseline;gap:8px;margin-bottom:4px;">
        <strong style="font-size:13px;color:var(--ink);">${{g.title}}</strong>
        <span style="font-family:var(--mono);font-size:10.5px;color:var(--embl-petrol);">${{g.cat}}</span>
      </div>
      <div style="font-size:12.3px;color:var(--ink-secondary);margin-bottom:6px;">${{g.meaning}}</div>
      <div style="font-family:var(--mono);font-size:11px;background:var(--surface);border:1px solid var(--line-light);border-radius:5px;padding:5px 8px;margin-bottom:5px;white-space:pre-wrap;word-break:break-word;">${{g.formula}}</div>
      <div style="font-size:11.5px;color:var(--ink-tertiary);"><strong>Interpretation:</strong> ${{g.thresholds}}</div>
    </div>
  `).join('');
}}

function renderAllTables() {{
  if (ATLAS.top_tumours && document.querySelector('#cohort-burden-table tbody')) {{
    document.querySelector('#cohort-burden-table tbody').innerHTML = ATLAS.top_tumours.map(r => {{
      const n = Number(r.N);
      const c10 = Number(r.patient_prop_core10_strict);
      const n10 = Math.round(c10 * n);
      return `<tr>
        <td>${{r.rank_core10_strict}}</td>
        <td><strong>${{r.tumour_acronym}}${{r.low_coverage_N_lt_100 === '1' ? '*' : ''}}</strong></td>
        <td>${{n.toLocaleString()}}</td>
        <td><strong>${{(c10 * 100).toFixed(1)}}%</strong> (${{n10}})</td>
        <td>${{(Number(r.patient_prop_extended12_strict) * 100).toFixed(1)}}%</td>
        <td>${{(Number(r.patient_prop_panel14_strict) * 100).toFixed(1)}}%</td>
        <td>${{(Number(r.patient_prop_core10_any_alter) * 100).toFixed(1)}}%</td>
        <td>${{(Number(r.ago1_3_4_1p34_concordant_cna_prop) * 100).toFixed(1)}}%</td>
      </tr>`;
    }}).join('');
  }}

  if (ATLAS.collateral && document.querySelector('#collat-table tbody')) {{
    document.querySelector('#collat-table tbody').innerHTML = ATLAS.collateral.map(c => `
      <tr>
        <td><code>${{c.test_id}}</code></td>
        <td><strong>${{c.crispr_target_gene}}</strong></td>
        <td style="white-space:normal;min-width:280px;">${{c.hypothesis}}</td>
        <td>${{c.n_altered_lines}} (${{Number(c.median_chronos_altered).toFixed(4)}})</td>
        <td>${{c.n_control_lines}} (${{Number(c.median_chronos_control).toFixed(4)}})</td>
        <td><strong>${{(Number(c.chronos_shift_altered_minus_control)>=0?'+':'') + Number(c.chronos_shift_altered_minus_control).toFixed(4)}}</strong></td>
        <td>${{c.p_value_two_sided}}</td>
      </tr>
    `).join('');
  }}

  if (ATLAS.excluded_cands && document.querySelector('#excluded-cands-table tbody')) {{
    document.querySelector('#excluded-cands-table tbody').innerHTML = ATLAS.excluded_cands.map(e => `
      <tr>
        <td><strong>${{e.candidateGene}}</strong></td>
        <td><code>${{e.cytoband}}</code></td>
        <td>${{e.molecularFunction}}</td>
        <td style="white-space:normal;">${{e.exclusionRationale}}</td>
      </tr>
    `).join('');
  }}

  if (ATLAS.consortia && document.querySelector('#consortia-sources-table tbody')) {{
    document.querySelector('#consortia-sources-table tbody').innerHTML = ATLAS.consortia.map(c => `
      <tr>
        <td><strong>${{c.consortium}}</strong></td>
        <td>${{c.scale}}</td>
        <td style="white-space:normal;min-width:320px;">${{c.layers}}</td>
        <td><code>${{c.join_key}}</code></td>
      </tr>
    `).join('');
  }}

  if (ATLAS.study_summary && document.querySelector('#study-decoder-table tbody')) {{
    document.querySelector('#study-decoder-table tbody').innerHTML = ATLAS.study_summary.map(s => `
      <tr>
        <td><strong>${{s.tumour_acronym}}${{s.low_coverage_flag_N_lt_100 === '1' ? '*' : ''}}</strong></td>
        <td>${{s.studyName.replace(' (TCGA, PanCancer Atlas)', '')}}</td>
        <td><strong>${{Number(s.intersection_sequenced_cna_n).toLocaleString()}}</strong></td>
        <td>${{Number(s.sampleCount_sequenced).toLocaleString()}}</td>
        <td>${{Number(s.sampleCount_cna).toLocaleString()}}</td>
        <td>${{Number(s.sampleCount_mrna_matched_n).toLocaleString()}}</td>
        <td>${{Number(s.sampleCount_mirnaseq_matched_n).toLocaleString()}}</td>
        <td>${{s.hypermutated_samples_n}}</td>
        <td>${{s.nonhypermutated_samples_n}}</td>
      </tr>
    `).join('');
  }}
  document.querySelector('#mir-bench-table tbody').innerHTML = ATLAS.mir_bench.map(b => `
    <tr>
      <td><strong>${{b.benchmark_id}}</strong></td>
      <td><code>${{b.metric}}</code></td>
      <td>${{b.group_A}} (n=${{b.n_group_A}})</td>
      <td>${{Number(b.median_group_A).toFixed(4)}}</td>
      <td>${{b.group_B}} (n=${{b.n_group_B}})</td>
      <td>${{Number(b.median_group_B).toFixed(4)}}</td>
      <td><strong>${{(Number(b.median_shift_A_minus_B)>=0?'+':'') + Number(b.median_shift_A_minus_B).toFixed(4)}}</strong></td>
      <td>${{b.p_value_two_sided}}</td>
    </tr>
  `).join('');

  const sortedDicer = [...ATLAS.dicer_pts].sort((a,b) => a.d53 - b.d53);
  document.querySelector('#dicer-mir-table tbody').innerHTML = sortedDicer.map(d => `
    <tr>
      <td><code>${{d.s}}</code></td>
      <td><strong>${{d.t}}</strong></td>
      <td>${{d.c === 'DICER1_RNaseIIIb_3D_Catalytic_Pocket' ? tagDot('dot-crimson','RNase IIIb 3D pocket') : d.c}}</td>
      <td><strong>${{d.d53.toFixed(4)}}</strong></td>
      <td>${{d.f5 !== null ? d.f5.toFixed(4) : '&mdash;'}}</td>
      <td>${{d.l7z !== null && d.l7z !== undefined ? (d.l7z>=0?'+':'') + d.l7z.toFixed(3) : '&mdash;'}}</td>
    </tr>
  `).join('');

  document.querySelector('#top-axes-table tbody').innerHTML = ATLAS.top_axes.map(ax => `
    <tr>
      <td>${{ax.rank}}</td>
      <td><code>${{ax.matureMiRNA}}</code></td>
      <td><strong>${{ax.targetGeneSymbol}}</strong></td>
      <td>${{ax.n_dual_tumours}}</td>
      <td><strong>${{Number(ax.pearson_r_log2).toFixed(4)}}</strong></td>
      <td>${{ax.t_statistic}}</td>
    </tr>
  `).join('');

  document.querySelector('#clinvar-roc-table tbody').innerHTML = ATLAS.clinvar_roc.map(r => `
    <tr>
      <td><code>${{r.benchmark_subset}}</code></td>
      <td><strong>${{r.predictor}}</strong></td>
      <td>${{r.n_pos}}</td>
      <td>${{r.n_neg}}</td>
      <td><strong>${{Number(r.roc_auc).toFixed(4)}}</strong> [${{Number(r.roc_ci95_low).toFixed(3)}}&ndash;${{Number(r.roc_ci95_high).toFixed(3)}}]</td>
      <td>${{Number(r.pr_auc).toFixed(4)}} (${{Number(r.pr_baseline).toFixed(3)}})</td>
      <td>${{r.mann_whitney_p}}</td>
    </tr>
  `).join('');

  document.querySelector('#msk-rep-table tbody').innerHTML = ATLAS.msk_rep.map(m => `
    <tr>
      <td><strong>${{m.geneSymbol}}</strong></td>
      <td>${{Number(m.msk_N_profiled).toLocaleString()}}</td>
      <td><strong>${{(Number(m.msk_mut_freq)*100).toFixed(2)}}%</strong> (${{m.msk_mut_samples_n}})</td>
      <td>${{(Number(m.tcga_mut_freq_N10195)*100).toFixed(2)}}%</td>
      <td>${{(Number(m.msk_nonhypermut_mut_freq)*100).toFixed(2)}}%</td>
      <td>${{(Number(m.msk_truncating_freq)*100).toFixed(2)}}% (${{m.msk_truncating_n}})</td>
      <td><strong>${{(Number(m.msk_3dpocket_or_hotspot_freq)*100).toFixed(2)}}%</strong> (${{m.msk_3dpocket_or_hotspot_n}})</td>
      <td>${{(Number(m.msk_primary_mut_freq)*100).toFixed(2)}}% / ${{(Number(m.msk_metastasis_mut_freq)*100).toFixed(2)}}%</td>
    </tr>
  `).join('');

  document.querySelector('#msk-dicer-table tbody').innerHTML = ATLAS.msk_dicer.map(d => `
    <tr>
      <td><strong>${{d.msk_cancer_type}}</strong></td>
      <td>${{d.N_samples}}</td>
      <td>${{d.dicer1_mut_n}} (${{(Number(d.dicer1_mut_freq)*100).toFixed(2)}}%)</td>
      <td><strong>${{d.dicer1_rnase3b_hotspot_n}} (${{(Number(d.dicer1_rnase3b_hotspot_freq)*100).toFixed(2)}}%)</strong></td>
      <td>${{d.dicer1_truncating_n}} (${{(Number(d.dicer1_truncating_freq)*100).toFixed(2)}}%)</td>
    </tr>
  `).join('');

  document.querySelector('#gene-ref-table tbody').innerHTML = ATLAS.panel.map(p => `
    <tr>
      <td><strong>${{p.geneSymbol}}</strong></td>
      <td>${{p.gnomad_pLI !== '' ? Number(p.gnomad_pLI).toFixed(3) : '&mdash;'}}</td>
      <td><strong>${{p.gnomad_LOEUF !== '' ? Number(p.gnomad_LOEUF).toFixed(3) : '&mdash;'}}</strong></td>
      <td>${{p.gnomad_mis_z !== '' ? Number(p.gnomad_mis_z).toFixed(2) : '&mdash;'}}</td>
      <td>${{p.hpa_main_subcellular_location || '&mdash;'}}</td>
      <td>${{p.gtex_pan54_median_tpm !== '' ? Number(p.gtex_pan54_median_tpm).toFixed(2) : '&mdash;'}}</td>
      <td>${{p.gtex_top_tissue}} (${{Number(p.gtex_top_tissue_tpm||0).toFixed(1)}})</td>
      <td>${{p.gtex_tau_tissue_specificity}}</td>
      <td>${{p.opentargets_sm_tractable === '1' ? tagDot('dot-emerald','Small-molecule') + ' ' : ''}}${{p.opentargets_protac_tractable === '1' ? tagDot('dot-petrol','PROTAC') : ''}}</td>
      <td>${{p.clinvar_pathogenic_or_likely_pathogenic_n}} / ${{p.clinvar_total_single_gene_n}}</td>
    </tr>
  `).join('');

  const cpMap = {{}}; ATLAS.cptac.forEach(c => {{ cpMap[c.geneSymbol] = c; }});
  const dmMap = {{}}; ATLAS.depmap.forEach(d => {{ dmMap[d.geneSymbol] = d; }});
  document.querySelector('#consortia-table tbody').innerHTML = ATLAS.panel.map(p => {{
    const g = p.geneSymbol, c = cpMap[g] || {{}}, d = dmMap[g] || {{}};
    return `<tr>
      <td><strong>${{g}}</strong></td>
      <td>${{p.panelRole}}</td>
      <td><code>${{p.pdbComplexIds}}</code></td>
      <td>${{c.spearman_CNA_to_mRNA || '&mdash;'}}</td>
      <td>${{c.spearman_CNA_to_Protein || '&mdash;'}}</td>
      <td>${{c.post_transcriptional_CNA_attenuation || '&mdash;'}}</td>
      <td><strong>${{d.median_chronos_gene_effect}}</strong></td>
      <td>${{(Number(d.prop_dependent_le_m05)*100).toFixed(1)}}%</td>
      <td>${{p.depmapIsCommonEssential === '1' ? tagDot('dot-crimson','Pan-essential') : tagDot('dot-petrol','Context / buffered')}}</td>
    </tr>`;
  }}).join('');

  document.querySelector('#surv-table tbody').innerHTML = ATLAS.surv_bench.map(s => `
    <tr>
      <td><strong>${{s.test_id}}</strong></td>
      <td>${{s.cohort}}</td>
      <td><code>${{s.endpoint}}</code></td>
      <td>${{s.altered_group}} (${{s.events_altered}} / ${{s.n_altered}})</td>
      <td>${{s.control_group}} (${{s.events_control}} / ${{s.n_control}})</td>
      <td><strong>${{Number(s.hazard_ratio_stratified).toFixed(3)}}</strong> [${{Number(s.hr_ci95_low).toFixed(3)}}&ndash;${{Number(s.hr_ci95_high).toFixed(3)}}]</td>
      <td>${{s.p_value_two_sided}}</td>
      <td>${{s.multivar_cox_hr !== '' && s.multivar_cox_hr !== undefined ? '<strong>' + Number(s.multivar_cox_hr).toFixed(3) + '</strong> [' + Number(s.multivar_ci95_low).toFixed(3) + '&ndash;' + Number(s.multivar_ci95_high).toFixed(3) + ']' : '&mdash;'}}</td>
      <td>${{s.multivar_p_value || '&mdash;'}}</td>
    </tr>
  `).join('');

  document.querySelector('#imm-table tbody').innerHTML = ATLAS.imm_bench.map(i => `
    <tr>
      <td><strong>${{i.benchmark_id}}</strong></td>
      <td><code>${{i.immune_or_instability_metric}}</code></td>
      <td>${{i.group_A}} (n=${{i.n_group_A}})</td>
      <td>${{Number(i.median_group_A).toFixed(4)}}</td>
      <td>${{i.group_B}} (n=${{i.n_group_B}})</td>
      <td>${{Number(i.median_group_B).toFixed(4)}}</td>
      <td><strong>${{(Number(i.median_shift_A_minus_B)>=0?'+':'') + Number(i.median_shift_A_minus_B).toFixed(4)}}</strong></td>
      <td>${{i.p_value_two_sided}}</td>
    </tr>
  `).join('');
}}

function downloadCSV(filename, rows) {{
  if (!rows.length) return;
  const keys = Object.keys(rows[0]);
  const csv = [keys.join(','), ...rows.map(r => keys.map(k => JSON.stringify(r[k] ?? '')).join(','))].join('\\n');
  const blob = new Blob([csv], {{type: 'text/csv;charset=utf-8;'}});
  const a = document.createElement('a');
  a.href = URL.createObjectURL(blob);
  a.download = filename;
  a.click();
}}
function exportLongCSV() {{ downloadCSV('miRBioGenAtlas_480_gene_study_summary.csv', ATLAS.long); }}
function exportFilteredVariantsCSV() {{ downloadCSV('miRBioGenAtlas_filtered_3D_variants.csv', currentFilteredVars); }}
function exportDicerMirnaCSV() {{ downloadCSV('miRBioGenAtlas_DICER1_mirnaseq_5p3p.csv', ATLAS.dicer_pts); }}

window.addEventListener('DOMContentLoaded', () => {{
  const selVar = document.getElementById('var-gene');
  const selLol = document.getElementById('lol-gene');
  const selGtex = document.getElementById('gtex-gene');
  ATLAS.panel.forEach(p => {{
    const g = p.geneSymbol;
    selVar.insertAdjacentHTML('beforeend', `<option value="${{g}}">${{g}} (${{p.uniprotAccession}})</option>`);
    selLol.insertAdjacentHTML('beforeend', `<option value="${{g}}" ${{g==='DICER1'?'selected':''}}>${{g}} · ${{p.uniprotAccession}} (${{p.proteinLengthAa}} aa)</option>`);
    selGtex.insertAdjacentHTML('beforeend', `<option value="${{g}}" ${{g==='LIN28B'?'selected':''}}>${{g}} (&tau; = ${{p.gtex_tau_tissue_specificity}})</option>`);
  }});
  document.getElementById('sat-query').value = 'E1705K';
  renderHeatmap();
  renderVariants();
  renderAllTables();
  renderGtexBars();
  renderFormulas();
  renderDatasetsCatalogue();
  renderGlossaryCards();
  renderLollipopAndLoadSat();
}});
</script>
</body>
</html>
"""
    out_html = DOCS_DIR / "index.html"
    out_html.write_text(html, encoding="utf-8")
    print(f"6. Built standalone interactive web explorer: {out_html.relative_to(REPO_ROOT)} ({out_html.stat().st_size / 1024:.1f} KB)")


if __name__ == "__main__":
    build_web_explorer()
