#!/usr/bin/env python3
"""Compile CSV/CSV.GZ tables and 3D atomic coordinates into docs/index.html, docs/sat/*.js, and docs/struct3d/*.js."""

import csv
import gzip
import json
import math
import pathlib
from collections import Counter, defaultdict

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
    }

    data_json = json.dumps(bundle, separators=(",", ":"))

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>miRBioGenAtlas v2.0.2 — Multi-Consortium Pan-Cancer miRNA Biogenesis &amp; 3D Structure Explorer</title>
<script src="vendor/3Dmol-min.js"></script>
<style>
  :root {{
    --bg: #f8fafc;
    --card: #ffffff;
    --ink: #0f172a;
    --muted: #475569;
    --border: #e2e8f0;
    --brand: #0f5499;
    --brand-light: #e0f2fe;
    --accent: #dc2626;
    --teal: #0d9488;
  }}
  * {{ box-sizing: border-box; }}
  body {{
    margin: 0;
    font-family: -apple-system, BlinkMacSystemFont, "Inter", "Segoe UI", Roboto, sans-serif;
    background: var(--bg);
    color: var(--ink);
    line-height: 1.5;
  }}
  header {{
    background: linear-gradient(135deg, #0b2545 0%, #0f5499 65%, #136f63 100%);
    color: #fff;
    padding: 26px 24px 20px;
    border-bottom: 1px solid #0b2545;
  }}
  .container {{ max-width: 1420px; margin: 0 auto; padding: 0 18px; }}
  header h1 {{ margin: 0 0 6px; font-size: 1.6rem; letter-spacing: -0.02em; }}
  header p {{ margin: 0; opacity: 0.93; font-size: 0.92rem; max-width: 1150px; }}
  .badges {{ display: flex; flex-wrap: wrap; gap: 8px; margin-top: 12px; }}
  .badge {{
    background: rgba(255,255,255,0.14);
    border: 1px solid rgba(255,255,255,0.25);
    color: #fff;
    padding: 4px 10px;
    border-radius: 999px;
    font-size: 0.77rem;
    font-weight: 600;
    text-decoration: none;
  }}
  .kpi-row {{
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(175px, 1fr));
    gap: 10px;
    margin: 16px 0;
  }}
  .kpi {{
    background: var(--card);
    border: 1px solid var(--border);
    border-radius: 10px;
    padding: 12px 14px;
    box-shadow: 0 1px 2px rgba(15,23,42,0.04);
  }}
  .kpi .num {{ font-size: 1.38rem; font-weight: 800; color: var(--brand); }}
  .kpi .lbl {{ font-size: 0.75rem; color: var(--muted); margin-top: 2px; }}
  .tabs {{
    display: flex;
    flex-wrap: wrap;
    gap: 4px;
    border-bottom: 2px solid var(--border);
    margin-bottom: 16px;
  }}
  .tab-btn {{
    border: none;
    background: transparent;
    padding: 9px 13px;
    font-size: 0.85rem;
    font-weight: 600;
    color: var(--muted);
    cursor: pointer;
    border-bottom: 3px solid transparent;
    margin-bottom: -2px;
  }}
  .tab-btn.active {{
    color: var(--brand);
    border-bottom-color: var(--brand);
    background: #fff;
    border-radius: 8px 8px 0 0;
  }}
  .panel {{ display: none; }}
  .panel.active {{ display: block; }}
  .card {{
    background: var(--card);
    border: 1px solid var(--border);
    border-radius: 10px;
    padding: 16px 18px;
    margin-bottom: 16px;
    box-shadow: 0 1px 3px rgba(15,23,42,0.04);
  }}
  .controls {{
    display: flex;
    flex-wrap: wrap;
    gap: 10px;
    align-items: center;
    margin-bottom: 12px;
  }}
  select, input[type="text"], input[type="number"] {{
    padding: 6px 10px;
    border: 1px solid #cbd5e1;
    border-radius: 6px;
    font-size: 0.84rem;
    background: #fff;
  }}
  button.action-btn {{
    background: var(--brand);
    color: #fff;
    border: none;
    padding: 7px 13px;
    border-radius: 6px;
    font-size: 0.82rem;
    font-weight: 600;
    cursor: pointer;
  }}
  button.action-btn:hover {{ opacity: 0.92; }}
  button.ghost-btn {{
    background: #f1f5f9;
    color: #0f172a;
    border: 1px solid #cbd5e1;
    padding: 5px 10px;
    border-radius: 6px;
    font-size: 0.78rem;
    font-weight: 600;
    cursor: pointer;
  }}
  button.ghost-btn:hover {{ background: #e2e8f0; }}
  .preset-chip {{
    background: #eff6ff;
    color: #1e3a8a;
    border: 1px solid #bfdbfe;
    padding: 4px 10px;
    border-radius: 999px;
    font-size: 0.76rem;
    font-weight: 600;
    cursor: pointer;
  }}
  .preset-chip:hover {{ background: #dbeafe; }}
  .aa-btn {{
    border: 1px solid #cbd5e1;
    border-radius: 6px;
    padding: 3px 7px;
    font-size: 0.76rem;
    font-weight: 700;
    cursor: pointer;
    font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
  }}
  .aa-btn.active {{
    outline: 2px solid #0f172a;
    transform: translateY(-1px);
  }}
  .struct-grid {{
    display: grid;
    grid-template-columns: 1.15fr 1fr;
    gap: 16px;
    margin-top: 14px;
  }}
  @media (max-width: 1024px) {{
    .struct-grid {{ grid-template-columns: 1fr; }}
  }}
  .mol3d-stage {{
    position: relative;
    width: 100%;
    height: 495px;
    background: #0b1120;
    border-radius: 10px;
    border: 1px solid #1e293b;
    overflow: hidden;
  }}
  .mol3d-legend {{
    position: absolute;
    bottom: 10px;
    left: 10px;
    right: 10px;
    background: rgba(15, 23, 42, 0.84);
    color: #e2e8f0;
    padding: 7px 11px;
    border-radius: 7px;
    font-size: 0.74rem;
    display: flex;
    flex-wrap: wrap;
    gap: 12px;
    pointer-events: none;
    z-index: 5;
  }}
  .mol3d-dot {{
    display: inline-block;
    width: 9px;
    height: 9px;
    border-radius: 50%;
    margin-right: 4px;
    vertical-align: middle;
  }}
  .prop-box {{
    background: #f8fafc;
    border: 1px solid #cbd5e1;
    border-radius: 10px;
    padding: 13px 15px;
  }}
  .mech-note {{
    background: #eff6ff;
    border-left: 4px solid #0f5499;
    padding: 9px 12px;
    border-radius: 0 6px 6px 0;
    font-size: 0.8rem;
    color: #1e293b;
    margin-top: 10px;
  }}
  .table-wrap {{ overflow-x: auto; max-height: 500px; overflow-y: auto; border: 1px solid var(--border); border-radius: 8px; }}
  table {{ width: 100%; border-collapse: collapse; font-size: 0.81rem; }}
  th, td {{ padding: 7px 9px; border-bottom: 1px solid var(--border); text-align: left; white-space: nowrap; }}
  th {{ background: #f1f5f9; position: sticky; top: 0; z-index: 2; font-weight: 700; color: #1e293b; }}
  tr:hover td {{ background: #f8fafc; }}
  .pill {{
    display: inline-block;
    padding: 2px 8px;
    border-radius: 999px;
    font-size: 0.72rem;
    font-weight: 700;
  }}
  .pill-red {{ background: #fee2e2; color: #991b1b; }}
  .pill-orange {{ background: #ffedd5; color: #9a3412; }}
  .pill-dark {{ background: #e2e8f0; color: #0f172a; }}
  .pill-blue {{ background: #dbeafe; color: #1e40af; }}
  .pill-amber {{ background: #fef3c7; color: #92400e; }}
  .pill-green {{ background: #dcfce7; color: #166534; }}
  .heatmap-grid {{ overflow-x: auto; }}
  .hm-cell {{ text-align: center; font-size: 0.72rem; font-variant-numeric: tabular-nums; padding: 5px 4px; cursor: default; }}
  .svg-box {{ width: 100%; overflow-x: auto; border: 1px solid var(--border); border-radius: 8px; background: #fcfdff; padding: 8px; }}
  footer {{
    text-align: center;
    padding: 22px 16px 34px;
    color: var(--muted);
    font-size: 0.81rem;
  }}
</style>
</head>
<body>
<header>
  <div class="container">
    <h1>miRBioGenAtlas v2.0.2 — Multi-Consortium Pan-Cancer miRNA Biogenesis &amp; 3D Saturation Atlas</h1>
    <p>Harmonised across 21,140 clinical &amp; primary human tumours (10,195 TCGA PanCanAtlas + 10,945 MSK-IMPACT), 9,578 mature miRNA-seq profiles, 13.37M miRNA&times;mRNA couplings (9,353 dual tumours), 213,940 saturation-mutagenesis substitutions with interactive 3D wild-type vs mutated side-chain &amp; physicochemical delta modelling (16 Cryo-EM/X-ray PDB complexes, AlphaFold v6, AlphaMissense, ESM-2 3B, ESM-1v, ESM-IF1 3D), 16,569 ClinVar clinical records, gnomAD v4 constraint, GTEx 54-tissue baselines, CPTAC proteomes (N=1,002), and DepMap CRISPR screens (N=1,208).</p>
    <div class="badges">
      <a class="badge" href="https://doi.org/10.5281/zenodo.24455067" target="_blank" rel="noopener">Zenodo Dataset DOI: 10.5281/zenodo.24455067</a>
      <a class="badge" href="https://doi.org/10.5281/zenodo.24337286" target="_blank" rel="noopener">Concept DOI: 10.5281/zenodo.24337286</a>
      <a class="badge" href="https://github.com/siddugayam/miRBioGenAtlas" target="_blank" rel="noopener">GitHub: siddugayam/miRBioGenAtlas</a>
      <span class="badge">Licence: CC BY 4.0 (Data) / MIT (Code)</span>
    </div>
  </div>
</header>

<main class="container">
  <div class="kpi-row">
    <div class="kpi"><div class="num">21,140</div><div class="lbl">Human Tumours (10,195 TCGA + 10,945 MSK-IMPACT)</div></div>
    <div class="kpi"><div class="num">213,940</div><div class="lbl">3D Mutated Side-Chain &amp; Saturation Variants</div></div>
    <div class="kpi"><div class="num">16,569</div><div class="lbl">ClinVar Records (Missense ROC-AUC = 0.978)</div></div>
    <div class="kpi"><div class="num">9,578</div><div class="lbl">Mature miRNA-seq &amp; 13.37M miRNA&times;mRNA Pairs</div></div>
    <div class="kpi"><div class="num">54 Tissues</div><div class="lbl">GTEx v10 Normal Adult TPM + HPA Subcellular</div></div>
    <div class="kpi"><div class="num">2,210</div><div class="lbl">Orthogonal Validation (1,002 CPTAC + 1,208 DepMap)</div></div>
  </div>

  <div class="tabs">
    <button class="tab-btn active" id="btn-tab-heatmap" onclick="switchTab('tab-heatmap', this)">1. Pan-Cancer Heatmap (15&times;32)</button>
    <button class="tab-btn" id="btn-tab-variants" onclick="switchTab('tab-variants', this)">2. 3D Mutated Structure, Lollipop &amp; 213,940 Saturation</button>
    <button class="tab-btn" id="btn-tab-mirna" onclick="switchTab('tab-mirna', this)">3. miRNA-seq &amp; 13.37M Target Couplings</button>
    <button class="tab-btn" id="btn-tab-clinvar" onclick="switchTab('tab-clinvar', this)">4. ClinVar ROC-AUC, gnomAD &amp; MSK-IMPACT</button>
    <button class="tab-btn" id="btn-tab-gtex" onclick="switchTab('tab-gtex', this)">5. GTEx 54-Tissue, HPA &amp; Tractability</button>
    <button class="tab-btn" id="btn-tab-consortia" onclick="switchTab('tab-consortia', this)">6. CPTAC, DepMap, Survival &amp; PanImmune</button>
  </div>

  <!-- TAB 1: HEATMAP -->
  <section id="tab-heatmap" class="panel active">
    <div class="card">
      <div class="controls">
        <label><strong>Alteration Metric:</strong>
          <select id="hm-metric" onchange="renderHeatmap()">
            <option value="hc">Strict High-Confidence (Mutation or GISTIC2 &plusmn;2)</option>
            <option value="fcna">High-Level Focal CNA (&plusmn;2 only)</option>
            <option value="mut">Somatic Nonsynonymous Mutation (All)</option>
            <option value="nhm">Somatic Mutation (Non-Hypermutated Only)</option>
            <option value="any">Broad Any Alteration (Mutation or GISTIC2 &ne;0)</option>
          </select>
        </label>
        <button class="action-btn" onclick="exportLongCSV()">Download 480-Row Summary CSV</button>
      </div>
      <div id="hm-container" class="heatmap-grid"></div>
      <p style="font-size:0.78rem;color:var(--muted);margin:10px 0 0;">Values represent exact cohort proportions in [0.00, 1.00]. Hover any cell for exact sample counts (N, mutated, high-level amplified +2, deep deleted -2). Asterisk (*) marks studies with N &lt; 100.</p>
    </div>
  </section>

  <!-- TAB 2: INTERACTIVE 3D MUTATED STRUCTURE, SVG LOLLIPOP & 213,940 SATURATION CALCULATOR -->
  <section id="tab-variants" class="panel">
    <div class="card" id="struct3d-card">
      <div style="display:flex;justify-content:space-between;align-items:center;flex-wrap:wrap;gap:8px;margin-bottom:8px;">
        <h3 style="margin:0;">Interactive 3D Mutated Structure, Side-Chain Rotamer &amp; Physicochemical Delta Inspector</h3>
        <span id="sat-status" style="font-size:0.8rem;color:var(--muted);"></span>
      </div>
      <div class="controls" style="margin-bottom:8px;">
        <label><strong>Protein:</strong>
          <select id="lol-gene" onchange="renderLollipopAndLoadSat()"></select>
        </label>
        <label><strong>Codon or Mutation:</strong>
          <input type="text" id="sat-query" placeholder="e.g. E1705K, D1810Y, 1705..." oninput="inspectSaturationCodon()" style="width:165px;font-weight:600;">
        </label>
        <span style="font-size:0.78rem;color:var(--muted);font-weight:600;">Presets:</span>
        <button class="preset-chip" onclick="loadPresetMutation('DICER1','E1705K')">DICER1 E1705K (RNase IIIb Mg&sup2;&plus;)</button>
        <button class="preset-chip" onclick="loadPresetMutation('DICER1','D1810Y')">DICER1 D1810Y</button>
        <button class="preset-chip" onclick="loadPresetMutation('DROSHA','E1147K')">DROSHA E1147K</button>
        <button class="preset-chip" onclick="loadPresetMutation('TP53','R175H')">TP53 R175H (Core Fold)</button>
        <button class="preset-chip" onclick="loadPresetMutation('TP53','R248Q')">TP53 R248Q (DNA Contact)</button>
        <button class="preset-chip" onclick="loadPresetMutation('SMAD4','R361H')">SMAD4 R361H (Interface)</button>
        <button class="preset-chip" onclick="loadPresetMutation('AGO2','M364T')">AGO2 M364T (PAZ)</button>
        <button class="preset-chip" onclick="loadPresetMutation('DICER1','R944Q')">DICER1 R944Q (Benign Control)</button>
      </div>

      <div id="lollipop-svg-wrap" class="svg-box"></div>

      <div class="struct-grid">
        <div>
          <div class="controls" style="margin-bottom:8px;gap:6px;">
            <select id="mol3d-mode" onchange="refresh3DScene()" title="3D Display Mode">
              <option value="overlay">Superimposed: WT (Cyan) + Mutated Side-Chain (Crimson)</option>
              <option value="mutant">Mutated 3D Structure Only (In-Situ Side-Chain)</option>
              <option value="wt">Wild-Type 3D Structure Only</option>
              <option value="rcsb">Experimental RCSB Cryo-EM/X-ray Complex (RNA/DNA + Metal)</option>
            </select>
            <select id="mol3d-color" onchange="refresh3DScene()" title="Backbone Coloring">
              <option value="plddt">Ribbon: AlphaFold v6 pLDDT</option>
              <option value="am">Ribbon: Codon AlphaMissense</option>
              <option value="pocket">Ribbon: 3D Active Pocket (&le;6&Aring;)</option>
            </select>
            <button class="ghost-btn" onclick="zoomToMutatedResidue()">Zoom 8&Aring; Pocket</button>
            <button class="ghost-btn" onclick="resetFullProtein3D()">Full Protein</button>
            <button class="ghost-btn" id="btn-nbrs" onclick="toggleNeighbors3D()">8&Aring; Sidechains: ON</button>
            <button class="ghost-btn" id="btn-spin" onclick="toggleSpin3D()">Spin</button>
            <button class="ghost-btn" onclick="downloadMutatedPDB()">Export PDB</button>
          </div>
          <div class="mol3d-stage">
            <div id="mol3d-canvas" style="width:100%;height:100%;position:relative;"></div>
            <div class="mol3d-legend" id="mol3d-caption">
              <span><i class="mol3d-dot" style="background:#06b6d4;"></i>Wild-Type Side-Chain</span>
              <span><i class="mol3d-dot" style="background:#ef4444;"></i>Mutated 3D Rotamer</span>
              <span><i class="mol3d-dot" style="background:#94a3b8;"></i>8&Aring; Pocket Neighbors</span>
              <span><i class="mol3d-dot" style="background:#f59e0b;"></i>Polar/Pocket Contact (&le;3.8&Aring;)</span>
              <span><i class="mol3d-dot" style="background:#ff2e4c;"></i>Steric Clash (&lt;2.4&Aring;)</span>
            </div>
          </div>
          <p style="font-size:0.74rem;color:var(--muted);margin:6px 0 0;">Click any residue on the SVG track, 3D ribbon, or tables below to inspect it in 3D. Wild-type coordinates: AlphaFold v6 monomeric heavy-atom models (or live experimental RCSB PDB complexes); mutant side chains: clash-minimised AMBER-relaxed rotamer placed onto the local N&ndash;C&alpha;&ndash;C backbone frame.</p>
        </div>

        <div class="prop-box" id="prop-delta-panel">
          <div style="color:var(--muted);font-size:0.84rem;">Loading 3D structural &amp; physicochemical delta profile...</div>
        </div>
      </div>

      <div id="sat-codon-detail" style="margin-top:14px;"></div>
    </div>

    <div class="card">
      <h3 style="margin:0 0 10px;">Observed TCGA PanCanAtlas Somatic Coding Variants (6,292 Calls Across 10,195 Tumours)</h3>
      <div class="controls">
        <label><strong>Gene:</strong>
          <select id="var-gene" onchange="renderVariants()">
            <option value="ALL">All 15 Genes</option>
          </select>
        </label>
        <label><strong>V2 Functional / 3D Tier:</strong>
          <select id="var-tier" onchange="renderVariants()">
            <option value="ALL">All Tiers</option>
            <option value="3D_Catalytic_or_Substrate_Pocket_Driver">3D Pocket &le;6 &Aring; to RNA/Metal/Partner</option>
            <option value="Pathogenic_Fold_or_Interface_Missense">Pathogenic Fold / Interface (AlphaMissense)</option>
            <option value="Truncating_LoF">Truncating LoF</option>
            <option value="Ambiguous_Missense_VUS">Ambiguous VUS</option>
            <option value="Likely_Benign_Passenger_Missense">Likely Benign Passenger VUS</option>
          </select>
        </label>
        <label><strong>ESM-IF1 3D Mechanism:</strong>
          <select id="var-mech" onchange="renderVariants()">
            <option value="ALL">All Mechanisms</option>
            <option value="Fold_Preserving_Catalytic_or_Interface">Fold-Preserving Active / Interface</option>
            <option value="Thermodynamic_Fold_Destabilizing">Predicted Fold-Destabilising</option>
            <option value="Tolerated_Neutral">Tolerated Neutral</option>
          </select>
        </label>
        <label><strong>Filter:</strong>
          <input type="text" id="var-q" placeholder="e.g. E1705K, R248Q, Pathogenic..." oninput="renderVariants()">
        </label>
        <button class="action-btn" onclick="exportFilteredVariantsCSV()">Export Filtered Variants CSV</button>
      </div>
      <div id="var-count" style="font-size:0.83rem;color:var(--muted);margin-bottom:8px;"></div>
      <div class="table-wrap">
        <table id="var-table">
          <thead>
            <tr>
              <th>3D</th>
              <th>Gene</th>
              <th>Change</th>
              <th>TCGA Calls</th>
              <th>ClinVar Clinical</th>
              <th>V2 Functional / 3D Tier</th>
              <th>Biophysical Mechanism</th>
              <th>PDB</th>
              <th>RNA/DNA (&Aring;)</th>
              <th>Metal (&Aring;)</th>
              <th>AlphaMissense</th>
              <th>ESM-2 3B</th>
              <th>ESM-1v</th>
              <th>ESM-IF1 &Delta;logP</th>
              <th>Top TCGA Studies</th>
            </tr>
          </thead>
          <tbody></tbody>
        </table>
      </div>
    </div>
  </section>

  <!-- TAB 3: MIRNA-SEQ & TARGET COUPLINGS -->
  <section id="tab-mirna" class="panel">
    <div class="card">
      <h3 style="margin-top:0;">Mature miRNA-seq Enzymatic &amp; Downstream Target mRNA Benchmarks (B1&ndash;B12)</h3>
      <p style="font-size:0.85rem;color:var(--muted);">Evaluated across 9,578 matched mature miRNA-seq profiles (192 canonical <code>-5p / -3p</code> duplexes) and 9,353 dual-profiled TCGA tumours (13,369,752 mature miRNA &times; mRNA pairs).</p>
      <div class="table-wrap">
        <table id="mir-bench-table">
          <thead>
            <tr>
              <th>Benchmark</th>
              <th>Readout Metric</th>
              <th>Group A (n)</th>
              <th>Median A</th>
              <th>Group B (n)</th>
              <th>Median B</th>
              <th>Shift</th>
              <th>Two-Sided p</th>
            </tr>
          </thead>
          <tbody></tbody>
        </table>
      </div>
    </div>
    <div class="card">
      <div style="display:flex;justify-content:space-between;align-items:center;flex-wrap:wrap;gap:10px;margin-bottom:10px;">
        <h3 style="margin:0;">Individual DICER1-Altered TCGA Tumours: 5p:3p Strand Balance &amp; 5p let-7 Target mRNA Z-Score</h3>
        <button class="action-btn" onclick="exportDicerMirnaCSV()">Download DICER1 miRNA-seq CSV</button>
      </div>
      <div class="table-wrap" style="max-height:320px;">
        <table id="dicer-mir-table">
          <thead>
            <tr>
              <th>Sample ID</th>
              <th>Study</th>
              <th>DICER1 3D / AlphaMissense Class</th>
              <th>Paired Duplex log2((5p+1)/(3p+1))</th>
              <th>Linear 5p Fraction</th>
              <th>let-7 (5p) Target mRNA Z</th>
            </tr>
          </thead>
          <tbody></tbody>
        </table>
      </div>
    </div>
    <div class="card">
      <h3 style="margin-top:0;">Top Pan-Cancer Anti-Correlated miRNA &times; mRNA Axes (13,369,752 Pairs, N = 9,353 Tumours)</h3>
      <div class="table-wrap" style="max-height:300px;">
        <table id="top-axes-table">
          <thead>
            <tr>
              <th>Rank</th>
              <th>Mature miRNA</th>
              <th>Target mRNA Gene</th>
              <th>Dual Tumours (n)</th>
              <th>Pearson r (log2)</th>
              <th>t-Statistic</th>
            </tr>
          </thead>
          <tbody></tbody>
        </table>
      </div>
    </div>
  </section>

  <!-- TAB 4: CLINVAR ROC-AUC, GNOMAD & MSK-IMPACT -->
  <section id="tab-clinvar" class="panel">
    <div class="card">
      <h3 style="margin-top:0;">Non-Circular Human Clinical Benchmark: NCBI ClinVar Pathogenic vs Benign Missense ROC-AUC &amp; PR-AUC</h3>
      <p style="font-size:0.85rem;color:var(--muted);">Evaluated across unambiguous ClinVar coding missense assertions joined against the 213,940-variant saturation atlas.</p>
      <div class="table-wrap">
        <table id="clinvar-roc-table">
          <thead>
            <tr>
              <th>ClinVar Benchmark Subset</th>
              <th>Predictor</th>
              <th>Pathogenic (n)</th>
              <th>Benign (n)</th>
              <th>ROC-AUC (95% CI)</th>
              <th>PR-AUC (vs Baseline)</th>
              <th>Mann&ndash;Whitney p</th>
            </tr>
          </thead>
          <tbody></tbody>
        </table>
      </div>
    </div>
    <div class="card">
      <h3 style="margin-top:0;">Independent Clinical Sequencing Replication: MSK-IMPACT Cohort (N = 10,945 Patients) vs TCGA (N = 10,195)</h3>
      <div class="table-wrap">
        <table id="msk-rep-table">
          <thead>
            <tr>
              <th>Gene</th>
              <th>MSK Profiled N</th>
              <th>MSK Mut % (n)</th>
              <th>TCGA Mut %</th>
              <th>MSK Non-Hypermut %</th>
              <th>MSK Truncating % (n)</th>
              <th>MSK 3D Pocket / Hotspot % (n)</th>
              <th>MSK Primary vs Metastasis Mut %</th>
            </tr>
          </thead>
          <tbody></tbody>
        </table>
      </div>
    </div>
    <div class="card">
      <h3 style="margin-top:0;">MSK-IMPACT Clinical Cancer-Type Spectrum of Somatic DICER1 Mutations &amp; RNase IIIb Hotspots</h3>
      <div class="table-wrap" style="max-height:300px;">
        <table id="msk-dicer-table">
          <thead>
            <tr>
              <th>Clinical Cancer Type</th>
              <th>Samples (N)</th>
              <th>DICER1 Mutated (n, %)</th>
              <th>RNase IIIb 3D Hotspot (n, %)</th>
              <th>Truncating LoF (n, %)</th>
            </tr>
          </thead>
          <tbody></tbody>
        </table>
      </div>
    </div>
  </section>

  <!-- TAB 5: GTEX 54-TISSUE, HPA & TRACTABILITY -->
  <section id="tab-gtex" class="panel">
    <div class="card">
      <h3 style="margin-top:0;">Germline Constraint (gnomAD v4), Subcellular Localisation (HPA), GTEx Tissue Specificity &amp; Open Targets Tractability</h3>
      <div class="table-wrap">
        <table id="gene-ref-table">
          <thead>
            <tr>
              <th>Gene</th>
              <th>gnomAD pLI</th>
              <th>LOEUF</th>
              <th>Missense Z</th>
              <th>HPA Subcellular Compartment</th>
              <th>GTEx Pan-54 Median TPM</th>
              <th>GTEx Peak Tissue (TPM)</th>
              <th>Specificity (&tau;)</th>
              <th>Small-Mol / PROTAC Tractable</th>
              <th>ClinVar Pathogenic / Total</th>
            </tr>
          </thead>
          <tbody></tbody>
        </table>
      </div>
    </div>
    <div class="card">
      <div class="controls">
        <h3 style="margin:0;">GTEx v10 Healthy Human Baseline Expression Across All 54 Adult Tissues:</h3>
        <select id="gtex-gene" onchange="renderGtexBars()"></select>
      </div>
      <div id="gtex-bars" class="table-wrap" style="max-height:420px;padding:10px;"></div>
    </div>
  </section>

  <!-- TAB 6: CPTAC, DEPMAP, SURVIVAL & PANIMMUNE -->
  <section id="tab-consortia" class="panel">
    <div class="card">
      <h3 style="margin-top:0;">CPTAC Mass-Spectrometry Proteogenomics (N = 1,002) &amp; Broad DepMap CRISPR-Cas9 Essentiality (N = 1,208)</h3>
      <div class="table-wrap">
        <table id="consortia-table">
          <thead>
            <tr>
              <th>Gene</th>
              <th>Role</th>
              <th>PDB Complexes</th>
              <th>CPTAC CNA&rarr;mRNA &rho;</th>
              <th>CPTAC CNA&rarr;Protein &rho;</th>
              <th>Attenuation (&Delta;&rho;)</th>
              <th>DepMap Median Chronos</th>
              <th>% Dependent (&le;&minus;0.5)</th>
              <th>Essentiality Class</th>
            </tr>
          </thead>
          <tbody></tbody>
        </table>
      </div>
    </div>
    <div class="card">
      <h3 style="margin-top:0;">Cancer-Type-Stratified Log-Rank &amp; Multivariable Breslow Cox Survival Benchmarks (TCGA-CDR &amp; MSK-IMPACT)</h3>
      <p style="font-size:0.84rem;color:var(--muted);">Multivariable TCGA Cox models adjust simultaneously for age at diagnosis, biological sex, AJCC pathologic stage (III/IV vs I/II), hypermutator status, and cancer-type baseline hazard strata; MSK-IMPACT models adjust for metastatic specimen status, TMB hypermutator status, panel version, and cancer-type strata.</p>
      <div class="table-wrap">
        <table id="surv-table">
          <thead>
            <tr>
              <th>Benchmark ID</th>
              <th>Cohort</th>
              <th>Endpoint</th>
              <th>Altered Group (Events / n)</th>
              <th>Control Group (Events / n)</th>
              <th>Stratified Log-Rank HR (95% CI)</th>
              <th>Log-Rank p</th>
              <th>Multivariable Cox HR (95% CI)</th>
              <th>Multivariable p</th>
            </tr>
          </thead>
          <tbody></tbody>
        </table>
      </div>
    </div>
    <div class="card">
      <h3 style="margin-top:0;">Thorsson Pan-Cancer Immune Landscape &amp; Homologous Recombination Deficiency (HRD) Benchmarks</h3>
      <div class="table-wrap">
        <table id="imm-table">
          <thead>
            <tr>
              <th>Benchmark ID</th>
              <th>Metric</th>
              <th>Group A (n)</th>
              <th>Median A</th>
              <th>Group B (n)</th>
              <th>Median B</th>
              <th>Median Shift</th>
              <th>Two-Sided p</th>
            </tr>
          </thead>
          <tbody></tbody>
        </table>
      </div>
    </div>
  </section>
</main>

<footer>
  miRBioGenAtlas v2.0.2 &bull; Manipal Academy of Higher Education (MAHE) &bull; Data: CC BY 4.0 &bull; Code: MIT &bull; Deterministically compiled from 46 SHA-256-verified snapshots.
</footer>

<script>
window.SAT_CACHE = window.SAT_CACHE || {{}};
window.STRUCT3D_CACHE = window.STRUCT3D_CACHE || {{}};
window.RCSB_CACHE = window.RCSB_CACHE || {{}};
const ATLAS = {data_json};
const LOW_N = new Set(["ACC","CHOL","DLBC","KICH","MESO","UCS","UVM"]);

const AA_PROPS = {{
  A: {{three:'ALA', name:'Alanine', cls:'Aliphatic Nonpolar', q:0.0, vol:88.6, mw:89.09, kd:1.8, hd:0, ha:0, arom:0, metal:0}},
  R: {{three:'ARG', name:'Arginine', cls:'Basic Guanidinium (+)', q:1.0, vol:173.4, mw:174.20, kd:-4.5, hd:5, ha:0, arom:0, metal:0}},
  N: {{three:'ASN', name:'Asparagine', cls:'Polar Carboxamide', q:0.0, vol:114.1, mw:132.12, kd:-3.5, hd:2, ha:1, arom:0, metal:1}},
  D: {{three:'ASP', name:'Aspartate', cls:'Acidic Carboxylate (-)', q:-1.0, vol:111.1, mw:133.10, kd:-3.5, hd:0, ha:2, arom:0, metal:1}},
  C: {{three:'CYS', name:'Cysteine', cls:'Thiol / Redox Sulfhydryl', q:0.0, vol:108.5, mw:121.16, kd:2.5, hd:1, ha:1, arom:0, metal:1}},
  Q: {{three:'GLN', name:'Glutamine', cls:'Polar Carboxamide', q:0.0, vol:143.8, mw:146.15, kd:-3.5, hd:2, ha:1, arom:0, metal:1}},
  E: {{three:'GLU', name:'Glutamate', cls:'Acidic Carboxylate (-)', q:-1.0, vol:138.4, mw:147.13, kd:-3.5, hd:0, ha:2, arom:0, metal:1}},
  G: {{three:'GLY', name:'Glycine', cls:'Achiral Flexible H-only', q:0.0, vol:60.1, mw:75.07, kd:-0.4, hd:0, ha:0, arom:0, metal:0}},
  H: {{three:'HIS', name:'Histidine', cls:'Basic Imidazole / Zn-Ligand', q:0.1, vol:153.2, mw:155.16, kd:-3.2, hd:1, ha:1, arom:1, metal:1}},
  I: {{three:'ILE', name:'Isoleucine', cls:'Branched Hydrophobic', q:0.0, vol:166.7, mw:131.17, kd:4.5, hd:0, ha:0, arom:0, metal:0}},
  L: {{three:'LEU', name:'Leucine', cls:'Branched Hydrophobic', q:0.0, vol:166.7, mw:131.17, kd:3.8, hd:0, ha:0, arom:0, metal:0}},
  K: {{three:'LYS', name:'Lysine', cls:'Basic Primary Amine (+)', q:1.0, vol:168.6, mw:146.19, kd:-3.9, hd:3, ha:0, arom:0, metal:0}},
  M: {{three:'MET', name:'Methionine', cls:'Thioether Nonpolar', q:0.0, vol:162.9, mw:149.21, kd:1.9, hd:0, ha:1, arom:0, metal:0}},
  F: {{three:'PHE', name:'Phenylalanine', cls:'Bulky Aromatic Nonpolar', q:0.0, vol:189.9, mw:165.19, kd:2.8, hd:0, ha:0, arom:1, metal:0}},
  P: {{three:'PRO', name:'Proline', cls:'Cyclic Imino Helix-Breaker', q:0.0, vol:112.7, mw:115.13, kd:-1.6, hd:0, ha:0, arom:0, metal:0}},
  S: {{three:'SER', name:'Serine', cls:'Small Polar Hydroxyl', q:0.0, vol:89.0, mw:105.09, kd:-0.8, hd:1, ha:1, arom:0, metal:0}},
  T: {{three:'THR', name:'Threonine', cls:'Beta-Branched Polar Hydroxyl', q:0.0, vol:116.1, mw:119.12, kd:-0.7, hd:1, ha:1, arom:0, metal:0}},
  W: {{three:'TRP', name:'Tryptophan', cls:'Bulky Indole Aromatic', q:0.0, vol:227.8, mw:204.23, kd:-0.9, hd:1, ha:0, arom:1, metal:0}},
  Y: {{three:'TYR', name:'Tyrosine', cls:'Phenolic Polar Aromatic', q:0.0, vol:193.6, mw:181.19, kd:-1.3, hd:1, ha:1, arom:1, metal:0}},
  V: {{three:'VAL', name:'Valine', cls:'Beta-Branched Hydrophobic', q:0.0, vol:140.0, mw:117.15, kd:4.2, hd:0, ha:0, arom:0, metal:0}}
}};

// Grantham (1974) Science 185:862 composition-polarity-volume parameters (c, p, v)
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
      refresh3DScene();
    }}, 40);
  }}
}}

function loadPresetMutation(gene, mutStr) {{
  const selLol = document.getElementById('lol-gene');
  selLol.value = gene;
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

function heatColor(val, maxScale) {{
  const t = Math.min(1, Math.max(0, val / maxScale));
  const r = Math.round(255 - t * (255 - 15));
  const g = Math.round(255 - t * (255 - 84));
  const b = Math.round(255 - t * (255 - 153));
  const txt = t > 0.55 ? '#ffffff' : '#0f172a';
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
      const tip = r ? `${{g}} in ${{s}} (N=${{r.N}}): value=${{v.toFixed(3)}} | mut_n=${{r.mut_n}}, amp_n=${{r.amp_n}}, deep_del_n=${{r.dd_n}}` : '';
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
    document.getElementById('sat-status').textContent = `Loaded ${{(nCod * 19).toLocaleString()}} substitutions (${{nCod}} codons) + 3D heavy-atom model for ${{gene}}.`;
    if (cb) cb();
    return;
  }}
  document.getElementById('sat-status').textContent = `Loading 19x saturation matrix & 3D structure for ${{gene}}...`;
  const sc = document.createElement('script');
  sc.src = `sat/${{gene}}.js`;
  sc.onload = () => {{
    const nCod = window.SAT_CACHE[gene] ? Object.keys(window.SAT_CACHE[gene]).length : 0;
    document.getElementById('sat-status').textContent = `Loaded ${{(nCod * 19).toLocaleString()}} substitutions (${{nCod}} codons) for ${{gene}}.`;
    if (cb) cb();
  }};
  sc.onerror = () => {{
    document.getElementById('sat-status').textContent = `Codon track active for ${{gene}}.`;
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
  sc.onerror = () => {{
    if (cb) cb(null);
  }};
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
  const W = 1260, H = 195, padL = 52, padR = 24;
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
      pockets += `<rect x="${{x-1.2}}" y="122" width="2.6" height="16" fill="#dc2626" opacity="0.68" style="cursor:pointer" onclick="selectCodonFromTrack(${{pos}})"><title>${{gene}} ${{wt}}${{pos}} (3D pocket &le;6&Aring; in ${{pdb}})</title></rect>`;
    }} else if (pdb) {{
      pockets += `<rect x="${{x-0.8}}" y="125" width="1.8" height="10" fill="#0ea5e9" opacity="0.35"></rect>`;
    }}
    if (am !== null) {{
      const yAm = 176 - am * 28;
      amPath += (amPath ? ' L' : 'M') + x.toFixed(1) + ',' + yAm.toFixed(1);
    }}
    if (tcgaN > 0 || cvN > 0) {{
      const h = 15 + (Math.log2(tcgaN + cvN + 1) / Math.log2(maxCnt + 1)) * 84;
      const yTop = 118 - h;
      const col = pkt === 1 ? '#dc2626' : (cvN > 0 ? '#7c3aed' : '#0f5499');
      const r = Math.min(6.5, 2.8 + Math.log2(tcgaN + cvN + 1));
      stems += `<g style="cursor:pointer" onclick="selectCodonFromTrack(${{pos}})">
        <line x1="${{x.toFixed(1)}}" y1="118" x2="${{x.toFixed(1)}}" y2="${{yTop.toFixed(1)}}" stroke="${{col}}" stroke-width="1.3" opacity="0.75"/>
        <circle cx="${{x.toFixed(1)}}" cy="${{yTop.toFixed(1)}}" r="${{r.toFixed(1)}}" fill="${{col}}" stroke="#fff" stroke-width="0.7">
          <title>${{gene}} ${{wt}}${{pos}} | TCGA somatic calls: ${{tcgaN}} | ClinVar Pathogenic: ${{cvN}} | Mean AlphaMissense: ${{am ?? 'NA'}} | Click to view in 3D</title>
        </circle>
      </g>`;
    }}
  }});

  let ticks = '';
  const step = L > 1200 ? 200 : (L > 500 ? 100 : 50);
  for (let p = 1; p <= L; p += step) {{
    const tx = xOf(p);
    ticks += `<text x="${{tx.toFixed(1)}}" y="190" font-size="10" fill="#475569" text-anchor="middle">${{p}}</text>`;
  }}
  ticks += `<text x="${{xOf(L).toFixed(1)}}" y="190" font-size="10" fill="#475569" text-anchor="end">${{L}} aa</text>`;

  const svg = `<svg viewBox="0 0 ${{W}} ${{H}}" width="100%" height="195">
    <text x="6" y="24" font-size="10.5" font-weight="700" fill="#0f172a">Somatic / ClinVar</text>
    <text x="6" y="134" font-size="10.5" font-weight="700" fill="#0f172a">3D PDB / Pocket</text>
    <text x="6" y="166" font-size="10.5" font-weight="700" fill="#0f172a">AlphaMissense</text>
    <rect x="${{padL}}" y="122" width="${{plotW}}" height="16" rx="4" fill="#e2e8f0"/>
    ${{pockets}}
    ${{stems}}
    <path d="${{amPath}}" fill="none" stroke="#ea580c" stroke-width="1.15" opacity="0.85"/>
    <line x1="${{padL}}" y1="${{176 - 0.564*28}}" x2="${{W-padR}}" y2="${{176 - 0.564*28}}" stroke="#dc2626" stroke-dasharray="3,3" stroke-width="0.8"/>
    ${{ticks}}
    <g transform="translate(${{W - 430}}, 14)" font-size="10.5" fill="#334155">
      <circle cx="0" cy="-3" r="4" fill="#dc2626"/><text x="8" y="0">3D Pocket (&le;6&Aring;)</text>
      <circle cx="115" cy="-3" r="4" fill="#7c3aed"/><text x="123" y="0">ClinVar Pathogenic</text>
      <circle cx="240" cy="-3" r="4" fill="#0f5499"/><text x="248" y="0">TCGA Somatic</text>
      <rect x="335" y="-7" width="8" height="8" fill="#0ea5e9" opacity="0.5"/><text x="347" y="0">PDB Resolved</text>
    </g>
  </svg>`;
  document.getElementById('lollipop-svg-wrap').innerHTML = svg;

  ensureSatLoaded(gene, () => {{
    const rawQ = document.getElementById('sat-query').value.trim().toUpperCase();
    const gCache = window.SAT_CACHE[gene] || {{}};
    const m = rawQ.match(/^([A-Z])?(\\d+)([A-Z])?$/);
    if (!m || !gCache[m[2]] || (m[1] && gCache[m[2]].m[0] !== m[1])) {{
      // Pick top mutated / pocket codon in this protein automatically when switching genes
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

  // Collect neighbor heavy atoms within 10 A (excluding current residue sidechain)
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

  // Wild-type closest heavy-atom contact
  let wtMinDist = 99.9;
  resRows.forEach(r => {{
    if (['N','CA','C','O'].includes(r[1])) return;
    envAtoms.forEach(ea => {{
      const d = Math.hypot(r[2]-ea[2], r[3]-ea[3], r[4]-ea[4]);
      if (d < wtMinDist) wtMinDist = d;
    }});
  }});
  if (wtMinDist > 50) wtMinDist = null;

  // Orthonormal frame at CA
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
      // Slight preference for native template angle (deg === 0)
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
  viewer3d = $3Dmol.createViewer(el, {{ backgroundColor: '#0b1120', antialias: true }});
  return viewer3d;
}}

function refresh3DScene(keepCamera) {{
  const v = ensureViewerInitialized();
  if (!v) return;
  const {{ gene, pos, wt, mut }} = active3DState;
  const mode = document.getElementById('mol3d-mode').value;
  const colorMode = document.getElementById('mol3d-color').value;

  if (mode === 'rcsb') {{
    renderExperimentalRCSB(gene, pos, wt, mut);
    return;
  }}

  ensureStruct3DLoaded(gene, st => {{
    if (!st) return;
    const geom = computeMutant3DGeometry(st, pos, wt, mut);
    const pdbs = buildPdbModels(gene, st, pos, wt, mut, geom, colorMode);
    active3DState.lastPdbText = pdbs.wtPdb + (mode !== 'wt' ? ('\\n' + pdbs.mutPdb) : '');

    v.removeAllModels();
    v.removeAllShapes();
    v.removeAllLabels();

    const m0 = v.addModel(pdbs.wtPdb, 'pdb');
    const colorFunc = atom => {{
      const b = atom.b;
      if (colorMode === 'plddt') {{
        if (b >= 90) return '#1d4ed8';
        if (b >= 70) return '#38bdf8';
        if (b >= 50) return '#facc15';
        return '#f97316';
      }} else if (colorMode === 'am') {{
        if (b >= 56.4) return '#dc2626';
        if (b >= 34.0) return '#f59e0b';
        return '#3b82f6';
      }} else {{
        return b >= 80 ? '#ef4444' : '#475569';
      }}
    }};
    v.setStyle({{model: 0}}, {{cartoon: {{colorfunc: colorFunc, opacity: 0.86}}}});

    // 8 A neighbor sidechains
    if (showNeighbors && geom && geom.bb.CA) {{
      const [cx, cy, cz] = geom.bb.CA;
      const nbrRes = [];
      st._caList.forEach(c => {{
        if (c[0] !== pos && Math.hypot(c[1]-cx, c[2]-cy, c[3]-cz) <= 8.5) nbrRes.push(c[0]);
      }});
      if (nbrRes.length) {{
        v.addStyle({{model: 0, resi: nbrRes}}, {{stick: {{radius: 0.11, colorscheme: 'Jmol', opacity: 0.55}}}});
      }}
    }}

    // Wild-type side chain highlight
    if (mode === 'overlay' || mode === 'wt') {{
      v.addStyle(
        {{model: 0, resi: pos}},
        {{
          stick: {{radius: 0.22, color: '#06b6d4', opacity: mode === 'overlay' ? 0.72 : 1.0}},
          sphere: {{scale: 0.26, color: '#06b6d4', opacity: mode === 'overlay' ? 0.72 : 1.0}}
        }}
      );
    }}

    // Mutated 3D rotamer model + contacts/clashes
    if ((mode === 'overlay' || mode === 'mutant') && pdbs.mutPdb && wt !== mut) {{
      v.addModel(pdbs.mutPdb, 'pdb');
      v.setStyle(
        {{model: 1}},
        {{
          stick: {{radius: 0.25, color: '#ef4444'}},
          sphere: {{scale: 0.30, color: '#ef4444'}}
        }}
      );
      if (geom) {{
        geom.contacts.forEach(ct => {{
          v.addCylinder({{
            start: {{x:ct.from[0], y:ct.from[1], z:ct.from[2]}},
            end: {{x:ct.to[0], y:ct.to[1], z:ct.to[2]}},
            radius: 0.05, dashed: true, dashLength: 0.22, gapLength: 0.16, color: '#f59e0b'
          }});
        }});
        geom.clashes.forEach(cl => {{
          v.addCylinder({{
            start: {{x:cl.from[0], y:cl.from[1], z:cl.from[2]}},
            end: {{x:cl.to[0], y:cl.to[1], z:cl.to[2]}},
            radius: 0.08, dashed: true, dashLength: 0.18, gapLength: 0.12, color: '#ff2e4c'
          }});
        }});
      }}
    }}

    if (geom && geom.bb.CA) {{
      const [lx, ly, lz] = geom.bb.CA;
      const lblTxt = wt === mut ? `${{gene}} ${{wt}}${{pos}} (WT)` : `${{gene}} ${{wt}}${{pos}}${{mut}}`;
      v.addLabel(lblTxt, {{
        position: {{x: lx, y: ly + 1.8, z: lz}},
        backgroundColor: wt === mut ? '#0e7490' : '#b91c1c',
        fontColor: '#ffffff',
        fontSize: 12,
        backgroundOpacity: 0.88,
        borderThickness: 0
      }});
    }}

    // Click any residue on 3D ribbon to inspect
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

  const drawCif = cifTxt => {{
    v.removeAllModels();
    v.removeAllShapes();
    v.removeAllLabels();
    v.addModel(cifTxt, 'cif');
    v.setStyle({{}}, {{cartoon: {{color: '#64748b', opacity: 0.82}}}});
    v.addStyle({{resn: ['A','U','G','C','DA','DT','DG','DC']}}, {{cartoon: {{color: '#f59e0b'}}, stick: {{radius: 0.14, color: '#fbbf24'}}}});
    v.addStyle({{elem: ['MG','ZN','CA','MN']}}, {{sphere: {{radius: 1.35, color: '#22c55e'}}}});
    v.addStyle({{resi: pos}}, {{stick: {{radius: 0.28, color: '#ef4444'}}, sphere: {{scale: 0.32, color: '#ef4444'}}}});
    v.addLabel(`Experimental PDB ${{pdbRaw}} • ${{gene}} ${{wt}}${{pos}}`, {{
      resi: pos, backgroundColor: '#991b1b', fontColor: '#fff', fontSize: 12, backgroundOpacity: 0.9
    }});
    v.zoomTo({{resi: pos}}, 450);
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
  viewer3d.zoomTo({{model: 0, resi: [Math.max(1, pos-4), pos, pos+4]}}, 400);
  viewer3d.zoom(0.82, 200);
}}

function resetFullProtein3D() {{
  if (!viewer3d) return;
  viewer3d.zoomTo({{model: 0}}, 450);
}}

function toggleNeighbors3D() {{
  showNeighbors = !showNeighbors;
  document.getElementById('btn-nbrs').textContent = `8Å Sidechains: ${{showNeighbors ? 'ON' : 'OFF'}}`;
  refresh3DScene(true);
}}

function toggleSpin3D() {{
  if (!viewer3d) return;
  spinActive = !spinActive;
  viewer3d.spin(spinActive ? 'y' : false);
  document.getElementById('btn-spin').textContent = spinActive ? 'Stop Spin' : 'Spin';
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
  const gLabel = gDist === 0 ? 'Identical (WT)' : (gDist <= 50 ? 'Conservative' : (gDist <= 100 ? 'Moderately Conservative' : (gDist <= 150 ? 'Moderately Radical' : 'Radical Shift')));
  const gPill = gDist <= 50 ? 'pill-blue' : (gDist <= 100 ? 'pill-amber' : 'pill-red');

  let chargeNote = 'Neutral / Preserved';
  let chargePill = 'pill-blue';
  if (wp.q * mp.q < -0.2) {{ chargeNote = 'Charge Inversion!'; chargePill = 'pill-red'; }}
  else if (Math.abs(dq) >= 0.8) {{ chargeNote = dq > 0 ? 'Gain of (+) Charge' : 'Gain of (-) / Loss of (+)'; chargePill = 'pill-orange'; }}

  const alerts = [];
  if (mut === 'P' && wt !== 'P') alerts.push('<span class="pill pill-red">Proline Helix/Sheet Breaker</span>');
  if (wt === 'G' && mut !== 'G') alerts.push('<span class="pill pill-orange">Loss of Glycine Torsional Freedom</span>');
  if ((wt === 'C') !== (mut === 'C')) alerts.push('<span class="pill pill-amber">Thiol / Disulfide Change</span>');
  if (wp.metal && !mp.metal && dmet !== null && dmet <= 6.0) alerts.push('<span class="pill pill-red">Catalytic Metal-Chelator Loss</span>');
  if (!wp.arom && mp.arom) alerts.push('<span class="pill pill-amber">Bulky Aromatic Ring Introduced</span>');

  const constrained = (am !== null && am >= 0.564) || (e2 !== null && e2 <= -5.0);
  const regimeBadge = wt === mut
    ? '<span class="pill pill-green">Wild-Type Reference</span>'
    : (constrained
        ? ((if1 !== null && if1 > -1.5) ? '<span class="pill pill-red">Fold-Preserving Catalytic / Interface</span>' : '<span class="pill pill-orange">Thermodynamic Fold-Destabilising</span>')
        : '<span class="pill pill-blue">Tolerated / Neutral</span>');

  const burialStr = geom
    ? (geom.caNeighbors8A >= 13 ? `Buried Core (${{geom.caNeighbors8A}} C&alpha; &le;8&Aring;)` : (geom.caNeighbors8A >= 8 ? `Partially Buried / Cleft (${{geom.caNeighbors8A}} C&alpha;)` : `Solvent-Exposed (${{geom.caNeighbors8A}} C&alpha;)`))
    : '&mdash;';
  const clashCount = geom ? geom.clashes.length : 0;
  const clashBadge = clashCount > 0
    ? `<span class="pill pill-red">${{clashCount}} Steric Clash${{clashCount>1?'es':''}} (&lt;2.4&Aring;; min ${{geom.mutMinDist ? geom.mutMinDist.toFixed(2) : '?'}}&Aring;)</span>`
    : `<span class="pill pill-green">0 Clashes (min ${{geom && geom.mutMinDist ? geom.mutMinDist.toFixed(2) : '&gt;2.5'}}&Aring;)</span>`;

  // 19-AA quick switcher buttons at this codon
  const aaBtns = entry.s.map(r => {{
    const a = r[0], rAm = r[1], rE2 = r[2], rIf1 = r[4];
    const cns = (rAm !== null && rAm >= 0.564) || (rE2 !== null && rE2 <= -5.0);
    const bg = cns ? ((rIf1 !== null && rIf1 > -1.5) ? '#fee2e2' : '#ffedd5') : '#eff6ff';
    const fg = cns ? ((rIf1 !== null && rIf1 > -1.5) ? '#991b1b' : '#9a3412') : '#1e40af';
    return `<button class="aa-btn ${{a===mut?'active':''}}" style="background:${{bg}};color:${{fg}};" onclick="selectMutAA(${{pos}},'${{wt}}','${{a}}')" title="Mutate ${{wt}}${{pos}} to ${{a}} (AM=${{rAm??'NA'}}, IF1=${{rIf1??'NA'}})">${{a}}</button>`;
  }}).join(' ');

  // Mechanistic plain-English synthesis
  let mechTxt = '';
  if (wt === mut) {{
    mechTxt = `Wild-type <strong>${{wp.name}} (${{wt}}${{pos}})</strong> resides in a ${{burialStr.toLowerCase().replace(/<[^>]+>/g,'')}} environment (AFDB pLDDT ${{plddt ?? 'NA'}}). Select any of the 19 amino acid buttons above to model the 3D side-chain substitution and physicochemical shift.`;
  }} else {{
    const pocketParts = [];
    if (dmet !== null && dmet <= 6.0) pocketParts.push(`${{dmet}} &Aring; from catalytic divalent metal`);
    if (dna !== null && dna <= 6.0) pocketParts.push(`${{dna}} &Aring; from bound RNA/DNA`);
    if (dprt !== null && dprt <= 6.0) pocketParts.push(`${{dprt}} &Aring; from partner subunit`);
    const locPhrase = pocketParts.length ? `within the experimental 3D active site (${{pocketParts.join(', ')}} in PDB ${{pdb}})` : `in a ${{burialStr.toLowerCase().replace(/<[^>]+>/g,'')}} region`;
    const stabPhrase = (if1 !== null && if1 <= -1.5)
      ? `destabilises local 3D backbone fold thermodynamics (ESM-IF1 &Delta;logP = ${{if1.toFixed(2)}}${{clashCount>0 ? ', ' + clashCount + ' steric clash(es) <2.4 Å' : ''}})`
      : (constrained
          ? `preserves global 3D backbone fold stability (ESM-IF1 &Delta;logP = ${{if1 !== null ? if1.toFixed(2) : 'NA'}}, ${{clashCount}} steric clashes) while perturbing active-site/interface chemistry`
          : `is biophysically tolerated on both evolutionary conservation (AlphaMissense ${{am !== null ? am.toFixed(3) : 'NA'}}) and 3D fold stability`);
    mechTxt = `<strong>${{gene}} ${{wt}}${{pos}}${{mut}} (${{wp.three}}&rarr;${{mp.three}})</strong> replaces ${{wp.cls.toLowerCase()}} <strong>${{wp.name}}</strong> with ${{mp.cls.toLowerCase()}} <strong>${{mp.name}}</strong> ${{locPhrase}} (&Delta;q = ${{dq>=0?'+':''}}${{dq.toFixed(1)}} e, &Delta;Volume = ${{dvol>=0?'+':''}}${{dvol.toFixed(1)}} &Aring;&sup3;, Grantham D = ${{gDist}}), which ${{stabPhrase}}.`;
  }}

  box.innerHTML = `
    <div style="display:flex;justify-content:space-between;align-items:center;flex-wrap:wrap;gap:6px;margin-bottom:8px;">
      <div>
        <span style="font-size:1.05rem;font-weight:800;color:#0f172a;">${{gene}} <code>p.${{wp.three}}${{pos}}${{mp.three}}</code> (${{wt}}${{pos}}${{mut}})</span>
        ${{tn > 0 ? `<span class="pill pill-dark" style="margin-left:4px;">TCGA n=${{tn}}</span>` : ''}}
        ${{cv ? `<span class="pill pill-amber" style="margin-left:4px;">ClinVar: ${{cv}}</span>` : ''}}
      </div>
      <div>${{regimeBadge}}</div>
    </div>

    <div style="margin-bottom:9px;">
      <div style="font-size:0.74rem;color:var(--muted);margin-bottom:4px;font-weight:600;">Click any amino acid to mutate codon ${{wt}}${{pos}} in 3D (<span style="color:#991b1b;">red</span>=pocket/interface, <span style="color:#9a3412;">orange</span>=fold-destabilising, <span style="color:#1e40af;">blue</span>=neutral):</div>
      <div style="display:flex;flex-wrap:wrap;gap:4px;">${{aaBtns}}</div>
    </div>

    <div class="table-wrap" style="max-height:none;">
      <table style="font-size:0.78rem;">
        <thead>
          <tr><th>Protein Property</th><th>Wild-Type (${{wt}})</th><th>Mutant (${{mut}})</th><th>&Delta; Change &amp; Structural Effect</th></tr>
        </thead>
        <tbody>
          <tr>
            <td><strong>Side-Chain Class</strong></td>
            <td>${{wp.cls}}</td>
            <td><strong>${{mp.cls}}</strong></td>
            <td>${{wp.cls === mp.cls ? '<span class="pill pill-blue">Same Class</span>' : '<span class="pill pill-orange">Class Switch</span>'}}</td>
          </tr>
          <tr>
            <td><strong>Net Charge (pH 7.4)</strong></td>
            <td>${{wp.q >= 0 ? '+' : ''}}${{wp.q.toFixed(1)}} e</td>
            <td><strong>${{mp.q >= 0 ? '+' : ''}}${{mp.q.toFixed(1)}} e</strong></td>
            <td><strong>${{dq >= 0 ? '+' : ''}}${{dq.toFixed(1)}} e</strong> &bull; <span class="pill ${{chargePill}}">${{chargeNote}}</span></td>
          </tr>
          <tr>
            <td><strong>Steric Volume</strong></td>
            <td>${{wp.vol.toFixed(1)}} &Aring;&sup3;</td>
            <td><strong>${{mp.vol.toFixed(1)}} &Aring;&sup3;</strong></td>
            <td><strong>${{dvol >= 0 ? '+' : ''}}${{dvol.toFixed(1)}} &Aring;&sup3;</strong> (${{dvol >= 0 ? '+' : ''}}${{((dvol/wp.vol)*100).toFixed(0)}}% ${{Math.abs(dvol) < 15 ? 'isosteric' : (dvol > 0 ? 'bulkier' : 'cavity shift')}})</td>
          </tr>
          <tr>
            <td><strong>Hydropathy (KD)</strong></td>
            <td>${{wp.kd.toFixed(1)}}</td>
            <td><strong>${{mp.kd.toFixed(1)}}</strong></td>
            <td><strong>${{dkd >= 0 ? '+' : ''}}${{dkd.toFixed(1)}}</strong> (${{dkd > 1.5 ? 'More Hydrophobic' : (dkd < -1.5 ? 'More Polar/Charged' : 'Similar Polarity')}})</td>
          </tr>
          <tr>
            <td><strong>H-Bonds / Metal</strong></td>
            <td>${{wp.hd}}D / ${{wp.ha}}A${{wp.metal ? ' (Chelator)' : ''}}</td>
            <td><strong>${{mp.hd}}D / ${{mp.ha}}A${{mp.metal ? ' (Chelator)' : ''}}</strong></td>
            <td>&Delta;Donors=${{mp.hd-wp.hd>=0?'+':''}}${{mp.hd-wp.hd}}, &Delta;Acceptors=${{mp.ha-wp.ha>=0?'+':''}}${{mp.ha-wp.ha}}</td>
          </tr>
          <tr>
            <td><strong>Grantham &amp; Alerts</strong></td>
            <td colspan="2">Grantham D = <strong>${{gDist}}</strong> <span class="pill ${{gPill}}">${{gLabel}}</span></td>
            <td>${{alerts.length ? alerts.join(' ') : '<span style="color:var(--muted);">No special backbone alert</span>'}}</td>
          </tr>
          <tr>
            <td><strong>3D Packing &amp; Clash</strong></td>
            <td>${{burialStr}}</td>
            <td>WT min: ${{geom && geom.wtMinDist ? geom.wtMinDist.toFixed(2)+'&Aring;' : '&mdash;'}}</td>
            <td>${{clashBadge}}</td>
          </tr>
          <tr>
            <td><strong>PLM &amp; 3D Scores</strong></td>
            <td colspan="2">AM: <strong>${{am !== null ? am.toFixed(3) : '&mdash;'}}</strong> &bull; ESM-2: <strong>${{e2 !== null ? e2.toFixed(2) : '&mdash;'}}</strong></td>
            <td>ESM-IF1 3D &Delta;logP: <strong>${{if1 !== null ? if1.toFixed(2) : '&mdash;'}}</strong></td>
          </tr>
        </tbody>
      </table>
    </div>
    <div class="mech-note">${{mechTxt}}</div>
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
    box.innerHTML = `<p style="font-size:0.82rem;color:#b91c1c;">Enter a codon number (e.g. <code>1705</code>) or single-AA substitution (e.g. <code>E1705K</code>).</p>`;
    return;
  }}
  const posStr = mVar[2];
  const pos = parseInt(posStr, 10);
  const entry = gCache[posStr];
  if (!entry) {{
    box.innerHTML = `<p style="font-size:0.82rem;color:#b91c1c;">Codon ${{pos}} is outside ${{gene}} canonical length (1..${{Object.keys(gCache).length}} aa).</p>`;
    return;
  }}
  const [wt, pdb, dna, dmet, dprt, plddt] = entry.m;
  let targetAlt = (mVar[3] && mVar[3] !== '*') ? mVar[3] : '';

  // Pick active mutant AA for 3D viewer (explicitly typed alt, or highest-impact/observed substitution at this codon)
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

  let html = `<div style="background:#f8fafc;border:1px solid #cbd5e1;border-radius:8px;padding:10px 12px;">
    <div style="font-size:0.85rem;margin-bottom:8px;display:flex;justify-content:space-between;flex-wrap:wrap;gap:8px;">
      <span>
        <strong>All Substitutions at ${{gene}} Codon ${{wt}}${{pos}}</strong> &bull; AFDB pLDDT: <strong>${{plddt ?? 'NA'}}</strong>
        &bull; Experimental PDB: <strong>${{pdb || 'Unresolved (AFDB)'}}</strong>
        &bull; Dist RNA/DNA: <strong>${{dna !== null ? dna + ' &Aring;' : '&mdash;'}}</strong>
        &bull; Dist Metal: <strong>${{dmet !== null ? dmet + ' &Aring;' : '&mdash;'}}</strong>
        &bull; Dist Partner: <strong>${{dprt !== null ? dprt + ' &Aring;' : '&mdash;'}}</strong>
      </span>
      ${{targetAlt ? `<button class="ghost-btn" onclick="document.getElementById('sat-query').value='${{pos}}';inspectSaturationCodon();">Show All 19 Substitutions at ${{wt}}${{pos}}</button>` : ''}}
    </div>
    <div class="table-wrap" style="max-height:245px;">
      <table>
        <thead><tr><th>3D</th><th>Substitution</th><th>Grantham D</th><th>&Delta;Charge</th><th>&Delta;Vol (&Aring;&sup3;)</th><th>AlphaMissense</th><th>ESM-2 (3B)</th><th>ESM-1v</th><th>ESM-IF1 3D</th><th>Biophysical Regime</th><th>TCGA</th><th>ClinVar Assertion</th></tr></thead>
        <tbody>`;
  rows.forEach(r => {{
    const [alt, am, e2, e1v, if1, tn, cv] = r;
    const wp = AA_PROPS[wt] || AA_PROPS.A, mp = AA_PROPS[alt] || wp;
    const dq = mp.q - wp.q, dv = mp.vol - wp.vol, gd = granthamDist(wt, alt);
    const constrained = (am !== null && am >= 0.564) || (e2 !== null && e2 <= -5.0);
    const regime = constrained ? ((if1 !== null && if1 > -1.5) ? '<span class="pill pill-red">Fold-Preserving Active/Interface</span>' : '<span class="pill pill-orange">Predicted Fold-Destabilising</span>') : '<span class="pill pill-blue">Tolerated / Neutral</span>';
    html += `<tr>
      <td><button class="ghost-btn" style="padding:2px 7px;font-size:0.73rem;" onclick="selectMutAA(${{pos}},'${{wt}}','${{alt}}')">3D</button></td>
      <td><strong><code>${{wt}}${{pos}}${{alt}}</code></strong></td>
      <td>${{gd}}</td>
      <td>${{dq>=0?'+':''}}${{dq.toFixed(1)}}</td>
      <td>${{dv>=0?'+':''}}${{dv.toFixed(0)}}</td>
      <td>${{am !== null ? am.toFixed(3) : '&mdash;'}}</td>
      <td>${{e2 !== null ? e2.toFixed(2) : '&mdash;'}}</td>
      <td>${{e1v !== null ? e1v.toFixed(2) : '&mdash;'}}</td>
      <td>${{if1 !== null ? if1.toFixed(2) : '&mdash;'}}</td>
      <td>${{regime}}</td>
      <td>${{tn > 0 ? '<strong>' + tn + '</strong>' : '0'}}</td>
      <td>${{cv ? '<span class="pill pill-amber">' + cv + '</span>' : '&mdash;'}}</td>
    </tr>`;
  }});
  html += `</tbody></table></div></div>`;
  box.innerHTML = html;
}}

function tierPill(t) {{
  if (t === '3D_Catalytic_or_Substrate_Pocket_Driver') return '<span class="pill pill-red">3D Pocket &le;6&Aring;</span>';
  if (t === 'Pathogenic_Fold_or_Interface_Missense') return '<span class="pill pill-orange">Pathogenic Fold</span>';
  if (t === 'Truncating_LoF') return '<span class="pill pill-dark">Truncating LoF</span>';
  if (t === 'Ambiguous_Missense_VUS') return '<span class="pill pill-amber">Ambiguous VUS</span>';
  if (t === 'Likely_Benign_Passenger_Missense') return '<span class="pill pill-blue">Likely Benign VUS</span>';
  return `<span class="pill pill-dark">${{t}}</span>`;
}}

function mechPill(m) {{
  if (m === 'Fold_Preserving_Catalytic_or_Interface') return '<span class="pill pill-red">Fold-Preserving Active/Interface</span>';
  if (m === 'Thermodynamic_Fold_Destabilizing') return '<span class="pill pill-orange">Fold-Destabilising Core</span>';
  if (m === 'Tolerated_Neutral') return '<span class="pill pill-blue">Tolerated Neutral</span>';
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
  document.getElementById('var-count').textContent = `Showing ${{Math.min(250, currentFilteredVars.length)}} of ${{currentFilteredVars.length}} unique protein variants (${{currentFilteredVars.reduce((a,b)=>a+b.n,0)}} TCGA somatic calls).`;
  document.querySelector('#var-table tbody').innerHTML = currentFilteredVars.slice(0, 250).map(v => {{
    const cleanP = v.p.replace(/^p\\./i, '');
    const isMissense = /^[A-Z]\\d+[A-Z]$/.test(cleanP);
    const jumpArg = isMissense ? cleanP : (v.pos > 0 ? String(v.pos) : '');
    const btn3d = jumpArg ? `<button class="ghost-btn" style="padding:2px 7px;font-size:0.73rem;" onclick="openVariantIn3D('${{v.g}}','${{jumpArg}}')">3D</button>` : '&mdash;';
    return `
    <tr>
      <td>${{btn3d}}</td>
      <td><strong>${{v.g}}</strong></td>
      <td><code>${{v.p}}</code></td>
      <td>${{v.n}}</td>
      <td>${{v.cv ? '<span class="pill pill-amber">' + v.cv + '</span>' : '&mdash;'}}</td>
      <td>${{tierPill(v.tier)}}</td>
      <td>${{mechPill(v.mech)}}</td>
      <td>${{v.pdb || '&mdash;'}}</td>
      <td>${{v.dna !== null ? v.dna.toFixed(2) : '&mdash;'}}</td>
      <td>${{v.dmet !== null ? v.dmet.toFixed(2) : '&mdash;'}}</td>
      <td>${{v.am !== null ? v.am.toFixed(3) + ' (' + v.amc + ')' : '&mdash;'}}</td>
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
  document.getElementById('gtex-bars').innerHTML = '<table><thead><tr><th>GTEx v10 Adult Tissue (54 Sites)</th><th>Median TPM</th><th>Relative Abundance</th></tr></thead><tbody>' +
    pairs.map(([t, v]) => `<tr>
      <td>${{t.replace(/_/g, ' ')}}</td>
      <td><strong>${{v.toFixed(2)}}</strong></td>
      <td><div style="width:${{Math.max(2, (v/mx)*280).toFixed(0)}}px;height:10px;background:#0f5499;border-radius:3px;"></div></td>
    </tr>`).join('') + '</tbody></table>';
}}

function renderAllTables() {{
  document.querySelector('#mir-bench-table tbody').innerHTML = ATLAS.mir_bench.map(b => `
    <tr>
      <td><strong>${{b.benchmark_id}}</strong></td>
      <td><code>${{b.metric}}</code></td>
      <td>${{b.group_A}} (n=${{b.n_group_A}})</td>
      <td>${{Number(b.median_group_A).toFixed(4)}}</td>
      <td>${{b.group_B}} (n=${{b.n_group_B}})</td>
      <td>${{Number(b.median_group_B).toFixed(4)}}</td>
      <td><strong>${{(Number(b.median_shift_A_minus_B)>=0?'+':'') + Number(b.median_shift_A_minus_B).toFixed(4)}}</strong></td>
      <td><strong>${{b.p_value_two_sided}}</strong></td>
    </tr>
  `).join('');

  const sortedDicer = [...ATLAS.dicer_pts].sort((a,b) => a.d53 - b.d53);
  document.querySelector('#dicer-mir-table tbody').innerHTML = sortedDicer.map(d => `
    <tr>
      <td><code>${{d.s}}</code></td>
      <td><strong>${{d.t}}</strong></td>
      <td>${{d.c === 'DICER1_RNaseIIIb_3D_Catalytic_Pocket' ? '<span class="pill pill-red">RNase IIIb 3D Pocket (&le;8&Aring;)</span>' : d.c}}</td>
      <td><strong>${{d.d53.toFixed(4)}}</strong></td>
      <td>${{d.f5 !== null ? d.f5.toFixed(4) : '&mdash;'}}</td>
      <td>${{d.l7z !== null && d.l7z !== undefined ? (d.l7z>=0?'+':'') + d.l7z.toFixed(3) : '&mdash;'}}</td>
    </tr>
  `).join('');

  document.querySelector('#top-axes-table tbody').innerHTML = ATLAS.top_axes.map(ax => `
    <tr>
      <td>#${{ax.rank}}</td>
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
      <td><strong>${{Number(r.pr_auc).toFixed(4)}}</strong> (base ${{Number(r.pr_baseline).toFixed(3)}})</td>
      <td><strong>${{r.mann_whitney_p}}</strong></td>
    </tr>
  `).join('');

  document.querySelector('#msk-rep-table tbody').innerHTML = ATLAS.msk_rep.map(m => `
    <tr>
      <td><strong>${{m.geneSymbol}}</strong></td>
      <td>${{Number(m.msk_N_profiled).toLocaleString()}}</td>
      <td><strong>${{(Number(m.msk_mut_freq)*100).toFixed(2)}}%</strong> (n=${{m.msk_mut_samples_n}})</td>
      <td>${{(Number(m.tcga_mut_freq_N10195)*100).toFixed(2)}}%</td>
      <td>${{(Number(m.msk_nonhypermut_mut_freq)*100).toFixed(2)}}%</td>
      <td>${{(Number(m.msk_truncating_freq)*100).toFixed(2)}}% (n=${{m.msk_truncating_n}})</td>
      <td><strong>${{(Number(m.msk_3dpocket_or_hotspot_freq)*100).toFixed(2)}}%</strong> (n=${{m.msk_3dpocket_or_hotspot_n}})</td>
      <td>${{(Number(m.msk_primary_mut_freq)*100).toFixed(2)}}% vs ${{(Number(m.msk_metastasis_mut_freq)*100).toFixed(2)}}%</td>
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
      <td>${{p.opentargets_sm_tractable === '1' ? '<span class="pill pill-green">SM Ligand/Pocket</span> ' : ''}}${{p.opentargets_protac_tractable === '1' ? '<span class="pill pill-blue">PROTAC</span>' : ''}}</td>
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
      <td>${{p.depmapIsCommonEssential === '1' ? '<span class="pill pill-red">Essential</span>' : '<span class="pill pill-blue">Buffered / Context</span>'}}</td>
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
      <td><strong>${{s.multivar_p_value || '&mdash;'}}</strong></td>
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
      <td><strong>${{i.p_value_two_sided}}</strong></td>
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
    selLol.insertAdjacentHTML('beforeend', `<option value="${{g}}" ${{g==='DICER1'?'selected':''}}>${{g}} (${{p.proteinLengthAa}} aa, ${{p.pdbComplexIds}})</option>`);
    selGtex.insertAdjacentHTML('beforeend', `<option value="${{g}}" ${{g==='LIN28B'?'selected':''}}>${{g}} (median ${{Number(p.gtex_pan54_median_tpm||0).toFixed(2)}} TPM, &tau;=${{p.gtex_tau_tissue_specificity}})</option>`);
  }});
  document.getElementById('sat-query').value = 'E1705K';
  renderHeatmap();
  renderVariants();
  renderAllTables();
  renderGtexBars();
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
