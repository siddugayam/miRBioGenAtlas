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
        <p>Sample-, variant-, and atomic-level reference of the human microRNA biogenesis machinery across TCGA, MSK-IMPACT, GDC, PDB, AlphaFold, ClinVar, gnomAD, GTEx, CPTAC, and DepMap.</p>
      </div>
      <div class="omni-search">
        <input type="text" id="global-search" placeholder="Jump to gene, mutation, or study (e.g. DICER1 E1705K, TP53 R175H)..." onkeydown="if(event.key==='Enter') runGlobalSearch(this.value)">
        <span class="omni-kbd">/</span>
      </div>
    </div>

    <div class="spec-strip">
      <div class="spec-item">
        <div class="spec-val">21,140</div>
        <div class="spec-lbl">Sequenced Tumours (TCGA &amp; MSK)</div>
      </div>
      <div class="spec-item">
        <div class="spec-val">15 Proteins</div>
        <div class="spec-lbl">11,260 Codons &middot; 89,132 Heavy Atoms</div>
      </div>
      <div class="spec-item">
        <div class="spec-val">213,940</div>
        <div class="spec-lbl">Missense Substitutions &amp; 3D Rotamers</div>
      </div>
      <div class="spec-item">
        <div class="spec-val">16 Complexes</div>
        <div class="spec-lbl">Experimental Cryo-EM &amp; X-ray PDBs</div>
      </div>
      <div class="spec-item">
        <div class="spec-val">9,578</div>
        <div class="spec-lbl">Mature miRNA-seq &amp; 13.37M Pairs</div>
      </div>
      <div class="spec-item">
        <div class="spec-val">0.978 AUC</div>
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
    </div>
  </div>
</nav>

<main class="container workspace">

  <!-- SECTION 1: STRUCTURE, LOLLIPOP & 213,940 SATURATION WORKBENCH -->
  <section id="tab-variants" class="panel active">
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
                <span><i class="mol3d-swatch" style="background:#e11d48;"></i>Mutant 3D rotamer</span>
                <span><i class="mol3d-swatch" style="background:#94a3b8;"></i>8&Aring; shell</span>
                <span><i class="mol3d-swatch" style="background:#f59e0b;"></i>Contact (&le;3.8&Aring;)</span>
                <span><i class="mol3d-swatch" style="background:#ff3b30;"></i>Steric clash (&lt;2.4&Aring;)</span>
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
              <th>TCGA</th>
              <th>ClinVar</th>
              <th>Structural Tier</th>
              <th>3D Mechanism</th>
              <th>PDB</th>
              <th>RNA/DNA (&Aring;)</th>
              <th>Metal (&Aring;)</th>
              <th>AlphaMissense</th>
              <th>ESM-2 3B</th>
              <th>ESM-1v</th>
              <th>ESM-IF1 &Delta;logP</th>
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
            <option value="nhm">Somatic nonsynonymous mutation (non-hypermutated, n = 9,360)</option>
            <option value="any">Broad alteration (mutation or any GISTIC2 &ne; 0)</option>
          </select>
        </div>
        <button class="btn-primary" onclick="exportLongCSV()">Download 480-Row CSV</button>
      </div>
      <div id="hm-container" class="heatmap-grid"></div>
      <p class="vp-footnote">Cells report exact study proportions in [0.00, 1.00]. Hover any cell for denominator N, somatic mutation count, +2 amplification count, and &minus;2 deep-deletion count. Studies marked with * have N &lt; 100.</p>
    </div>
  </section>

  <!-- SECTION 3: MIRNA-SEQ & TARGET COUPLINGS -->
  <section id="tab-mirna" class="panel">
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
              <th>Readout Metric</th>
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
              <th>DICER1 Class</th>
              <th>Paired Duplex log2((5p+1)/(3p+1))</th>
              <th>Linear 5p Fraction</th>
              <th>let-7 Target mRNA Z</th>
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
              <th>Pearson r (log2)</th>
              <th>t-Statistic</th>
            </tr>
          </thead>
          <tbody></tbody>
        </table>
      </div>
    </div>
  </section>

  <!-- SECTION 4: CLINVAR, GNOMAD & MSK-IMPACT -->
  <section id="tab-clinvar" class="panel">
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
              <th>ROC-AUC (95% CI)</th>
              <th>PR-AUC (Baseline)</th>
              <th>Mann&ndash;Whitney p</th>
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
              <th>MSK Non-Hypermut %</th>
              <th>MSK Truncating % (n)</th>
              <th>MSK 3D Pocket / Hotspot % (n)</th>
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
              <th>RNase IIIb Hotspot (n, %)</th>
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
              <th>gnomAD pLI</th>
              <th>LOEUF</th>
              <th>Missense Z</th>
              <th>HPA Subcellular Compartment</th>
              <th>GTEx Pan-54 Median TPM</th>
              <th>GTEx Peak Tissue (TPM)</th>
              <th>Specificity (&tau;)</th>
              <th>Tractability</th>
              <th>ClinVar Pathogenic / Total</th>
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
              <th>PDB Complexes</th>
              <th>CPTAC CNA&rarr;mRNA &rho;</th>
              <th>CPTAC CNA&rarr;Protein &rho;</th>
              <th>Attenuation (&Delta;&rho;)</th>
              <th>DepMap Median Chronos</th>
              <th>% Dependent (&le;&minus;0.5)</th>
              <th>Dependency Class</th>
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
              <th>Endpoint</th>
              <th>Altered (Events / n)</th>
              <th>Control (Events / n)</th>
              <th>Stratified HR (95% CI)</th>
              <th>Log-Rank p</th>
              <th>Multivariable Cox HR (95% CI)</th>
              <th>Multivariable p</th>
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
  if (/^[A-Z]?\\d+[A-Z]?$/.test(q)) {{
    document.getElementById('sat-query').value = q;
    switchTab('tab-variants', document.getElementById('btn-tab-variants'));
    inspectSaturationCodon();
    return;
  }}
  document.getElementById('var-q').value = raw.trim();
  switchTab('tab-variants', document.getElementById('btn-tab-variants'));
  renderVariants();
  document.getElementById('var-table').scrollIntoView({{behavior:'smooth'}});
}}

document.addEventListener('keydown', e => {{
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
            <td><strong>Net charge (pH 7.4)</strong></td>
            <td>${{wp.q >= 0 ? '+' : ''}}${{wp.q.toFixed(1)}} e</td>
            <td>${{mp.q >= 0 ? '+' : ''}}${{mp.q.toFixed(1)}} e</td>
            <td>${{dq >= 0 ? '+' : ''}}${{dq.toFixed(1)}} e &middot; ${{tagDot(chargeDot, chargeNote)}}</td>
          </tr>
          <tr>
            <td><strong>Side-chain volume</strong></td>
            <td>${{wp.vol.toFixed(1)}} &Aring;&sup3;</td>
            <td>${{mp.vol.toFixed(1)}} &Aring;&sup3;</td>
            <td>${{dvol >= 0 ? '+' : ''}}${{dvol.toFixed(1)}} &Aring;&sup3; (${{dvol >= 0 ? '+' : ''}}${{((dvol/wp.vol)*100).toFixed(0)}}%)</td>
          </tr>
          <tr>
            <td><strong>Hydropathy (KD)</strong></td>
            <td>${{wp.kd.toFixed(1)}}</td>
            <td>${{mp.kd.toFixed(1)}}</td>
            <td>${{dkd >= 0 ? '+' : ''}}${{dkd.toFixed(1)}} (${{dkd > 1.5 ? 'More hydrophobic' : (dkd < -1.5 ? 'More polar' : 'Comparable')}})</td>
          </tr>
          <tr>
            <td><strong>H-bond / metal</strong></td>
            <td>${{wp.hd}}D / ${{wp.ha}}A${{wp.metal ? ' · Metal' : ''}}</td>
            <td>${{mp.hd}}D / ${{mp.ha}}A${{mp.metal ? ' · Metal' : ''}}</td>
            <td>&Delta;D=${{mp.hd-wp.hd>=0?'+':''}}${{mp.hd-wp.hd}}, &Delta;A=${{mp.ha-wp.ha>=0?'+':''}}${{mp.ha-wp.ha}}</td>
          </tr>
          <tr>
            <td><strong>Grantham / alerts</strong></td>
            <td colspan="2">D = <strong>${{gDist}}</strong> &middot; ${{tagDot(gDot, gLabel)}}</td>
            <td>${{alerts.length ? alerts.join('; ') : 'None'}}</td>
          </tr>
          <tr>
            <td><strong>3D environment</strong></td>
            <td>${{burialStr}}</td>
            <td>WT ${{geom && geom.wtMinDist ? geom.wtMinDist.toFixed(2)+'&Aring;' : '&mdash;'}}</td>
            <td>${{clashTag}}</td>
          </tr>
          <tr>
            <td><strong>Model scores</strong></td>
            <td colspan="2">AM ${{am !== null ? am.toFixed(3) : '&mdash;'}} &middot; ESM-2 ${{e2 !== null ? e2.toFixed(2) : '&mdash;'}}</td>
            <td>ESM-IF1 &Delta;logP: <strong>${{if1 !== null ? if1.toFixed(2) : '&mdash;'}}</strong></td>
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
        <thead><tr><th>3D</th><th>Variant</th><th>Grantham</th><th>&Delta;q</th><th>&Delta;Vol (&Aring;&sup3;)</th><th>AlphaMissense</th><th>ESM-2 3B</th><th>ESM-1v</th><th>ESM-IF1 3D</th><th>Regime</th><th>TCGA</th><th>ClinVar</th></tr></thead>
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
  document.getElementById('gtex-bars').innerHTML = '<table><thead><tr><th>GTEx v10 Adult Tissue (54 Sites)</th><th>Median TPM</th><th>Relative Expression</th></tr></thead><tbody>' +
    pairs.map(([t, v]) => `<tr>
      <td>${{t.replace(/_/g, ' ')}}</td>
      <td><strong>${{v.toFixed(2)}}</strong></td>
      <td><div style="width:${{Math.max(2, (v/mx)*280).toFixed(0)}}px;height:8px;background:var(--embl-petrol);border-radius:2px;"></div></td>
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
