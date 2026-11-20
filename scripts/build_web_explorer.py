#!/usr/bin/env python3
"""Compile CSV/CSV.GZ tables into docs/index.html and docs/sat/<GENE>.js."""

import csv
import gzip
import json
import pathlib
from collections import Counter, defaultdict

REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent
DATA_DIR = REPO_ROOT / "data"
DOCS_DIR = REPO_ROOT / "docs"
SAT_DIR = DOCS_DIR / "sat"


def read_csv_rows(path: pathlib.Path):
    if not path.exists():
        return []
    opener = gzip.open if path.name.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def build_saturation_js_chunks(clinvar_rows):
    """Exports 15 compact `docs/sat/<GENE>.js` files covering all 213,940 variants + codon tracks."""
    sat_gz = DATA_DIR / "miRNA_biogenesis_saturation_mutagenesis_213940_atlas.csv.gz"
    res_gz = DATA_DIR / "miRNA_biogenesis_residue_biophysical_11260_profile.csv.gz"
    if not (sat_gz.exists() and res_gz.exists()):
        return {}

    SAT_DIR.mkdir(parents=True, exist_ok=True)

    # Count ClinVar pathogenic/benign per (gene, pos) and lookup per (gene, variant)
    cv_pos_path = Counter()
    cv_var_sig = {}
    for r in clinvar_rows:
        g = r["geneSymbol"]
        pos = int(r["codonPos"]) if str(r.get("codonPos", "")).isdigit() else 0
        if r.get("significanceCategory") == "Pathogenic_or_Likely_Pathogenic" and pos > 0:
            cv_pos_path[(g, pos)] += 1
        cv_var_sig[(g, r["proteinChange"])] = r["clinicalSignificance"]

    # Build compact codon biophysical tracks for SVG lollipop viewer
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

    # Group all 213,940 substitutions by gene -> codon -> list of [altAA, am, e2, e1v, if1, tcga_n, clinvar_sig]
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

    # Compact GTEx matrix: tissues list + {gene: [54 floats]}
    gtex_tissues = sorted({r["tissueSiteDetailId"] for r in gtex_rows})
    t_idx = {t: i for i, t in enumerate(gtex_tissues)}
    gtex_by_gene = {p["geneSymbol"]: [0.0] * len(gtex_tissues) for p in panel}
    for r in gtex_rows:
        if r["geneSymbol"] in gtex_by_gene:
            gtex_by_gene[r["geneSymbol"]][t_idx[r["tissueSiteDetailId"]]] = round(float(r["medianTPM"]), 2)

    bundle = {
        "meta": {
            "version": "2.0.1",
            "doi": "10.5281/zenodo.24393421",
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
<title>miRBioGenAtlas v2.0 — Multi-Consortium Pan-Cancer miRNA Biogenesis Explorer</title>
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
  .container {{ max-width: 1400px; margin: 0 auto; padding: 0 18px; }}
  header h1 {{ margin: 0 0 6px; font-size: 1.6rem; letter-spacing: -0.02em; }}
  header p {{ margin: 0; opacity: 0.93; font-size: 0.92rem; max-width: 1100px; }}
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
    <h1>miRBioGenAtlas v2.0 — Multi-Consortium Pan-Cancer miRNA Biogenesis &amp; Saturation Atlas</h1>
    <p>Harmonised across 21,140 clinical &amp; primary human tumours (10,195 TCGA PanCanAtlas + 10,945 MSK-IMPACT), 9,578 mature miRNA-seq profiles, 13.37M miRNA&times;mRNA couplings (9,353 dual tumours), 213,940 saturation-mutagenesis substitutions (16 Cryo-EM/X-ray PDB complexes, AlphaFold v6, AlphaMissense, ESM-2 3B, ESM-1v, ESM-IF1 3D), 16,569 ClinVar clinical records, gnomAD v4 constraint, GTEx 54-tissue baselines, HPA subcellular localisation, CPTAC proteomes (N=1,002), DepMap CRISPR screens (N=1,208), TCGA-CDR survival, and Thorsson PanImmune profiles.</p>
    <div class="badges">
      <a class="badge" href="https://doi.org/10.5281/zenodo.24393421" target="_blank" rel="noopener">Zenodo Dataset DOI: 10.5281/zenodo.24393421</a>
      <a class="badge" href="https://github.com/siddugayam/miRBioGenAtlas" target="_blank" rel="noopener">GitHub: siddugayam/miRBioGenAtlas</a>
      <span class="badge">Licence: CC BY 4.0 (Data) / MIT (Code)</span>
    </div>
  </div>
</header>

<main class="container">
  <div class="kpi-row">
    <div class="kpi"><div class="num">21,140</div><div class="lbl">Human Tumours (10,195 TCGA + 10,945 MSK-IMPACT)</div></div>
    <div class="kpi"><div class="num">213,940</div><div class="lbl">Saturation Variants (AM, ESM-2 3B, ESM-1v, ESM-IF1)</div></div>
    <div class="kpi"><div class="num">16,569</div><div class="lbl">ClinVar Records (Missense ROC-AUC = 0.978)</div></div>
    <div class="kpi"><div class="num">9,578</div><div class="lbl">Mature miRNA-seq &amp; 13.37M miRNA&times;mRNA Pairs</div></div>
    <div class="kpi"><div class="num">54 Tissues</div><div class="lbl">GTEx v10 Normal Adult TPM + HPA Subcellular</div></div>
    <div class="kpi"><div class="num">2,210</div><div class="lbl">Orthogonal Validation (1,002 CPTAC + 1,208 DepMap)</div></div>
  </div>

  <div class="tabs">
    <button class="tab-btn active" onclick="switchTab('tab-heatmap', this)">1. Pan-Cancer Heatmap (15&times;32)</button>
    <button class="tab-btn" onclick="switchTab('tab-variants', this)">2. SVG Lollipop &amp; 213,940 Saturation Lookup</button>
    <button class="tab-btn" onclick="switchTab('tab-mirna', this)">3. miRNA-seq &amp; 13.37M Target Couplings</button>
    <button class="tab-btn" onclick="switchTab('tab-clinvar', this)">4. ClinVar ROC-AUC, gnomAD &amp; MSK-IMPACT</button>
    <button class="tab-btn" onclick="switchTab('tab-gtex', this)">5. GTEx 54-Tissue, HPA &amp; Tractability</button>
    <button class="tab-btn" onclick="switchTab('tab-consortia', this)">6. CPTAC, DepMap, Survival &amp; PanImmune</button>
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

  <!-- TAB 2: INTERACTIVE SVG LOLLIPOP & 213,940 SATURATION CALCULATOR -->
  <section id="tab-variants" class="panel">
    <div class="card">
      <h3 style="margin:0 0 8px;">Interactive 3D Pocket, ClinVar Pathogenic &amp; TCGA Somatic Codon Architecture + 213,940-Variant Lookup</h3>
      <div class="controls">
        <label><strong>Protein Track &amp; Saturation Gene:</strong>
          <select id="lol-gene" onchange="renderLollipopAndLoadSat()"></select>
        </label>
        <label><strong>Inspect Any Codon or Variant (even unseen in TCGA):</strong>
          <input type="text" id="sat-query" placeholder="e.g. 1705, E1705K, E1705V, D597A..." oninput="inspectSaturationCodon()">
        </label>
        <span id="sat-status" style="font-size:0.8rem;color:var(--muted);"></span>
      </div>
      <div id="lollipop-svg-wrap" class="svg-box"></div>
      <div id="sat-codon-detail" style="margin-top:12px;"></div>
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
  miRBioGenAtlas v2.0 &bull; Manipal Academy of Higher Education (MAHE) &bull; Data: CC BY 4.0 &bull; Code: MIT &bull; Deterministically compiled from 45 SHA-256-verified snapshots.
</footer>

<script>
window.SAT_CACHE = window.SAT_CACHE || {{}};
const ATLAS = {data_json};
const LOW_N = new Set(["ACC","CHOL","DLBC","KICH","MESO","UCS","UVM"]);

function switchTab(id, btn) {{
  document.querySelectorAll('.panel').forEach(p => p.classList.remove('active'));
  document.querySelectorAll('.tab-btn').forEach(b => b.classList.remove('active'));
  document.getElementById(id).classList.add('active');
  btn.classList.add('active');
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
    document.getElementById('sat-status').textContent = `Loaded all ${{(Object.keys(window.SAT_CACHE[gene]).length * 19).toLocaleString()}} single-AA substitutions for ${{gene}}.`;
    if (cb) cb();
    return;
  }}
  document.getElementById('sat-status').textContent = `Loading 19x saturation matrix for ${{gene}}...`;
  const sc = document.createElement('script');
  sc.src = `sat/${{gene}}.js`;
  sc.onload = () => {{
    const nCod = window.SAT_CACHE[gene] ? Object.keys(window.SAT_CACHE[gene]).length : 0;
    document.getElementById('sat-status').textContent = `Loaded all ${{(nCod * 19).toLocaleString()}} single-AA substitutions (${{nCod}} codons) for ${{gene}}. Click any codon on the SVG track or type above!`;
    if (cb) cb();
  }};
  sc.onerror = () => {{
    document.getElementById('sat-status').textContent = `Codon track active (place docs/sat/${{gene}}.js alongside index.html for 19x per-residue lookup).`;
  }};
  document.head.appendChild(sc);
}}

function selectCodonFromTrack(pos) {{
  document.getElementById('sat-query').value = String(pos);
  inspectSaturationCodon();
}}

function renderLollipopAndLoadSat() {{
  const gene = document.getElementById('lol-gene').value;
  const tr = ATLAS.tracks[gene] || [];
  if (!tr.length) return;
  const L = tr.length;
  const W = 1260, H = 205, padL = 52, padR = 24;
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
      pockets += `<rect x="${{x-1.2}}" y="128" width="2.6" height="18" fill="#dc2626" opacity="0.65"><title>${{gene}} ${{wt}}${{pos}} (3D pocket &le;6&Aring; in ${{pdb}})</title></rect>`;
    }} else if (pdb) {{
      pockets += `<rect x="${{x-0.8}}" y="132" width="1.8" height="10" fill="#0ea5e9" opacity="0.35"></rect>`;
    }}
    if (am !== null) {{
      const yAm = 186 - am * 30;
      amPath += (amPath ? ' L' : 'M') + x.toFixed(1) + ',' + yAm.toFixed(1);
    }}
    if (tcgaN > 0 || cvN > 0) {{
      const h = 16 + (Math.log2(tcgaN + cvN + 1) / Math.log2(maxCnt + 1)) * 88;
      const yTop = 124 - h;
      const col = pkt === 1 ? '#dc2626' : (cvN > 0 ? '#7c3aed' : '#0f5499');
      const r = Math.min(6.5, 2.8 + Math.log2(tcgaN + cvN + 1));
      stems += `<g style="cursor:pointer" onclick="selectCodonFromTrack(${{pos}})">
        <line x1="${{x.toFixed(1)}}" y1="124" x2="${{x.toFixed(1)}}" y2="${{yTop.toFixed(1)}}" stroke="${{col}}" stroke-width="1.3" opacity="0.75"/>
        <circle cx="${{x.toFixed(1)}}" cy="${{yTop.toFixed(1)}}" r="${{r.toFixed(1)}}" fill="${{col}}" stroke="#fff" stroke-width="0.7">
          <title>${{gene}} ${{wt}}${{pos}} | TCGA somatic calls: ${{tcgaN}} | ClinVar Pathogenic: ${{cvN}} | Mean AlphaMissense: ${{am ?? 'NA'}} | Click to inspect all 19 substitutions</title>
        </circle>
      </g>`;
    }}
  }});

  let ticks = '';
  const step = L > 1200 ? 200 : (L > 500 ? 100 : 50);
  for (let p = 1; p <= L; p += step) {{
    const tx = xOf(p);
    ticks += `<text x="${{tx.toFixed(1)}}" y="200" font-size="10" fill="#475569" text-anchor="middle">${{p}}</text>`;
  }}
  ticks += `<text x="${{xOf(L).toFixed(1)}}" y="200" font-size="10" fill="#475569" text-anchor="end">${{L}} aa</text>`;

  const svg = `<svg viewBox="0 0 ${{W}} ${{H}}" width="100%" height="205">
    <text x="6" y="24" font-size="10.5" font-weight="700" fill="#0f172a">Somatic / ClinVar</text>
    <text x="6" y="140" font-size="10.5" font-weight="700" fill="#0f172a">3D PDB / Pocket</text>
    <text x="6" y="174" font-size="10.5" font-weight="700" fill="#0f172a">AlphaMissense</text>
    <rect x="${{padL}}" y="128" width="${{plotW}}" height="18" rx="4" fill="#e2e8f0"/>
    ${{pockets}}
    ${{stems}}
    <path d="${{amPath}}" fill="none" stroke="#ea580c" stroke-width="1.15" opacity="0.85"/>
    <line x1="${{padL}}" y1="${{186 - 0.564*30}}" x2="${{W-padR}}" y2="${{186 - 0.564*30}}" stroke="#dc2626" stroke-dasharray="3,3" stroke-width="0.8"/>
    ${{ticks}}
    <g transform="translate(${{W - 430}}, 14)" font-size="10.5" fill="#334155">
      <circle cx="0" cy="-3" r="4" fill="#dc2626"/><text x="8" y="0">3D Pocket (&le;6&Aring;)</text>
      <circle cx="115" cy="-3" r="4" fill="#7c3aed"/><text x="123" y="0">ClinVar Pathogenic</text>
      <circle cx="240" cy="-3" r="4" fill="#0f5499"/><text x="248" y="0">TCGA Somatic</text>
      <rect x="335" y="-7" width="8" height="8" fill="#0ea5e9" opacity="0.5"/><text x="347" y="0">PDB Resolved</text>
    </g>
  </svg>`;
  document.getElementById('lollipop-svg-wrap').innerHTML = svg;
  ensureSatLoaded(gene, inspectSaturationCodon);
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
  const mVar = rawQ.match(/^([A-Z])?(\\d+)([A-Z])?$/);
  if (!mVar) {{
    box.innerHTML = `<p style="font-size:0.82rem;color:#b91c1c;">Enter a codon number (e.g. <code>1705</code>) or single-AA substitution (e.g. <code>E1705K</code>).</p>`;
    return;
  }}
  const pos = mVar[2];
  const targetAlt = mVar[3] || '';
  const entry = gCache[pos];
  if (!entry) {{
    box.innerHTML = `<p style="font-size:0.82rem;color:#b91c1c;">Codon ${{pos}} is outside ${{gene}} canonical length (1..${{Object.keys(gCache).length}} aa).</p>`;
    return;
  }}
  const [wt, pdb, dna, dmet, dprt, plddt] = entry.m;
  let rows = entry.s;
  if (targetAlt) rows = rows.filter(r => r[0] === targetAlt);

  let html = `<div style="background:#f8fafc;border:1px solid #cbd5e1;border-radius:8px;padding:10px 12px;">
    <div style="font-size:0.85rem;margin-bottom:8px;">
      <strong>${{gene}} Codon ${{wt}}${{pos}}</strong> &bull; AFDB pLDDT: <strong>${{plddt ?? 'NA'}}</strong>
      &bull; Experimental PDB: <strong>${{pdb || 'Unresolved (AFDB)'}}</strong>
      &bull; Dist RNA/DNA: <strong>${{dna !== null ? dna + ' &Aring;' : '&mdash;'}}</strong>
      &bull; Dist Metal: <strong>${{dmet !== null ? dmet + ' &Aring;' : '&mdash;'}}</strong>
      &bull; Dist Partner: <strong>${{dprt !== null ? dprt + ' &Aring;' : '&mdash;'}}</strong>
    </div>
    <div class="table-wrap" style="max-height:260px;">
      <table>
        <thead><tr><th>Substitution</th><th>AlphaMissense</th><th>ESM-2 (3B) LLR</th><th>ESM-1v LLR</th><th>ESM-IF1 3D LLR</th><th>Biophysical Regime</th><th>TCGA Tumours</th><th>ClinVar Assertion</th></tr></thead>
        <tbody>`;
  rows.forEach(r => {{
    const [alt, am, e2, e1v, if1, tn, cv] = r;
    const constrained = (am !== null && am >= 0.564) || (e2 !== null && e2 <= -5.0);
    const regime = constrained ? ((if1 !== null && if1 > -1.5) ? '<span class="pill pill-red">Fold-Preserving Active/Interface</span>' : '<span class="pill pill-orange">Predicted Fold-Destabilising</span>') : '<span class="pill pill-blue">Tolerated / Neutral</span>';
    html += `<tr>
      <td><strong><code>${{wt}}${{pos}}${{alt}}</code></strong></td>
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
  document.querySelector('#var-table tbody').innerHTML = currentFilteredVars.slice(0, 250).map(v => `
    <tr>
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
    </tr>
  `).join('');
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
