#!/usr/bin/env python3
"""Integrate 3D structural distances, AlphaMissense, PLM/IF1 scores, mature miRNA-seq, CPTAC, and DepMap layers."""

import csv
import gzip
import json
import math
import pathlib
import re
import statistics
from collections import Counter, defaultdict

MISSENSE_PARSE_RE = re.compile(r"^p?\.?([A-Z]\d+[A-Z])$")


def mann_whitney_z_and_p(xs, ys):
    """Tie-corrected two-sided Mann-Whitney U normal approximation."""
    nx, ny = len(xs), len(ys)
    if nx < 2 or ny < 2:
        return float("nan"), float("nan"), float("nan")
    combined = [(v, 0) for v in xs] + [(v, 1) for v in ys]
    combined.sort(key=lambda t: t[0])
    n_tot = nx + ny
    ranks = [0.0] * n_tot
    tie_sum = 0.0
    i = 0
    while i < n_tot:
        j = i
        while j + 1 < n_tot and combined[j + 1][0] == combined[i][0]:
            j += 1
        t = j - i + 1
        if t > 1:
            tie_sum += t ** 3 - t
        avg_r = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            ranks[k] = avg_r
        i = j + 1
    r_x = sum(ranks[idx] for idx in range(n_tot) if combined[idx][1] == 0)
    u_x = r_x - nx * (nx + 1) / 2.0
    mu_u = nx * ny / 2.0
    var_u = (nx * ny / 12.0) * ((n_tot + 1) - tie_sum / (n_tot * (n_tot - 1)))
    if var_u <= 0:
        return u_x, 0.0, 1.0
    z = (u_x - mu_u) / math.sqrt(var_u)
    p_two = math.erfc(abs(z) / math.sqrt(2.0))
    r_rb = 2.0 * (u_x / (nx * ny)) - 1.0
    return round(r_rb, 4), round(z, 4), p_two


def load_gz_json(path: pathlib.Path):
    with gzip.open(path, "rt", encoding="utf-8") as f:
        return json.load(f)


def run_multimodal_enrichment(
    repo_root: pathlib.Path,
    data_dir: pathlib.Path,
    gene_panel: list,
    gene_order: list,
    panel_rows: list,
    annotated_mut_rows: list,
    patient_burden_rows: list,
    study_summary_rows: list,
    sample_level_rows: list,
    spearman_rho_fn,
    write_csv_fn,
):
    ext_dir = data_dir / "raw_snapshot" / "external"
    if not (ext_dir / "structural_3d_and_uniprot_15genes.json.gz").exists():
        return

    print("\n4. Integrating 3D structures, AlphaMissense, PLM/IF1, TCGA miRNA-seq, CPTAC, and DepMap...")
    struct_snap = load_gz_json(ext_dir / "structural_3d_and_uniprot_15genes.json.gz")
    am_snap = load_gz_json(ext_dir / "alphamissense_15genes_lookup.json.gz")
    mir_snap = load_gz_json(ext_dir / "tcga_pancan_mirnaseq_per_sample.json.gz")
    cptac_snap = load_gz_json(ext_dir / "cptac_pancan_15genes.json.gz")
    depmap_snap = load_gz_json(ext_dir / "depmap_crispr_and_ccle_15genes.json.gz")

    plm_var_map = {}
    mir_pca_map = {}
    risc_derep_map = {}
    top_axes = []
    plm_snap_gz = ext_dir / "plm_if1_saturation_and_mirna_svd.json.gz"
    if plm_snap_gz.exists():
        psnap = load_gz_json(plm_snap_gz)
        for r in psnap.get("variant_scores", []):
            plm_var_map[(r["geneSymbol"], r["proteinChange"])] = r
        for r in psnap.get("sample_mirna_svd", []):
            mir_pca_map[r["sampleId"]] = r
        for r in psnap.get("sample_risc_derepression", []):
            risc_derep_map[r["sampleId"]] = r
        top_axes = psnap.get("top250_repressive_axes", [])

    var_csv = data_dir / "tcga_somatic_variant_plm_if1_scores.csv"
    if var_csv.exists():
        with open(var_csv, "r", encoding="utf-8") as f:
            for r in csv.DictReader(f):
                plm_var_map[(r["geneSymbol"], r["proteinChange"])] = r
    mir_pca_csv = data_dir / "tcga_pancan_sample_mirna_pca_scores.csv"
    if mir_pca_csv.exists():
        with open(mir_pca_csv, "r", encoding="utf-8") as f:
            for r in csv.DictReader(f):
                mir_pca_map[r["sampleId"]] = r
    derep_csv = data_dir / "tcga_pancan_sample_risc_target_derepression.csv"
    if derep_csv.exists():
        with open(derep_csv, "r", encoding="utf-8") as f:
            for r in csv.DictReader(f):
                risc_derep_map[r["sampleId"]] = r
    if top_axes:
        write_csv_fn(
            data_dir / "tcga_pancan_mirna_mrna_top250_repressive_axes.csv",
            list(top_axes[0].keys()),
            top_axes,
        )

    ess_flags = depmap_snap.get("gene_is_common_essential_opentargets", {})
    for row in panel_rows:
        sym = row["geneSymbol"]
        gm = struct_snap["gene_metadata"][sym]
        row["uniprotAccession"] = gm["uniprotAccession"]
        row["afdbEntryId"] = gm["afdbEntryId"]
        row["afdbGlobalPlddt"] = gm["afdbGlobalPlddt"]
        row["afdbDistinctDomains"] = gm["afdbDistinctDomains"]
        row["pdbComplexIds"] = ";".join(gm["pdbStructuresUsed"])
        row["pdbResolvedResidues"] = gm["pdbResolvedResidues"]
        row["pdbCoverageFraction"] = gm["pdbCoverageFraction"]
        row["depmapIsCommonEssential"] = int(bool(ess_flags.get(sym, False)))
    write_csv_fn(
        data_dir / "miRNA_biogenesis_gene_panel_metadata.csv",
        list(panel_rows[0].keys()),
        panel_rows,
    )

    am_full = {"P": "likely_pathogenic", "A": "ambiguous", "B": "likely_benign"}
    v2_tier_counts = Counter()
    sample_has_3d_driver = defaultdict(set)
    sample_dicer1_class = {}

    for mrow in annotated_mut_rows:
        sym = mrow["geneSymbol"]
        gm = struct_snap["gene_metadata"][sym]
        per_res = struct_snap["per_residue"][sym]
        am_gene = am_snap["genes"][sym]

        pos_s = str(mrow.get("proteinPosStart") or "").strip()
        res_info = per_res.get(pos_s)
        plddt = res_info[0] if res_info else ""
        pdb_id = res_info[1] if res_info else ""
        d_na = res_info[2] if (res_info and res_info[2] >= 0) else ""
        d_met = res_info[3] if (res_info and res_info[3] >= 0) else ""
        d_prt = res_info[4] if (res_info and res_info[4] >= 0) else ""

        is_3d_pocket = 0
        for dval in (d_na, d_met, d_prt):
            if dval != "" and float(dval) <= 6.0:
                is_3d_pocket = 1

        pch = mrow["proteinChange"].strip()
        m_sub = MISSENSE_PARSE_RE.match(pch)
        am_score = ""
        am_cls = "non_missense"
        if m_sub:
            hit = am_gene["subs"].get(m_sub.group(1))
            if hit:
                am_score = hit[0]
                am_cls = am_full.get(hit[1], hit[1])
        am_res_mean = am_gene["res_mean"].get(pos_s, "")

        base_tier = mrow["consequenceTier"]
        if base_tier == "Truncating_LoF":
            v2_tier = "Truncating_LoF"
        elif is_3d_pocket == 1 and (mrow["isCuratedHotspot"] == 1 or (am_score != "" and am_score >= 0.564)):
            v2_tier = "3D_Catalytic_or_Substrate_Pocket_Driver"
        elif mrow["isCuratedHotspot"] == 1 or am_cls == "likely_pathogenic":
            v2_tier = "Pathogenic_Fold_or_Interface_Missense"
        elif am_cls == "likely_benign":
            v2_tier = "Likely_Benign_Passenger_Missense"
        elif am_cls == "ambiguous":
            v2_tier = "Ambiguous_Missense_VUS"
        else:
            v2_tier = base_tier

        v2_tier_counts[v2_tier] += 1
        if v2_tier in ("Truncating_LoF", "3D_Catalytic_or_Substrate_Pocket_Driver", "Pathogenic_Fold_or_Interface_Missense"):
            sample_has_3d_driver[mrow["sampleId"]].add(sym)

        if sym == "DICER1":
            sid = mrow["sampleId"]
            pos_int = int(pos_s) if pos_s.isdigit() else 0
            is_r3b_cleft = ((1690 <= pos_int <= 1830) or pos_int == 1344) and (
                mrow["isCuratedHotspot"] == 1
                or (d_met != "" and float(d_met) <= 8.0)
                or (d_na != "" and float(d_na) <= 6.0 and am_cls == "likely_pathogenic")
            )
            if is_r3b_cleft:
                sample_dicer1_class[sid] = "DICER1_RNaseIIIb_3D_Catalytic_Pocket"
            elif sid not in sample_dicer1_class or sample_dicer1_class[sid] != "DICER1_RNaseIIIb_3D_Catalytic_Pocket":
                if base_tier == "Truncating_LoF":
                    sample_dicer1_class[sid] = "DICER1_Truncating_LoF"
                elif am_cls == "likely_pathogenic" or mrow["isCuratedHotspot"] == 1:
                    sample_dicer1_class[sid] = "DICER1_Other_Domain_Pathogenic_Missense"
                elif am_cls == "likely_benign" and sample_dicer1_class.get(sid) is None:
                    sample_dicer1_class[sid] = "DICER1_Likely_Benign_Passenger_VUS"
                elif sample_dicer1_class.get(sid) is None:
                    sample_dicer1_class[sid] = "DICER1_Ambiguous_VUS"

        mrow["uniprotAccession"] = gm["uniprotAccession"]
        mrow["afdbResiduePlddt"] = plddt
        mrow["pdbComplexId"] = pdb_id
        mrow["distToSubstrateA"] = d_na
        mrow["distToCatalyticMetalA"] = d_met
        mrow["distToPartnerInterfaceA"] = d_prt
        mrow["is3dPocketOrInterface_le6A"] = is_3d_pocket
        mrow["alphamissenseScore"] = am_score
        mrow["alphamissenseClass"] = am_cls
        mrow["alphamissenseResidueMean"] = am_res_mean
        mrow["functionalVariantTierV2"] = v2_tier

        if plm_var_map:
            prec = plm_var_map.get((sym, pch), {})
            mrow["esm2_3b_llr"] = prec.get("esm2_3b_llr", "")
            mrow["esm1v_650m_llr"] = prec.get("esm1v_650m_llr", "")
            mrow["esm2_3b_effect_tier"] = prec.get("esm2_3b_effect_tier", "")
            mrow["esm_if1_3d_dlogp"] = prec.get("esm_if1_3d_dlogp", "")
            mrow["esm_if1_3d_marginal_llr"] = prec.get("esm_if1_3d_marginal_llr", "")
            mrow["biophysicalMechanismClass"] = prec.get("biophysicalMechanismClass", "")

    write_csv_fn(
        data_dir / "tcga_pancan_miRNA_biogenesis_somatic_mutations_annotated.csv",
        list(annotated_mut_rows[0].keys()),
        annotated_mut_rows,
    )

    if plm_var_map:
        uniq_by_gene = defaultdict(dict)
        for m in annotated_mut_rows:
            if m["alphamissenseScore"] != "" and m.get("esm2_3b_llr", "") != "" and m.get("esm1v_650m_llr", "") != "":
                uniq_by_gene[m["geneSymbol"]][m["proteinChange"]] = m
        plm_audit_rows = []
        all_u = [m for gmap in uniq_by_gene.values() for m in gmap.values()]
        for sym in list(gene_order) + ["PAN_15_POOLED"]:
            sub = all_u if sym == "PAN_15_POOLED" else list(uniq_by_gene[sym].values())
            if len(sub) < 5:
                continue
            am_v = [float(x["alphamissenseScore"]) for x in sub]
            e2_v = [-float(x["esm2_3b_llr"]) for x in sub]
            e1_v = [-float(x["esm1v_650m_llr"]) for x in sub]
            if1_sub = [x for x in sub if x.get("esm_if1_3d_dlogp", "") != ""]
            pock_e2 = [float(x["esm2_3b_llr"]) for x in sub if x["functionalVariantTierV2"] == "3D_Catalytic_or_Substrate_Pocket_Driver"]
            ben_e2 = [float(x["esm2_3b_llr"]) for x in sub if x["functionalVariantTierV2"] == "Likely_Benign_Passenger_Missense"]
            pock_if1 = [float(x["esm_if1_3d_dlogp"]) for x in if1_sub if x["functionalVariantTierV2"] == "3D_Catalytic_or_Substrate_Pocket_Driver"]
            fold_if1 = [float(x["esm_if1_3d_dlogp"]) for x in if1_sub if x["functionalVariantTierV2"] == "Pathogenic_Fold_or_Interface_Missense"]
            ben_if1 = [float(x["esm_if1_3d_dlogp"]) for x in if1_sub if x["functionalVariantTierV2"] == "Likely_Benign_Passenger_Missense"]
            plm_audit_rows.append({
                "geneSymbol": sym,
                "n_unique_missense_scored": len(sub),
                "spearman_AlphaMissense_vs_negESM2_3B": round(spearman_rho_fn(am_v, e2_v), 4),
                "spearman_AlphaMissense_vs_negESM1v_650M": round(spearman_rho_fn(am_v, e1_v), 4),
                "spearman_ESM2_3B_vs_ESM1v_650M": round(spearman_rho_fn(e2_v, e1_v), 4),
                "spearman_AlphaMissense_vs_negESM_IF1_3D": round(
                    spearman_rho_fn([float(x["alphamissenseScore"]) for x in if1_sub], [-float(x["esm_if1_3d_dlogp"]) for x in if1_sub]), 4
                ) if len(if1_sub) >= 5 else "",
                "median_ESM2_3B_LLR_3DPocket": round(statistics.median(pock_e2), 4) if pock_e2 else "",
                "median_ESM2_3B_LLR_LikelyBenign": round(statistics.median(ben_e2), 4) if ben_e2 else "",
                "median_ESM_IF1_3D_dlogP_3DPocket": round(statistics.median(pock_if1), 4) if pock_if1 else "",
                "median_ESM_IF1_3D_dlogP_FoldPathogenic": round(statistics.median(fold_if1), 4) if fold_if1 else "",
                "median_ESM_IF1_3D_dlogP_LikelyBenign": round(statistics.median(ben_if1), 4) if ben_if1 else "",
            })
        write_csv_fn(
            data_dir / "validation_plm_alphamissense_concordance.csv",
            list(plm_audit_rows[0].keys()),
            plm_audit_rows,
        )

    mir_by_sid = mir_snap["sample_readouts"]
    mirna_table_rows = []
    mir_matched_by_study = Counter()
    core10_syms = ("DROSHA", "DGCR8", "XPO5", "DICER1", "TARBP2", "PRKRA", "AGO1", "AGO2", "AGO3", "AGO4")

    for pb in patient_burden_rows:
        sid = pb["sampleId"]
        acro = pb["tumour_acronym"]
        mrec = mir_by_sid.get(sid)
        drec = risc_derep_map.get(sid, {})
        pb["core10_has_3d_or_am_pathogenic_mut"] = int(any(g in sample_has_3d_driver.get(sid, set()) for g in core10_syms))
        if mrec:
            mir_matched_by_study[acro] += 1
            pb["hasMatchedMirnaSeq"] = 1
            pb["mirna_5p_arm_fraction"] = mrec["arm_5p_frac"] if mrec["arm_5p_frac"] is not None else ""
            pb["mirna_duplex_5p_3p_log2_ratio"] = mrec["duplex_5p_3p_log2"] if mrec["duplex_5p_3p_log2"] is not None else ""
            pb["mirna_let7_family_log2_rpm"] = mrec["let7_log2_rpm"] if mrec["let7_log2_rpm"] is not None else ""
            pb["mirna_mir200_family_log2_rpm"] = mrec["mir200_log2_rpm"] if mrec["mir200_log2_rpm"] is not None else ""
            pb["mirna_mir34_family_log2_rpm"] = mrec["mir34_log2_rpm"] if mrec["mir34_log2_rpm"] is not None else ""
            pb["mirna_mir451a_minus_144_5p_log2"] = mrec["mir451a_minus_144_5p"] if mrec["mir451a_minus_144_5p"] is not None else ""
            pb["mirna_global_median_log2_rpm"] = mrec["global_median_log2_rpm"]
            pb["mirna_shannon_entropy_bits"] = mrec["shannon_entropy_bits"]

            row_out = {
                "tumour_acronym": acro,
                "studyId": pb["studyId"],
                "sampleId": sid,
                "patientId": pb["patientId"],
                "dicer1FunctionalClass": sample_dicer1_class.get(sid, "DICER1_WildType"),
                "mirna_5p_arm_fraction": pb["mirna_5p_arm_fraction"],
                "mirna_duplex_5p_3p_log2_ratio": pb["mirna_duplex_5p_3p_log2_ratio"],
                "n_expressed_canonical_duplexes": mrec["n_expr_duplexes"],
                "let7_family_mean_log2_rpm": pb["mirna_let7_family_log2_rpm"],
                "mir200_family_mean_log2_rpm": pb["mirna_mir200_family_log2_rpm"],
                "mir34_family_mean_log2_rpm": pb["mirna_mir34_family_log2_rpm"],
                "mir451a_log2_rpm": mrec["mir451a_log2_rpm"] if mrec["mir451a_log2_rpm"] is not None else "",
                "mir451a_minus_mir144_5p_log2": pb["mirna_mir451a_minus_144_5p_log2"],
                "global_canonical_median_log2_rpm": mrec["global_median_log2_rpm"],
                "n_detected_mature_mirnas_rpm_ge1": mrec["n_detected_rpm_ge1"],
                "mirnaome_shannon_entropy_bits": mrec["shannon_entropy_bits"],
            }
            if mir_pca_map:
                pc1 = mir_pca_map.get(sid, {}).get("mirnaome_pc1_score", "")
                row_out["mirnaome_pc1_score"] = pc1
                pb["mirnaome_pc1_score"] = pc1
            if risc_derep_map:
                for dcol in ("let7_target_oncofetal_mrna_z", "mir200_target_emt_mrna_z", "empirical_5p_target_derepression_z", "global_risc_target_derepression_z"):
                    row_out[dcol] = drec.get(dcol, "")
                    pb[dcol] = drec.get(dcol, "")
            mirna_table_rows.append(row_out)
        else:
            pb["hasMatchedMirnaSeq"] = 0
            pb["mirna_5p_arm_fraction"] = ""
            pb["mirna_duplex_5p_3p_log2_ratio"] = ""
            pb["mirna_let7_family_log2_rpm"] = ""
            pb["mirna_mir200_family_log2_rpm"] = ""
            pb["mirna_mir34_family_log2_rpm"] = ""
            pb["mirna_mir451a_minus_144_5p_log2"] = ""
            pb["mirna_global_median_log2_rpm"] = ""
            pb["mirna_shannon_entropy_bits"] = ""
            if mir_pca_map:
                pb["mirnaome_pc1_score"] = ""
            if risc_derep_map:
                for dcol in ("let7_target_oncofetal_mrna_z", "mir200_target_emt_mrna_z", "empirical_5p_target_derepression_z", "global_risc_target_derepression_z"):
                    pb[dcol] = ""

    write_csv_fn(
        data_dir / "tcga_pancan_miRNA_biogenesis_patient_level_burden.csv",
        list(patient_burden_rows[0].keys()),
        patient_burden_rows,
    )
    write_csv_fn(
        data_dir / "tcga_pancan_miRNA_biogenesis_mirnaseq_functional_readouts.csv",
        list(mirna_table_rows[0].keys()),
        mirna_table_rows,
    )

    for srow in study_summary_rows:
        srow["sampleCount_mirnaseq_matched_n"] = mir_matched_by_study[srow["tumour_acronym"]]
    write_csv_fn(
        data_dir / "tcga_pancan_atlas_sample_set_summary.csv",
        list(study_summary_rows[0].keys()),
        study_summary_rows,
    )

    sl_lookup = {(r["sampleId"], r["geneSymbol"]): r for r in sample_level_rows}
    bench_rows = []

    def add_benchmark(test_id, group_a_label, group_b_label, metric_name, vals_a, vals_b):
        if len(vals_a) < 2 or len(vals_b) < 2:
            return
        r_rb, z_val, p_val = mann_whitney_z_and_p(vals_a, vals_b)
        bench_rows.append({
            "benchmark_id": test_id,
            "metric": metric_name,
            "group_A": group_a_label,
            "n_group_A": len(vals_a),
            "median_group_A": round(statistics.median(vals_a), 5),
            "mean_group_A": round(statistics.mean(vals_a), 5),
            "group_B": group_b_label,
            "n_group_B": len(vals_b),
            "median_group_B": round(statistics.median(vals_b), 5),
            "mean_group_B": round(statistics.mean(vals_b), 5),
            "median_shift_A_minus_B": round(statistics.median(vals_a) - statistics.median(vals_b), 5),
            "rank_biserial_r": r_rb,
            "mann_whitney_z": z_val,
            "p_value_two_sided": f"{p_val:.3e}",
        })

    wt_5p_frac, r3b_5p_frac, ben_5p_frac, trunc_5p_frac = [], [], [], []
    wt_dup, r3b_dup, ben_dup, trunc_dup = [], [], [], []
    wt_451_144, r3b_451_144 = [], []

    for r in mirna_table_rows:
        cls = r["dicer1FunctionalClass"]
        f5 = r["mirna_5p_arm_fraction"]
        d53 = r["mirna_duplex_5p_3p_log2_ratio"]
        m451 = r["mir451a_minus_mir144_5p_log2"]
        if f5 != "" and d53 != "":
            if cls == "DICER1_WildType":
                wt_5p_frac.append(float(f5))
                wt_dup.append(float(d53))
                if m451 != "":
                    wt_451_144.append(float(m451))
            elif cls == "DICER1_RNaseIIIb_3D_Catalytic_Pocket":
                r3b_5p_frac.append(float(f5))
                r3b_dup.append(float(d53))
                if m451 != "":
                    r3b_451_144.append(float(m451))
            elif cls == "DICER1_Likely_Benign_Passenger_VUS":
                ben_5p_frac.append(float(f5))
                ben_dup.append(float(d53))
            elif cls == "DICER1_Truncating_LoF":
                trunc_5p_frac.append(float(f5))
                trunc_dup.append(float(d53))

    add_benchmark("B1_DICER1_RNaseIIIb_vs_WT_DuplexLog2", "DICER1_RNaseIIIb_3D_Catalytic_Pocket", "DICER1_WildType", "mirna_duplex_5p_3p_log2_ratio", r3b_dup, wt_dup)
    add_benchmark("B2_DICER1_RNaseIIIb_vs_BenignVUS_DuplexLog2", "DICER1_RNaseIIIb_3D_Catalytic_Pocket", "DICER1_Likely_Benign_Passenger_VUS", "mirna_duplex_5p_3p_log2_ratio", r3b_dup, ben_dup)
    add_benchmark("B3_DICER1_BenignVUS_vs_WT_NegativeControl", "DICER1_Likely_Benign_Passenger_VUS", "DICER1_WildType", "mirna_duplex_5p_3p_log2_ratio", ben_dup, wt_dup)
    add_benchmark("B4_DICER1_RNaseIIIb_vs_WT_5pLinearFraction", "DICER1_RNaseIIIb_3D_Catalytic_Pocket", "DICER1_WildType", "mirna_5p_arm_fraction", r3b_5p_frac, wt_5p_frac)
    add_benchmark("B5_DICER1_RNaseIIIb_AGO2_Bypass_miR451a_vs_miR144", "DICER1_RNaseIIIb_3D_Catalytic_Pocket", "DICER1_WildType", "mir451a_minus_mir144_5p_log2", r3b_451_144, wt_451_144)

    lin28_hi_let7, lin28_lo_let7 = [], []
    for r in mirna_table_rows:
        sid = r["sampleId"]
        l7 = r["let7_family_mean_log2_rpm"]
        if l7 == "":
            continue
        za = sl_lookup.get((sid, "LIN28A"), {}).get("mrnaZscore", "")
        zb = sl_lookup.get((sid, "LIN28B"), {}).get("mrnaZscore", "")
        z_vals = [float(z) for z in (za, zb) if z != ""]
        if not z_vals:
            continue
        max_z = max(z_vals)
        if max_z >= 2.0:
            lin28_hi_let7.append(float(l7))
        elif max_z <= 0.0:
            lin28_lo_let7.append(float(l7))
    add_benchmark("B6_LIN28AB_Reactivated_vs_Silent_Let7Family", "LIN28A_or_LIN28B_mRNA_Z_ge_2", "LIN28AB_Silent_Z_le_0", "let7_family_mean_log2_rpm", lin28_hi_let7, lin28_lo_let7)

    tp53_alt_m34, tp53_wt_m34 = [], []
    for r in mirna_table_rows:
        sid = r["sampleId"]
        m34 = r["mir34_family_mean_log2_rpm"]
        if m34 == "":
            continue
        t_rec = sl_lookup.get((sid, "TP53"), {})
        if t_rec.get("hasHighConfAlteration") == 1:
            tp53_alt_m34.append(float(m34))
        elif t_rec.get("hasAnyAlteration") == 0:
            tp53_wt_m34.append(float(m34))
    add_benchmark("B7_TP53_Inactivated_vs_WT_miR34Family", "TP53_Mutated_or_DeepDel", "TP53_WildType_Diploid", "mir34_family_mean_log2_rpm", tp53_alt_m34, tp53_wt_m34)

    ago2_amp_glob, ago2_del_glob = [], []
    for r in mirna_table_rows:
        sid = r["sampleId"]
        gmed = r["global_canonical_median_log2_rpm"]
        gcna = sl_lookup.get((sid, "AGO2"), {}).get("gistic2DiscreteCNA", 0)
        if gcna == 2:
            ago2_amp_glob.append(float(gmed))
        elif gcna in (-1, -2):
            ago2_del_glob.append(float(gmed))
    add_benchmark("B8_AGO2_HighAmp_p2_vs_Loss_GlobalMatureMiRNA", "AGO2_GISTIC_Amp_p2", "AGO2_GISTIC_Del_m1_m2", "global_canonical_median_log2_rpm", ago2_amp_glob, ago2_del_glob)

    if risc_derep_map:
        r3b_l7t, wt_l7t = [], []
        r3b_m200t, wt_m200t = [], []
        lin28_hi_l7t, lin28_lo_l7t = [], []
        m200_q_pairs = []
        for r in mirna_table_rows:
            sid = r["sampleId"]
            if sid not in risc_derep_map:
                continue
            l7_tz = float(risc_derep_map[sid]["let7_target_oncofetal_mrna_z"])
            m200_tz = float(risc_derep_map[sid]["mir200_target_emt_mrna_z"])
            cls = r["dicer1FunctionalClass"]
            if cls == "DICER1_RNaseIIIb_3D_Catalytic_Pocket":
                r3b_l7t.append(l7_tz)
                r3b_m200t.append(m200_tz)
            elif cls == "DICER1_WildType":
                wt_l7t.append(l7_tz)
                wt_m200t.append(m200_tz)

            za = sl_lookup.get((sid, "LIN28A"), {}).get("mrnaZscore", "")
            zb = sl_lookup.get((sid, "LIN28B"), {}).get("mrnaZscore", "")
            z_vals = [float(z) for z in (za, zb) if z != ""]
            if z_vals:
                mz = max(z_vals)
                if mz >= 2.0:
                    lin28_hi_l7t.append(l7_tz)
                elif mz <= 0.0:
                    lin28_lo_l7t.append(l7_tz)

            if r["mir200_family_mean_log2_rpm"] != "":
                m200_q_pairs.append((float(r["mir200_family_mean_log2_rpm"]), m200_tz))

        add_benchmark("B9_DICER1_RNaseIIIb_vs_WT_Let7_5pTargetMRNA_DeRepression", "DICER1_RNaseIIIb_3D_Catalytic_Pocket", "DICER1_WildType", "let7_target_oncofetal_mrna_z", r3b_l7t, wt_l7t)
        add_benchmark("B10_DICER1_RNaseIIIb_vs_WT_miR200_3pTargetMRNA_SparedControl", "DICER1_RNaseIIIb_3D_Catalytic_Pocket", "DICER1_WildType", "mir200_target_emt_mrna_z", r3b_m200t, wt_m200t)
        add_benchmark("B11_LIN28AB_Reactivated_vs_Silent_Let7TargetMRNA_DeRepression", "LIN28A_or_LIN28B_mRNA_Z_ge_2", "LIN28AB_Silent_Z_le_0", "let7_target_oncofetal_mrna_z", lin28_hi_l7t, lin28_lo_l7t)
        if len(m200_q_pairs) >= 100:
            m200_q_pairs.sort(key=lambda t: t[0])
            q = len(m200_q_pairs) // 4
            add_benchmark("B12_miR200_LowQ1_vs_HighQ4_EMTTargetMRNA_DeRepression", "miR200_BottomQuartile_Q1", "miR200_TopQuartile_Q4", "mir200_target_emt_mrna_z", [t[1] for t in m200_q_pairs[:q]], [t[1] for t in m200_q_pairs[-q:]])

    write_csv_fn(
        data_dir / "validation_mirnaseq_enzymatic_benchmarks.csv",
        list(bench_rows[0].keys()),
        bench_rows,
    )

    ent_to_sym = {g["entrezGeneId"]: g["geneSymbol"] for g in gene_panel}
    cptac_prot_mrna_pairs = {g: ([], []) for g in gene_order}
    cptac_cna_mrna_pairs = {g: ([], []) for g in gene_order}
    cptac_cna_prot_pairs = {g: ([], []) for g in gene_order}
    cptac_complex_pairs = {
        "DROSHA_vs_DGCR8": ([], []),
        "DICER1_vs_TARBP2": ([], []),
        "DICER1_vs_PRKRA": ([], []),
        "SMAD2_vs_SMAD4": ([], []),
    }
    cptac_cohort_rows = []

    def zscore_dict(val_map):
        if len(val_map) < 5:
            return {}
        vs = list(val_map.values())
        mu = statistics.mean(vs)
        sd = statistics.stdev(vs)
        if sd <= 1e-8:
            return {}
        return {k: (v - mu) / sd for k, v in val_map.items()}

    for st_entry in cptac_snap["studies"]:
        sid = st_entry["study"]["studyId"]
        n_samp = len(st_entry["sample_ids"])

        prot_raw = defaultdict(dict)
        for r in st_entry["protein_data"]:
            sym = ent_to_sym.get(r.get("entrezGeneId"))
            if sym and r.get("value") is not None:
                try:
                    v = float(r["value"])
                    if not math.isnan(v):
                        prot_raw[sym][r["sampleId"]] = v
                except Exception:
                    pass

        mrna_raw = defaultdict(dict)
        for r in st_entry["mrna_data"]:
            sym = ent_to_sym.get(r.get("entrezGeneId"))
            if sym and r.get("value") is not None:
                try:
                    v = float(r["value"])
                    if not math.isnan(v):
                        mrna_raw[sym][r["sampleId"]] = math.log2(v + 1.0) if v > 50 else v
                except Exception:
                    pass

        cna_raw = defaultdict(dict)
        cna_src = st_entry["cna_linear_data"] if st_entry["cna_linear_data"] else st_entry["cna_discrete_data"]
        for r in cna_src:
            sym = ent_to_sym.get(r.get("entrezGeneId"))
            if sym and r.get("value") is not None:
                try:
                    v = float(r["value"])
                    if not math.isnan(v):
                        cna_raw[sym][r["sampleId"]] = v
                except Exception:
                    pass

        prot_z = {sym: zscore_dict(prot_raw[sym]) for sym in gene_order}
        mrna_z = {sym: zscore_dict(mrna_raw[sym]) for sym in gene_order}
        cna_z = {sym: zscore_dict(cna_raw[sym]) for sym in gene_order}

        for sym in gene_order:
            for s in sorted(set(prot_z[sym].keys()) & set(mrna_z[sym].keys())):
                cptac_prot_mrna_pairs[sym][0].append(mrna_z[sym][s])
                cptac_prot_mrna_pairs[sym][1].append(prot_z[sym][s])
            for s in sorted(set(cna_z[sym].keys()) & set(mrna_z[sym].keys())):
                cptac_cna_mrna_pairs[sym][0].append(cna_z[sym][s])
                cptac_cna_mrna_pairs[sym][1].append(mrna_z[sym][s])
            for s in sorted(set(cna_z[sym].keys()) & set(prot_z[sym].keys())):
                cptac_cna_prot_pairs[sym][0].append(cna_z[sym][s])
                cptac_cna_prot_pairs[sym][1].append(prot_z[sym][s])

        study_drosha_dgcr8_s = sorted(set(prot_z["DROSHA"].keys()) & set(prot_z["DGCR8"].keys()))
        rho_mp = spearman_rho_fn(
            [prot_z["DROSHA"][s] for s in study_drosha_dgcr8_s],
            [prot_z["DGCR8"][s] for s in study_drosha_dgcr8_s],
        ) if len(study_drosha_dgcr8_s) >= 10 else float("nan")

        study_dicer_tarbp2_s = sorted(set(prot_z["DICER1"].keys()) & set(prot_z["TARBP2"].keys()))
        rho_rlc = spearman_rho_fn(
            [prot_z["DICER1"][s] for s in study_dicer_tarbp2_s],
            [prot_z["TARBP2"][s] for s in study_dicer_tarbp2_s],
        ) if len(study_dicer_tarbp2_s) >= 10 else float("nan")

        for label, (g1, g2) in [
            ("DROSHA_vs_DGCR8", ("DROSHA", "DGCR8")),
            ("DICER1_vs_TARBP2", ("DICER1", "TARBP2")),
            ("DICER1_vs_PRKRA", ("DICER1", "PRKRA")),
            ("SMAD2_vs_SMAD4", ("SMAD2", "SMAD4")),
        ]:
            for s in sorted(set(prot_z[g1].keys()) & set(prot_z[g2].keys())):
                cptac_complex_pairs[label][0].append(prot_z[g1][s])
                cptac_complex_pairs[label][1].append(prot_z[g2][s])

        cptac_cohort_rows.append({
            "cptacStudyId": sid,
            "studyName": st_entry["study"].get("name", sid),
            "n_total_samples": n_samp,
            "n_drosha_dgcr8_protein_paired": len(study_drosha_dgcr8_s),
            "spearman_drosha_vs_dgcr8_protein": round(rho_mp, 4) if not math.isnan(rho_mp) else "",
            "n_dicer1_tarbp2_protein_paired": len(study_dicer_tarbp2_s),
            "spearman_dicer1_vs_tarbp2_protein": round(rho_rlc, 4) if not math.isnan(rho_rlc) else "",
        })

    cptac_gene_rows = []
    for g in gene_panel:
        sym = g["geneSymbol"]
        cm_x, cm_y = cptac_cna_mrna_pairs[sym]
        cp_x, cp_y = cptac_cna_prot_pairs[sym]
        mp_x, mp_y = cptac_prot_mrna_pairs[sym]
        rho_cm = spearman_rho_fn(cm_x, cm_y)
        rho_cp = spearman_rho_fn(cp_x, cp_y)
        rho_mp = spearman_rho_fn(mp_x, mp_y)
        attenuation = round(rho_cm - rho_cp, 4) if (not math.isnan(rho_cm) and not math.isnan(rho_cp)) else ""
        cptac_gene_rows.append({
            "geneSymbol": sym,
            "panelRole": g["panelRole"],
            "n_cptac_mrna_protein_paired": len(mp_x),
            "n_cptac_cna_protein_paired": len(cp_x),
            "spearman_CNA_to_mRNA": round(rho_cm, 4) if not math.isnan(rho_cm) else "",
            "spearman_CNA_to_Protein": round(rho_cp, 4) if not math.isnan(rho_cp) else "",
            "spearman_mRNA_to_Protein": round(rho_mp, 4) if not math.isnan(rho_mp) else "",
            "post_transcriptional_CNA_attenuation": attenuation,
        })

    write_csv_fn(
        data_dir / "cptac_pancan_miRNA_biogenesis_proteogenomic_validation.csv",
        list(cptac_gene_rows[0].keys()),
        cptac_gene_rows,
    )
    write_csv_fn(
        data_dir / "cptac_cohort_complex_stoichiometry.csv",
        list(cptac_cohort_rows[0].keys()),
        cptac_cohort_rows,
    )

    cl_map = depmap_snap["cell_lines"]
    depmap_gene_rows = []
    for g in gene_panel:
        sym = g["geneSymbol"]
        effs = []
        expr_pairs = ([], [])
        for cl in cl_map.values():
            if sym in cl["geneEffect"]:
                ev = cl["geneEffect"][sym]
                effs.append(ev)
                if sym in cl.get("expression", {}):
                    expr_pairs[0].append(cl["expression"][sym])
                    expr_pairs[1].append(ev)
        n_eff = len(effs)
        n_dep_05 = sum(1 for v in effs if v <= -0.5)
        n_dep_10 = sum(1 for v in effs if v <= -1.0)
        rho_ex = spearman_rho_fn(expr_pairs[0], expr_pairs[1])
        depmap_gene_rows.append({
            "geneSymbol": sym,
            "panelRole": g["panelRole"],
            "opentargetsIsCommonEssential": int(bool(ess_flags.get(sym, False))),
            "n_depmap_crispr_cell_lines": n_eff,
            "median_chronos_gene_effect": round(statistics.median(effs), 4) if effs else "",
            "mean_chronos_gene_effect": round(statistics.mean(effs), 4) if effs else "",
            "iqr_chronos_gene_effect": round(statistics.quantiles(effs, n=4)[2] - statistics.quantiles(effs, n=4)[0], 4) if len(effs) >= 10 else "",
            "n_dependent_le_m05": n_dep_05,
            "prop_dependent_le_m05": round(n_dep_05 / n_eff, 4) if n_eff else 0.0,
            "n_strongly_essential_le_m10": n_dep_10,
            "prop_strongly_essential_le_m10": round(n_dep_10 / n_eff, 4) if n_eff else 0.0,
            "spearman_expression_vs_chronos": round(rho_ex, 4) if not math.isnan(rho_ex) else "",
        })

    write_csv_fn(
        data_dir / "depmap_crispr_miRNA_biogenesis_essentiality.csv",
        list(depmap_gene_rows[0].keys()),
        depmap_gene_rows,
    )

    collat_rows = []
    def add_collateral_test(tid, desc, target_gene, grp_del_label, grp_wt_label, del_vals, wt_vals):
        if len(del_vals) < 5 or len(wt_vals) < 5:
            return
        r_rb, z_val, p_val = mann_whitney_z_and_p(del_vals, wt_vals)
        collat_rows.append({
            "test_id": tid,
            "hypothesis": desc,
            "crispr_target_gene": target_gene,
            "altered_group": grp_del_label,
            "n_altered_lines": len(del_vals),
            "median_chronos_altered": round(statistics.median(del_vals), 4),
            "mean_chronos_altered": round(statistics.mean(del_vals), 4),
            "control_group": grp_wt_label,
            "n_control_lines": len(wt_vals),
            "median_chronos_control": round(statistics.median(wt_vals), 4),
            "mean_chronos_control": round(statistics.mean(wt_vals), 4),
            "chronos_shift_altered_minus_control": round(statistics.median(del_vals) - statistics.median(wt_vals), 4),
            "rank_biserial_r": r_rb,
            "mann_whitney_z": z_val,
            "p_value_two_sided": f"{p_val:.3e}",
        })

    ago2_in_1p_del, ago2_in_1p_intact = [], []
    for cl in cl_map.values():
        gmap = cl.get("gistic")
        if not gmap:
            continue
        if "AGO2" in cl["geneEffect"]:
            ago134 = [gmap.get(k, 0) for k in ("AGO1", "AGO3", "AGO4")]
            if all(v <= -1 for v in ago134):
                ago2_in_1p_del.append(cl["geneEffect"]["AGO2"])
            elif all(v >= 0 for v in ago134):
                ago2_in_1p_intact.append(cl["geneEffect"]["AGO2"])

    add_collateral_test(
        "C1_1p34_AGO1_3_4_Del_to_AGO2_Dependency",
        "Paralogous collateral vulnerability: 1p34.3 AGO1/3/4 co-deletion heightens dependency on AGO2",
        "AGO2",
        "1p34.3_AGO1_3_4_CoDeleted_GISTIC_le_m1",
        "1p34.3_AGO1_3_4_CopyNeutral_or_Gain",
        ago2_in_1p_del,
        ago2_in_1p_intact,
    )

    for sym in ("DROSHA", "DGCR8", "XPO5", "DICER1", "AGO1", "AGO2"):
        loss_v, neut_v = [], []
        for cl in cl_map.values():
            gmap = cl.get("gistic")
            if not gmap or sym not in cl["geneEffect"] or sym not in gmap:
                continue
            if gmap[sym] <= -1:
                loss_v.append(cl["geneEffect"][sym])
            elif gmap[sym] == 0:
                neut_v.append(cl["geneEffect"][sym])
        add_collateral_test(
            f"C2_CYCLOPS_HemizygousLoss_{sym}",
            f"Copy-loss vulnerability (CYCLOPS): {sym} hemizygous deletion vs copy-neutral diploid lines",
            sym,
            f"{sym}_GISTIC_Del_le_m1",
            f"{sym}_GISTIC_Diploid_0",
            loss_v,
            neut_v,
        )

    ago134_scored = []
    for cl in cl_map.values():
        ex = cl.get("expression", {})
        ge = cl.get("geneEffect", {})
        if all(k in ex for k in ("AGO1", "AGO3", "AGO4")) and "AGO2" in ge:
            ago134_scored.append((ex["AGO1"] + ex["AGO3"] + ex["AGO4"], ge["AGO2"]))
    ago134_scored.sort(key=lambda t: t[0])
    if len(ago134_scored) >= 40:
        q = len(ago134_scored) // 4
        add_collateral_test(
            "C3_1p34_Cluster_AGO1_3_4_LowExpr_to_AGO2_Dependency",
            "Paralogous RISC buffering: low combined 1p34.3 AGO1+AGO3+AGO4 expression (Q1 vs Q4) sensitizes cells to AGO2 knockout",
            "AGO2",
            "1p34_AGO1_3_4_CombinedExpr_Q1_Low",
            "1p34_AGO1_3_4_CombinedExpr_Q4_High",
            [x[1] for x in ago134_scored[:q]],
            [x[1] for x in ago134_scored[-q:]],
        )

    l28b_on = [cl["geneEffect"]["LIN28B"] for cl in cl_map.values() if cl.get("expression", {}).get("LIN28B", -1) >= 2.0 and "LIN28B" in cl["geneEffect"]]
    l28b_off = [cl["geneEffect"]["LIN28B"] for cl in cl_map.values() if 0 <= cl.get("expression", {}).get("LIN28B", -1) < 0.2 and "LIN28B" in cl["geneEffect"]]
    add_collateral_test(
        "C4_LIN28B_Reactivated_Selective_Dependency",
        "Context-specific oncogene addiction: LIN28B-reactivated lines (log2TPM >= 2.0) depend selectively on LIN28B vs silent lines",
        "LIN28B",
        "LIN28B_Reactivated_log2TPM_ge_2",
        "LIN28B_Silent_log2TPM_lt_0p2",
        l28b_on,
        l28b_off,
    )

    tp53_wt = [cl["geneEffect"]["TP53"] for cl in cl_map.values() if not cl.get("mutation", {}).get("TP53") and "TP53" in cl["geneEffect"]]
    tp53_mut = [cl["geneEffect"]["TP53"] for cl in cl_map.values() if cl.get("mutation", {}).get("TP53") and "TP53" in cl["geneEffect"]]
    add_collateral_test(
        "C5_TP53_WT_vs_Mutant_CRISPR_Positive_Control",
        "Orthogonal screen QC benchmark: TP53 knockout confers antiproliferative release selectively in TP53-WT vs TP53-mutant lines",
        "TP53",
        "TP53_WildType_CellLines",
        "TP53_Mutated_CellLines",
        tp53_wt,
        tp53_mut,
    )

    write_csv_fn(
        data_dir / "validation_depmap_collateral_vulnerability.csv",
        list(collat_rows[0].keys()),
        collat_rows,
    )

    print(f"   Matched TCGA mature miRNA-seq tumours: {len(mirna_table_rows)} / {len(patient_burden_rows)}")
    print(f"   Enzymatic benchmarks: {len(bench_rows)}; CPTAC cohorts: {len(cptac_cohort_rows)}; DepMap tests: {len(collat_rows)}")

    from clinical_and_external import integrate_clinical_and_external
    integrate_clinical_and_external(
        repo_root=repo_root,
        data_dir=data_dir,
        gene_panel=gene_panel,
        gene_order=gene_order,
        panel_rows=panel_rows,
        annotated_mut_rows=annotated_mut_rows,
        patient_burden_rows=patient_burden_rows,
        sample_level_rows=sample_level_rows,
        spearman_rho_fn=spearman_rho_fn,
        mann_whitney_fn=mann_whitney_z_and_p,
        write_csv_fn=write_csv_fn,
    )
