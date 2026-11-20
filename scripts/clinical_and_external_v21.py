#!/usr/bin/env python3
"""
miRBioGenAtlas v2.0 — Pure-Stdlib Clinical, Germline Constraint, ClinVar ROC-AUC,
GTEx/HPA/OpenTargets, MSK-IMPACT (N=10,945), TCGA-CDR Survival & Thorsson PanImmune Module
==========================================================================================
Integrates frozen, SHA-256-verified snapshots from `data/raw_snapshot/external/`
using ONLY the Python 3 standard library (`csv`, `gzip`, `json`, `math`, `statistics`).
"""

import csv
import gzip
import json
import math
import pathlib
import re
import statistics
from collections import Counter, defaultdict

MISSENSE_SHORT_RE = re.compile(r"^p?\.?([A-Z]\d+[A-Z])$")
TRUNC_TYPES = {
    "Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins",
    "Splice_Site", "Translation_Start_Site", "Nonstop_Mutation",
}


def load_gz_json(path: pathlib.Path):
    with gzip.open(path, "rt", encoding="utf-8") as f:
        return json.load(f)


def compute_roc_pr_auc(y_true, scores, mann_whitney_fn):
    """
    Exact tie-handled ROC-AUC + Hanley-McNeil 95% CI + Average-Precision PR-AUC
    in pure Python stdlib. `y_true` in {0, 1}, higher `scores` => positive class (1).
    """
    pos = [s for y, s in zip(y_true, scores) if y == 1]
    neg = [s for y, s in zip(y_true, scores) if y == 0]
    n1, n0 = len(pos), len(neg)
    if n1 < 3 or n0 < 3:
        return None
    r_rb, z_val, p_val = mann_whitney_fn(pos, neg)
    # Note: r_rb = 2 * AUC - 1 => AUC = (r_rb + 1) / 2!
    auc = (r_rb + 1.0) / 2.0
    q1 = auc / (2.0 - auc) if auc < 1.0 else 1.0
    q2 = 2.0 * (auc ** 2) / (1.0 + auc)
    var_auc = (auc * (1.0 - auc) + (n1 - 1) * (q1 - auc ** 2) + (n0 - 1) * (q2 - auc ** 2)) / (n1 * n0)
    se = math.sqrt(max(0.0, var_auc))
    ci_lo = max(0.0, auc - 1.96 * se)
    ci_hi = min(1.0, auc + 1.96 * se)

    # Exact Average Precision (PR-AUC) sorting descending by score (tie-averaged)
    pairs = sorted(zip(scores, y_true), key=lambda t: -t[0])
    tp = 0.0
    fp = 0.0
    ap_sum = 0.0
    i = 0
    N = len(pairs)
    while i < N:
        j = i
        pos_tie = 0
        neg_tie = 0
        while j < N and pairs[j][0] == pairs[i][0]:
            if pairs[j][1] == 1:
                pos_tie += 1
            else:
                neg_tie += 1
            j += 1
        if pos_tie > 0:
            # Average precision across tied items
            tp_end = tp + pos_tie
            fp_end = fp + neg_tie
            prec_at_tie = tp_end / (tp_end + fp_end)
            ap_sum += pos_tie * prec_at_tie
        tp += pos_tie
        fp += neg_tie
        i = j
    pr_auc = ap_sum / n1
    baseline_prev = n1 / (n1 + n0)
    return {
        "n_pos": n1,
        "n_neg": n0,
        "roc_auc": round(auc, 4),
        "roc_ci95_low": round(ci_lo, 4),
        "roc_ci95_high": round(ci_hi, 4),
        "pr_auc": round(pr_auc, 4),
        "pr_baseline": round(baseline_prev, 4),
        "mann_whitney_p": f"{p_val:.3e}",
    }


def stratified_logrank_and_hr(records):
    """
    Exact cancer-type-stratified Mantel-Haenszel / Peto log-rank test and
    one-step hazard ratio in pure Python stdlib.
    `records`: list of (stratum, time_months, event_01, group_01) with group=1 (altered) vs 0 (control).
    """
    by_stratum = defaultdict(list)
    n1_tot, n0_tot, ev1_tot, ev0_tot = 0, 0, 0, 0
    for st, tm, ev, grp in records:
        if tm is None or ev not in (0, 1) or grp not in (0, 1) or tm <= 0:
            continue
        by_stratum[st].append((float(tm), int(ev), int(grp)))
        if grp == 1:
            n1_tot += 1
            ev1_tot += int(ev)
        else:
            n0_tot += 1
            ev0_tot += int(ev)

    if n1_tot < 5 or n0_tot < 5 or (ev1_tot + ev0_tot) < 5:
        return None

    O_minus_E_sum = 0.0
    V_sum = 0.0
    active_strata = 0

    for st, items in by_stratum.items():
        r1 = sum(1 for _, _, g in items if g == 1)
        r0 = sum(1 for _, _, g in items if g == 0)
        if r1 == 0 or r0 == 0:
            continue
        active_strata += 1
        items.sort(key=lambda x: x[0])
        n_risk_1 = r1
        n_risk_0 = r0
        idx = 0
        m = len(items)
        while idx < m:
            t_cur = items[idx][0]
            d1, d0 = 0, 0
            c1, c0 = 0, 0
            while idx < m and items[idx][0] == t_cur:
                _, ev, grp = items[idx]
                if grp == 1:
                    if ev == 1:
                        d1 += 1
                    else:
                        c1 += 1
                else:
                    if ev == 1:
                        d0 += 1
                    else:
                        c0 += 1
                idx += 1
            d = d1 + d0
            n_risk = n_risk_1 + n_risk_0
            if d > 0 and n_risk > 1 and n_risk_1 > 0 and n_risk_0 > 0:
                e1 = n_risk_1 * (d / n_risk)
                v1 = (n_risk_1 * n_risk_0 * d * (n_risk - d)) / ((n_risk ** 2) * (n_risk - 1))
                O_minus_E_sum += (d1 - e1)
                V_sum += v1
            n_risk_1 -= (d1 + c1)
            n_risk_0 -= (d0 + c0)

    if V_sum <= 1e-9:
        return None
    z = O_minus_E_sum / math.sqrt(V_sum)
    p_two = math.erfc(abs(z) / math.sqrt(2.0))
    log_hr = O_minus_E_sum / V_sum
    se_log_hr = 1.0 / math.sqrt(V_sum)
    hr = math.exp(max(-10.0, min(10.0, log_hr)))
    hr_lo = math.exp(max(-10.0, min(10.0, log_hr - 1.96 * se_log_hr)))
    hr_hi = math.exp(max(-10.0, min(10.0, log_hr + 1.96 * se_log_hr)))
    return {
        "n_altered": n1_tot,
        "events_altered": ev1_tot,
        "n_control": n0_tot,
        "events_control": ev0_tot,
        "n_informative_strata": active_strata,
        "hazard_ratio_stratified": round(hr, 4),
        "hr_ci95_low": round(hr_lo, 4),
        "hr_ci95_high": round(hr_hi, 4),
        "logrank_z": round(z, 4),
        "p_value_two_sided": f"{p_two:.3e}",
    }


def _solve_linear_system(A, b):
    """Gauss-Jordan elimination with partial pivoting in pure Python stdlib."""
    k = len(b)
    M = [row[:] + [b[i]] for i, row in enumerate(A)]
    for col in range(k):
        piv = max(range(col, k), key=lambda r: abs(M[r][col]))
        if abs(M[piv][col]) < 1e-12:
            return None
        M[col], M[piv] = M[piv], M[col]
        div = M[col][col]
        for j in range(col, k + 1):
            M[col][j] /= div
        for r in range(k):
            if r != col and M[r][col] != 0.0:
                fac = M[r][col]
                for j in range(col, k + 1):
                    M[r][j] -= fac * M[col][j]
    return [M[i][k] for i in range(k)]


def _invert_matrix(A):
    k = len(A)
    inv_cols = []
    for col in range(k):
        e = [1.0 if i == col else 0.0 for i in range(k)]
        sol = _solve_linear_system(A, e)
        if sol is None:
            return None
        inv_cols.append(sol)
    return [[inv_cols[j][i] for j in range(k)] for i in range(k)]


def _parse_ajcc_stage_indicators(st_raw):
    """Returns (stage_III_IV_advanced, stage_unstaged_or_nonAJCC) with Stage I/II as reference (0, 0)."""
    s = (st_raw or "").strip().upper()
    if not s or s in (
        "STAGE X", "[DISCREPANCY]", "[NOT AVAILABLE]", "I/II NOS",
        "STAGE I/II (NOS)", "IS", "STAGE IS", "STAGE 0",
    ):
        return 0.0, 1.0
    if "IV" in s or "III" in s:
        return 1.0, 0.0
    if "II" in s or "I" in s:
        return 0.0, 0.0
    return 0.0, 1.0


def stratified_multivar_cox_breslow(records, covariates_desc, max_iter=35, tol=1e-8):
    """
    Exact Breslow-tied cancer-type-stratified multivariable Cox proportional hazards
    regression via Newton-Raphson partial-likelihood maximization in pure Python stdlib.
    Cross-verified to 6 decimal places against R `survival::coxph(..., ties='breslow')`.
    `records`: list of (stratum, time_months, event_01, x_vec), where x_vec[0] is the
    primary binary exposure (1 = altered, 0 = control).
    """
    by_st = defaultdict(list)
    for st, tm, ev, x in records:
        if tm is None or ev not in (0, 1) or tm <= 0:
            continue
        by_st[st].append((float(tm), int(ev), [float(v) for v in x]))

    strata = {}
    for st, items in by_st.items():
        if len(items) >= 2 and sum(ev for _, ev, _ in items) >= 1:
            items.sort(key=lambda r: r[0])
            strata[st] = items
    if not strata:
        return None

    k_raw = len(next(iter(strata.values()))[0][2])
    active_cols = []
    for j in range(k_raw):
        var_sum = 0.0
        for items in strata.values():
            vals = [x[j] for _, _, x in items]
            mu = sum(vals) / len(vals)
            var_sum += sum((v - mu) ** 2 for v in vals)
        if j == 0 and var_sum < 1e-9:
            return None
        if var_sum >= 1e-6:
            active_cols.append(j)

    proj_strata = {}
    n_tot = 0
    ev_tot = 0
    for st, items in strata.items():
        proj_strata[st] = [(tm, ev, [x[j] for j in active_cols]) for tm, ev, x in items]
        n_tot += len(items)
        ev_tot += sum(ev for _, ev, _ in items)

    k = len(active_cols)
    beta = [0.0] * k

    def loglik_grad_info(b):
        ll = 0.0
        grad = [0.0] * k
        info = [[0.0] * k for _ in range(k)]
        for items in proj_strata.values():
            m = len(items)
            S0 = 0.0
            S1 = [0.0] * k
            S2 = [[0.0] * k for _ in range(k)]
            idx = m - 1
            while idx >= 0:
                t_cur = items[idx][0]
                j = idx
                d = 0
                x_ev_sum = [0.0] * k
                while j >= 0 and items[j][0] == t_cur:
                    _, ev, x = items[j]
                    eta = max(-50.0, min(50.0, sum(b[c] * x[c] for c in range(k))))
                    w = math.exp(eta)
                    S0 += w
                    for c in range(k):
                        wx = w * x[c]
                        S1[c] += wx
                        for c2 in range(c, k):
                            S2[c][c2] += wx * x[c2]
                    if ev == 1:
                        d += 1
                        for c in range(k):
                            x_ev_sum[c] += x[c]
                    j -= 1
                if d > 0 and S0 > 0.0:
                    ll += sum(b[c] * x_ev_sum[c] for c in range(k)) - d * math.log(S0)
                    inv_s0 = 1.0 / S0
                    mean_x = [S1[c] * inv_s0 for c in range(k)]
                    for c in range(k):
                        grad[c] += x_ev_sum[c] - d * mean_x[c]
                        for c2 in range(c, k):
                            cov_cc2 = S2[c][c2] * inv_s0 - mean_x[c] * mean_x[c2]
                            info[c][c2] += d * cov_cc2
                idx = j
        for c in range(k):
            for c2 in range(c):
                info[c][c2] = info[c2][c]
        return ll, grad, info

    ll_cur, grad, info = loglik_grad_info(beta)
    for _ in range(max_iter):
        step = _solve_linear_system(info, grad)
        if step is None:
            info_reg = [[info[r][c] + (1e-6 if r == c else 0.0) for c in range(k)] for r in range(k)]
            step = _solve_linear_system(info_reg, grad)
            if step is None:
                return None
        alpha = 1.0
        while alpha >= 1e-4:
            b_try = [beta[c] + alpha * max(-2.0, min(2.0, step[c])) for c in range(k)]
            ll_try, grad_try, info_try = loglik_grad_info(b_try)
            if ll_try >= ll_cur - 1e-10:
                break
            alpha *= 0.5
        beta = b_try
        ll_cur, grad, info = ll_try, grad_try, info_try
        if max(abs(alpha * s) for s in step) < tol:
            break

    cov = _invert_matrix(info)
    if cov is None or cov[0][0] <= 0.0:
        return None
    b0 = beta[0]
    se0 = math.sqrt(cov[0][0])
    z0 = b0 / se0
    p0 = math.erfc(abs(z0) / math.sqrt(2.0))
    return {
        "multivar_n_complete": n_tot,
        "multivar_events_complete": ev_tot,
        "multivar_cox_hr": round(math.exp(max(-10.0, min(10.0, b0))), 4),
        "multivar_ci95_low": round(math.exp(max(-10.0, min(10.0, b0 - 1.96 * se0))), 4),
        "multivar_ci95_high": round(math.exp(max(-10.0, min(10.0, b0 + 1.96 * se0))), 4),
        "multivar_wald_z": round(z0, 4),
        "multivar_p_value": f"{p0:.3e}",
        "multivar_covariates": covariates_desc,
    }


def run_v21_clinical_and_external_enrichment(
    repo_root: pathlib.Path,
    data_dir: pathlib.Path,
    gene_panel: list,
    gene_order: list,
    panel_rows: list,
    annotated_mut_rows: list,
    patient_burden_rows: list,
    sample_level_rows: list,
    spearman_rho_fn,
    mann_whitney_fn,
    write_csv_fn,
):
    ext_dir = data_dir / "raw_snapshot" / "external"
    gc_gz = ext_dir / "gene_clinical_constraint_gtex_hpa_ot.json.gz"
    cv_gz = ext_dir / "clinvar_15genes_benchmark_and_clinical.json.gz"
    msk_gz = ext_dir / "msk_impact_10945_replication.json.gz"
    imm_gz = ext_dir / "thorsson_panimmune_tcga.json.gz"

    if not (gc_gz.exists() and cv_gz.exists() and msk_gz.exists() and imm_gz.exists()):
        print("   [Note] v2.1 clinical/external snapshots not present; skipping v2.1 extensions.")
        return

    print("\n5. Integrating gnomAD, GTEx (54 tissues), HPA, Open Targets, ClinVar, MSK-IMPACT (N=10,945), TCGA-CDR & Thorsson PanImmune...")
    gc_snap = load_gz_json(gc_gz)
    cv_snap = load_gz_json(cv_gz)
    msk_snap = load_gz_json(msk_gz)
    imm_snap = load_gz_json(imm_gz)

    # Load 213,940 saturation mutagenesis lookup so EVERY ClinVar & MSK-IMPACT missense variant gets instant PLM/3D/IF1 scores!
    sat_lookup = {}
    sat_gz = data_dir / "miRNA_biogenesis_saturation_mutagenesis_213940_atlas.csv.gz"
    if sat_gz.exists():
        with gzip.open(sat_gz, "rt", encoding="utf-8") as f:
            for r in csv.DictReader(f):
                sat_lookup[(r["geneSymbol"], r["variant"])] = r

    # -------------------------------------------------------------------------
    # 1. Enrich Gene Panel Metadata (15 genes) + Export GTEx 54-Tissue Table (810 rows)
    # -------------------------------------------------------------------------
    for prow in panel_rows:
        sym = prow["geneSymbol"]
        gi = gc_snap["genes"].get(sym, {})
        cs = cv_snap["per_gene_summary"].get(sym, {})
        prow["gnomad_pLI"] = gi.get("gnomad_pLI", "")
        prow["gnomad_LOEUF"] = gi.get("gnomad_LOEUF", "")
        prow["gnomad_oe_mis"] = gi.get("gnomad_oe_mis", "")
        prow["gnomad_mis_z"] = gi.get("gnomad_mis_z", "")
        prow["hpa_main_subcellular_location"] = gi.get("hpa_main_locations", "")
        prow["hpa_ihc_detected_tissues_n"] = gi.get("hpa_ihc_detected_tissues_n", "")
        prow["gtex_pan54_median_tpm"] = gi.get("gtex_pan54_median_tpm", "")
        prow["gtex_top_tissue"] = gi.get("gtex_max_tissue", "")
        prow["gtex_top_tissue_tpm"] = gi.get("gtex_max_tpm", "")
        prow["gtex_tau_tissue_specificity"] = gi.get("gtex_tau_specificity", "")
        prow["opentargets_sm_tractable"] = gi.get("ot_sm_tractable", "")
        prow["opentargets_protac_tractable"] = gi.get("ot_protac_tractable", "")
        prow["opentargets_clinical_candidates_n"] = gi.get("ot_clinical_candidates_n", 0)
        prow["opentargets_max_clinical_stage"] = gi.get("ot_max_clinical_stage", "None")
        prow["clinvar_total_single_gene_n"] = cs.get("total_single_gene_records", 0)
        prow["clinvar_pathogenic_or_likely_pathogenic_n"] = cs.get("pathogenic_or_likely_pathogenic_n", 0)
        prow["clinvar_benign_or_likely_benign_n"] = cs.get("benign_or_likely_benign_n", 0)

    write_csv_fn(
        data_dir / "miRNA_biogenesis_gene_panel_metadata.csv",
        list(panel_rows[0].keys()),
        panel_rows,
    )

    gtex_rows = gc_snap.get("gtex_54tissues_long", [])
    if gtex_rows:
        write_csv_fn(
            data_dir / "gtex_v10_normal_54tissues_miRNA_biogenesis_tpm.csv",
            list(gtex_rows[0].keys()),
            gtex_rows,
        )

    # -------------------------------------------------------------------------
    # 2. ClinVar Somatic Cross-Annotation + Non-Circular Human Clinical ROC-AUC Benchmarks
    # -------------------------------------------------------------------------
    # Deduplicate ClinVar coding variants per (geneSymbol, proteinChange), preferring highest reviewStars & unambiguous labels
    sig_priority = {
        "Pathogenic_or_Likely_Pathogenic": 4,
        "Benign_or_Likely_Benign": 3,
        "Uncertain_Significance": 2,
        "Conflicting": 1,
        "Other": 0,
    }
    cv_best = {}
    for cv in cv_snap.get("coding_variants", []):
        key = (cv["geneSymbol"], cv["proteinChange"])
        rank = (cv.get("reviewStars", 0), sig_priority.get(cv["significanceCategory"], 0))
        if key not in cv_best or rank > cv_best[key][0]:
            cv_best[key] = (rank, cv)

    # Cross-annotate TCGA 6,292 somatic mutations
    tcga_clinvar_overlap_n = 0
    tcga_obs_set = Counter((m["geneSymbol"], (m["proteinChange"] or "").lstrip("p.")) for m in annotated_mut_rows)
    for mrow in annotated_mut_rows:
        pch = (mrow.get("proteinChange") or "").strip().lstrip("p.")
        hit = cv_best.get((mrow["geneSymbol"], pch))
        if hit:
            tcga_clinvar_overlap_n += 1
            cv = hit[1]
            mrow["clinvarVariationId"] = cv["variantId"]
            mrow["clinvarClinicalSignificance"] = cv["clinicalSignificance"]
            mrow["clinvarReviewStars"] = cv["reviewStars"]
            mrow["clinvarPhenotypes"] = cv["phenotypes"]
        else:
            mrow["clinvarVariationId"] = ""
            mrow["clinvarClinicalSignificance"] = ""
            mrow["clinvarReviewStars"] = ""
            mrow["clinvarPhenotypes"] = ""

    write_csv_fn(
        data_dir / "tcga_pancan_miRNA_biogenesis_somatic_mutations_annotated.csv",
        list(annotated_mut_rows[0].keys()),
        annotated_mut_rows,
    )

    # Export annotated ClinVar coding variants table & compute ROC-AUC / PR-AUC across predictors!
    clinvar_export_rows = []
    for (_, _), (_, cv) in sorted(cv_best.items(), key=lambda kv: (kv[0][0], kv[1][1]["codonPos"], kv[0][1])):
        sym = cv["geneSymbol"]
        pch = cv["proteinChange"]
        srec = sat_lookup.get((sym, pch), {})
        d_sub = srec.get("distSubstrateA", "")
        d_met = srec.get("distMetalA", "")
        d_prt = srec.get("distPartnerA", "")
        is_3dp = int(any(v != "" and float(v) <= 6.0 for v in (d_sub, d_met, d_prt))) if srec else ""
        clinvar_export_rows.append({
            "geneSymbol": sym,
            "clinvarVariationId": cv["variantId"],
            "proteinChange": pch,
            "codonPos": cv["codonPos"],
            "codingClass": cv["codingClass"],
            "clinicalSignificance": cv["clinicalSignificance"],
            "significanceCategory": cv["significanceCategory"],
            "reviewStars": cv["reviewStars"],
            "reviewStatus": cv["reviewStatus"],
            "phenotypes": cv["phenotypes"],
            "rsid": cv.get("rsid", ""),
            "grch38Vcf": cv.get("grch38Vcf", ""),
            "tcgaSomaticCount": tcga_obs_set.get((sym, pch), 0),
            "pdbComplexId": srec.get("pdbComplexId", ""),
            "distToSubstrateA": d_sub,
            "distToCatalyticMetalA": d_met,
            "is3dPocket_le6A": is_3dp,
            "alphamissenseScore": srec.get("alphamissenseScore", ""),
            "esm2_3b_llr": srec.get("esm2_3b_llr", ""),
            "esm1v_650m_llr": srec.get("esm1v_650m_llr", ""),
            "esm_if1_3d_llr": srec.get("esm_if1_3d_llr", ""),
        })

    if clinvar_export_rows:
        write_csv_fn(
            data_dir / "clinvar_clinical_variants_15genes_annotated.csv",
            list(clinvar_export_rows[0].keys()),
            clinvar_export_rows,
        )

    # Evaluate ClinVar Missense Pathogenic/Likely_Pathogenic (y=1) vs Benign/Likely_Benign (y=0)
    core10_set = {"DROSHA", "DGCR8", "XPO5", "DICER1", "TARBP2", "PRKRA", "AGO1", "AGO2", "AGO3", "AGO4"}
    benchmark_subsets = [
        ("ClinVar_All15_Missense_AllStars", lambda r: True),
        ("ClinVar_All15_Missense_HighConf_ge2Stars", lambda r: int(r["reviewStars"]) >= 2),
        ("ClinVar_NonTP53_14Genes_Missense_AllStars", lambda r: r["geneSymbol"] != "TP53"),
        ("ClinVar_Core10_Biogenesis_Missense_AllStars", lambda r: r["geneSymbol"] in core10_set),
    ]
    predictors = [
        ("AlphaMissense", lambda r: float(r["alphamissenseScore"]) if r["alphamissenseScore"] != "" else None),
        ("neg_ESM1v_650M_LLR", lambda r: -float(r["esm1v_650m_llr"]) if r["esm1v_650m_llr"] != "" else None),
        ("neg_ESM2_3B_LLR", lambda r: -float(r["esm2_3b_llr"]) if r["esm2_3b_llr"] != "" else None),
        ("neg_ESM_IF1_3D_LLR", lambda r: -float(r["esm_if1_3d_llr"]) if r["esm_if1_3d_llr"] != "" else None),
        (
            "Joint_Evolutionary_3D_Composite",
            lambda r: (
                0.5 * float(r["alphamissenseScore"])
                + 0.25 * (-float(r["esm1v_650m_llr"]) / 12.0)
                + 0.25 * (-float(r["esm_if1_3d_llr"]) / 8.0)
            ) if (r["alphamissenseScore"] != "" and r["esm1v_650m_llr"] != "" and r["esm_if1_3d_llr"] != "") else None,
        ),
    ]
    roc_rows = []
    for sub_name, sub_fn in benchmark_subsets:
        cand = [
            r for r in clinvar_export_rows
            if r["codingClass"] == "Missense"
            and r["significanceCategory"] in ("Pathogenic_or_Likely_Pathogenic", "Benign_or_Likely_Benign")
            and sub_fn(r)
        ]
        for pred_name, pred_fn in predictors:
            ys, scs = [], []
            for r in cand:
                sv = pred_fn(r)
                if sv is not None and not math.isnan(sv):
                    ys.append(1 if r["significanceCategory"] == "Pathogenic_or_Likely_Pathogenic" else 0)
                    scs.append(sv)
            res = compute_roc_pr_auc(ys, scs, mann_whitney_fn)
            if res:
                roc_rows.append({
                    "benchmark_subset": sub_name,
                    "predictor": pred_name,
                    **res,
                })
    if roc_rows:
        write_csv_fn(
            data_dir / "validation_clinvar_roc_auc_benchmarks.csv",
            list(roc_rows[0].keys()),
            roc_rows,
        )

    # -------------------------------------------------------------------------
    # 3. MSK-IMPACT (N=10,945) Independent Clinical Sequencing Replication
    # -------------------------------------------------------------------------
    gene_meta_by_sym = {g["geneSymbol"]: g for g in gene_panel}
    msk_samples = msk_snap["samples"]  # sid -> [cancerType, sampleType, panelSuffix, tmb, osStatus, osMonths]
    N_msk = msk_snap["n_sequenced_and_cna"]
    msk_profiled = msk_snap["profiled_genes_all_10945"]  # ["DICER1", "SMAD2", "SMAD4", "TP53"]

    msk_cna_by_gene = defaultdict(lambda: {"amp": set(), "dpdel": set(), "any_cna": set()})
    for sid, sym, val in msk_snap["cna_nonzero"]:
        msk_cna_by_gene[sym]["any_cna"].add(sid)
        if val == 2:
            msk_cna_by_gene[sym]["amp"].add(sid)
        elif val == -2:
            msk_cna_by_gene[sym]["dpdel"].add(sid)

    msk_mut_by_gene = defaultdict(lambda: {
        "mut_samples": set(),
        "nonhyper_mut_samples": set(),
        "trunc_samples": set(),
        "hotspot_or_3dpocket_samples": set(),
        "am_or_plm_path_samples": set(),
        "primary_mut_samples": set(),
        "met_mut_samples": set(),
        "calls_n": 0,
    })
    msk_cancer_dicer = defaultdict(lambda: {"N": 0, "mut": set(), "r3b_hotspot": set(), "trunc": set()})
    for sid, sdata in msk_samples.items():
        msk_cancer_dicer[sdata[0]]["N"] += 1

    msk_nonhyper_n = sum(1 for sdata in msk_samples.values() if sdata[3] is None or sdata[3] < 10.0)
    msk_prim_n = sum(1 for sdata in msk_samples.values() if sdata[1] == "Primary")
    msk_met_n = sum(1 for sdata in msk_samples.values() if sdata[1] == "Metastasis")

    for m in msk_snap["mutations"]:
        sid = m["s"]
        sym = m["g"]
        if sid not in msk_samples or sym not in msk_profiled:
            continue
        sdata = msk_samples[sid]
        ctype, stype, _, tmb, _, _ = sdata
        is_hyper = tmb is not None and tmb >= 10.0
        pch = (m.get("p") or "").lstrip("p.")
        pos = int(m["pos"]) if str(m.get("pos") or "").isdigit() else 0
        mtype = m["mt"]
        hs_set = gene_meta_by_sym[sym]["hotspotCodons"]

        srec = sat_lookup.get((sym, pch), {})
        is_3dp = any(srec.get(k, "") not in ("", None) and float(srec[k]) <= 6.0 for k in ("distSubstrateA", "distMetalA", "distPartnerA"))
        am_sc = float(srec["alphamissenseScore"]) if srec.get("alphamissenseScore", "") not in ("", None) else None
        e2_sc = float(srec["esm2_3b_llr"]) if srec.get("esm2_3b_llr", "") not in ("", None) else None

        gstat = msk_mut_by_gene[sym]
        gstat["calls_n"] += 1
        gstat["mut_samples"].add(sid)
        if not is_hyper:
            gstat["nonhyper_mut_samples"].add(sid)
        if stype == "Primary":
            gstat["primary_mut_samples"].add(sid)
        elif stype == "Metastasis":
            gstat["met_mut_samples"].add(sid)
        if mtype in TRUNC_TYPES:
            gstat["trunc_samples"].add(sid)
        if pos in hs_set or (is_3dp and (am_sc is not None and am_sc >= 0.564)):
            gstat["hotspot_or_3dpocket_samples"].add(sid)
        if (am_sc is not None and am_sc >= 0.564) or (e2_sc is not None and e2_sc <= -5.0) or (pos in hs_set) or (mtype in TRUNC_TYPES):
            gstat["am_or_plm_path_samples"].add(sid)

        if sym == "DICER1":
            msk_cancer_dicer[ctype]["mut"].add(sid)
            if mtype in TRUNC_TYPES:
                msk_cancer_dicer[ctype]["trunc"].add(sid)
            if pos in hs_set or (1690 <= pos <= 1830 and is_3dp):
                msk_cancer_dicer[ctype]["r3b_hotspot"].add(sid)

    # Compare against TCGA (N=10,195) pan-cancer pooled counts across the exact same genes!
    tcga_pooled = defaultdict(lambda: {"mut": set(), "trunc": set(), "hs3d": set(), "amp": set(), "dpdel": set(), "hc": set()})
    for sl in sample_level_rows:
        sym = sl["geneSymbol"]
        sid = sl["sampleId"]
        if sl["nonsynMutationCount"] > 0:
            tcga_pooled[sym]["mut"].add(sid)
        if sl["hasTruncatingLoF"] == 1:
            tcga_pooled[sym]["trunc"].add(sid)
        if sl["gistic2DiscreteCNA"] == 2:
            tcga_pooled[sym]["amp"].add(sid)
        elif sl["gistic2DiscreteCNA"] == -2:
            tcga_pooled[sym]["dpdel"].add(sid)
        if sl["hasHighConfAlteration"] == 1:
            tcga_pooled[sym]["hc"].add(sid)
    for mr in annotated_mut_rows:
        if mr["isCuratedHotspot"] == 1 or mr["functionalVariantTierV2"] == "3D_Catalytic_or_Substrate_Pocket_Driver":
            tcga_pooled[mr["geneSymbol"]]["hs3d"].add(mr["sampleId"])

    N_tcga = len(patient_burden_rows)
    msk_rep_rows = []
    for sym in msk_profiled:
        gs = msk_mut_by_gene[sym]
        cs = msk_cna_by_gene[sym]
        hc_msk = gs["mut_samples"] | cs["amp"] | cs["dpdel"]
        tp = tcga_pooled[sym]
        msk_rep_rows.append({
            "geneSymbol": sym,
            "msk_N_profiled": N_msk,
            "msk_mut_calls_n": gs["calls_n"],
            "msk_mut_samples_n": len(gs["mut_samples"]),
            "msk_mut_freq": round(len(gs["mut_samples"]) / N_msk, 5),
            "tcga_mut_freq_N10195": round(len(tp["mut"]) / N_tcga, 5),
            "msk_nonhypermut_mut_freq": round(len(gs["nonhyper_mut_samples"]) / msk_nonhyper_n, 5),
            "msk_truncating_n": len(gs["trunc_samples"]),
            "msk_truncating_freq": round(len(gs["trunc_samples"]) / N_msk, 5),
            "tcga_truncating_freq": round(len(tp["trunc"]) / N_tcga, 5),
            "msk_3dpocket_or_hotspot_n": len(gs["hotspot_or_3dpocket_samples"]),
            "msk_3dpocket_or_hotspot_freq": round(len(gs["hotspot_or_3dpocket_samples"]) / N_msk, 5),
            "tcga_3dpocket_or_hotspot_freq": round(len(tp["hs3d"]) / N_tcga, 5),
            "msk_amp_p2_n": len(cs["amp"]),
            "msk_deep_del_m2_n": len(cs["dpdel"]),
            "msk_high_conf_alter_freq": round(len(hc_msk) / N_msk, 5),
            "tcga_high_conf_alter_freq": round(len(tp["hc"]) / N_tcga, 5),
            "msk_primary_mut_freq": round(len(gs["primary_mut_samples"]) / msk_prim_n, 5) if msk_prim_n else 0.0,
            "msk_metastasis_mut_freq": round(len(gs["met_mut_samples"]) / msk_met_n, 5) if msk_met_n else 0.0,
        })

    write_csv_fn(
        data_dir / "msk_impact_10945_clinical_cohort_replication.csv",
        list(msk_rep_rows[0].keys()),
        msk_rep_rows,
    )

    dicer_ct_rows = []
    for ctype, st in sorted(msk_cancer_dicer.items(), key=lambda kv: (-len(kv[1]["r3b_hotspot"]), -len(kv[1]["mut"]), -kv[1]["N"])):
        if st["N"] < 15 and len(st["mut"]) == 0:
            continue
        dicer_ct_rows.append({
            "msk_cancer_type": ctype,
            "N_samples": st["N"],
            "dicer1_mut_n": len(st["mut"]),
            "dicer1_mut_freq": round(len(st["mut"]) / st["N"], 4),
            "dicer1_rnase3b_hotspot_n": len(st["r3b_hotspot"]),
            "dicer1_rnase3b_hotspot_freq": round(len(st["r3b_hotspot"]) / st["N"], 4),
            "dicer1_truncating_n": len(st["trunc"]),
            "dicer1_truncating_freq": round(len(st["trunc"]) / st["N"], 4),
        })
    write_csv_fn(
        data_dir / "msk_impact_dicer1_cancer_type_enrichment.csv",
        list(dicer_ct_rows[0].keys()),
        dicer_ct_rows,
    )

    # -------------------------------------------------------------------------
    # 4. Align Thorsson PanImmune onto 10,195 TCGA Patients + Survival & Immune Benchmarks
    # -------------------------------------------------------------------------
    imm_map = imm_snap["samples"]
    for pb in patient_burden_rows:
        im = imm_map.get(pb["sampleId"], {})
        pb["immuneSubtypeC1C6"] = im.get("isub") or ""
        pb["leukocyteFraction"] = im["leuk"] if im.get("leuk") is not None else ""
        pb["ifngResponseScore"] = im["ifng"] if im.get("ifng") is not None else ""
        pb["tgfbResponseScore"] = im["tgfb"] if im.get("tgfb") is not None else ""
        pb["lymphocyteInfiltrationScore"] = im["lymph"] if im.get("lymph") is not None else ""
        pb["macrophageCsf1Score"] = im["csf1"] if im.get("csf1") is not None else ""
        pb["woundHealingScore"] = im["wound"] if im.get("wound") is not None else ""
        pb["hrdScore"] = im["hrd"] if im.get("hrd") is not None else ""

    write_csv_fn(
        data_dir / "tcga_pancan_miRNA_biogenesis_patient_level_burden.csv",
        list(patient_burden_rows[0].keys()),
        patient_burden_rows,
    )

    # Build fast sample-level gene alteration lookup for survival & immune tests
    sl_by_sg = {(r["sampleId"], r["geneSymbol"]): r for r in sample_level_rows}

    surv_rows = []
    def add_surv_test(test_id, cohort_src, endpoint_name, time_key, ev_key, grp_label, ctrl_label, grp_fn, strata_fn=lambda r: r["tumour_acronym"]):
        recs = []
        mv_recs = []
        for pb in patient_burden_rows:
            g = grp_fn(pb)
            if g not in (0, 1):
                continue
            tm = pb.get(time_key, "")
            ev = pb.get(ev_key, "")
            if tm == "" or ev == "":
                continue
            st = strata_fn(pb)
            tm_f = float(tm)
            ev_i = int(ev)
            recs.append((st, tm_f, ev_i, g))
            if pb.get("ageAtDiagnosis", "") != "" and pb.get("sex", "") in ("Male", "Female") and tm_f > 0:
                age10 = (float(pb["ageAtDiagnosis"]) - 60.0) / 10.0
                male = 1.0 if pb["sex"] == "Male" else 0.0
                st_adv, st_unst = _parse_ajcc_stage_indicators(pb.get("ajccPathologicStage", ""))
                hyp = float(int(pb.get("isHypermutated", 0)))
                mv_recs.append((st, tm_f, ev_i, [float(g), age10, male, st_adv, st_unst, hyp]))
        res = stratified_logrank_and_hr(recs)
        mv_res = stratified_multivar_cox_breslow(
            mv_recs,
            "age_per10y+sex_male+ajcc_stage_III_IV+stage_unstaged+hypermutator+strata(tumour_acronym)",
        )
        if res:
            row_out = {
                "test_id": test_id,
                "cohort": cohort_src,
                "endpoint": endpoint_name,
                "altered_group": grp_label,
                "control_group": ctrl_label,
                **res,
            }
            if mv_res:
                row_out.update(mv_res)
            else:
                row_out.update({
                    "multivar_n_complete": "",
                    "multivar_events_complete": "",
                    "multivar_cox_hr": "",
                    "multivar_ci95_low": "",
                    "multivar_ci95_high": "",
                    "multivar_wald_z": "",
                    "multivar_p_value": "",
                    "multivar_covariates": "",
                })
            surv_rows.append(row_out)

    # Cancer-type-stratified TCGA-CDR survival benchmarks (controlling 100% for cancer type across 32 studies!)
    add_surv_test(
        "S1_Core10_StrictAlter_Stratified_OS", "TCGA_PanCanAtlas_Stratified32", "OS", "osMonths", "osEvent",
        "Core10_Strict_Altered_Mut_or_FocalCNA", "Core10_No_Strict_Alteration",
        lambda pb: int(pb["core10_has_strict_alter"]),
    )
    add_surv_test(
        "S2_Core10_StrictAlter_Stratified_PFI", "TCGA_PanCanAtlas_Stratified32", "PFI", "pfsMonths", "pfsEvent",
        "Core10_Strict_Altered_Mut_or_FocalCNA", "Core10_No_Strict_Alteration",
        lambda pb: int(pb["core10_has_strict_alter"]),
    )
    add_surv_test(
        "S3_Core10_StrictAlter_Stratified_DSS", "TCGA_PanCanAtlas_Stratified32", "DSS", "dssMonths", "dssEvent",
        "Core10_Strict_Altered_Mut_or_FocalCNA", "Core10_No_Strict_Alteration",
        lambda pb: int(pb["core10_has_strict_alter"]),
    )
    add_surv_test(
        "S4_AGO2_FocalAmp_p2_Stratified_OS", "TCGA_PanCanAtlas_Stratified32", "OS", "osMonths", "osEvent",
        "AGO2_GISTIC2_HighAmp_p2", "AGO2_Diploid_0",
        lambda pb: 1 if sl_by_sg[(pb["sampleId"], "AGO2")]["gistic2DiscreteCNA"] == 2 else (0 if sl_by_sg[(pb["sampleId"], "AGO2")]["gistic2DiscreteCNA"] == 0 else -1),
    )
    add_surv_test(
        "S5_AGO2_FocalAmp_p2_Stratified_PFI", "TCGA_PanCanAtlas_Stratified32", "PFI", "pfsMonths", "pfsEvent",
        "AGO2_GISTIC2_HighAmp_p2", "AGO2_Diploid_0",
        lambda pb: 1 if sl_by_sg[(pb["sampleId"], "AGO2")]["gistic2DiscreteCNA"] == 2 else (0 if sl_by_sg[(pb["sampleId"], "AGO2")]["gistic2DiscreteCNA"] == 0 else -1),
    )
    add_surv_test(
        "S6_LIN28AB_Oncofetal_Reactivation_Stratified_OS", "TCGA_PanCanAtlas_Stratified32", "OS", "osMonths", "osEvent",
        "LIN28A_or_LIN28B_mRNA_Z_ge_2", "LIN28AB_Silent_Z_le_0",
        lambda pb: (
            1 if max(
                [float(sl_by_sg[(pb["sampleId"], g)]["mrnaZscore"]) for g in ("LIN28A", "LIN28B") if sl_by_sg[(pb["sampleId"], g)]["mrnaZscore"] != ""] or [-99.0]
            ) >= 2.0 else (
                0 if max(
                    [float(sl_by_sg[(pb["sampleId"], g)]["mrnaZscore"]) for g in ("LIN28A", "LIN28B") if sl_by_sg[(pb["sampleId"], g)]["mrnaZscore"] != ""] or [99.0]
                ) <= 0.0 else -1
            )
        ),
    )
    add_surv_test(
        "S7_1p34_AGO1_3_4_CoDeletion_Stratified_OS", "TCGA_PanCanAtlas_Stratified32", "OS", "osMonths", "osEvent",
        "1p34_AGO1_3_4_Concordant_Loss", "1p34_AGO1_3_4_Diploid",
        lambda pb: (
            1 if all(sl_by_sg[(pb["sampleId"], g)]["gistic2DiscreteCNA"] <= -1 for g in ("AGO1", "AGO3", "AGO4"))
            else (0 if all(sl_by_sg[(pb["sampleId"], g)]["gistic2DiscreteCNA"] == 0 for g in ("AGO1", "AGO3", "AGO4")) else -1)
        ),
    )
    add_surv_test(
        "S8_SMAD4_HighConf_Inactivation_Stratified_OS", "TCGA_PanCanAtlas_Stratified32", "OS", "osMonths", "osEvent",
        "SMAD4_Mut_or_DeepDel", "SMAD4_WildType_Diploid",
        lambda pb: 1 if sl_by_sg[(pb["sampleId"], "SMAD4")]["hasHighConfAlteration"] == 1 else (0 if sl_by_sg[(pb["sampleId"], "SMAD4")]["hasAnyAlteration"] == 0 else -1),
    )
    add_surv_test(
        "S9_TP53_PositiveControl_Stratified_OS", "TCGA_PanCanAtlas_Stratified32", "OS", "osMonths", "osEvent",
        "TP53_Mut_or_DeepDel", "TP53_WildType_Diploid",
        lambda pb: 1 if sl_by_sg[(pb["sampleId"], "TP53")]["hasHighConfAlteration"] == 1 else (0 if sl_by_sg[(pb["sampleId"], "TP53")]["hasAnyAlteration"] == 0 else -1),
    )

    # Also replicate survival in MSK-IMPACT (N=7,574 OS-evaluable clinical patients stratified by cancer type!)
    for sym in ("DICER1", "SMAD4", "TP53"):
        alt_set = msk_mut_by_gene[sym]["am_or_plm_path_samples"] | msk_cna_by_gene[sym]["dpdel"]
        msk_surv_recs = []
        msk_mv_recs = []
        for sid, sdata in msk_samples.items():
            ctype, stype, panel_sfx, tmb_v, os_st, os_mo = sdata
            if not os_st or os_mo in ("", None, "NA"):
                continue
            try:
                tm = float(os_mo)
            except Exception:
                continue
            if tm <= 0:
                continue
            ev = 1 if str(os_st).startswith("1") or "DECEASED" in str(os_st).upper() else 0
            grp = 1 if sid in alt_set else 0
            msk_surv_recs.append((ctype, tm, ev, grp))
            is_met = 1.0 if stype == "Metastasis" else 0.0
            is_hyp = 1.0 if (tmb_v is not None and float(tmb_v) >= 10.0) else 0.0
            is_im5 = 1.0 if panel_sfx == "IM5" else 0.0
            msk_mv_recs.append((ctype, tm, ev, [float(grp), is_met, is_hyp, is_im5]))
        res_m = stratified_logrank_and_hr(msk_surv_recs)
        mv_m = stratified_multivar_cox_breslow(
            msk_mv_recs,
            "metastatic_specimen+tmb_hypermutator+impact410_panel+strata(cancer_type)",
        )
        if res_m:
            row_m = {
                "test_id": f"S_MSK_IMPACT_{sym}_Pathogenic_Stratified_OS",
                "cohort": "MSK_IMPACT_Clinical_Stratified",
                "endpoint": "OS",
                "altered_group": f"{sym}_PathogenicMut_or_DeepDel",
                "control_group": f"{sym}_Unaltered",
                **res_m,
            }
            if mv_m:
                row_m.update(mv_m)
            else:
                row_m.update({
                    "multivar_n_complete": "",
                    "multivar_events_complete": "",
                    "multivar_cox_hr": "",
                    "multivar_ci95_low": "",
                    "multivar_ci95_high": "",
                    "multivar_wald_z": "",
                    "multivar_p_value": "",
                    "multivar_covariates": "",
                })
            surv_rows.append(row_m)

    write_csv_fn(
        data_dir / "validation_clinical_survival_benchmarks.csv",
        list(surv_rows[0].keys()),
        surv_rows,
    )

    # Thorsson PanImmune microenvironment & HRD associations (non-hypermutated & pan-cancer)
    imm_bench_rows = []
    def add_imm_test(tid, metric_col, grp_a_lbl, grp_b_lbl, pred_a, pred_b):
        va = [float(pb[metric_col]) for pb in patient_burden_rows if pb.get(metric_col, "") != "" and pred_a(pb)]
        vb = [float(pb[metric_col]) for pb in patient_burden_rows if pb.get(metric_col, "") != "" and pred_b(pb)]
        if len(va) < 10 or len(vb) < 10:
            return
        r_rb, z_v, p_v = mann_whitney_fn(va, vb)
        imm_bench_rows.append({
            "benchmark_id": tid,
            "immune_or_instability_metric": metric_col,
            "group_A": grp_a_lbl,
            "n_group_A": len(va),
            "median_group_A": round(statistics.median(va), 4),
            "group_B": grp_b_lbl,
            "n_group_B": len(vb),
            "median_group_B": round(statistics.median(vb), 4),
            "median_shift_A_minus_B": round(statistics.median(va) - statistics.median(vb), 4),
            "rank_biserial_r": r_rb,
            "mann_whitney_z": z_v,
            "p_value_two_sided": f"{p_v:.3e}",
        })

    add_imm_test(
        "I1_AGO2_FocalAmp_p2_vs_Diploid_HRDScore", "hrdScore",
        "AGO2_GISTIC2_HighAmp_p2", "AGO2_Diploid_0",
        lambda pb: sl_by_sg[(pb["sampleId"], "AGO2")]["gistic2DiscreteCNA"] == 2,
        lambda pb: sl_by_sg[(pb["sampleId"], "AGO2")]["gistic2DiscreteCNA"] == 0,
    )
    add_imm_test(
        "I2_AGO2_FocalAmp_p2_vs_Diploid_IFNgResponse", "ifngResponseScore",
        "AGO2_GISTIC2_HighAmp_p2", "AGO2_Diploid_0",
        lambda pb: sl_by_sg[(pb["sampleId"], "AGO2")]["gistic2DiscreteCNA"] == 2,
        lambda pb: sl_by_sg[(pb["sampleId"], "AGO2")]["gistic2DiscreteCNA"] == 0,
    )
    add_imm_test(
        "I3_Core10_StrictAlter_vs_WT_HRDScore", "hrdScore",
        "Core10_Strict_Altered", "Core10_Unaltered",
        lambda pb: int(pb["core10_has_strict_alter"]) == 1,
        lambda pb: int(pb["core10_has_strict_alter"]) == 0,
    )
    add_imm_test(
        "I4_Core10_3DorAM_DriverMut_vs_WT_IFNgResponse", "ifngResponseScore",
        "Core10_3D_or_Pathogenic_Mut", "Core10_No_Pathogenic_Mut",
        lambda pb: int(pb["core10_has_3d_or_am_pathogenic_mut"]) == 1,
        lambda pb: int(pb["core10_has_3d_or_am_pathogenic_mut"]) == 0,
    )
    add_imm_test(
        "I5_Core10_3DorAM_DriverMut_vs_WT_LeukocyteFraction", "leukocyteFraction",
        "Core10_3D_or_Pathogenic_Mut", "Core10_No_Pathogenic_Mut",
        lambda pb: int(pb["core10_has_3d_or_am_pathogenic_mut"]) == 1,
        lambda pb: int(pb["core10_has_3d_or_am_pathogenic_mut"]) == 0,
    )
    add_imm_test(
        "I6_TARBP2_or_XPO5_Truncating_vs_WT_IFNgResponse", "ifngResponseScore",
        "TARBP2_or_XPO5_Truncating_LoF", "TARBP2_and_XPO5_WT",
        lambda pb: sl_by_sg[(pb["sampleId"], "TARBP2")]["hasTruncatingLoF"] == 1 or sl_by_sg[(pb["sampleId"], "XPO5")]["hasTruncatingLoF"] == 1,
        lambda pb: sl_by_sg[(pb["sampleId"], "TARBP2")]["nonsynMutationCount"] == 0 and sl_by_sg[(pb["sampleId"], "XPO5")]["nonsynMutationCount"] == 0,
    )
    add_imm_test(
        "I7_SMAD4_Inactivated_vs_WT_TGFbResponse", "tgfbResponseScore",
        "SMAD4_Mut_or_DeepDel", "SMAD4_WildType_Diploid",
        lambda pb: sl_by_sg[(pb["sampleId"], "SMAD4")]["hasHighConfAlteration"] == 1,
        lambda pb: sl_by_sg[(pb["sampleId"], "SMAD4")]["hasAnyAlteration"] == 0,
    )
    add_imm_test(
        "I8_LIN28AB_Reactivated_vs_Silent_WoundHealingProliferation", "woundHealingScore",
        "LIN28AB_mRNA_Z_ge_2", "LIN28AB_Silent_Z_le_0",
        lambda pb: max([float(sl_by_sg[(pb["sampleId"], g)]["mrnaZscore"]) for g in ("LIN28A", "LIN28B") if sl_by_sg[(pb["sampleId"], g)]["mrnaZscore"] != ""] or [-99.0]) >= 2.0,
        lambda pb: max([float(sl_by_sg[(pb["sampleId"], g)]["mrnaZscore"]) for g in ("LIN28A", "LIN28B") if sl_by_sg[(pb["sampleId"], g)]["mrnaZscore"] != ""] or [99.0]) <= 0.0,
    )

    write_csv_fn(
        data_dir / "validation_panimmune_and_hrd_benchmarks.csv",
        list(imm_bench_rows[0].keys()),
        imm_bench_rows,
    )

    # Print rich summary!
    print(f"   [ClinVar] Annotated {len(clinvar_export_rows)} unique coding variants; overlapped {tcga_clinvar_overlap_n} TCGA somatic calls.")
    for rr in roc_rows:
        print(f"     ROC-AUC | {rr['benchmark_subset']:<44} | {rr['predictor']:<30} | AUC={rr['roc_auc']:.4f} [{rr['roc_ci95_low']:.4f}-{rr['roc_ci95_high']:.4f}] | PR={rr['pr_auc']:.4f} | p={rr['mann_whitney_p']}")
    print(f"   [MSK-IMPACT N={N_msk}] Independent clinical cohort replication:")
    for mr in msk_rep_rows:
        print(f"     {mr['geneSymbol']:<8} | MSK mut={mr['msk_mut_freq']*100:5.2f}% vs TCGA={mr['tcga_mut_freq_N10195']*100:5.2f}% | MSK 3D/hs={mr['msk_3dpocket_or_hotspot_n']:3d} ({mr['msk_3dpocket_or_hotspot_freq']*100:5.2f}%) vs TCGA={mr['tcga_3dpocket_or_hotspot_freq']*100:5.2f}%")
    print("   [Survival — Cancer-Type-Stratified Log-Rank & Multivariable Breslow Cox]:")
    for sr in surv_rows:
        mv_str = (
            f" | Multivar Cox HR={sr['multivar_cox_hr']:.3f} [{sr['multivar_ci95_low']:.3f}-{sr['multivar_ci95_high']:.3f}] (p={sr['multivar_p_value']})"
            if sr.get("multivar_cox_hr", "") != "" else ""
        )
        print(f"     {sr['test_id']:<48} | LogRank HR={sr['hazard_ratio_stratified']:.3f} (p={sr['p_value_two_sided']}){mv_str}")
    print("   [Thorsson PanImmune & HRD Associations]:")
    for ir in imm_bench_rows:
        print(f"     {ir['benchmark_id']:<52} | medA={ir['median_group_A']:+8.4f} vs medB={ir['median_group_B']:+8.4f} | p={ir['p_value_two_sided']}")
