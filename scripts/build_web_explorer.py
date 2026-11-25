#!/usr/bin/env python3
"""Compile CSV/CSV.GZ tables, 3D atomic coordinates, and KaTeX LaTeX equations into a clean minimal docs/index.html."""

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
    SECTION_GUIDES,
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
    snap_gz = DATA_DIR / "raw_snapshot" / "external" / "afdb_v6_15genes_heavy_atoms.json.gz"
    if not snap_gz.exists():
        return {}

    STRUCT3D_DIR.mkdir(parents=True, exist_ok=True)
    with gzip.open(snap_gz, "rt", encoding="utf-8") as f:
        snap = json.load(f)

    candidates_by_aa = defaultdict(list)
    total_kb = 0.0

    for gene, gd in sorted(snap.get("genes", {}).items()):
        out_obj = {"u": gd["uniprot"], "s": gd["sequence"], "a": gd["atoms"]}
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
        out_obj = {str(pos): {"m": per_gene_meta[(g, pos)], "s": cmap[pos]} for pos in sorted(cmap.keys())}
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
        pk, desc = DATASET_DESCRIPTIONS.get(fpath.name, ("Key columns", "Derived table in data/."))
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

    cds_json_path = DATA_DIR / "canonical_cds_genomic_15genes.json"
    cds_genomic = json.loads(cds_json_path.read_text(encoding="utf-8")).get("genes", {}) if cds_json_path.exists() else {}

    id_lookup = {}
    clinvar_by_gene = defaultdict(list)
    for r in clinvar_rows:
        g = r["geneSymbol"]
        pc = r["proteinChange"]
        if r.get("rsid"):
            id_lookup[r["rsid"].lower()] = f"{g}:{pc}"
        if r.get("clinvarVariationId"):
            id_lookup[f"vcv{int(r['clinvarVariationId']):09d}"] = f"{g}:{pc}"
            id_lookup[f"cv:{r['clinvarVariationId']}"] = f"{g}:{pc}"
        if r.get("significanceCategory") == "Pathogenic_or_Likely_Pathogenic":
            clinvar_by_gene[g].append({
                "clinvarVariationId": r["clinvarVariationId"],
                "proteinChange": pc,
                "codonPos": int(r["codonPos"]) if str(r.get("codonPos", "")).isdigit() else None,
                "clinicalSignificance": r["clinicalSignificance"],
                "reviewStars": int(r["reviewStars"]) if str(r.get("reviewStars", "")).isdigit() else 0,
                "rsid": r.get("rsid", ""),
                "grch38Vcf": r.get("grch38Vcf", ""),
                "phenotypes": r.get("phenotypes", ""),
            })

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
        "guides": SECTION_GUIDES,
        "glossary": TERM_GLOSSARY,
        "cds_genomic": cds_genomic,
        "id_lookup": id_lookup,
    }

    # Build static REST-style JSON API v2 under docs/api/v2/
    api_dir = DOCS_DIR / "api" / "v2"
    api_genes_dir = api_dir / "genes"
    api_genes_dir.mkdir(parents=True, exist_ok=True)

    (api_dir / "index.json").write_text(json.dumps({
        "resource": "miRBioGenAtlas",
        "version": "2.0.2",
        "doi": "10.5281/zenodo.24455067",
        "license_data": "CC-BY-4.0",
        "license_code": "MIT",
        "assemblies_supported": ["GRCh38", "GRCh37"],
        "counts": bundle["meta"],
        "endpoints": {
            "manifest": "api/v2/index.json",
            "panel": "api/v2/panel.json",
            "pancancer": "api/v2/pancancer.json",
            "gene_template": "api/v2/genes/{GENE}.json",
            "genes": [p["geneSymbol"] for p in panel],
        },
    }, indent=2), encoding="utf-8")

    (api_dir / "panel.json").write_text(json.dumps({
        "version": "2.0.2",
        "n_genes": len(panel),
        "genes": panel,
    }, indent=2), encoding="utf-8")

    (api_dir / "pancancer.json").write_text(json.dumps({
        "version": "2.0.2",
        "n_studies": len(study_summary),
        "n_gene_study_pairs": len(compact_long),
        "pathway_burden_ranking_excluding_TP53": top_tumours,
        "study_denominators": study_summary,
        "alteration_frequencies": compact_long,
    }, indent=2), encoding="utf-8")

    vars_by_gene = defaultdict(list)
    for v in compact_vars:
        vars_by_gene[v["g"]].append(v)

    for p in panel:
        g = p["geneSymbol"]
        g_cds = cds_genomic.get(g, {})
        gene_api_obj = {
            "geneSymbol": g,
            "version": "2.0.2",
            "metadata": p,
            "genomic_coordinates": {
                "ensembl_gene_id": g_cds.get("ensg", p["ensemblGeneId"]),
                "ensembl_canonical_transcript": g_cds.get("enst", ""),
                "refseq_canonical": g_cds.get("refseq", p["canonicalRefSeq"]),
                "chromosome": g_cds.get("chr", ""),
                "strand": g_cds.get("strand", 1),
                "cds_exons_grch38": g_cds.get("cds_exons_grch38", []),
                "cds_exons_grch37": g_cds.get("cds_exons_grch37", []),
            },
            "tcga_somatic_variants": vars_by_gene.get(g, []),
            "clinvar_pathogenic_or_likely_pathogenic": clinvar_by_gene.get(g, []),
            "residue_track_columns": ["codonPos", "wtAA", "afdbPlddt", "is3dPocket_le6A", "mean_alphamissense", "mean_esm2_3b_llr", "mean_esm_if1_3d_llr", "observedTcgaCodonMutCount", "clinvarPathogenicCount", "pdbComplexId"],
            "residue_track": residue_tracks.get(g, []),
        }
        (api_genes_dir / f"{g}.json").write_text(json.dumps(gene_api_obj, separators=(",", ":")), encoding="utf-8")
    print("   Exported 18 static REST JSON API endpoints in docs/api/v2/")

    data_json = json.dumps(bundle, separators=(",", ":"))

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>miRBioGenAtlas — MicroRNA Biogenesis Structural &amp; Genomic Reference</title>
<meta name="description" content="Open-access structural, genomic, and small-RNA reference of the human microRNA biogenesis machinery across 21,140 sequenced tumours and 213,940 3D missense variants.">
<meta property="og:title" content="miRBioGenAtlas (v2.0.2) — Structural &amp; Pan-Cancer miRNA Biogenesis Reference">
<meta property="og:description" content="Interactive 3D mutant rotamer viewer, 2D DMS saturation heatmap, GRCh38/GRCh37 VCF translator, 5p:3p small-RNA benchmarks, and static JSON API across 21,140 tumours.">
<meta property="og:type" content="website">
<meta property="og:url" content="https://siddugayam.github.io/miRBioGenAtlas/">
<meta name="twitter:card" content="summary_large_image">
<link rel="canonical" href="https://siddugayam.github.io/miRBioGenAtlas/">
<link rel="icon" type="image/svg+xml" href="data:image/svg+xml,%3Csvg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 32 32'%3E%3Crect width='32' height='32' rx='8' fill='%23111315'/%3E%3Ccircle cx='11.5' cy='16' r='5' fill='%23009f4d'/%3E%3Ccircle cx='20.5' cy='16' r='5' fill='%2300709b'/%3E%3C/svg%3E">
<script type="application/ld+json">
{{
  "@context": "https://schema.org",
  "@type": "Dataset",
  "name": "miRBioGenAtlas: Pan-cancer alteration and 3D saturation mutagenesis atlas of the human miRNA biogenesis machinery",
  "version": "2.0.2",
  "identifier": "https://doi.org/10.5281/zenodo.24455067",
  "url": "https://siddugayam.github.io/miRBioGenAtlas/",
  "sameAs": "https://github.com/siddugayam/miRBioGenAtlas",
  "description": "Multi-consortium structural, genomic, and functional small-RNA reference covering 15 human miRNA biogenesis and regulatory proteins across 21,140 sequenced tumours (TCGA PanCancer Atlas and MSK-IMPACT), 213,940 single-amino-acid substitutions, 16 experimental PDB complexes, and 9,578 mature miRNA-seq transcriptomes.",
  "license": "https://creativecommons.org/licenses/by/4.0/",
  "isAccessibleForFree": true,
  "keywords": ["microRNA biogenesis", "DICER1", "DROSHA", "AGO2", "saturation mutagenesis", "AlphaMissense", "ESM-IF1", "TCGA PanCancer Atlas", "5p:3p strand ratio"],
  "creator": [
    {{"@type": "Person", "name": "Gayam Prasanna Kumar Reddy", "identifier": "https://orcid.org/0000-0001-5032-364X"}},
    {{"@type": "Person", "name": "Venkat Abhiram Earny", "identifier": "https://orcid.org/0009-0008-6926-553X"}},
    {{"@type": "Person", "name": "Jesil Mathew Aranjani", "identifier": "https://orcid.org/0000-0002-7666-3252"}},
    {{"@type": "Person", "name": "Fayaz Shaik Mahammad", "identifier": "https://orcid.org/0000-0002-1212-7952"}}
  ],
  "distribution": [
    {{"@type": "DataDownload", "encodingFormat": "application/zip", "contentUrl": "https://doi.org/10.5281/zenodo.24455067"}},
    {{"@type": "DataDownload", "encodingFormat": "application/json", "contentUrl": "https://siddugayam.github.io/miRBioGenAtlas/api/v2/index.json"}}
  ]
}}
</script>
<link rel="stylesheet" href="vendor/katex/katex.min.css">
<script src="vendor/katex/katex.min.js"></script>
<script src="vendor/3Dmol-min.js"></script>
<style>
  :root {{
    --bg: #ffffff;
    --canvas: #f8f9fa;
    --surface: #ffffff;
    --subtle: #f3f4f6;
    --ink: #202124;
    --ink-2: #5f6368;
    --ink-3: #80868b;
    --line: #dadce0;
    --line-soft: #ebeef2;
    --accent: #00709b;
    --accent-soft: rgba(0, 112, 155, 0.08);
    --crimson: #c5221f;
    --amber: #b06000;
    --emerald: #137333;
    --sans: -apple-system, BlinkMacSystemFont, "SF Pro Text", "Inter", "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
    --mono: "SF Mono", "JetBrains Mono", ui-monospace, Menlo, Consolas, monospace;
  }}
  * {{ box-sizing: border-box; }}
  html {{ -webkit-font-smoothing: antialiased; -moz-osx-font-smoothing: grayscale; }}
  body {{
    margin: 0;
    font-family: var(--sans);
    background: var(--bg);
    color: var(--ink);
    font-size: 14px;
    line-height: 1.52;
  }}
  .wrap {{ max-width: 1520px; margin: 0 auto; padding: 0 28px; }}

  /* Minimal Google / Apple Top Navigation Bar */
  .topbar {{
    position: sticky;
    top: 0;
    z-index: 40;
    background: rgba(255, 255, 255, 0.92);
    backdrop-filter: saturate(180%) blur(16px);
    -webkit-backdrop-filter: saturate(180%) blur(16px);
    border-bottom: 1px solid var(--line-soft);
    height: 60px;
  }}
  .topbar-inner {{
    height: 100%;
    display: flex;
    align-items: center;
    justify-content: space-between;
    gap: 16px;
  }}
  .brand {{
    display: flex;
    align-items: center;
    gap: 9px;
    font-size: 16.5px;
    font-weight: 600;
    letter-spacing: -0.02em;
    color: var(--ink);
    cursor: pointer;
    white-space: nowrap;
  }}
  .brand-dot {{
    width: 9px;
    height: 9px;
    border-radius: 50%;
    background: #009f4d;
    box-shadow: 11px 0 0 var(--accent);
    margin-right: 11px;
  }}
  .ver-pill {{
    font-family: var(--mono);
    font-size: 11px;
    font-weight: 500;
    color: var(--ink-2);
    background: var(--subtle);
    padding: 1px 7px;
    border-radius: 999px;
  }}
  .nav-tabs {{
    display: flex;
    align-items: center;
    gap: 5px;
    overflow-x: auto;
    scrollbar-width: none;
  }}
  .nav-tabs::-webkit-scrollbar {{ display: none; }}
  .nav-tab {{
    border: none;
    background: transparent;
    padding: 8px 14px;
    border-radius: 999px;
    font-family: var(--sans);
    font-size: 13.5px;
    font-weight: 500;
    color: var(--ink-2);
    cursor: pointer;
    white-space: nowrap;
    transition: background 0.15s ease, color 0.15s ease;
  }}
  .nav-tab:hover {{ background: var(--subtle); color: var(--ink); }}
  .nav-tab.active {{
    background: var(--accent-soft);
    color: var(--accent);
    font-weight: 600;
  }}
  .top-links {{
    display: flex;
    align-items: center;
    gap: 14px;
    font-size: 12.5px;
    white-space: nowrap;
  }}
  .top-links a {{
    color: var(--ink-2);
    text-decoration: none;
  }}
  .top-links a:hover {{ color: var(--accent); }}

  /* Google-Simple Centered Search Hero */
  .search-hero {{
    padding: 30px 0 22px;
    text-align: center;
    border-bottom: 1px solid var(--line-soft);
    background: var(--bg);
  }}
  .search-shell {{
    position: relative;
    max-width: 720px;
    margin: 0 auto;
  }}
  .search-box {{
    display: flex;
    align-items: center;
    height: 50px;
    padding: 0 16px 0 18px;
    background: #ffffff;
    border: 1px solid var(--line);
    border-radius: 999px;
    box-shadow: 0 1px 6px rgba(32, 33, 36, 0.08);
    transition: box-shadow 0.18s ease, border-color 0.18s ease;
  }}
  .search-box:focus-within {{
    border-color: rgba(0, 112, 155, 0.55);
    box-shadow: 0 2px 14px rgba(0, 112, 155, 0.15);
  }}
  .search-icon {{
    width: 18px;
    height: 18px;
    color: var(--ink-3);
    margin-right: 11px;
    flex-shrink: 0;
  }}
  .search-box input {{
    flex: 1;
    border: none;
    outline: none;
    font-size: 15px;
    font-family: var(--sans);
    color: var(--ink);
    background: transparent;
  }}
  .search-kbd {{
    font-family: var(--mono);
    font-size: 11px;
    color: var(--ink-3);
    border: 1px solid var(--line);
    border-radius: 5px;
    padding: 1px 6px;
    margin-left: 8px;
  }}
  .suggest-menu {{
    position: absolute;
    top: 56px;
    left: 0;
    right: 0;
    background: #ffffff;
    border: 1px solid var(--line);
    border-radius: 14px;
    box-shadow: 0 10px 28px rgba(32, 33, 36, 0.14);
    text-align: left;
    padding: 6px;
    z-index: 50;
    display: none;
  }}
  .suggest-menu.open {{ display: block; }}
  .sug-row {{
    display: flex;
    justify-content: space-between;
    align-items: center;
    padding: 9px 13px;
    border-radius: 8px;
    cursor: pointer;
    font-size: 13.5px;
  }}
  .sug-row:hover {{ background: var(--canvas); }}
  .sug-cat {{
    font-family: var(--mono);
    font-size: 10.5px;
    color: var(--ink-3);
  }}
  .quick-Row {{
    display: flex;
    justify-content: center;
    align-items: center;
    flex-wrap: wrap;
    gap: 8px;
    margin-top: 13px;
  }}
  .chip {{
    border: 1px solid var(--line-soft);
    background: var(--canvas);
    color: var(--ink-2);
    font-size: 12.5px;
    padding: 5px 12px;
    border-radius: 999px;
    cursor: pointer;
    transition: all 0.14s ease;
  }}
  .chip:hover {{
    border-color: var(--accent);
    color: var(--accent);
    background: #ffffff;
  }}
  .stat-ribbon {{
    margin-top: 13px;
    display: flex;
    justify-content: center;
    flex-wrap: wrap;
    gap: 20px;
    font-size: 12.5px;
    color: var(--ink-3);
  }}
  .stat-ribbon span[data-term] {{
    border-bottom: 1px dotted var(--line);
    padding-bottom: 1px;
  }}
  .stat-ribbon strong {{ color: var(--ink-2); font-weight: 600; }}

  /* Main Workspace & Clean Surface Cards */
  main.workspace {{
    background: var(--canvas);
    min-height: calc(100vh - 210px);
    padding: 26px 0 52px;
  }}
  .panel {{ display: none; }}
  .panel.active {{ display: block; }}
  .subpanel {{ display: none; }}
  .subpanel.active {{ display: block; }}

  .view-bar {{
    display: flex;
    justify-content: space-between;
    align-items: center;
    flex-wrap: wrap;
    gap: 12px;
    margin-bottom: 16px;
  }}
  /* Apple-style Segmented Control */
  .seg-group {{
    display: inline-flex;
    flex-wrap: wrap;
    gap: 2px;
    background: #e9ecef;
    padding: 3.5px;
    border-radius: 10px;
  }}
  .seg-btn {{
    border: none;
    background: transparent;
    padding: 7px 15px;
    border-radius: 7px;
    font-family: var(--sans);
    font-size: 13px;
    font-weight: 500;
    color: var(--ink-2);
    cursor: pointer;
    transition: all 0.14s ease;
  }}
  .seg-btn:hover {{ color: var(--ink); }}
  .seg-btn.active {{
    background: #ffffff;
    color: var(--ink);
    font-weight: 600;
    box-shadow: 0 1px 3px rgba(0,0,0,0.08);
  }}
  .guide-btn {{
    display: inline-flex;
    align-items: center;
    gap: 6px;
    height: 34px;
    padding: 0 13px;
    border-radius: 999px;
    border: 1px solid var(--line);
    background: #ffffff;
    color: var(--ink-2);
    font-size: 12.5px;
    font-weight: 500;
    cursor: pointer;
    transition: all 0.14s ease;
  }}
  .guide-btn:hover {{ border-color: var(--accent); color: var(--accent); }}

  .card {{
    background: var(--surface);
    border: 1px solid var(--line-soft);
    border-radius: 15px;
    padding: 22px 24px;
    box-shadow: 0 1px 3px rgba(32, 33, 36, 0.03);
  }}
  .card-head {{
    display: flex;
    justify-content: space-between;
    align-items: center;
    flex-wrap: wrap;
    gap: 10px;
    margin-bottom: 15px;
  }}
  .card-head h2 {{
    margin: 0;
    font-size: 16.5px;
    font-weight: 600;
    letter-spacing: -0.015em;
  }}
  .meta-note {{ font-size: 12.5px; color: var(--ink-3); }}

  /* Minimal Controls */
  .controls {{
    display: flex;
    flex-wrap: wrap;
    gap: 10px;
    align-items: center;
    margin-bottom: 14px;
  }}
  select, input[type="text"].ctrl-input {{
    height: 36px;
    padding: 0 12px;
    border: 1px solid var(--line);
    border-radius: 8px;
    font-size: 13px;
    font-family: var(--sans);
    color: var(--ink);
    background: #ffffff;
  }}
  select:focus, input[type="text"].ctrl-input:focus, textarea.batch-area:focus {{
    outline: none;
    border-color: var(--accent);
  }}
  textarea.batch-area {{
    width: 100%;
    min-height: 92px;
    padding: 10px 12px;
    border: 1px solid var(--line);
    border-radius: 10px;
    font-size: 13px;
    font-family: var(--mono);
    color: var(--ink);
    background: #ffffff;
    resize: vertical;
    line-height: 1.5;
  }}
  .btn-action {{
    height: 36px;
    background: var(--ink);
    color: #ffffff;
    border: none;
    padding: 0 15px;
    border-radius: 8px;
    font-size: 13px;
    font-weight: 500;
    cursor: pointer;
  }}
  .btn-action:hover {{ opacity: 0.88; }}
  .btn-subtle {{
    height: 30px;
    background: #ffffff;
    color: var(--ink-2);
    border: 1px solid var(--line);
    padding: 0 11px;
    border-radius: 7px;
    font-size: 12px;
    cursor: pointer;
    transition: background 0.12s ease, color 0.12s ease, border-color 0.12s ease;
  }}
  .btn-subtle:hover {{ color: var(--ink); background: var(--canvas); border-color: #bdc1c6; }}
  .step-group {{
    display: inline-flex;
    align-items: center;
    border: 1px solid var(--line);
    border-radius: 8px;
    overflow: hidden;
    background: #ffffff;
  }}
  .step-btn {{
    height: 34px;
    width: 30px;
    border: none;
    background: #ffffff;
    color: var(--ink-2);
    font-size: 16px;
    line-height: 1;
    cursor: pointer;
    padding: 0;
  }}
  .step-btn:hover {{ background: var(--canvas); color: var(--accent); }}
  .step-group input.ctrl-input {{
    border: none;
    border-left: 1px solid var(--line-soft);
    border-right: 1px solid var(--line-soft);
    border-radius: 0;
    height: 34px;
  }}

  /* Lollipop & Enlarged 3D Workbench */
  .svg-track {{
    width: 100%;
    overflow-x: auto;
    border: 1px solid var(--line-soft);
    border-radius: 11px;
    background: #fbfcfd;
    padding: 8px 12px 4px;
    margin-bottom: 18px;
  }}
  .workbench {{
    display: grid;
    grid-template-columns: 1.52fr 1fr;
    gap: 22px;
    align-items: stretch;
  }}
  .workbench.theatre {{
    grid-template-columns: 1fr;
  }}
  @media (max-width: 1100px) {{
    .workbench {{ grid-template-columns: 1fr; }}
  }}
  .vp-shell {{
    border-radius: 13px;
    overflow: hidden;
    background: #111315;
    border: 1px solid #26292e;
    box-shadow: 0 4px 16px rgba(17, 19, 21, 0.10);
  }}
  .vp-bar {{
    background: #1a1d21;
    padding: 9px 12px;
    display: flex;
    flex-wrap: wrap;
    gap: 7px;
    justify-content: space-between;
    align-items: center;
  }}
  .vp-bar select {{
    height: 30px;
    background: #272b31;
    color: #f1f3f4;
    border: 1px solid #3c4043;
    font-size: 12px;
    border-radius: 6px;
    padding: 0 9px;
  }}
  .vp-btn {{
    height: 30px;
    background: #272b31;
    color: #e8eaed;
    border: 1px solid #3c4043;
    border-radius: 6px;
    padding: 0 10px;
    font-size: 12px;
    cursor: pointer;
    transition: background 0.12s ease;
  }}
  .vp-btn:hover {{ background: #353a42; }}
  .vp-btn.active-vp {{ background: var(--accent); border-color: var(--accent); color: #ffffff; }}
  .mol3d-stage {{
    position: relative;
    width: 100%;
    height: 600px;
    background: #111315;
    transition: height 0.18s ease;
  }}
  .workbench.theatre .mol3d-stage {{
    height: 760px;
  }}
  @media (max-width: 720px) {{
    .mol3d-stage {{ height: 440px; }}
  }}
  .mol3d-hud {{
    position: absolute;
    top: 12px;
    left: 14px;
    color: #e8eaed;
    font-family: var(--mono);
    font-size: 11.5px;
    background: rgba(26, 29, 33, 0.80);
    padding: 5px 10px;
    border-radius: 7px;
    pointer-events: none;
    z-index: 5;
  }}
  .mol3d-legend {{
    position: absolute;
    bottom: 12px;
    left: 14px;
    right: 14px;
    background: rgba(20, 22, 26, 0.84);
    color: #dadce0;
    padding: 7px 12px;
    border-radius: 7px;
    font-size: 11.5px;
    display: flex;
    flex-wrap: wrap;
    gap: 15px;
    z-index: 5;
  }}
  .swatch {{ display: inline-block; width: 8.5px; height: 8.5px; border-radius: 50%; margin-right: 5px; }}

  /* Clean Visual Residue Scorecard */
  .scorecard {{
    border: 1px solid var(--line-soft);
    border-radius: 13px;
    padding: 20px 22px;
    background: #ffffff;
    display: flex;
    flex-direction: column;
    justify-content: space-between;
    gap: 14px;
  }}
  .aa-strip {{
    display: grid;
    grid-template-columns: repeat(19, minmax(22px, 1fr));
    gap: 3.5px;
    margin: 12px 0 16px;
  }}
  .aa-cell {{
    background: var(--canvas);
    border: 1px solid var(--line);
    border-top-width: 3px;
    border-radius: 5px;
    height: 32px;
    font-family: var(--mono);
    font-size: 12px;
    font-weight: 600;
    color: var(--ink);
    cursor: pointer;
    padding: 0;
  }}
  .aa-cell:hover {{ background: #e8eaed; }}
  .aa-cell.active {{ background: var(--ink); color: #ffffff; border-color: var(--ink); }}

  .metric-tiles {{
    display: grid;
    grid-template-columns: repeat(3, 1fr);
    gap: 10px;
    margin-bottom: 12px;
  }}
  @media (max-width: 560px) {{
    .metric-tiles {{ grid-template-columns: repeat(2, 1fr); }}
  }}
  .m-tile {{
    background: var(--canvas);
    border: 1px solid var(--line-soft);
    border-radius: 10px;
    padding: 11px 13px;
  }}
  .m-lbl {{
    font-size: 11.5px;
    color: var(--ink-3);
    display: block;
    margin-bottom: 3px;
  }}
  .m-val {{
    font-size: 17px;
    font-weight: 600;
    color: var(--ink);
    font-variant-numeric: tabular-nums;
  }}
  .m-sub {{
    font-size: 11.5px;
    color: var(--ink-2);
    margin-top: 2px;
  }}
  .summary-pill {{
    padding: 12px 14px;
    border-radius: 10px;
    background: var(--canvas);
    font-size: 13px;
    color: var(--ink-2);
    line-height: 1.5;
  }}

  /* Clean Minimal Tables with Universal One-Click Sorting */
  .table-wrap {{
    overflow-x: auto;
    max-height: 620px;
    overflow-y: auto;
    border: 1px solid var(--line-soft);
    border-radius: 11px;
  }}
  table {{
    width: 100%;
    border-collapse: collapse;
    font-size: 13px;
    font-variant-numeric: tabular-nums;
  }}
  th, td {{
    padding: 10px 13px;
    border-bottom: 1px solid var(--line-soft);
    text-align: left;
    white-space: nowrap;
  }}
  th {{
    background: #fafafa;
    position: sticky;
    top: 0;
    z-index: 2;
    font-size: 11.5px;
    font-weight: 600;
    color: var(--ink-2);
    border-bottom: 1px solid var(--line);
  }}
  th.sortable {{
    cursor: pointer;
    user-select: none;
    transition: background 0.12s ease, color 0.12s ease;
  }}
  th.sortable:hover {{
    background: #f1f3f5;
    color: var(--ink);
  }}
  th.sortable::after {{
    content: " ↕";
    font-size: 10px;
    color: #bdc1c6;
    font-weight: 400;
  }}
  th.sortable.sort-asc {{ color: var(--accent); }}
  th.sortable.sort-asc::after {{ content: " ▲"; color: var(--accent); font-size: 9.5px; }}
  th.sortable.sort-desc {{ color: var(--accent); }}
  th.sortable.sort-desc::after {{ content: " ▼"; color: var(--accent); font-size: 9.5px; }}
  tbody tr:hover td {{ background: #f8fafd; }}
  tbody tr.row-flash td {{ background: rgba(0, 112, 155, 0.11); transition: background 0.6s ease; }}
  code {{ font-family: var(--mono); font-size: 12.5px; }}

  .tag {{
    display: inline-flex;
    align-items: center;
    gap: 5px;
    font-size: 11.5px;
    font-weight: 500;
    padding: 2px 8px;
    border-radius: 999px;
    background: var(--subtle);
    color: var(--ink);
  }}
  .dot {{ width: 6px; height: 6px; border-radius: 50%; display: inline-block; }}
  .dot-crimson {{ background: var(--crimson); }}
  .dot-amber {{ background: var(--amber); }}
  .dot-petrol {{ background: var(--accent); }}
  .dot-emerald {{ background: var(--emerald); }}
  .dot-slate {{ background: var(--ink-3); }}

  .hm-cell {{
    text-align: center;
    font-family: var(--mono);
    font-size: 11.5px;
    padding: 7px 5px;
    cursor: pointer;
    transition: outline 0.08s ease;
  }}
  .hm-cell:hover {{
    outline: 2px solid var(--ink);
    outline-offset: -2px;
    position: relative;
    z-index: 1;
  }}
  .hm-ramp {{
    width: 110px;
    height: 9px;
    border-radius: 999px;
    border: 1px solid var(--line);
    background: linear-gradient(90deg, #ffffff 0%, #7fb8cd 50%, #00587c 100%);
    display: inline-block;
  }}

  /* Dotted Interactive Terms & Light Floating KaTeX Popover */
  [data-term] {{
    border-bottom: 1px dotted rgba(95, 99, 104, 0.55);
    cursor: help;
  }}
  th.sortable[data-term] {{
    cursor: pointer;
  }}
  [data-term]:hover {{
    color: var(--accent);
    border-bottom-color: var(--accent);
  }}
  .popover {{
    position: fixed;
    z-index: 9999;
    width: 416px;
    max-width: calc(100vw - 24px);
    background: #ffffff;
    color: var(--ink);
    border: 1px solid var(--line);
    border-radius: 14px;
    padding: 14px 16px;
    box-shadow: 0 12px 32px rgba(32, 33, 36, 0.16);
    display: none;
  }}
  .popover.visible {{ display: block; }}
  .pop-cat {{
    font-family: var(--mono);
    font-size: 10px;
    text-transform: uppercase;
    letter-spacing: 0.05em;
    color: var(--accent);
  }}
  .pop-title {{ font-size: 14px; font-weight: 600; margin: 2px 0 6px; }}
  .pop-desc {{ font-size: 12.5px; color: var(--ink-2); margin-bottom: 10px; }}
  .pop-math {{
    background: var(--canvas);
    border: 1px solid var(--line-soft);
    border-radius: 8px;
    padding: 8px 10px;
    margin-bottom: 8px;
    overflow-x: auto;
    text-align: center;
    font-size: 13px;
  }}
  .pop-thresh {{ font-size: 11.5px; color: var(--ink-3); }}
  .pop-cite {{
    margin-top: 7px;
    padding-top: 7px;
    border-top: 1px solid var(--line-soft);
    font-size: 11px;
    color: var(--ink-3);
    display: flex;
    flex-wrap: wrap;
    justify-content: space-between;
    align-items: center;
    gap: 6px;
  }}
  .pop-cite a {{
    color: var(--accent);
    text-decoration: none;
    font-weight: 500;
  }}
  .pop-cite a:hover {{ text-decoration: underline; }}

  /* Slide-Over Contextual Guide Drawer */
  .drawer-backdrop {{
    position: fixed;
    inset: 0;
    background: rgba(32, 33, 36, 0.25);
    z-index: 90;
    display: none;
  }}
  .drawer-backdrop.open {{ display: block; }}
  .drawer {{
    position: fixed;
    top: 0;
    right: 0;
    bottom: 0;
    width: 450px;
    max-width: 92vw;
    background: #ffffff;
    border-left: 1px solid var(--line);
    box-shadow: -8px 0 28px rgba(32, 33, 36, 0.12);
    z-index: 95;
    padding: 24px;
    overflow-y: auto;
    transform: translateX(100%);
    transition: transform 0.2s ease;
  }}
  .drawer.open {{ transform: translateX(0); }}

  /* Spacious LaTeX Equation Gallery Cards */
  .eq-grid {{
    display: grid;
    grid-template-columns: repeat(2, 1fr);
    gap: 16px;
  }}
  @media (max-width: 1024px) {{
    .eq-grid {{ grid-template-columns: 1fr; }}
  }}
  .eq-card {{
    background: #ffffff;
    border: 1px solid var(--line-soft);
    border-radius: 12px;
    padding: 18px 20px;
    display: flex;
    flex-direction: column;
    justify-content: space-between;
    gap: 12px;
  }}
  .eq-math-box {{
    background: var(--canvas);
    border-radius: 9px;
    padding: 14px 12px;
    overflow-x: auto;
    text-align: center;
    min-height: 60px;
    display: flex;
    align-items: center;
    justify-content: center;
  }}
  .eq-vars {{
    margin-top: 9px;
    padding: 8px 11px;
    background: #fafbfc;
    border: 1px solid var(--line-soft);
    border-radius: 8px;
    display: flex;
    flex-direction: column;
    gap: 5px;
  }}
  .eq-vars-hdr {{
    font-family: var(--mono);
    font-size: 10px;
    font-weight: 600;
    text-transform: uppercase;
    letter-spacing: 0.04em;
    color: var(--ink-3);
    margin-bottom: 1px;
  }}
  .eq-var-row {{
    display: flex;
    flex-wrap: wrap;
    align-items: baseline;
    gap: 4px 9px;
    font-size: 11.5px;
    line-height: 1.42;
  }}
  .eq-var-sym {{
    display: inline-block;
    background: #ffffff;
    border: 1px solid var(--line-soft);
    border-radius: 5px;
    padding: 1px 7px;
    font-size: 11.5px;
    color: var(--ink);
    flex-shrink: 0;
  }}
  .eq-var-desc {{
    color: var(--ink-2);
    flex: 1;
    min-width: 170px;
  }}
  .eq-cite {{
    margin-top: 4px;
    padding-top: 9px;
    border-top: 1px solid var(--line-soft);
    font-size: 11.5px;
    color: var(--ink-3);
    display: flex;
    flex-wrap: wrap;
    align-items: center;
    gap: 5px 8px;
    line-height: 1.45;
  }}
  .eq-cite a {{
    color: var(--accent);
    text-decoration: none;
    font-weight: 500;
  }}
  .eq-cite a:hover {{ text-decoration: underline; }}
  .pathway-steps {{
    display: grid;
    grid-template-columns: repeat(5, 1fr);
    gap: 12px;
    margin-bottom: 18px;
  }}
  @media (max-width: 980px) {{
    .pathway-steps {{ grid-template-columns: repeat(2, 1fr); }}
  }}
  .p-card {{
    background: var(--canvas);
    border: 1px solid var(--line-soft);
    border-radius: 10px;
    padding: 14px;
  }}

  footer {{
    border-top: 1px solid var(--line-soft);
    background: #ffffff;
    padding: 20px 0;
    color: var(--ink-3);
    font-size: 12.5px;
  }}

  /* Clean Home-Page About, Authors & Citation Section */
  .about-section {{
    background: #ffffff;
    border-top: 1px solid var(--line-soft);
    padding: 40px 0 36px;
  }}
  .about-grid {{
    display: grid;
    grid-template-columns: 1.25fr 1fr;
    gap: 24px;
    margin-bottom: 24px;
  }}
  @media (max-width: 960px) {{
    .about-grid {{ grid-template-columns: 1fr; }}
  }}
  .pillars-row {{
    display: grid;
    grid-template-columns: repeat(3, 1fr);
    gap: 10px;
    margin-top: 14px;
  }}
  @media (max-width: 680px) {{
    .pillars-row {{ grid-template-columns: 1fr; }}
  }}
  .pillar-tile {{
    background: var(--canvas);
    border: 1px solid var(--line-soft);
    border-radius: 10px;
    padding: 11px 13px;
  }}
  .authors-grid {{
    display: grid;
    grid-template-columns: repeat(4, 1fr);
    gap: 12px;
    margin-bottom: 20px;
  }}
  @media (max-width: 1024px) {{
    .authors-grid {{ grid-template-columns: repeat(2, 1fr); }}
  }}
  @media (max-width: 560px) {{
    .authors-grid {{ grid-template-columns: 1fr; }}
  }}
  .author-card {{
    border: 1px solid var(--line-soft);
    border-radius: 11px;
    padding: 13px 15px;
    background: var(--canvas);
    display: flex;
    flex-direction: column;
    justify-content: space-between;
    gap: 8px;
  }}
  .orcid-link {{
    display: inline-flex;
    align-items: center;
    gap: 5px;
    font-family: var(--mono);
    font-size: 11px;
    color: var(--emerald);
    text-decoration: none;
  }}
  .orcid-link:hover {{ text-decoration: underline; }}
  .cite-box {{
    background: var(--canvas);
    border: 1px solid var(--line-soft);
    border-radius: 11px;
    padding: 14px 16px;
    font-size: 12.5px;
    color: var(--ink-2);
  }}

  /* Dedicated Sub-Header Public Datasets Row, Visual Legends & Table Guides */
  .public-bar {{
    background: #f8fafc;
    border-bottom: 1px solid var(--line-soft);
    padding: 7px 0;
    font-size: 11.5px;
  }}
  .public-bar-inner {{
    display: flex;
    align-items: center;
    justify-content: space-between;
    flex-wrap: wrap;
    gap: 7px 12px;
  }}
  .pub-left {{
    display: inline-flex;
    align-items: center;
    gap: 7px;
    color: var(--ink-2);
    font-weight: 600;
    white-space: nowrap;
  }}
  .pub-gene-sel {{
    height: 25px !important;
    padding: 0 8px !important;
    font-size: 11.5px !important;
    font-weight: 600;
    border-radius: 999px !important;
    background: #ffffff;
    color: var(--accent) !important;
    border: 1px solid rgba(0, 112, 155, 0.30) !important;
    cursor: pointer;
  }}
  .pub-links {{
    display: flex;
    flex-wrap: wrap;
    align-items: center;
    gap: 5px;
  }}
  .xref-pill {{
    display: inline-flex;
    align-items: center;
    gap: 4px;
    padding: 2px 9px;
    border-radius: 999px;
    border: 1px solid var(--line-soft);
    background: #ffffff;
    color: var(--ink-2);
    text-decoration: none;
    font-size: 11px;
    transition: all 0.12s ease;
    white-space: nowrap;
  }}
  .xref-pill:hover {{
    border-color: var(--accent);
    color: var(--accent);
    background: var(--accent-soft);
  }}
  .loeuf-pill {{
    display: inline-flex;
    align-items: center;
    gap: 5px;
    padding: 2px 9px;
    border-radius: 999px;
    font-family: var(--mono);
    font-size: 11px;
    font-weight: 600;
  }}
  .plot-legend {{
    display: flex;
    flex-wrap: wrap;
    align-items: center;
    gap: 8px 15px;
    padding: 7px 11px;
    margin-top: 4px;
    background: #f4f6f9;
    border-top: 1px solid var(--line-soft);
    border-radius: 7px;
    font-size: 11.5px;
    color: var(--ink-2);
  }}
  .table-guide {{
    display: flex;
    flex-wrap: wrap;
    align-items: center;
    justify-content: space-between;
    gap: 6px 14px;
    padding: 8px 12px;
    margin-bottom: 11px;
    background: var(--canvas);
    border: 1px solid var(--line-soft);
    border-left: 3px solid var(--accent);
    border-radius: 8px;
    font-size: 12px;
    color: var(--ink-2);
    line-height: 1.45;
  }}
  .leg-item {{
    display: inline-flex;
    align-items: center;
    gap: 5px;
    white-space: nowrap;
  }}
  .leg-line {{
    display: inline-block;
    width: 14px;
    height: 2.5px;
    border-radius: 1px;
  }}
  .leg-box {{
    display: inline-block;
    width: 10px;
    height: 10px;
    border-radius: 2.5px;
  }}
  .contact-hud {{
    position: absolute;
    top: 12px;
    right: 14px;
    width: 214px;
    background: rgba(20, 22, 26, 0.92);
    border: 1px solid #3c4043;
    border-radius: 10px;
    padding: 9px 10px;
    z-index: 8;
    display: none;
    box-shadow: 0 8px 24px rgba(0,0,0,0.38);
  }}
  .contact-hud.open {{ display: block; }}
  .score-bottom {{
    display: grid;
    grid-template-columns: 1fr auto;
    gap: 12px;
    align-items: center;
  }}
  @media (max-width: 620px) {{
    .score-bottom {{ grid-template-columns: 1fr; }}
  }}
  .radar-box {{
    width: 132px;
    min-height: 128px;
    padding: 4px 4px 5px;
    border: 1px solid var(--line-soft);
    border-radius: 10px;
    background: var(--canvas);
    display: flex;
    flex-direction: column;
    align-items: center;
    justify-content: center;
    flex-shrink: 0;
  }}

  /* 1-Page Printable Variant Report Mode */
  #print-report-banner, #print-3d-img {{ display: none; }}
  @media print {{
    body {{ background: #ffffff; font-size: 11.5px; }}
    .topbar, .public-bar, .search-hero, .view-bar, #exp-dms, #exp-tcga, #exp-batch, .table-guide, .plot-legend, .contact-hud, .about-section, footer, .vp-bar, .guide-btn, .no-print, #lollipop-svg-wrap, .popover, .drawer, .drawer-backdrop {{
      display: none !important;
    }}
    main.workspace {{ background: #ffffff; padding: 0; min-height: auto; }}
    .wrap {{ max-width: 100%; padding: 0; }}
    #print-report-banner {{
      display: flex !important;
      justify-content: space-between;
      align-items: baseline;
      border-bottom: 2px solid #202124;
      padding-bottom: 8px;
      margin-bottom: 12px;
    }}
    .card {{ border: none; padding: 0; box-shadow: none; }}
    #struct3d-card > .card-head {{ display: none !important; }}
    .workbench {{ grid-template-columns: 1fr 1fr !important; gap: 14px; }}
    .mol3d-stage {{ display: none !important; }}
    #print-3d-img {{
      display: block !important;
      width: 100%;
      max-height: 310px;
      object-fit: contain;
      background: #111315;
    }}
    .aa-strip {{ display: none !important; }}
    #exp-codon {{ display: block !important; margin-top: 14px; }}
    #exp-codon .table-wrap {{ max-height: none !important; overflow: visible !important; }}
  }}
</style>
</head>
<body>

<!-- Minimal Top Header -->
<header class="topbar">
  <div class="wrap topbar-inner">
    <div class="brand" onclick="switchTab('tab-variants', document.getElementById('btn-tab-variants'))">
      <span class="brand-dot"></span>
      <span>miRBioGenAtlas</span>
      <span class="ver-pill">v2.0.2</span>
    </div>

    <nav class="nav-tabs" aria-label="Sections">
      <button class="nav-tab active" id="btn-tab-variants" onclick="switchTab('tab-variants', this)">Explorer</button>
      <button class="nav-tab" id="btn-tab-heatmap" onclick="switchTab('tab-heatmap', this)">Pan-Cancer</button>
      <button class="nav-tab" id="btn-tab-mirna" onclick="switchTab('tab-mirna', this)">Small-RNA</button>
      <button class="nav-tab" id="btn-tab-clinical" onclick="switchTab('tab-clinical', this)">Clinical &amp; Biology</button>
      <button class="nav-tab" id="btn-tab-docs" onclick="switchTab('tab-docs', this)">Formulas &amp; Data</button>
      <button class="nav-tab" onclick="document.getElementById('about-atlas').scrollIntoView({{behavior:'smooth'}})">About &amp; Cite</button>
    </nav>

    <div class="top-links">
      <a href="https://doi.org/10.5281/zenodo.24455067" target="_blank" rel="noopener">DOI</a>
      <a href="https://github.com/siddugayam/miRBioGenAtlas" target="_blank" rel="noopener">GitHub</a>
    </div>
  </div>
</header>

<!-- Dedicated Sub-Header Row: Access from Public Datasets -->
<div class="public-bar no-print" role="region" aria-label="Access from Public Datasets">
  <div class="wrap public-bar-inner">
    <div class="pub-left">
      <span>Access from Public Datasets:</span>
      <select id="pubbar-gene" class="pub-gene-sel" aria-label="Select gene for external public datasets" onchange="syncGeneFromPublicBar(this.value)"></select>
    </div>
    <div class="pub-links" id="public-datasets-links"></div>
  </div>
</div>

<!-- Google-Simple Omni-Search Hero -->
<section class="search-hero">
  <div class="wrap">
    <div class="search-shell">
      <div class="search-box">
        <svg class="search-icon" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
          <circle cx="11" cy="11" r="7"></circle>
          <line x1="20" y1="20" x2="16.35" y2="16.35"></line>
        </svg>
        <input type="text" id="global-search" autocomplete="off" placeholder="Search gene, AA/cDNA/GRCh38 VCF, rsID, or formula (e.g. DICER1 E1705K, chr14:95094139:C>T, c.5113G>A)..." oninput="updateSearchSuggestions(this.value)" onfocus="showRecentOrSuggestions(this.value)" onkeydown="handleSearchKey(event)">
        <span class="search-kbd">&#8984;K</span>
      </div>
      <div id="suggest-menu" class="suggest-menu"></div>
    </div>

    <div class="quick-Row">
      <button class="chip" onclick="loadPresetMutation('DICER1','E1705K')">DICER1 E1705K</button>
      <button class="chip" onclick="runGlobalSearch('chr14:95094139:C>T')" title="Resolve GRCh38 VCF coordinate directly to 3D protein">chr14:95094139:C&gt;T</button>
      <button class="chip" onclick="loadPresetMutation('DROSHA','E1147K')">DROSHA E1147K</button>
      <button class="chip" onclick="loadPresetMutation('TP53','R175H')">TP53 R175H</button>
      <button class="chip" onclick="openDmsLandscape()">2D DMS Heatmap</button>
      <button class="chip" onclick="openBatchTriage()">.VCF Triage</button>
      <button class="chip" onclick="openApiDocs()">REST API &amp; CLI</button>
      <button class="chip" onclick="openFormulaTab()">LaTeX Equations &rarr;</button>
    </div>

    <div class="stat-ribbon">
      <span data-term="kpi_tumours"><strong>21,140</strong> tumours</span>
      <span data-term="kpi_proteins"><strong>15</strong> proteins</span>
      <span data-term="kpi_sat"><strong>213,940</strong> 3D variants</span>
      <span data-term="kpi_pdb"><strong>16</strong> PDB complexes</span>
      <span data-term="kpi_mirna"><strong>9,578</strong> miRNA-seq</span>
      <span data-term="kpi_auc"><strong>0.978</strong> ClinVar AUC</span>
    </div>
  </div>
</section>

<main class="workspace">
  <div class="wrap">

    <!-- TAB 1: 3D VARIANT & SATURATION EXPLORER -->
    <section id="tab-variants" class="panel active">
      <div id="print-report-banner">
        <div>
          <strong style="font-size:15px;">miRBioGenAtlas (v2.0.2) — 3D Structural &amp; Biophysical Variant Report</strong>
          <div style="font-size:11px;color:#5f6368;">DOI: 10.5281/zenodo.24455067 &middot; https://siddugayam.github.io/miRBioGenAtlas/</div>
        </div>
        <span id="print-report-date" style="font-family:var(--mono);font-size:11px;color:#5f6368;"></span>
      </div>

      <div class="card" id="struct3d-card">
        <div class="card-head">
          <div style="display:flex;align-items:center;gap:10px;flex-wrap:wrap;">
            <select id="lol-gene" onchange="renderLollipopAndLoadSat()" style="font-weight:600;"></select>
            <div class="step-group" title="Step along residues with ‹ / › or Left / Right arrow keys">
              <button class="step-btn" onclick="stepCodon(-1)" aria-label="Previous residue">&lsaquo;</button>
              <input type="text" class="ctrl-input" id="sat-query" placeholder="Codon or variant (E1705K)" oninput="inspectSaturationCodon()" style="width:150px;font-family:var(--mono);font-weight:500;">
              <button class="step-btn" onclick="stepCodon(1)" aria-label="Next residue">&rsaquo;</button>
            </div>
            <span class="meta-note" id="sat-status"></span>
            <span id="gene-bio-badges" style="display:inline-flex;gap:6px;align-items:center;flex-wrap:wrap;"></span>
          </div>
          <div style="display:flex;gap:7px;align-items:center;">
            <button class="btn-subtle" onclick="downloadLollipopSVG()" title="Export publication-ready vector SVG of domain &amp; lollipop track">Track SVG</button>
            <button class="guide-btn" onclick="openGuideDrawer('explorer')">&#9432; How this works</button>
          </div>
        </div>

        <div id="lollipop-svg-wrap" class="svg-track"></div>

        <div class="workbench" id="explorer-workbench">
          <div>
            <div class="vp-shell">
              <div class="vp-bar">
                <div style="display:flex;gap:6px;flex-wrap:wrap;">
                  <select id="mol3d-mode" onchange="refresh3DScene()" aria-label="3D Mode">
                    <option value="overlay">Overlay: WT + Mutant</option>
                    <option value="mutant">Mutant Only</option>
                    <option value="wt">Wild-Type Only</option>
                    <option value="rcsb">Experimental PDB</option>
                  </select>
                  <select id="mol3d-color" onchange="refresh3DScene()" aria-label="Coloring">
                    <option value="plddt">Color: pLDDT</option>
                    <option value="domain">Color: Domains</option>
                    <option value="am">Color: AlphaMissense</option>
                    <option value="pocket">Color: Active Pocket</option>
                  </select>
                </div>
                <div style="display:flex;gap:5px;flex-wrap:wrap;">
                  <button class="vp-btn" onclick="zoomToMutatedResidue()">8&Aring; Focus</button>
                  <button class="vp-btn" onclick="resetFullProtein3D()">Full</button>
                  <button class="vp-btn" id="btn-nbrs" onclick="toggleNeighbors3D()">Shell: On</button>
                  <button class="vp-btn" id="btn-surf" onclick="toggleSurface3D()" title="Render translucent 8.5Å pocket Van der Waals surface">Surface</button>
                  <button class="vp-btn" id="btn-2dmap" onclick="toggleContactMap2D()" title="AlphaFold 2D Cα domain packing &amp; contact matrix">2D Map</button>
                  <button class="vp-btn" id="btn-spin" onclick="toggleSpin3D()">Spin</button>
                  <button class="vp-btn" onclick="download3DSnapshotPNG()" title="Save 3D PNG image">PNG</button>
                  <button class="vp-btn" onclick="downloadMutatedPDB()" title="Download atomic PDB coordinates">PDB</button>
                  <button class="vp-btn" id="btn-expand3d" onclick="toggleTheatre3D()" title="Expand 3D viewer to full width">Expand</button>
                </div>
              </div>
              <img id="print-3d-img" alt="3D Structural Snapshot">
              <div class="mol3d-stage">
                <div class="mol3d-hud" id="mol3d-hud">DICER1</div>
                <div class="contact-hud" id="contact-hud">
                  <div style="display:flex;justify-content:space-between;align-items:center;margin-bottom:5px;color:#e8eaed;font-size:11px;font-family:var(--mono);">
                    <span>2D C&alpha; Packing (&le;24&Aring;)</span>
                    <button onclick="toggleContactMap2D()" style="border:none;background:none;color:#9aa0a6;cursor:pointer;font-size:14px;line-height:1;">&times;</button>
                  </div>
                  <canvas id="contact-canvas" width="192" height="192" style="width:192px;height:192px;border-radius:6px;cursor:crosshair;display:block;background:#111315;" onclick="handleContactMapClick(event)" title="Click anywhere on 2D matrix to jump to that residue in 3D"></canvas>
                  <div style="display:flex;align-items:center;justify-content:space-between;gap:5px;margin-top:5px;font-size:9.5px;color:#cbd5e1;font-family:var(--mono);">
                    <span>0&Aring;</span>
                    <span style="flex:1;height:5px;border-radius:999px;background:linear-gradient(90deg,#12bfa9 0%,#135664 55%,#16191e 100%);border:1px solid #3c4043;"></span>
                    <span>24&Aring;</span>
                  </div>
                  <div style="display:flex;justify-content:space-between;margin-top:3px;font-size:9.5px;color:#9aa0a6;font-family:var(--mono);">
                    <span>1 aa</span><span style="color:#ff4d6d;">&#10010; Active aa</span><span id="contact-max-lbl">L aa</span>
                  </div>
                </div>
                <div id="mol3d-canvas" style="width:100%;height:100%;position:relative;"></div>
                <div class="mol3d-legend" id="mol3d-legend">
                  <span><i class="swatch" style="background:#0095c8;"></i>WT</span>
                  <span data-term="rotamer_clash"><i class="swatch" style="background:#e11d48;"></i>Mutant rotamer</span>
                  <span data-term="burial_8a"><i class="swatch" style="background:#94a3b8;"></i>8&Aring; shell</span>
                  <span data-term="dist_na"><i class="swatch" style="background:#f59e0b;"></i>Contact</span>
                  <span data-term="rotamer_clash"><i class="swatch" style="background:#ff3b30;"></i>Clash (&lt;2.4&Aring;)</span>
                </div>
              </div>
            </div>
          </div>

          <div class="scorecard" id="prop-delta-panel"></div>
        </div>
      </div>

      <!-- Segmented Sub-View Below 3D Workbench -->
      <div class="view-bar" style="margin-top:20px;">
        <div class="seg-group">
          <button class="seg-btn active" id="seg-exp-codon" onclick="switchSub('exp-sub', 'exp-codon', this)">Codon 19-AA Scan</button>
          <button class="seg-btn" id="seg-exp-dms" onclick="switchSub('exp-sub', 'exp-dms', this); renderDmsHeatmap();">2D DMS Landscape</button>
          <button class="seg-btn" id="seg-exp-tcga" onclick="switchSub('exp-sub', 'exp-tcga', this)">Observed Tumour Variants (6,292)</button>
          <button class="seg-btn" id="seg-exp-batch" onclick="switchSub('exp-sub', 'exp-batch', this)">.VCF Variant Triage</button>
        </div>
      </div>

      <div id="exp-codon" class="subpanel exp-sub active">
        <div class="card" id="sat-codon-detail"></div>
      </div>

      <div id="exp-dms" class="subpanel exp-sub">
        <div class="card">
          <div class="card-head">
            <div>
              <h2>2D Deep Mutational Scanning (DMS) Saturation Heatmap</h2>
              <div class="meta-note" style="margin-top:2px;">20 amino acids &times; sliding sequence window &middot; Click any cell to model its 3D side-chain rotamer above</div>
            </div>
            <div style="display:flex;gap:8px;align-items:center;flex-wrap:wrap;">
              <select id="dms-metric" onchange="renderDmsHeatmap()">
                <option value="am">Metric: AlphaMissense Pathogenicity</option>
                <option value="if1">Metric: ESM-IF1 3D Fold &Delta;logP</option>
                <option value="e2">Metric: ESM-2 Evolutionary &Delta;&ell;</option>
              </select>
              <select id="dms-span" onchange="renderDmsHeatmap()">
                <option value="25">Window: &plusmn;25 aa (51 codons)</option>
                <option value="40">Window: &plusmn;40 aa (81 codons)</option>
              </select>
              <select id="dms-domain-sel" onchange="jumpDmsToDomain(this.value)">
                <option value="">Jump to domain...</option>
              </select>
              <div class="step-group">
                <button class="step-btn" onclick="stepDmsWindow(-25)" title="Shift window N-terminal">&lsaquo;</button>
                <button class="step-btn" onclick="stepDmsWindow(25)" title="Shift window C-terminal">&rsaquo;</button>
              </div>
            </div>
          </div>
          <div id="dms-legend-bar" class="table-guide"></div>
          <div id="dms-svg-wrap" class="svg-track" style="margin-bottom:0;"></div>
        </div>
      </div>

      <div id="exp-tcga" class="subpanel exp-sub">
        <div class="card">
          <div class="card-head">
            <h2>Observed TCGA Somatic Coding Variants</h2>
            <span class="meta-note" id="var-count"></span>
          </div>
          <div class="controls">
            <select id="var-gene" onchange="renderVariants()"><option value="ALL">All 15 genes</option></select>
            <select id="var-tier" onchange="renderVariants()">
              <option value="ALL">All structural tiers</option>
              <option value="3D_Catalytic_or_Substrate_Pocket_Driver">3D pocket (&le;6 &Aring;)</option>
              <option value="Pathogenic_Fold_or_Interface_Missense">Pathogenic fold / interface</option>
              <option value="Truncating_LoF">Truncating LoF</option>
              <option value="Ambiguous_Missense_VUS">Ambiguous VUS</option>
              <option value="Likely_Benign_Passenger_Missense">Likely benign</option>
            </select>
            <select id="var-mech" onchange="renderVariants()">
              <option value="ALL">All 3D mechanisms</option>
              <option value="Fold_Preserving_Catalytic_or_Interface">Active / interface</option>
              <option value="Thermodynamic_Fold_Destabilizing">Fold-destabilising</option>
              <option value="Tolerated_Neutral">Neutral</option>
            </select>
            <input type="text" class="ctrl-input" id="var-q" placeholder="Filter variant, study, ClinVar..." oninput="renderVariants()">
            <button class="btn-action" onclick="exportFilteredVariantsCSV()">Export CSV</button>
          </div>
          <div class="table-guide">
            <span><strong>Reading Guide:</strong> <code>RNA / Metal (&Aring;)</code> = heavy-atom distance in Cryo-EM/X-ray PDB (<strong>&le;6.0 &Aring;</strong> = active pocket) &middot; <code>AlphaMissense &ge; 0.564</code> pathogenic &middot; <code>ESM-IF1 &gt; &minus;1.5</code> fold intact vs <code>&le; &minus;1.5</code> unfolding.</span>
            <span class="meta-note">Click <strong>View</strong> to load in 3D &middot; Click <code>&#8597;</code> headers to sort &middot; Hover dotted headers for LaTeX</span>
          </div>
          <div class="table-wrap">
            <table id="var-table">
              <thead>
                <tr>
                  <th>3D</th><th>Gene</th><th>Variant</th>
                  <th data-term="kpi_tumours">TCGA</th>
                  <th data-term="kpi_auc">ClinVar</th>
                  <th data-term="struct_tier">Tier</th>
                  <th data-term="mech_regime">Mechanism</th>
                  <th data-term="kpi_pdb">PDB</th>
                  <th data-term="dist_na">RNA (&Aring;)</th>
                  <th data-term="dist_metal">Metal (&Aring;)</th>
                  <th data-term="alphamissense">AlphaMissense</th>
                  <th data-term="esm2_llr">ESM-2</th>
                  <th data-term="esm1v_llr">ESM-1v</th>
                  <th data-term="esm_if1">ESM-IF1</th>
                  <th>Cohorts</th>
                </tr>
              </thead>
              <tbody></tbody>
            </table>
          </div>
        </div>
      </div>

      <div id="exp-batch" class="subpanel exp-sub">
        <div class="card" ondragover="event.preventDefault()" ondrop="handleBatchFileDrop(event)">
          <div class="card-head">
            <div>
              <h2>.VCF Clinical &amp; Somatic Variant Triage (GRCh38 / GRCh37)</h2>
              <div class="meta-note" style="margin-top:2px;">Upload or drop a standard uncompressed <code>.vcf</code> file, or paste VCF/HGVS variant calls below (100% client-side in browser; zero server upload)</div>
            </div>
            <span class="meta-note" id="batch-summary"></span>
          </div>
          <textarea id="batch-input" class="batch-area" placeholder="#CHROM	POS	ID	REF	ALT&#10;chr14	95094139	.	C	T&#10;chr5	31409348	.	C	T&#10;chr17	7675088	.	C	T"></textarea>
          <div class="controls" style="margin-top:10px;">
            <button class="btn-action" onclick="runBatchTriage()">Score .VCF Variants</button>
            <label class="btn-subtle" style="height:34px;display:inline-flex;align-items:center;cursor:pointer;font-weight:600;color:var(--accent);border-color:rgba(0,112,155,0.35);">
              Upload .VCF
              <input type="file" id="vcf-file-input" accept=".vcf" style="display:none;" onchange="handleBatchFileSelect(this.files)">
            </label>
            <button class="btn-subtle" style="height:34px;" onclick="loadBatchExample()">Load Sample .VCF</button>
            <button class="btn-subtle" style="height:34px;" onclick="clearBatchTriage()">Clear</button>
            <button class="btn-subtle" style="height:34px;" onclick="exportBatchTriageCSV()">Export Triage CSV</button>
          </div>
          <div class="table-guide">
            <span><strong>.VCF Triage Key:</strong> Auto-detects <strong>GRCh38</strong> &amp; <strong>GRCh37/hg19</strong> canonical MANE Select CDS exons across all 15 genes (non-panel chromosomes are skipped automatically) &middot; Classifies <span class="leg-item"><i class="dot dot-crimson"></i>Active site</span> vs <span class="leg-item"><i class="dot dot-amber"></i>Unfolding</span> vs <span class="leg-item"><i class="dot dot-petrol"></i>Neutral</span>.</span>
            <span class="meta-note">File upload accepts <code>.vcf</code> strictly</span>
          </div>
          <div class="table-wrap">
            <table id="batch-table">
              <thead>
                <tr>
                  <th>3D</th><th>Input / Genomic</th><th>Gene</th><th>Protein &amp; cDNA</th>
                  <th data-term="mech_regime">3D Mechanism</th>
                  <th data-term="alphamissense">AlphaMissense</th>
                  <th data-term="esm_if1">ESM-IF1 3D &Delta;logP</th>
                  <th data-term="esm2_llr">ESM-2</th>
                  <th data-term="grantham">Grantham D</th>
                  <th data-term="dist_na">Pocket (&Aring;)</th>
                  <th>TCGA</th><th>ClinVar</th>
                </tr>
              </thead>
              <tbody></tbody>
            </table>
          </div>
        </div>
      </div>
    </section>

    <!-- TAB 2: PAN-CANCER LANDSCAPE -->
    <section id="tab-heatmap" class="panel">
      <div class="view-bar">
        <div class="seg-group">
          <button class="seg-btn active" onclick="switchSub('pc-sub', 'pc-matrix', this)">15 &times; 32 Gene Matrix</button>
          <button class="seg-btn" onclick="switchSub('pc-sub', 'pc-ranking', this)">32-Study Pathway Burden</button>
          <button class="seg-btn" onclick="switchSub('pc-sub', 'pc-codes', this)">Study Denominators (N = 10,195)</button>
        </div>
        <button class="guide-btn" onclick="openGuideDrawer('pancancer')">&#9432; How this works</button>
      </div>

      <div id="pc-matrix" class="subpanel pc-sub active">
        <div class="card">
          <div class="card-head">
            <div>
              <h2>Pan-Cancer Alteration Matrix (15 Genes &times; 32 TCGA Studies)</h2>
              <div class="meta-note" style="margin-top:2px;">Matched exome + CNA cohorts (N = 10,195 tumours) &middot; Hover any cell for exact counts or click to inspect variants</div>
            </div>
            <div style="display:flex;align-items:center;gap:12px;flex-wrap:wrap;">
              <div id="hm-legend" style="display:flex;align-items:center;gap:7px;font-size:11.5px;color:var(--ink-2);font-family:var(--mono);"></div>
              <select id="hm-metric" onchange="renderHeatmap()">
                <option value="hc">Strict: Mutation or Focal GISTIC2 &plusmn;2</option>
                <option value="fcna">Focal High-Level CNA Only (&plusmn;2)</option>
                <option value="mut">Somatic Coding Mutation (All)</option>
                <option value="nhm">Somatic Mutation (Non-Hypermutated)</option>
                <option value="any">Broad: Mutation or Any GISTIC2 &ne; 0</option>
              </select>
              <button class="btn-action" onclick="exportLongCSV()">Export CSV</button>
            </div>
          </div>
          <div id="hm-explainer" class="table-guide"></div>
          <div id="hm-container" class="table-wrap" style="max-height:none;"></div>
        </div>
      </div>

      <div id="pc-ranking" class="subpanel pc-sub">
        <div class="card">
          <div class="card-head">
            <h2>Synteny-Aware Patient Pathway Burden (Excluding TP53)</h2>
            <span class="meta-note">Collapses tandem 1p34.3 (AGO1/3/4) &amp; 18q21.1 (SMAD2/4) clusters</span>
          </div>
          <div class="table-guide">
            <span><strong>Reading Guide:</strong> <strong>Core-10</strong> = dedicated maturation/RISC enzymes &middot; <strong>Extended-12</strong> = +<code>LIN28A/B</code> &middot; <strong>Panel-14</strong> = +<code>SMAD2/4</code> (<code>TP53</code> excluded so ubiquitous p53 loss does not mask miRNA-pathway specificity). Co-deleted <code>1p34.3</code> or <code>18q21.1</code> genes count once per patient.</span>
            <span class="meta-note"><code>*</code> = Low-N study (N &lt; 100) &middot; Click <code>&#8597;</code> to sort</span>
          </div>
          <div class="table-wrap">
            <table id="cohort-burden-table">
              <thead>
                <tr>
                  <th>Rank</th><th>Study</th><th>Matched N</th>
                  <th data-term="tier_core10">Core-10 Strict % (n)</th>
                  <th data-term="tier_ext12">Extended-12 %</th>
                  <th data-term="tier_panel14">Panel-14 %</th>
                  <th data-term="broad_any">Core-10 Broad %</th>
                  <th data-term="synteny">1p34.3 Co-CNA %</th>
                </tr>
              </thead>
              <tbody></tbody>
            </table>
          </div>
        </div>
      </div>

      <div id="pc-codes" class="subpanel pc-sub">
        <div class="card">
          <div class="card-head">
            <h2>TCGA 32-Study Acronym Decoder &amp; Sample Counts</h2>
            <span class="meta-note">Exact S_seq &cap; S_CNA denominators</span>
          </div>
          <div class="table-guide">
            <span><strong>Reading Guide:</strong> <strong>Matched N</strong> = exact intersection of exome-sequenced (<code>Exome</code>) and SNP6 copy-number profiled (<code>CNA</code>) primary tumours audited against Hoadley et al. (<em>Cell</em> 2018) &middot; <strong>Hypermut</strong> = POLE/MSI tumours (&gt;10 mut/Mb or &ge;500 coding mutations).</span>
            <span class="meta-note"><code>*</code> = Study N &lt; 100 &middot; Click <code>&#8597;</code> to sort</span>
          </div>
          <div class="table-wrap">
            <table id="study-decoder-table">
              <thead>
                <tr>
                  <th>Code</th><th>Cancer Type</th><th>Matched N</th><th>Exome</th><th>CNA</th><th>mRNA</th><th>miRNA</th>
                  <th data-term="hypermutator">Hypermut</th><th data-term="hypermutator">Non-Hypermut</th>
                </tr>
              </thead>
              <tbody></tbody>
            </table>
          </div>
        </div>
      </div>
    </section>

    <!-- TAB 3: SMALL-RNA & TARGETS -->
    <section id="tab-mirna" class="panel">
      <div class="view-bar">
        <div class="seg-group">
          <button class="seg-btn active" onclick="switchSub('mir-sub', 'mir-bench', this)">Enzymatic Benchmarks</button>
          <button class="seg-btn" onclick="switchSub('mir-sub', 'mir-dicer', this)">DICER1 5p:3p Tumours</button>
          <button class="seg-btn" onclick="switchSub('mir-sub', 'mir-pairs', this)">Top miRNA &times; mRNA Pairs</button>
        </div>
        <button class="guide-btn" onclick="openGuideDrawer('mirna')">&#9432; How this works</button>
      </div>

      <div id="mir-bench" class="subpanel mir-sub active">
        <div class="card">
          <div class="card-head">
            <h2>Mature miRNA-seq &amp; RISC Target Benchmarks</h2>
            <span class="meta-note">9,578 miRNA-seq &middot; 192 canonical duplexes</span>
          </div>
          <div class="table-guide">
            <span><strong>Reading Guide:</strong> <strong>Shift</strong> = <code>Median(Group A) &minus; Median(Group B)</code> evaluated by tie-corrected two-sided Wilcoxon&ndash;Mann&ndash;Whitney test &middot; Negative <code>Duplex 5p:3p</code> shift = selective 5p cleavage loss; positive <code>Target Z</code> shift = target mRNA derepression; <code>hsa-miR-451a</code> = DICER1-independent AGO2-Slicer control.</span>
            <span class="meta-note">Hover dotted headers for LaTeX equations &amp; citations</span>
          </div>
          <div class="table-wrap">
            <table id="mir-bench-table">
              <thead>
                <tr>
                  <th>Benchmark</th><th data-term="duplex_5p3p">Metric</th><th>Group A (n)</th><th>Median A</th><th>Group B (n)</th><th>Median B</th>
                  <th data-term="rank_shift_p">Shift</th><th data-term="rank_shift_p">Two-Sided p</th>
                </tr>
              </thead>
              <tbody></tbody>
            </table>
          </div>
        </div>
      </div>

      <div id="mir-dicer" class="subpanel mir-sub">
        <div class="card">
          <div class="card-head">
            <div>
              <h2>DICER1-Altered Tumours: Endogenous 5p:3p Strand Collapse &amp; let-7 Derepression</h2>
              <div class="meta-note" style="margin-top:2px;">RNase IIIb cleaves 5p strands while RNase IIIa cleaves 3p strands &middot; Within-hairpin 5p:3p ratios cancel transcriptional noise</div>
            </div>
            <div style="display:flex;gap:8px;align-items:center;">
              <button class="btn-subtle" id="dicer-reset-btn" style="display:none;" onclick="resetDicerSelection()">Show All Tumours</button>
              <button class="btn-action" onclick="exportDicerMirnaCSV()">Export CSV</button>
            </div>
          </div>
          <div id="dicer-svg-plot" class="svg-track" style="margin-bottom:14px;"></div>
          <div class="table-guide">
            <span><strong>Table Key:</strong> <code>Duplex log2 &lt; 0</code> = inversion from normal 5p dominance (+0.475 WT median) to 3p dominance &middot; <code>5p Fraction</code> = global proportion of 5p reads (0&ndash;1) &middot; <code>let-7 Target Z &gt; 0</code> = loss of <em>let-7</em> mRNA silencing in that tumour.</span>
            <span class="meta-note">Click any circle or class label above to filter/jump</span>
          </div>
          <div class="table-wrap">
            <table id="dicer-mir-table">
              <thead>
                <tr>
                  <th>Sample</th><th>Study</th><th data-term="struct_tier">DICER1 Class</th>
                  <th data-term="duplex_5p3p">Duplex log2((5p+1)/(3p+1))</th>
                  <th data-term="frac_5p">5p Fraction</th>
                  <th data-term="let7_target_z">let-7 Target Z</th>
                </tr>
              </thead>
              <tbody></tbody>
            </table>
          </div>
        </div>
      </div>

      <div id="mir-pairs" class="subpanel mir-sub">
        <div class="card">
          <div class="card-head">
            <h2>Top Repressive miRNA &times; mRNA Couplings</h2>
            <span class="meta-note">13,369,752 pairs across N = 9,353 dual tumours</span>
          </div>
          <div class="table-guide">
            <span><strong>Reading Guide:</strong> Ranked by strongest inverse expression coupling (<code>Pearson r &lt; 0</code> on <code>log2(RPM+1)</code> vs <code>log2(RSEM+1)</code>) across 9,353 dual-profiled TCGA tumours, capturing canonical RISC-mediated target mRNA destabilisation and lineage switches.</span>
            <span class="meta-note">Click <code>&#8597;</code> headers to sort</span>
          </div>
          <div class="table-wrap">
            <table id="top-axes-table">
              <thead>
                <tr>
                  <th>Rank</th><th>Mature miRNA</th><th>Target mRNA</th><th>Tumours (n)</th>
                  <th data-term="pearson_couplings">Pearson r</th><th data-term="pearson_couplings">t-Stat</th>
                </tr>
              </thead>
              <tbody></tbody>
            </table>
          </div>
        </div>
      </div>
    </section>

    <!-- TAB 4: CLINICAL & BIOLOGY (CONSOLIDATED CLEANLY) -->
    <section id="tab-clinical" class="panel">
      <div class="view-bar">
        <div class="seg-group">
          <button class="seg-btn active" onclick="switchSub('clin-sub', 'clin-roc', this)">ClinVar ROC-AUC</button>
          <button class="seg-btn" onclick="switchSub('clin-sub', 'clin-msk', this)">MSK-IMPACT (10,945)</button>
          <button class="seg-btn" onclick="switchSub('clin-sub', 'clin-gtex', this)">Germline &amp; GTEx Tissue</button>
          <button class="seg-btn" onclick="switchSub('clin-sub', 'clin-cptac', this)">Proteomics &amp; CRISPR</button>
          <button class="seg-btn" onclick="switchSub('clin-sub', 'clin-surv', this)">Survival &amp; HRD</button>
        </div>
        <button class="guide-btn" onclick="openGuideDrawer('clinical')">&#9432; How this works</button>
      </div>

      <div id="clin-roc" class="subpanel clin-sub active">
        <div class="card">
          <div class="card-head">
            <h2>ClinVar Pathogenic vs Benign Missense Discrimination</h2>
            <span class="meta-note">Non-circular clinical benchmark</span>
          </div>
          <div class="table-guide">
            <span><strong>Reading Guide:</strong> Evaluates zero-shot predictors on independent NCBI ClinVar Pathogenic/LP vs Benign/LB missense variants &middot; <strong>ROC-AUC</strong>: <code>0.50</code> = random, <code>1.00</code> = ideal separation [Hanley&ndash;McNeil 95% CI] &middot; <strong>PR-AUC (Base)</strong>: precision-recall area vs pathogenic prevalence baseline.</span>
            <span class="meta-note"><code>1star+</code> restricts to multi-submitter / expert-reviewed stars &ge;1</span>
          </div>
          <div class="table-wrap">
            <table id="clinvar-roc-table">
              <thead>
                <tr>
                  <th>Subset</th><th>Predictor</th><th>Pathogenic</th><th>Benign</th>
                  <th data-term="roc_pr_auc">ROC-AUC (95% CI)</th>
                  <th data-term="roc_pr_auc">PR-AUC (Base)</th>
                  <th data-term="rank_shift_p">Mann&ndash;Whitney p</th>
                </tr>
              </thead>
              <tbody></tbody>
            </table>
          </div>
        </div>
      </div>

      <div id="clin-msk" class="subpanel clin-sub">
        <div class="card" style="margin-bottom:16px;">
          <div class="card-head">
            <h2>MSK-IMPACT (N = 10,945) vs TCGA (N = 10,195) Replication</h2>
          </div>
          <div class="table-guide">
            <span><strong>Reading Guide:</strong> Compares somatic coding mutation prevalence in prospective clinical sequencing (<strong>MSK-IMPACT</strong>, primary + metastatic) against <strong>TCGA</strong> primary tumours &middot; <strong>3D Pocket %</strong> = mutations hitting active-site residues (&le;6 &Aring;).</span>
            <span class="meta-note">Genes absent from earlier MSK-IMPACT bait sets show lower profiled N</span>
          </div>
          <div class="table-wrap">
            <table id="msk-rep-table">
              <thead>
                <tr>
                  <th>Gene</th><th>MSK N</th><th>MSK Mut %</th><th>TCGA Mut %</th>
                  <th data-term="hypermutator">Non-Hypermut %</th><th>Truncating %</th>
                  <th data-term="struct_tier">3D Pocket %</th><th>Primary / Met %</th>
                </tr>
              </thead>
              <tbody></tbody>
            </table>
          </div>
        </div>
        <div class="card">
          <div class="card-head">
            <h2>MSK-IMPACT DICER1 Cancer-Type Breakdown</h2>
            <span class="meta-note">RNase IIIb metal-binding hotspots: E1705, D1709, G1809, D1810, E1813</span>
          </div>
          <div class="table-wrap" style="max-height:300px;">
            <table id="msk-dicer-table">
              <thead>
                <tr><th>Cancer Type</th><th>N</th><th>DICER1 Mut</th><th data-term="struct_tier">RNase IIIb Hotspot</th><th>Truncating</th></tr>
              </thead>
              <tbody></tbody>
            </table>
          </div>
        </div>
      </div>

      <div id="clin-gtex" class="subpanel clin-sub">
        <div class="card" style="margin-bottom:16px;">
          <div class="card-head">
            <h2>Germline Constraint (gnomAD v4), Localisation (HPA) &amp; Tractability</h2>
          </div>
          <div class="table-guide">
            <span><strong>Reading Guide:</strong> <code>LOEUF &lt; 0.35</code> &amp; <code>pLI &ge; 0.90</code> = severe human germline LoF intolerance (&gt;800k exomes/genomes) &middot; <code>Missense Z &ge; 3.09</code> = missense-depleted &middot; <strong>Specificity (&tau;)</strong>: <code>0</code> = ubiquitous housekeeping, <code>1</code> = tissue-restricted (e.g. testis <code>LIN28A/B &gt; 0.97</code>).</span>
            <span class="meta-note">Open Targets <code>Ligandable</code> = small-molecule binding pocket</span>
          </div>
          <div class="table-wrap">
            <table id="gene-ref-table">
              <thead>
                <tr>
                  <th>Gene</th><th data-term="gnomad_pli">pLI</th><th data-term="gnomad_loeuf">LOEUF</th>
                  <th data-term="gnomad_misz">Missense Z</th><th>Compartment</th>
                  <th data-term="gtex_tpm">Median TPM</th><th>Peak Tissue</th>
                  <th data-term="tau_spec">Specificity (&tau;)</th><th data-term="tractability">Tractability</th><th>ClinVar P / Total</th>
                </tr>
              </thead>
              <tbody></tbody>
            </table>
          </div>
        </div>
        <div class="card">
          <div class="card-head">
            <div>
              <h2>GTEx v10 Healthy Adult Expression (54 Tissues)</h2>
              <div class="meta-note" style="margin-top:2px;">Switch between ubiquitous Core-10 enzymes (&tau; &lt; 0.38) and oncofetal LIN28A/B (&tau; &gt; 0.97)</div>
            </div>
            <select id="gtex-gene" onchange="renderGtexBars()"></select>
          </div>
          <div id="gtex-bars" class="table-wrap" style="max-height:360px;"></div>
        </div>
      </div>

      <div id="clin-cptac" class="subpanel clin-sub">
        <div class="card" style="margin-bottom:16px;">
          <div class="card-head">
            <h2>CPTAC Proteogenomics (N = 1,002) &amp; DepMap CRISPR Essentiality (N = 1,208)</h2>
          </div>
          <div class="table-guide">
            <span><strong>Reading Guide:</strong> <strong>Buffering &Delta;&rho;</strong> = <code>&rho;(CNA,mRNA) &minus; &rho;(CNA,Protein)</code> (&gt;+0.15 indicates post-transcriptional degradation of unbound complex partners) &middot; <strong>Chronos</strong>: <code>0</code> = neutral, <code>&le;&minus;0.5</code> = cell-line dependency, <code>&le;&minus;1.0</code> = pan-essential.</span>
            <span class="meta-note"><span class="leg-item"><i class="dot dot-crimson"></i>Essential</span> vs <span class="leg-item"><i class="dot dot-petrol"></i>Buffered paralog</span></span>
          </div>
          <div class="table-wrap">
            <table id="consortia-table">
              <thead>
                <tr>
                  <th>Gene</th><th>Role</th><th>PDB</th>
                  <th data-term="cptac_rho">CNA&rarr;mRNA &rho;</th><th data-term="cptac_rho">CNA&rarr;Protein &rho;</th>
                  <th data-term="cptac_atten">Buffering &Delta;&rho;</th><th data-term="depmap_chronos">Chronos</th>
                  <th data-term="depmap_chronos">% Dep (&le;&minus;0.5)</th><th>Class</th>
                </tr>
              </thead>
              <tbody></tbody>
            </table>
          </div>
        </div>
        <div class="card">
          <div class="card-head"><h2>CRISPR Collateral Synthetic Lethality &amp; Oncofetal Addiction</h2></div>
          <div class="table-guide">
            <span><strong>Reading Guide:</strong> Negative <strong>&Delta;Chronos</strong> (<code>Median_Altered &minus; Median_Control &lt; 0</code>) indicates stronger CRISPR knockout lethality in the altered cell-line group (e.g. <code>1p34.3 AGO1/3/4</code> loss sensitises cells to <code>AGO2</code> knockout).</span>
          </div>
          <div class="table-wrap">
            <table id="collat-table">
              <thead>
                <tr>
                  <th>Test</th><th>Target</th><th>Hypothesis</th><th>Altered (n, Med)</th><th>Control (n, Med)</th>
                  <th data-term="depmap_chronos">&Delta;Chronos</th><th data-term="rank_shift_p">Two-Sided p</th>
                </tr>
              </thead>
              <tbody></tbody>
            </table>
          </div>
        </div>
      </div>

      <div id="clin-surv" class="subpanel clin-sub">
        <div class="card" style="margin-bottom:16px;">
          <div class="card-head">
            <h2>Stratified Log-Rank &amp; Multivariable Breslow Cox Survival</h2>
          </div>
          <div class="table-guide">
            <span><strong>Reading Guide:</strong> <strong>HR &gt; 1.0</strong> = higher hazard (shorter survival) compared strictly within histological cancer-type strata (<code>OS</code> Overall, <code>DSS</code> Disease-Specific, <code>PFI</code> Progression-Free) &middot; <strong>Multivariable Cox HR</strong> adjusts simultaneously for Age, Sex, AJCC Stage &amp; TMB.</span>
          </div>
          <div class="table-wrap">
            <table id="surv-table">
              <thead>
                <tr>
                  <th>Benchmark</th><th>Cohort</th><th data-term="surv_endpoints">Endpoint</th>
                  <th>Altered</th><th>Control</th>
                  <th data-term="strat_hr">Stratified HR (95% CI)</th><th>Log-Rank p</th>
                  <th data-term="multivar_cox">Multivariable Cox HR</th><th>Cox p</th>
                </tr>
              </thead>
              <tbody></tbody>
            </table>
          </div>
        </div>
        <div class="card">
          <div class="card-head"><h2>Homologous Recombination Deficiency (HRD) &amp; Immune Benchmarks</h2></div>
          <div class="table-guide">
            <span><strong>Reading Guide:</strong> <strong>HRD Score</strong> = <code>LOH + LST + TAI</code> genomic scar sum across 9,899 tumours (<code>&ge;42</code> clinical HRD threshold; focal <code>8q24.3 AGO2 +2</code> co-amplifies with chromosomal instability, median 33 vs 7).</span>
          </div>
          <div class="table-wrap">
            <table id="imm-table">
              <thead>
                <tr>
                  <th>Benchmark</th><th data-term="hrd_score">Metric</th><th>Group A</th><th>Med A</th><th>Group B</th><th>Med B</th>
                  <th data-term="rank_shift_p">Shift</th><th data-term="rank_shift_p">Two-Sided p</th>
                </tr>
              </thead>
              <tbody></tbody>
            </table>
          </div>
        </div>
      </div>
    </section>

    <!-- TAB 5: FORMULAS (LATEX), PATHWAY, DATASETS, API & GLOSSARY -->
    <section id="tab-docs" class="panel">
      <div class="view-bar">
        <div class="seg-group">
          <button class="seg-btn active" id="seg-eq" onclick="switchSub('doc-sub', 'doc-formulas', this)">Equations (LaTeX)</button>
          <button class="seg-btn" id="seg-path" onclick="switchSub('doc-sub', 'doc-pathway', this)">Pathway &amp; Tiers</button>
          <button class="seg-btn" id="seg-ds" onclick="switchSub('doc-sub', 'doc-datasets', this)">Datasets (41)</button>
          <button class="seg-btn" id="seg-api" onclick="switchSub('doc-sub', 'doc-api', this); renderApiSnippet();">REST API &amp; CLI</button>
          <button class="seg-btn" id="seg-gl" onclick="switchSub('doc-sub', 'doc-glossary', this)">Glossary (50)</button>
        </div>
      </div>

      <div id="doc-formulas" class="subpanel doc-sub active">
        <div class="table-guide" style="margin-bottom:14px;">
          <span><strong>Mathematical Methods &amp; Primary Literature Citations:</strong> All 16 equations powering cohort denominators, GISTIC2 filtering, 3D rotamers, AI language models, 5p:3p enzymology, proteomics, and survival models, rendered in KaTeX with clickable peer-reviewed DOIs.</span>
          <span class="meta-note">Direct links: <code>#eq=F1</code> &ndash; <code>#eq=F16</code></span>
        </div>
        <div class="controls" style="margin-bottom:16px;">
          <select id="formula-cat" onchange="renderFormulas()">
            <option value="ALL">All 16 Equations</option>
            <option value="Genomics">Genomics &amp; Cohorts</option>
            <option value="3D Structure">3D Structure &amp; Rotamers</option>
            <option value="AI & Biophysics">AI Language Models &amp; Inverse Folding</option>
            <option value="Small-RNA">Small-RNA Enzymology &amp; Couplings</option>
            <option value="Proteomics & Tissue">Tissue Specificity &amp; Proteomics</option>
            <option value="Biostatistics">Biostatistics, ROC &amp; Survival</option>
          </select>
          <input type="text" class="ctrl-input" id="formula-q" placeholder="Filter by equation, author, or journal (e.g. Cox, Grantham, Cheng, GISTIC)..." oninput="renderFormulas()" style="width:320px;">
        </div>
        <div id="formula-cards-wrap" class="eq-grid"></div>
      </div>

      <div id="doc-pathway" class="subpanel doc-sub">
        <div class="card">
          <div class="card-head">
            <h2>Canonical miRNA Biogenesis Cascade &amp; Nested Gene Tiers</h2>
          </div>
          <div class="pathway-steps">
            <div class="p-card">
              <div style="font-size:11px;color:var(--accent);font-weight:600;">1 &middot; Nucleus</div>
              <div style="font-weight:600;margin:3px 0;">Microprocessor</div>
              <div style="font-size:12px;color:var(--ink-2);"><strong>DROSHA</strong> + <strong>DGCR8</strong> crop pri-miRNAs into ~70-nt pre-miRNAs (<code>9ASM</code>).</div>
            </div>
            <div class="p-card">
              <div style="font-size:11px;color:var(--accent);font-weight:600;">2 &middot; Pore</div>
              <div style="font-weight:600;margin:3px 0;">Nuclear Export</div>
              <div style="font-size:12px;color:var(--ink-2);"><strong>XPO5</strong> exports 2-nt 3&prime;-overhang hairpins with RanGTP (<code>3A6P</code>).</div>
            </div>
            <div class="p-card">
              <div style="font-size:11px;color:var(--accent);font-weight:600;">3 &middot; Cytoplasm</div>
              <div style="font-weight:600;margin:3px 0;">5p / 3p Dicing</div>
              <div style="font-size:12px;color:var(--ink-2);"><strong>DICER1</strong> + <strong>TARBP2 / PRKRA</strong> cleave 5p (IIIb) &amp; 3p (IIIa) strands (<code>7XW2</code>).</div>
            </div>
            <div class="p-card">
              <div style="font-size:11px;color:var(--accent);font-weight:600;">4 &middot; RISC</div>
              <div style="font-weight:600;margin:3px 0;">Argonaute Slicing</div>
              <div style="font-size:12px;color:var(--ink-2);"><strong>AGO2</strong> (8q24.3 Slicer) &amp; 1p34.3 <strong>AGO1/3/4</strong> guide target silencing (<code>4W5N</code>).</div>
            </div>
            <div class="p-card">
              <div style="font-size:11px;color:var(--accent);font-weight:600;">5 &middot; Regulation</div>
              <div style="font-weight:600;margin:3px 0;">Let-7 &amp; Accessories</div>
              <div style="font-size:12px;color:var(--ink-2);"><strong>LIN28A/B</strong> block <em>let-7</em>; <strong>SMAD2/4</strong> &amp; <strong>TP53</strong> modulate pri-miRNAs.</div>
            </div>
          </div>

          <div class="table-wrap" style="margin-bottom:16px;">
            <table>
              <thead><tr><th>Tier</th><th>Genes</th><th>Included Symbols</th><th>Scope</th></tr></thead>
              <tbody>
                <tr><td><strong>Core-10</strong></td><td>10</td><td><code>DROSHA, DGCR8, XPO5, DICER1, TARBP2, PRKRA, AGO1–4</code></td><td>Dedicated maturation &amp; RISC machinery (collapses 1p34.3).</td></tr>
                <tr><td><strong>Extended-12</strong></td><td>12</td><td>Core-10 + <code>LIN28A, LIN28B</code></td><td>Adds direct oncofetal <em>let-7</em> blockade.</td></tr>
                <tr><td><strong>Panel-14</strong></td><td>14</td><td>Extended-12 + <code>SMAD2, SMAD4</code></td><td>Adds 18q21.1 SMAD Microprocessor cofactors.</td></tr>
                <tr><td><strong>All-15</strong></td><td>15</td><td>Panel-14 + <code>TP53</code></td><td>Gene/3D comparison only (excluded from pathway union).</td></tr>
              </tbody>
            </table>
          </div>

          <h3 style="font-size:13.5px;margin:16px 0 8px;">Excluded Pleiotropic / Downstream Candidates</h3>
          <div class="table-wrap">
            <table id="excluded-cands-table">
              <thead><tr><th>Gene</th><th>Locus</th><th>Function</th><th>Exclusion Rationale</th></tr></thead>
              <tbody></tbody>
            </table>
          </div>
        </div>
      </div>

      <div id="doc-datasets" class="subpanel doc-sub">
        <div class="card" style="margin-bottom:16px;">
          <div class="card-head">
            <h2>10 Integrated Public Consortia</h2>
            <span class="meta-note">46 SHA-256-verified raw archives</span>
          </div>
          <div class="table-wrap">
            <table id="consortia-sources-table">
              <thead><tr><th>Consortium</th><th>Scale</th><th>Layers</th><th>Join Key</th></tr></thead>
              <tbody></tbody>
            </table>
          </div>
        </div>
        <div class="card">
          <div class="card-head">
            <h2>41 Derived CSV / CSV.GZ Datasets</h2>
            <input type="text" class="ctrl-input" id="ds-q" placeholder="Filter datasets..." oninput="renderDatasetsCatalogue()" style="width:220px;">
          </div>
          <div class="table-wrap">
            <table id="datasets-cat-table">
              <thead><tr><th>File (<code>data/</code>)</th><th>Rows &times; Cols</th><th>Size</th><th>Key</th><th>Description</th></tr></thead>
              <tbody></tbody>
            </table>
          </div>
        </div>
      </div>

      <div id="doc-api" class="subpanel doc-sub">
        <div class="card" style="margin-bottom:16px;">
          <div class="card-head">
            <div>
              <h2>Static JSON REST API (<code>/api/v2/</code>) &amp; Offline Python CLI</h2>
              <div class="meta-note" style="margin-top:2px;">CORS-open static JSON endpoints hosted directly on GitHub Pages — zero API key or rate limit</div>
            </div>
            <div style="display:flex;gap:8px;align-items:center;">
              <select id="api-gene-sel" onchange="renderApiSnippet()"></select>
              <a id="api-open-link" href="api/v2/genes/DICER1.json" target="_blank" rel="noopener" class="btn-action" style="text-decoration:none;display:inline-flex;align-items:center;height:34px;">Open Gene JSON &#8599;</a>
            </div>
          </div>
          <div class="table-wrap" style="margin-bottom:16px;">
            <table>
              <thead><tr><th>Endpoint</th><th>Scope</th><th>Action</th></tr></thead>
              <tbody>
                <tr><td><code>api/v2/index.json</code></td><td>Atlas manifest, version, DOI, assemblies &amp; endpoint index</td><td><a href="api/v2/index.json" target="_blank" rel="noopener" style="color:var(--accent);text-decoration:none;font-weight:600;">Open JSON &#8599;</a></td></tr>
                <tr><td><code>api/v2/panel.json</code></td><td>All 15 biogenesis proteins with gnomAD v4, GTEx v10, CPTAC, DepMap &amp; PDB metadata</td><td><a href="api/v2/panel.json" target="_blank" rel="noopener" style="color:var(--accent);text-decoration:none;font-weight:600;">Open JSON &#8599;</a></td></tr>
                <tr><td><code>api/v2/pancancer.json</code></td><td>480 gene &times; study alteration frequencies, 32-study pathway burden &amp; sample denominators</td><td><a href="api/v2/pancancer.json" target="_blank" rel="noopener" style="color:var(--accent);text-decoration:none;font-weight:600;">Open JSON &#8599;</a></td></tr>
                <tr><td><code id="api-gene-path">api/v2/genes/DICER1.json</code></td><td>Per-gene GRCh38/GRCh37 CDS exons, TCGA somatic variants, ClinVar P/LP &amp; residue 3D constraint track</td><td><a id="api-gene-row-link" href="api/v2/genes/DICER1.json" target="_blank" rel="noopener" style="color:var(--accent);text-decoration:none;font-weight:600;">Open JSON &#8599;</a></td></tr>
              </tbody>
            </table>
          </div>
          <div style="display:flex;justify-content:space-between;align-items:center;margin-bottom:6px;">
            <strong style="font-size:13px;">Programmatic Quick-Start (curl, Python &amp; Local CLI)</strong>
            <button class="btn-subtle" id="btn-copy-api" onclick="copyApiCode()">Copy Snippet</button>
          </div>
          <pre id="api-code-box" style="margin:0;padding:14px 16px;background:#111315;color:#e8eaed;border-radius:10px;font-family:var(--mono);font-size:12px;line-height:1.55;overflow-x:auto;"></pre>
        </div>
      </div>

      <div id="doc-glossary" class="subpanel doc-sub">
        <div class="table-guide" style="margin-bottom:14px;">
          <span><strong>Interactive Term Glossary:</strong> Every dotted metric header across the portal links to these 50 definitions, biological thresholds, LaTeX equations, and primary literature citations.</span>
        </div>
        <div class="controls" style="margin-bottom:16px;">
          <input type="text" class="ctrl-input" id="gloss-q" placeholder="Filter glossary terms (e.g. LOEUF, GISTIC2, ESM-IF1, Chronos, HRD)..." oninput="renderGlossaryCards()" style="width:360px;">
          <span class="meta-note" id="gloss-count"></span>
        </div>
        <div id="glossary-cards-wrap" class="eq-grid"></div>
      </div>
    </section>

  </div>
</main>

<!-- Light Floating KaTeX Popover -->
<div id="term-popover" class="popover" role="dialog" aria-live="polite">
  <div style="display:flex;justify-content:space-between;align-items:flex-start;">
    <span class="pop-cat" id="tp-cat"></span>
    <button onclick="hideTermPopover(true)" style="border:none;background:none;color:var(--ink-3);cursor:pointer;font-size:15px;line-height:1;">&times;</button>
  </div>
  <div class="pop-title" id="tp-title"></div>
  <div class="pop-desc" id="tp-meaning"></div>
  <div class="pop-math" id="tp-math"></div>
  <div id="tp-vars" class="eq-vars" style="margin-top:0;margin-bottom:8px;"></div>
  <div class="pop-thresh" id="tp-thresh"></div>
  <div class="pop-cite" id="tp-cite"></div>
</div>

<!-- Contextual Right Slide-Over Guide Drawer -->
<div id="drawer-backdrop" class="drawer-backdrop" onclick="closeGuideDrawer()"></div>
<aside id="guide-drawer" class="drawer" aria-label="Section guide">
  <div style="display:flex;justify-content:space-between;align-items:center;margin-bottom:14px;">
    <span class="pop-cat">View Guide</span>
    <button class="btn-subtle" onclick="closeGuideDrawer()">Close</button>
  </div>
  <h3 id="gd-title" style="margin:0 0 4px;font-size:17px;"></h3>
  <p id="gd-sub" style="margin:0 0 14px;font-size:12.5px;color:var(--ink-2);"></p>
  <div id="gd-math" class="eq-math-box" style="margin-bottom:8px;"></div>
  <div id="gd-vars" class="eq-vars" style="margin-top:0;margin-bottom:16px;"></div>
  <ul id="gd-points" style="padding-left:18px;margin:0 0 20px;color:var(--ink-2);font-size:13px;line-height:1.6;"></ul>
  <button class="btn-action" style="width:100%;" onclick="closeGuideDrawer();openFormulaTab();">Browse All 16 LaTeX Equations &rarr;</button>
</aside>

<!-- Clean Home-Page About, Authors, Citation & Access Section -->
<section id="about-atlas" class="about-section">
  <div class="wrap">
    <div class="about-grid">
      <div>
        <div style="font-family:var(--mono);font-size:11px;color:var(--accent);font-weight:600;text-transform:uppercase;letter-spacing:0.05em;margin-bottom:4px;">About the Resource</div>
        <h2 style="margin:0 0 8px;font-size:18px;font-weight:650;letter-spacing:-0.02em;">miRBioGenAtlas (v2.0.2)</h2>
        <p style="margin:0;color:var(--ink-2);font-size:13px;line-height:1.58;">
          <strong>miRBioGenAtlas</strong> is an open-access, multi-consortium reference of the human microRNA biogenesis machinery integrating somatic coding mutations, chromosome-arm-deconfounded copy-number alterations, 213,940-variant 3D structural and language-model saturation scores, and endogenous mature <code>5p:3p</code> small-RNA readouts across 21,140 sequenced human tumours. All 3D side-chain modelling and queries execute locally in your browser with zero tracking or data upload.
        </p>

        <div class="pillars-row">
          <div class="pillar-tile">
            <div style="font-weight:600;font-size:12.5px;margin-bottom:2px;">Reproducible by Design</div>
            <div style="font-size:11.5px;color:var(--ink-2);">Built deterministically from 46 SHA-256-verified raw snapshots with 7 automated integrity tests.</div>
          </div>
          <div class="pillar-tile">
            <div style="font-weight:600;font-size:12.5px;margin-bottom:2px;">Open Dual Licensing</div>
            <div style="font-size:11.5px;color:var(--ink-2);">Derived tables &amp; figures: <strong>CC BY 4.0</strong> · Pipeline, CLI &amp; 3D web explorer: <strong>MIT License</strong>.</div>
          </div>
          <div class="pillar-tile">
            <div style="font-weight:600;font-size:12.5px;margin-bottom:2px;">Instant CLI &amp; Bulk Data</div>
            <div style="font-size:11.5px;color:var(--ink-2);">41 CSV/CSV.GZ tables on Zenodo plus offline Python CLI: <code>query_variant.py DICER1 E1705K</code>.</div>
          </div>
        </div>
      </div>

      <div class="cite-box" style="display:flex;flex-direction:column;justify-content:space-between;gap:12px;">
        <div>
          <div style="display:flex;justify-content:space-between;align-items:center;gap:8px;margin-bottom:6px;">
            <strong style="color:var(--ink);font-size:13px;">How to Cite &amp; Data Access</strong>
            <div style="display:flex;gap:6px;">
              <button class="btn-subtle" id="btn-copy-cite" onclick="copyCitationText('text')">Copy Citation</button>
              <button class="btn-subtle" id="btn-copy-bib" onclick="copyCitationText('bibtex')">Copy BibTeX</button>
            </div>
          </div>
          <div style="background:#ffffff;border:1px solid var(--line-soft);border-radius:8px;padding:10px 12px;font-size:12px;line-height:1.5;color:var(--ink);">
            Gayam, P. K. R., Earny, V. A., Aranjani, J. M. &amp; Mahammad, F. S. <em>miRBioGenAtlas: Pan-cancer alteration landscape of the miRNA biogenesis machinery across 32 TCGA cancer studies</em> (v2.0.2). <strong>Zenodo</strong> <a href="https://doi.org/10.5281/zenodo.24455067" target="_blank" rel="noopener" style="color:var(--accent);text-decoration:none;">https://doi.org/10.5281/zenodo.24455067</a> (2026).
          </div>
        </div>

        <div style="display:flex;flex-wrap:wrap;gap:8px;align-items:center;justify-content:space-between;">
          <div style="display:flex;gap:8px;flex-wrap:wrap;">
            <a href="https://doi.org/10.5281/zenodo.24455067" target="_blank" rel="noopener" class="btn-action" style="text-decoration:none;display:inline-flex;align-items:center;height:31px;font-size:12px;">Download Zenodo Deposit (44 MB)</a>
            <a href="https://github.com/siddugayam/miRBioGenAtlas" target="_blank" rel="noopener" class="btn-subtle" style="text-decoration:none;display:inline-flex;align-items:center;height:31px;">GitHub &amp; CLI</a>
          </div>
          <span style="font-family:var(--mono);font-size:11px;color:var(--ink-3);">Concept DOI: 10.5281/zenodo.24337286</span>
        </div>
      </div>
    </div>

    <!-- Authors & Affiliations -->
    <div style="margin-bottom:8px;display:flex;justify-content:space-between;align-items:baseline;flex-wrap:wrap;gap:8px;">
      <div style="font-size:13px;font-weight:600;color:var(--ink);">Authors &amp; Institutional Affiliations — Manipal Academy of Higher Education (MAHE), Manipal, India</div>
      <a href="mailto:fayaz.shaik@manipal.edu" style="font-size:12px;color:var(--accent);text-decoration:none;">Correspondence: fayaz.shaik@manipal.edu</a>
    </div>

    <div class="authors-grid">
      <div class="author-card">
        <div>
          <div style="font-weight:600;font-size:13px;color:var(--ink);">Gayam Prasanna Kumar Reddy</div>
          <div style="font-size:11.5px;color:var(--ink-2);margin-top:3px;">Dept. of Pharmaceutical Biotechnology, Manipal College of Pharmaceutical Sciences (MCOPS), MAHE</div>
        </div>
        <a class="orcid-link" href="https://orcid.org/0000-0001-5032-364X" target="_blank" rel="noopener">ORCID &middot; 0000-0001-5032-364X</a>
      </div>

      <div class="author-card">
        <div>
          <div style="font-weight:600;font-size:13px;color:var(--ink);">Venkat Abhiram Earny</div>
          <div style="font-size:11.5px;color:var(--ink-2);margin-top:3px;">Dept. of Microbiology, Kasturba Medical College (KMC), MAHE</div>
        </div>
        <a class="orcid-link" href="https://orcid.org/0009-0008-6926-553X" target="_blank" rel="noopener">ORCID &middot; 0009-0008-6926-553X</a>
      </div>

      <div class="author-card">
        <div>
          <div style="font-weight:600;font-size:13px;color:var(--ink);">Jesil Mathew Aranjani</div>
          <div style="font-size:11.5px;color:var(--ink-2);margin-top:3px;">Dept. of Pharmaceutical Biotechnology, Manipal College of Pharmaceutical Sciences (MCOPS), MAHE</div>
        </div>
        <a class="orcid-link" href="https://orcid.org/0000-0002-7666-3252" target="_blank" rel="noopener">ORCID &middot; 0000-0002-7666-3252</a>
      </div>

      <div class="author-card">
        <div>
          <div style="display:flex;justify-content:space-between;align-items:center;gap:6px;">
            <span style="font-weight:600;font-size:13px;color:var(--ink);">Fayaz Shaik Mahammad*</span>
            <span class="ver-pill" style="color:var(--accent);background:var(--accent-soft);">Corresponding</span>
          </div>
          <div style="font-size:11.5px;color:var(--ink-2);margin-top:3px;">Manipal Institute of Technology (MIT), MAHE, Manipal 576104, Karnataka, India</div>
        </div>
        <a class="orcid-link" href="https://orcid.org/0000-0002-1212-7952" target="_blank" rel="noopener">ORCID &middot; 0000-0002-1212-7952</a>
      </div>
    </div>

    <div style="font-size:11.5px;color:var(--ink-3);line-height:1.5;">
      <strong>Consortium Attribution &amp; Research Disclaimer:</strong> Harmonises public open-access data from TCGA PanCancer Atlas / NCI GDC, cBioPortal, MSK-IMPACT, RCSB PDB, AlphaFold DB, AlphaMissense, Meta ESM, NCBI ClinVar, gnomAD v4, GTEx v10, Human Protein Atlas, CPTAC, Broad DepMap, and Open Targets. Provided strictly for scientific research and educational use; not intended for direct clinical diagnostic decision-making without validated laboratory confirmation.
    </div>
  </div>
</section>

<footer>
  <div class="wrap" style="display:flex;justify-content:space-between;flex-wrap:wrap;gap:10px;">
    <span>&copy; 2026 miRBioGenAtlas Authors &middot; Manipal Academy of Higher Education (MAHE), Manipal, India</span>
    <span><a href="https://doi.org/10.5281/zenodo.24455067" target="_blank" rel="noopener" style="color:var(--accent);text-decoration:none;">DOI: 10.5281/zenodo.24455067</a> &middot; <a href="https://github.com/siddugayam/miRBioGenAtlas" target="_blank" rel="noopener" style="color:var(--accent);text-decoration:none;">GitHub</a> &middot; CC BY 4.0 / MIT</span>
  </div>
</footer>

<script>
window.SAT_CACHE = window.SAT_CACHE || {{}};
window.STRUCT3D_CACHE = window.STRUCT3D_CACHE || {{}};
window.RCSB_CACHE = window.RCSB_CACHE || {{}};
const ATLAS = {data_json};
const LOW_N = new Set(["ACC","CHOL","DLBC","KICH","MESO","UCS","UVM"]);
const GENE_SET = new Set(ATLAS.panel.map(p => p.geneSymbol));

function renderLatex(tex, el, displayMode) {{
  if (!el) return;
  if (!tex) {{ el.style.display = 'none'; return; }}
  el.style.display = '';
  if (typeof katex !== 'undefined') {{
    try {{
      katex.render(tex, el, {{ displayMode: !!displayMode, throwOnError: false }});
      return;
    }} catch (e) {{}}
  }}
  el.textContent = tex;
}}

const AA_PROPS = {{
  A: {{three:'ALA', name:'Alanine', cls:'Aliphatic nonpolar', q:0.0, vol:88.6, kd:1.8, hd:0, ha:0, arom:0, metal:0}},
  R: {{three:'ARG', name:'Arginine', cls:'Basic (+)', q:1.0, vol:173.4, kd:-4.5, hd:5, ha:0, arom:0, metal:0}},
  N: {{three:'ASN', name:'Asparagine', cls:'Polar amide', q:0.0, vol:114.1, kd:-3.5, hd:2, ha:1, arom:0, metal:1}},
  D: {{three:'ASP', name:'Aspartate', cls:'Acidic (-)', q:-1.0, vol:111.1, kd:-3.5, hd:0, ha:2, arom:0, metal:1}},
  C: {{three:'CYS', name:'Cysteine', cls:'Thiol', q:0.0, vol:108.5, kd:2.5, hd:1, ha:1, arom:0, metal:1}},
  Q: {{three:'GLN', name:'Glutamine', cls:'Polar amide', q:0.0, vol:143.8, kd:-3.5, hd:2, ha:1, arom:0, metal:1}},
  E: {{three:'GLU', name:'Glutamate', cls:'Acidic (-)', q:-1.0, vol:138.4, kd:-3.5, hd:0, ha:2, arom:0, metal:1}},
  G: {{three:'GLY', name:'Glycine', cls:'Achiral', q:0.0, vol:60.1, kd:-0.4, hd:0, ha:0, arom:0, metal:0}},
  H: {{three:'HIS', name:'Histidine', cls:'Imidazole', q:0.1, vol:153.2, kd:-3.2, hd:1, ha:1, arom:1, metal:1}},
  I: {{three:'ILE', name:'Isoleucine', cls:'Hydrophobic', q:0.0, vol:166.7, kd:4.5, hd:0, ha:0, arom:0, metal:0}},
  L: {{three:'LEU', name:'Leucine', cls:'Hydrophobic', q:0.0, vol:166.7, kd:3.8, hd:0, ha:0, arom:0, metal:0}},
  K: {{three:'LYS', name:'Lysine', cls:'Basic (+)', q:1.0, vol:168.6, kd:-3.9, hd:3, ha:0, arom:0, metal:0}},
  M: {{three:'MET', name:'Methionine', cls:'Thioether', q:0.0, vol:162.9, kd:1.9, hd:0, ha:1, arom:0, metal:0}},
  F: {{three:'PHE', name:'Phenylalanine', cls:'Aromatic', q:0.0, vol:189.9, kd:2.8, hd:0, ha:0, arom:1, metal:0}},
  P: {{three:'PRO', name:'Proline', cls:'Cyclic rigid', q:0.0, vol:112.7, kd:-1.6, hd:0, ha:0, arom:0, metal:0}},
  S: {{three:'SER', name:'Serine', cls:'Polar OH', q:0.0, vol:89.0, kd:-0.8, hd:1, ha:1, arom:0, metal:0}},
  T: {{three:'THR', name:'Threonine', cls:'Polar OH', q:0.0, vol:116.1, kd:-0.7, hd:1, ha:1, arom:0, metal:0}},
  W: {{three:'TRP', name:'Tryptophan', cls:'Indole aromatic', q:0.0, vol:227.8, kd:-0.9, hd:1, ha:0, arom:1, metal:0}},
  Y: {{three:'TYR', name:'Tyrosine', cls:'Phenolic aromatic', q:0.0, vol:193.6, kd:-1.3, hd:1, ha:1, arom:1, metal:0}},
  V: {{three:'VAL', name:'Valine', cls:'Hydrophobic', q:0.0, vol:140.0, kd:4.2, hd:0, ha:0, arom:0, metal:0}}
}};

const THREE_TO_ONE = {{}};
Object.entries(AA_PROPS).forEach(([one, info]) => {{ THREE_TO_ONE[info.three] = one; }});

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

const PROTEIN_DOMAINS = {{
  DROSHA: [
    [1, 215, 'P-rich', '#64748b'], [216, 390, 'RS-rich', '#0284c7'],
    [391, 875, 'CED Central', '#0d9488'], [876, 1056, 'RNase IIIa', '#d97706'],
    [1107, 1233, 'RNase IIIb', '#c5221f'], [1260, 1334, 'dsRBD', '#4f46e5']
  ],
  DGCR8: [
    [276, 498, 'Rhed Heme-Dimer', '#0d9488'], [511, 576, 'dsRBD 1', '#00709b'],
    [620, 684, 'dsRBD 2', '#4f46e5'], [700, 750, 'CTT', '#d97706']
  ],
  XPO5: [
    [26, 106, 'IBN-N', '#00709b'], [120, 1190, 'HEAT Solenoid', '#0d9488']
  ],
  DICER1: [
    [51, 227, 'Helicase N', '#0284c7'], [235, 428, 'HEL2i', '#0d9488'],
    [433, 602, 'Helicase C', '#00709b'], [630, 722, 'DUF283', '#64748b'],
    [752, 890, 'Platform', '#6366f1'], [895, 1042, 'PAZ', '#4f46e5'],
    [1276, 1403, 'RNase IIIa', '#d97706'], [1666, 1824, 'RNase IIIb', '#c5221f'],
    [1849, 1914, 'dsRBD', '#0f766e']
  ],
  TARBP2: [
    [30, 96, 'dsRBD 1', '#00709b'], [159, 226, 'dsRBD 2', '#4f46e5'], [298, 365, 'Medipal / Dicer-BD', '#0d9488']
  ],
  PRKRA: [
    [35, 99, 'dsRBD 1', '#00709b'], [127, 192, 'dsRBD 2', '#4f46e5'], [245, 311, 'Medipal', '#d97706']
  ],
  AGO1: [[25, 165, 'N-term', '#64748b'], [225, 367, 'PAZ', '#4f46e5'], [430, 572, 'MID', '#0d9488'], [578, 820, 'PIWI', '#c5221f']],
  AGO2: [[28, 167, 'N-term', '#64748b'], [227, 369, 'PAZ', '#4f46e5'], [432, 574, 'MID', '#0d9488'], [580, 822, 'PIWI Slicer', '#c5221f']],
  AGO3: [[28, 168, 'N-term', '#64748b'], [228, 370, 'PAZ', '#4f46e5'], [433, 575, 'MID', '#0d9488'], [581, 823, 'PIWI', '#c5221f']],
  AGO4: [[24, 164, 'N-term', '#64748b'], [222, 365, 'PAZ', '#4f46e5'], [428, 573, 'MID', '#0d9488'], [579, 824, 'PIWI', '#00709b']],
  LIN28A: [[39, 112, 'CSD Cold-Shock', '#00709b'], [137, 178, 'CCHC Zn-Knuckles', '#c5221f']],
  LIN28B: [[29, 102, 'CSD Cold-Shock', '#00709b'], [127, 168, 'CCHC Zn-Knuckles', '#c5221f']],
  SMAD2: [[14, 176, 'MH1 DNA-BD', '#00709b'], [177, 271, 'Linker', '#94a3b8'], [272, 467, 'MH2 Interface', '#c5221f']],
  SMAD4: [[18, 142, 'MH1 DNA-BD', '#00709b'], [143, 318, 'Linker', '#94a3b8'], [319, 552, 'MH2 Interface', '#c5221f']],
  TP53: [[1, 61, 'TAD', '#0d9488'], [62, 94, 'PRD', '#64748b'], [95, 292, 'Core DNA-Binding', '#c5221f'], [325, 356, 'Tetramer', '#4f46e5'], [363, 393, 'CTD', '#d97706']]
}};

function domainAtPos(gene, pos) {{
  const doms = PROTEIN_DOMAINS[gene] || [];
  for (let i = 0; i < doms.length; i++) {{
    if (pos >= doms[i][0] && pos <= doms[i][1]) return {{ idx: i, start: doms[i][0], end: doms[i][1], name: doms[i][2], color: doms[i][3] }};
  }}
  return null;
}}

let viewer3d = null;
let spinActive = false;
let showNeighbors = true;
let showSurface3D = false;
let showContact2D = false;
let theatre3D = false;
let suppressHashSync = false;
let dmsCenterPos = 1705;
let active3DState = {{ gene: 'DICER1', pos: 1705, wt: 'E', mut: 'K', lastPdbText: '' }};

function recordRecentVariant(gene, label) {{
  try {{
    const key = 'mirbiogen_recent_v2';
    const cur = JSON.parse(sessionStorage.getItem(key) || '[]');
    const item = `${{gene}} ${{label}}`;
    const next = [item, ...cur.filter(x => x !== item)].slice(0, 5);
    sessionStorage.setItem(key, JSON.stringify(next));
  }} catch (e) {{}}
}}

function showRecentOrSuggestions(val) {{
  if ((val || '').trim()) {{ updateSearchSuggestions(val); return; }}
  try {{
    const rec = JSON.parse(sessionStorage.getItem('mirbiogen_recent_v2') || '[]');
    if (!rec.length) return;
    const menu = document.getElementById('suggest-menu');
    window._SUG_FNS = rec.map(r => () => runGlobalSearch(r));
    menu.innerHTML = rec.map((r, i) => `
      <div class="sug-row" onclick="window._SUG_FNS[${{i}}]();document.getElementById('suggest-menu').classList.remove('open');">
        <span>${{r}}</span>
        <span class="sug-cat">Recent</span>
      </div>
    `).join('');
    menu.classList.add('open');
  }} catch (e) {{}}
}}

function syncUrlHash(customHash) {{
  if (suppressHashSync || !window.history || !window.history.replaceState) return;
  try {{
    window.history.replaceState(null, '', customHash);
  }} catch (e) {{}}
}}

function switchTab(id, btn, skipHash) {{
  document.querySelectorAll('.panel').forEach(p => p.classList.remove('active'));
  document.querySelectorAll('.nav-tab').forEach(b => b.classList.remove('active'));
  const panelEl = document.getElementById(id);
  if (panelEl) panelEl.classList.add('active');
  if (btn) btn.classList.add('active');
  if (id === 'tab-variants') {{
    setTimeout(() => {{
      if (viewer3d) viewer3d.resize();
      refresh3DScene(true);
    }}, 40);
    if (!skipHash) {{
      const {{ gene, wt, pos, mut }} = active3DState;
      syncUrlHash(`#variant=${{gene}}:${{wt}}${{pos}}${{mut}}`);
    }}
  }} else if (!skipHash) {{
    const mapTab = {{ 'tab-heatmap':'pancancer', 'tab-mirna':'mirna', 'tab-clinical':'clinical', 'tab-docs':'docs' }};
    if (mapTab[id]) syncUrlHash(`#tab=${{mapTab[id]}}`);
  }}
}}

function switchSub(groupCls, targetId, btn) {{
  document.querySelectorAll('.' + groupCls).forEach(el => el.classList.remove('active'));
  document.getElementById(targetId).classList.add('active');
  if (btn && btn.parentElement) {{
    btn.parentElement.querySelectorAll('.seg-btn').forEach(b => b.classList.remove('active'));
    btn.classList.add('active');
  }}
}}

function openFormulaTab(eqId) {{
  switchTab('tab-docs', document.getElementById('btn-tab-docs'), true);
  switchSub('doc-sub', 'doc-formulas', document.getElementById('seg-eq'));
  if (eqId) {{
    document.getElementById('formula-cat').value = 'ALL';
    document.getElementById('formula-q').value = '';
    renderFormulas();
    syncUrlHash(`#eq=${{eqId}}`);
    setTimeout(() => {{
      const el = document.getElementById(eqId);
      if (el) el.scrollIntoView({{behavior: 'smooth', block: 'center'}});
    }}, 40);
  }} else {{
    syncUrlHash('#tab=docs');
  }}
}}

function renderFormulaVars(varsList, boxEl, hdrLabel) {{
  if (!boxEl) return;
  if (!Array.isArray(varsList) || !varsList.length) {{
    boxEl.innerHTML = '';
    boxEl.style.display = 'none';
    return;
  }}
  boxEl.style.display = 'flex';
  boxEl.innerHTML = `<div class="eq-vars-hdr">${{hdrLabel || 'Where (Symbols & Variables)'}}</div>` +
    varsList.map(pair => `<div class="eq-var-row"><span class="eq-var-sym"></span><span class="eq-var-desc">${{pair[1] || ''}}</span></div>`).join('');
  boxEl.querySelectorAll('.eq-var-sym').forEach((sp, k) => {{
    renderLatex((varsList[k] && varsList[k][0]) || '', sp, false);
  }});
}}

function openGuideDrawer(key) {{
  const g = ATLAS.guides && ATLAS.guides[key];
  if (!g) return;
  document.getElementById('gd-title').textContent = g.title;
  document.getElementById('gd-sub').textContent = g.subtitle;
  renderLatex(g.latex, document.getElementById('gd-math'), true);
  renderFormulaVars(g.vars, document.getElementById('gd-vars'), 'Symbols in Formula');
  document.getElementById('gd-points').innerHTML = g.points.map(p => `<li style="margin-bottom:8px;">${{p}}</li>`).join('');
  document.getElementById('drawer-backdrop').classList.add('open');
  document.getElementById('guide-drawer').classList.add('open');
}}

function closeGuideDrawer() {{
  document.getElementById('drawer-backdrop').classList.remove('open');
  document.getElementById('guide-drawer').classList.remove('open');
}}

const CODON_TO_AA = {{
  TTT:'F',TTC:'F',TTA:'L',TTG:'L',CTT:'L',CTC:'L',CTA:'L',CTG:'L',
  ATT:'I',ATC:'I',ATA:'I',ATG:'M',GTT:'V',GTC:'V',GTA:'V',GTG:'V',
  TCT:'S',TCC:'S',TCA:'S',TCG:'S',CCT:'P',CCC:'P',CCA:'P',CCG:'P',
  ACT:'T',ACC:'T',ACA:'T',ACG:'T',GCT:'A',GCC:'A',GCA:'A',GCG:'A',
  TAT:'Y',TAC:'Y',TAA:'*',TAG:'*',CAT:'H',CAC:'H',CAA:'Q',CAG:'Q',
  AAT:'N',AAC:'N',AAA:'K',AAG:'K',GAT:'D',GAC:'D',GAA:'E',GAG:'E',
  TGT:'C',TGC:'C',TGA:'*',TGG:'W',CGT:'R',CGC:'R',CGA:'R',CGG:'R',
  AGT:'S',AGC:'S',AGA:'R',AGG:'R',GGT:'G',GGC:'G',GGA:'G',GGG:'G'
}};
const DNA_COMP = {{ A:'T', T:'A', C:'G', G:'C' }};

function mapGenomicSnvToProtein(chromRaw, posNum, refRaw, altRaw, preferredAssembly) {{
  const chr = String(chromRaw).replace(/^chr/i, '').toUpperCase();
  const pos = Number(posNum);
  const ref = String(refRaw || '').toUpperCase();
  const alt = String(altRaw || '').toUpperCase();
  if (!pos || ref.length !== 1 || alt.length !== 1 || !DNA_COMP[ref] || !DNA_COMP[alt]) return null;
  const cdsMap = ATLAS.cds_genomic || {{}};
  const builds = preferredAssembly === 'GRCh37' ? ['GRCh37', 'GRCh38'] : ['GRCh38', 'GRCh37'];

  for (const bld of builds) {{
    const exKey = bld === 'GRCh38' ? 'cds_exons_grch38' : 'cds_exons_grch37';
    const stKey = bld === 'GRCh38' ? 'strand' : 'strand_grch37';
    for (const [gene, gd] of Object.entries(cdsMap)) {{
      if (String(gd.chr).toUpperCase() !== chr) continue;
      const exons = gd[exKey] || [];
      const strand = gd[stKey] || gd.strand || 1;
      let cum = 0;
      for (let i = 0; i < exons.length; i++) {{
        const [exNum, cs, ce] = exons[i];
        if (pos >= cs && pos <= ce) {{
          const off = strand === 1 ? (pos - cs) : (ce - pos);
          const cIdx = cum + off;
          const cRef = strand === 1 ? ref : DNA_COMP[ref];
          const cAlt = strand === 1 ? alt : DNA_COMP[alt];
          const seq = gd.cds_seq || '';
          if (cIdx < 0 || cIdx >= seq.length) break;
          if (seq[cIdx] !== cRef) {{
            return {{ err: `${{bld}} ref allele mismatch at chr${{chr}}:${{pos}} (expected ${{strand===1?seq[cIdx]:DNA_COMP[seq[cIdx]]}}, got ${{ref}})` }};
          }}
          const codStart = Math.floor(cIdx / 3) * 3;
          const wtCod = seq.slice(codStart, codStart + 3);
          const mutArr = wtCod.split('');
          mutArr[cIdx % 3] = cAlt;
          const mutCod = mutArr.join('');
          const wtAA = CODON_TO_AA[wtCod] || 'X';
          const mutAA = CODON_TO_AA[mutCod] || 'X';
          const codonPos = Math.floor(cIdx / 3) + 1;
          return {{
            gene, wt: wtAA, pos: codonPos, alt: mutAA,
            norm: `${{wtAA}}${{codonPos}}${{mutAA}}`,
            cdna: `c.${{cIdx + 1}}${{cRef}}>${{cAlt}}`,
            wtCod, mutCod, exon: exNum, totalExons: gd.total_exons,
            assembly: bld, vcf: `chr${{chr}}:${{pos}}:${{ref}}>${{alt}}`
          }};
        }}
        cum += (ce - cs + 1);
      }}
    }}
  }}
  return null;
}}

function mapCdnaHgvsToProtein(geneSym, cdnaToken) {{
  const gd = (ATLAS.cds_genomic || {{}})[geneSym];
  if (!gd || !gd.cds_seq) return null;
  const m = String(cdnaToken || '').trim().match(/^c\\.(\\d+)([ACGT])>([ACGT])$/i);
  if (!m) return null;
  const cIdx = parseInt(m[1], 10) - 1;
  const cRef = m[2].toUpperCase(), cAlt = m[3].toUpperCase();
  if (cIdx < 0 || cIdx >= gd.cds_seq.length) return {{ err: `c.${{cIdx+1}} outside ${{geneSym}} CDS` }};
  if (gd.cds_seq[cIdx] !== cRef) return {{ err: `${{geneSym}} c.${{cIdx+1}} reference base is ${{gd.cds_seq[cIdx]}} (not ${{cRef}})` }};
  const codStart = Math.floor(cIdx / 3) * 3;
  const wtCod = gd.cds_seq.slice(codStart, codStart + 3);
  const mutArr = wtCod.split(''); mutArr[cIdx % 3] = cAlt;
  const mutCod = mutArr.join('');
  const wtAA = CODON_TO_AA[wtCod] || 'X', mutAA = CODON_TO_AA[mutCod] || 'X';
  const codonPos = Math.floor(cIdx / 3) + 1;
  return {{ gene: geneSym, wt: wtAA, pos: codonPos, alt: mutAA, norm: `${{wtAA}}${{codonPos}}${{mutAA}}`, cdna: `c.${{cIdx+1}}${{cRef}}>${{cAlt}}`, wtCod, mutCod }};
}}

function reverseMapVariantToGenomic(gene, codonPos, wtAA, mutAA) {{
  const gd = (ATLAS.cds_genomic || {{}})[gene];
  if (!gd || !gd.cds_seq) return null;
  const codStart = (codonPos - 1) * 3;
  if (codStart < 0 || codStart + 3 > gd.cds_seq.length) return null;
  const wtCod = gd.cds_seq.slice(codStart, codStart + 3);

  const mapCdsIdxToGenomic = (cIdx, exList, strand) => {{
    let cum = 0;
    for (let i = 0; i < exList.length; i++) {{
      const [exNum, cs, ce] = exList[i];
      const len = ce - cs + 1;
      if (cIdx < cum + len) {{
        const off = cIdx - cum;
        const gPos = strand === 1 ? (cs + off) : (ce - off);
        return {{ exNum, gPos }};
      }}
      cum += len;
    }}
    return null;
  }};

  const gStart = mapCdsIdxToGenomic(codStart, gd.cds_exons_grch38 || [], gd.strand);
  if (!gStart) return null;
  const strandSym = gd.strand === 1 ? '+' : '−';

  if (wtAA === mutAA) {{
    return {{ wtCod, mutCod: wtCod, exon: gStart.exNum, totalExons: gd.total_exons, strandSym, chr: gd.chr, gPos: gStart.gPos, snvVcf: `chr${{gd.chr}}:${{gStart.gPos}}`, cdna: `c.${{codStart+1}}_${{codStart+3}}` }};
  }}

  let bestCod = '', bestDiffs = 99, diffIdx = -1;
  Object.entries(CODON_TO_AA).forEach(([cod, aa]) => {{
    if (aa !== mutAA) return;
    let d = 0, lastI = -1;
    for (let k = 0; k < 3; k++) {{ if (cod[k] !== wtCod[k]) {{ d++; lastI = k; }} }}
    if (d < bestDiffs) {{ bestDiffs = d; bestCod = cod; diffIdx = lastI; }}
  }});

  if (bestDiffs === 1 && diffIdx >= 0) {{
    const cIdx = codStart + diffIdx;
    const hit38 = mapCdsIdxToGenomic(cIdx, gd.cds_exons_grch38 || [], gd.strand);
    const cRef = wtCod[diffIdx], cAlt = bestCod[diffIdx];
    const gRef = gd.strand === 1 ? cRef : DNA_COMP[cRef];
    const gAlt = gd.strand === 1 ? cAlt : DNA_COMP[cAlt];
    return {{
      wtCod, mutCod: bestCod, nDiff: 1, exon: hit38.exNum, totalExons: gd.total_exons,
      strandSym, chr: gd.chr, gPos: hit38.gPos,
      snvVcf: `chr${{gd.chr}}:${{hit38.gPos}}:${{gRef}}>${{gAlt}}`,
      cdna: `c.${{cIdx + 1}}${{cRef}}>${{cAlt}}`
    }};
  }}
  return {{
    wtCod, mutCod: bestCod, nDiff: bestDiffs, exon: gStart.exNum, totalExons: gd.total_exons,
    strandSym, chr: gd.chr, gPos: gStart.gPos,
    snvVcf: `chr${{gd.chr}}:${{gStart.gPos}} (${{bestDiffs}}-nt codon edit)`,
    cdna: `c.${{codStart+1}}_${{codStart+3}}${{wtCod}}>${{bestCod}}`
  }};
}}

function resolveUniversalVariantQuery(rawLine) {{
  const s = (rawLine || '').trim();
  if (!s || s.startsWith('#')) return null;

  // 1. Check rsID or ClinVar VCV/ID
  const low = s.toLowerCase().replace(/^clinvar[:\\s]*/i, 'cv:');
  const idMap = ATLAS.id_lookup || {{}};
  if (idMap[low]) {{
    const [g, pc] = idMap[low].split(':');
    const m = pc.match(/^([A-Z])(\\d+)([A-Z*])$/);
    if (m) return {{ raw: s, gene: g, wt: m[1], pos: parseInt(m[2], 10), alt: m[3], norm: pc, note: s }};
  }}

  // 2. Check tab/space multi-column VCF data row: CHROM POS ID REF ALT ...
  const cols = s.split(/\\t+/);
  if (cols.length >= 5 && /^(?:chr)?(?:[1-9]|1\\d|2[0-2]|X|Y)$/i.test(cols[0]) && /^\\d+$/.test(cols[1])) {{
    const altFirst = cols[4].split(',')[0];
    const hit = mapGenomicSnvToProtein(cols[0], parseInt(cols[1], 10), cols[3], altFirst);
    if (!hit) return {{ skipVcfNonPanel: true }};
    if (hit.err) return {{ raw: `${{cols[0]}}:${{cols[1]}}:${{cols[3]}}>${{altFirst}}`, err: hit.err }};
    return {{ raw: hit.vcf, ...hit }};
  }}

  // 3. Check compact VCF / gnomAD / hg19 strings: [hg19:]chr14:95094139:C>T or 14-95094139-C-T or 14:95094139:C:T
  const vcfMatch = s.match(/^(?:(hg19|grch37|hg38|grch38)[:\\s-]+)?(?:chr)?([0-9XY]+)[:\\s-]+(\\d+)[:\\s-]+([ACGT])(?:>|:|-|\\s+)([ACGT])$/i);
  if (vcfMatch) {{
    const pref = (vcfMatch[1] && /19|37/i.test(vcfMatch[1])) ? 'GRCh37' : 'GRCh38';
    const hit = mapGenomicSnvToProtein(vcfMatch[2], parseInt(vcfMatch[3], 10), vcfMatch[4], vcfMatch[5], pref);
    if (!hit) return {{ raw: s, err: 'Coordinate outside canonical CDS exons of the 15-gene panel' }};
    if (hit.err) return {{ raw: s, err: hit.err }};
    return {{ raw: s, ...hit }};
  }}

  // 4. Check GENE + cDNA or Protein (e.g. DICER1 c.5113G>A or DICER1 p.Glu1705Lys or DICER1 E1705K)
  const parts = s.split(/[\\s:]+/);
  if (parts.length >= 2 && GENE_SET.has(parts[0].toUpperCase())) {{
    const gene = parts[0].toUpperCase();
    if (/^c\\./i.test(parts[1])) {{
      const ch = mapCdnaHgvsToProtein(gene, parts[1]);
      if (!ch) return {{ raw: s, gene, err: 'Unrecognised cDNA HGVS (use e.g. c.5113G>A)' }};
      if (ch.err) return {{ raw: s, gene, err: ch.err }};
      return {{ raw: s, ...ch }};
    }}
    const mut = normalizeMutToken(parts[1]);
    const m = mut.match(/^([A-Z])(\\d+)([A-Z*])$/);
    if (!m) return {{ raw: s, gene, err: 'Unrecognised substitution format' }};
    return {{ raw: s, gene, wt: m[1], pos: parseInt(m[2], 10), alt: m[3], norm: `${{m[1]}}${{m[2]}}${{m[3]}}` }};
  }}

  return {{ raw: s, err: 'Use GENE MUT, chr:pos:ref>alt, GENE c.posRef>Alt, or rsID' }};
}}

function normalizeMutToken(rawMut) {{
  let s = (rawMut || '').trim().replace(/^p\\./i, '');
  const m3 = s.match(/^([A-Za-z]{{3}})(\\d+)([A-Za-z]{{3}}|\\*)$/);
  if (m3) {{
    const w1 = THREE_TO_ONE[m3[1].toUpperCase()];
    const a1 = m3[3] === '*' ? '*' : THREE_TO_ONE[m3[3].toUpperCase()];
    if (w1 && a1) return w1 + m3[2] + a1;
  }}
  return s.toUpperCase();
}}

/* Google-style Instant Search Suggestions */
function updateSearchSuggestions(val) {{
  const menu = document.getElementById('suggest-menu');
  const q = (val || '').trim();
  if (!q) {{ menu.classList.remove('open'); return; }}
  const qUp = q.toUpperCase();
  const qLow = q.toLowerCase();
  const hits = [];

  const univ = resolveUniversalVariantQuery(q);
  if (univ && univ.gene && univ.norm && !univ.err) {{
    const extra = univ.cdna ? ` · ${{univ.cdna}}` : '';
    const bld = univ.assembly ? ` [${{univ.assembly}}]` : '';
    hits.push({{
      label: `${{univ.gene}} ${{univ.norm}}${{extra}}${{bld}}`,
      cat: 'Resolved 3D Variant',
      fn: () => openVariantIn3D(univ.gene, univ.norm)
    }});
  }}

  ATLAS.panel.forEach(p => {{
    if (p.geneSymbol.includes(qUp) || p.pathwayStep.toLowerCase().includes(qLow)) {{
      hits.push({{ label: `${{p.geneSymbol}} · ${{p.pathwayStep}} (${{p.cytoband}})`, cat: 'Gene 3D', fn: () => loadPresetMutation(p.geneSymbol, '') }});
    }}
  }});

  ATLAS.variants.slice(0, 400).forEach(v => {{
    if (hits.length < 8 && (`${{v.g}} ${{v.p}}`).toUpperCase().includes(qUp)) {{
      hits.push({{ label: `${{v.g}} ${{v.p}} (TCGA n=${{v.n}})`, cat: 'Variant', fn: () => openVariantIn3D(v.g, v.p.replace(/^p\\./i,'')) }});
    }}
  }});

  (ATLAS.formulas || []).forEach(f => {{
    if (hits.length < 9 && (`${{f.id}} ${{f.title}} ${{f.summary}}`).toLowerCase().includes(qLow)) {{
      hits.push({{ label: `${{f.id}}: ${{f.title}}`, cat: 'LaTeX Formula', fn: () => openFormulaTab(f.id) }});
    }}
  }});

  Object.values(ATLAS.glossary || {{}}).forEach(g => {{
    if (hits.length < 10 && (`${{g.title}} ${{g.meaning}}`).toLowerCase().includes(qLow)) {{
      hits.push({{
        label: g.title,
        cat: 'Glossary',
        fn: () => {{
          switchTab('tab-docs', document.getElementById('btn-tab-docs'));
          switchSub('doc-sub', 'doc-glossary', document.getElementById('seg-gl'));
          document.getElementById('gloss-q').value = q;
          renderGlossaryCards();
        }}
      }});
    }}
  }});

  if (!hits.length) {{ menu.classList.remove('open'); return; }}
  window._SUG_FNS = hits.map(h => h.fn);
  menu.innerHTML = hits.slice(0, 7).map((h, i) => `
    <div class="sug-row" onclick="window._SUG_FNS[${{i}}]();document.getElementById('suggest-menu').classList.remove('open');">
      <span>${{h.label}}</span>
      <span class="sug-cat">${{h.cat}}</span>
    </div>
  `).join('');
  menu.classList.add('open');
}}

function handleSearchKey(e) {{
  if (e.key === 'Escape') {{
    document.getElementById('suggest-menu').classList.remove('open');
    return;
  }}
  if (e.key === 'Enter') {{
    document.getElementById('suggest-menu').classList.remove('open');
    runGlobalSearch(e.target.value);
  }}
}}

function runGlobalSearch(raw) {{
  const trimmed = (raw || '').trim();
  if (!trimmed) return;
  const univ = resolveUniversalVariantQuery(trimmed);
  if (univ && univ.gene && univ.norm && !univ.err) {{
    loadPresetMutation(univ.gene, univ.norm);
    return;
  }}
  if (/^(api|rest|json|cli)$/i.test(trimmed)) {{
    openApiDocs();
    return;
  }}
  const parts = trimmed.split(/[\\s:]+/);
  if (parts.length >= 2 && GENE_SET.has(parts[0].toUpperCase())) {{
    loadPresetMutation(parts[0].toUpperCase(), normalizeMutToken(parts[1]));
    return;
  }}
  const q = normalizeMutToken(trimmed);
  if (GENE_SET.has(q)) {{
    loadPresetMutation(q, '');
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
    switchTab('tab-heatmap', document.getElementById('btn-tab-heatmap'));
    switchSub('pc-sub', 'pc-codes', document.querySelectorAll('#tab-heatmap .seg-btn')[2]);
    return;
  }}
  const qLow = trimmed.toLowerCase();
  const fHit = (ATLAS.formulas || []).find(f => (`${{f.id}} ${{f.title}}`).toLowerCase().includes(qLow));
  if (fHit) {{
    openFormulaTab(fHit.id);
    return;
  }}
  const gHit = Object.values(ATLAS.glossary || {{}}).find(g => (`${{g.title}} ${{g.cat}}`).toLowerCase().includes(qLow));
  if (gHit) {{
    switchTab('tab-docs', document.getElementById('btn-tab-docs'));
    switchSub('doc-sub', 'doc-glossary', document.getElementById('seg-gl'));
    document.getElementById('gloss-q').value = trimmed;
    renderGlossaryCards();
    return;
  }}
  switchTab('tab-variants', document.getElementById('btn-tab-variants'));
  switchSub('exp-sub', 'exp-tcga', document.getElementById('seg-exp-tcga'));
  document.getElementById('var-q').value = trimmed;
  renderVariants();
}}

/* Light KaTeX Popover Logic with Primary Citations */
const FORMULA_BY_ID = Object.fromEntries((ATLAS.formulas || []).map(f => [f.id, f]));
let tpPinned = false;
let tpTimer = null;
function showTermPopover(targetEl, pin) {{
  const key = targetEl.getAttribute('data-term');
  const item = ATLAS.glossary && ATLAS.glossary[key];
  if (!item) return;
  clearTimeout(tpTimer);
  if (pin) tpPinned = true;

  const pop = document.getElementById('term-popover');
  document.getElementById('tp-cat').textContent = item.cat || '';
  document.getElementById('tp-title').textContent = item.title || key;
  document.getElementById('tp-meaning').textContent = item.meaning || '';
  renderLatex(item.latex || '', document.getElementById('tp-math'), false);
  renderFormulaVars(item.vars, document.getElementById('tp-vars'), 'Symbols in Formula');
  document.getElementById('tp-thresh').textContent = item.thresholds || '';

  const citeEl = document.getElementById('tp-cite');
  const fRef = item.anchor ? FORMULA_BY_ID[item.anchor] : null;
  if (citeEl) {{
    if (fRef && fRef.refs && fRef.refs.length) {{
      const links = fRef.refs.slice(0, 2).map(r => `<a href="${{r.url}}" target="_blank" rel="noopener">${{r.cite}} &#8599;</a>`).join(' &middot; ');
      citeEl.innerHTML = `<span>Ref: ${{links}}</span><a href="#eq=${{fRef.id}}" onclick="event.preventDefault();hideTermPopover(true);openFormulaTab('${{fRef.id}}');">${{fRef.id}} &rarr;</a>`;
      citeEl.style.display = 'flex';
    }} else {{
      citeEl.style.display = 'none';
    }}
  }}

  pop.classList.add('visible');
  const rect = targetEl.getBoundingClientRect();
  const pw = pop.offsetWidth || 380;
  const ph = pop.offsetHeight || 200;
  let left = Math.min(Math.max(12, rect.left), window.innerWidth - pw - 14);
  let top = rect.bottom + 8;
  if (top + ph > window.innerHeight - 12 && rect.top - ph - 8 > 8) top = rect.top - ph - 8;
  pop.style.left = left + 'px';
  pop.style.top = top + 'px';
}}

function hideTermPopover(force) {{
  if (tpPinned && !force) return;
  tpPinned = false;
  document.getElementById('term-popover').classList.remove('visible');
}}

document.addEventListener('mouseover', e => {{
  const pop = document.getElementById('term-popover');
  if (pop && pop.contains(e.target)) {{ clearTimeout(tpTimer); return; }}
  const t = e.target.closest('[data-term]');
  if (t && !tpPinned) showTermPopover(t, false);
}});
document.addEventListener('mouseout', e => {{
  const pop = document.getElementById('term-popover');
  if (pop && pop.contains(e.target) && !pop.contains(e.relatedTarget)) {{
    if (!tpPinned) tpTimer = setTimeout(() => hideTermPopover(false), 220);
    return;
  }}
  const t = e.target.closest('[data-term]');
  if (t && !tpPinned) tpTimer = setTimeout(() => hideTermPopover(false), 240);
}});
document.addEventListener('click', e => {{
  const pop = document.getElementById('term-popover');
  const sug = document.getElementById('suggest-menu');
  if (sug && !e.target.closest('.search-shell')) sug.classList.remove('open');
  if (pop && pop.contains(e.target)) return;
  if (e.target.closest('th.sortable')) {{
    hideTermPopover(true);
    return;
  }}
  const t = e.target.closest('[data-term]');
  if (t) {{ e.stopPropagation(); showTermPopover(t, true); return; }}
  if (tpPinned) hideTermPopover(true);
}});
document.addEventListener('keydown', e => {{
  const tag = document.activeElement ? document.activeElement.tagName : '';
  const inField = (tag === 'INPUT' || tag === 'TEXTAREA' || tag === 'SELECT');
  if (e.key === 'Escape') {{ hideTermPopover(true); closeGuideDrawer(); document.getElementById('suggest-menu').classList.remove('open'); }}
  if ((e.metaKey || e.ctrlKey) && e.key.toLowerCase() === 'k') {{
    e.preventDefault();
    const gs = document.getElementById('global-search');
    gs.focus(); gs.select();
    showRecentOrSuggestions(gs.value);
    return;
  }}
  if (!inField && e.key === '/') {{
    e.preventDefault();
    const gs = document.getElementById('global-search');
    gs.focus();
    showRecentOrSuggestions(gs.value);
  }}
  if (!inField && document.getElementById('tab-variants').classList.contains('active') && !e.metaKey && !e.ctrlKey && !e.altKey) {{
    if (e.key === 'ArrowLeft') {{ e.preventDefault(); stepCodon(-1); }}
    else if (e.key === 'ArrowRight') {{ e.preventDefault(); stepCodon(1); }}
  }}
}});

function loadPresetMutation(gene, mutStr) {{
  document.getElementById('lol-gene').value = gene;
  document.getElementById('sat-query').value = normalizeMutToken(mutStr || '');
  switchTab('tab-variants', document.getElementById('btn-tab-variants'), true);
  renderLollipopAndLoadSat();
}}

function openVariantIn3D(gene, mutStr) {{
  loadPresetMutation(gene, mutStr);
  document.getElementById('struct3d-card').scrollIntoView({{behavior: 'smooth', block: 'start'}});
}}

function stepCodon(delta) {{
  const gene = document.getElementById('lol-gene').value;
  const gCache = window.SAT_CACHE ? window.SAT_CACHE[gene] : null;
  if (!gCache) return;
  const maxPos = Object.keys(gCache).length;
  const nextPos = Math.min(maxPos, Math.max(1, (active3DState.pos || 1) + delta));
  if (nextPos === active3DState.pos) return;
  selectCodonFromTrack(nextPos);
}}

function shareCurrentVariant(btnEl) {{
  const {{ gene, wt, pos, mut }} = active3DState;
  const url = `${{window.location.origin}}${{window.location.pathname}}#variant=${{gene}}:${{wt}}${{pos}}${{mut}}`;
  syncUrlHash(`#variant=${{gene}}:${{wt}}${{pos}}${{mut}}`);
  if (navigator.clipboard && navigator.clipboard.writeText) {{
    const orig = btnEl.textContent;
    navigator.clipboard.writeText(url).then(() => {{
      btnEl.textContent = 'Copied Link ✓';
      setTimeout(() => {{ btnEl.textContent = orig; }}, 1600);
    }});
  }}
}}

function heatColor(val, maxScale) {{
  const t = Math.min(1, Math.max(0, val / maxScale));
  const r = Math.round(255 - t * 255);
  const g = Math.round(255 - t * (255 - 88));
  const b = Math.round(255 - t * (255 - 124));
  return `background:rgb(${{r}},${{g}},${{b}});color:${{t > 0.52 ? '#fff' : '#202124'}};`;
}}

function drillDownHeatmapCell(gene, study) {{
  document.getElementById('lol-gene').value = gene;
  document.getElementById('sat-query').value = '';
  switchTab('tab-variants', document.getElementById('btn-tab-variants'));
  renderLollipopAndLoadSat();
  switchSub('exp-sub', 'exp-tcga', document.getElementById('seg-exp-tcga'));
  document.getElementById('var-gene').value = gene;
  document.getElementById('var-tier').value = 'ALL';
  document.getElementById('var-mech').value = 'ALL';
  document.getElementById('var-q').value = study;
  renderVariants();
  setTimeout(() => {{
    document.getElementById('exp-tcga').scrollIntoView({{behavior: 'smooth', block: 'start'}});
  }}, 60);
}}

function renderHeatmap() {{
  const metric = document.getElementById('hm-metric').value;
  const maxScale = (metric === 'mut' || metric === 'nhm') ? 0.12 : (metric === 'any' ? 0.85 : 0.40);
  const legEl = document.getElementById('hm-legend');
  if (legEl) {{
    legEl.innerHTML = `<span>0.00 (0%)</span><i class="hm-ramp"></i><span>${{maxScale.toFixed(2)}} (${{Math.round(maxScale * 100)}}%+)</span>`;
  }}
  const hmNotes = {{
    hc: '<strong>Strict High-Confidence (Default):</strong> Proportion (<code>0.00&ndash;1.00</code>) with nonsynonymous somatic mutation OR focal high-level <code>GISTIC2 &plusmn;2</code> CNA (strips shallow <code>&plusmn;1</code> whole-arm aneuploidy, which inflates raw calls 5.3-fold).',
    fcna: '<strong>Focal High-Level CNA Only (<code>GISTIC2 &plusmn;2</code>):</strong> Proportion of tumours with deep homozygous deletion (<code>&minus;2</code>) or high-level focal amplification (<code>+2</code>), excluding point mutations and arm-level <code>&plusmn;1</code> shifts.',
    mut: '<strong>Somatic Coding Mutation (All Matched):</strong> Proportion of tumours carrying &ge;1 nonsynonymous SNV, indel, or splice-site mutation across <code>N = 10,195</code> matched samples.',
    nhm: '<strong>Deconfounded Somatic Mutation (Non-Hypermutated):</strong> Restricts each study to non-POLE/MSI tumours (<code>TMB &le; 10 mut/Mb</code> &amp; <code>&lt;500</code> coding mutations) to isolate driver selection from hypermutator passenger noise.',
    any: '<strong>Broad Any Alteration:</strong> Nonsynonymous mutation OR any <code>GISTIC2 &ne; 0</code> call (includes shallow <code>&plusmn;1</code> chromosome-arm gains/losses &mdash; shown to illustrate aneuploidy inflation).'
  }};
  const expEl = document.getElementById('hm-explainer');
  if (expEl) {{
    expEl.innerHTML = `<span>${{hmNotes[metric] || hmNotes.hc}}</span><span class="meta-note"><code>0.10</code> = 10% of cohort &middot; <code>*</code> = Low-N cohort (N &lt; 100) &middot; Click any cell &rarr; inspect variants</span>`;
  }}
  const studies = [...new Set(ATLAS.long.map(r => r.t))].sort();
  const genes = ATLAS.panel.map(p => p.geneSymbol);
  const map = {{}};
  ATLAS.long.forEach(r => {{ map[r.g + '|' + r.t] = r; }});

  let html = '<table class="no-auto-sort"><thead><tr><th>Gene</th>';
  studies.forEach(s => {{ html += `<th style="text-align:center">${{s}}${{LOW_N.has(s)?'*':''}}</th>`; }});
  html += '</tr></thead><tbody>';
  genes.forEach(g => {{
    html += `<tr><td><strong>${{g}}</strong></td>`;
    studies.forEach(s => {{
      const r = map[g + '|' + s];
      const v = r ? r[metric] : 0;
      const tip = r ? `${{g}} in ${{s}} (N=${{r.N}}): ${{(v*100).toFixed(1)}}% (${{v.toFixed(3)}}) | mut=${{r.mut_n}}, +2=${{r.amp_n}}, -2=${{r.dd_n}} — Click to inspect variants` : '';
      html += `<td class="hm-cell" style="${{heatColor(v, maxScale)}}" title="${{tip}}" onclick="drillDownHeatmapCell('${{g}}','${{s}}')">${{v.toFixed(2)}}</td>`;
    }});
    html += '</tr>';
  }});
  html += '</tbody></table>';
  document.getElementById('hm-container').innerHTML = html;
}}

function ensureSatLoaded(gene, cb) {{
  if (window.SAT_CACHE && window.SAT_CACHE[gene]) {{
    const nCod = Object.keys(window.SAT_CACHE[gene]).length;
    document.getElementById('sat-status').textContent = `${{nCod.toLocaleString()}} aa · ${{(nCod * 19).toLocaleString()}} variants`;
    if (cb) cb();
    return;
  }}
  document.getElementById('sat-status').textContent = `Loading ${{gene}}...`;
  const sc = document.createElement('script');
  sc.src = `sat/${{gene}}.js`;
  sc.onload = () => {{
    const nCod = window.SAT_CACHE[gene] ? Object.keys(window.SAT_CACHE[gene]).length : 0;
    document.getElementById('sat-status').textContent = `${{nCod.toLocaleString()}} aa · ${{(nCod * 19).toLocaleString()}} variants`;
    if (cb) cb();
  }};
  sc.onerror = () => {{ if (cb) cb(); }};
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
  const byRes = {{}}, caList = [];
  st.a.forEach(row => {{
    const [pos, aname, x, y, z, plddt] = row;
    if (!byRes[pos]) byRes[pos] = [];
    byRes[pos].push(row);
    if (aname === 'CA') caList.push([pos, x, y, z, plddt]);
  }});
  st._byRes = byRes; st._caList = caList; st._indexed = true;
}}

function selectCodonFromTrack(pos) {{
  const gene = document.getElementById('lol-gene').value;
  const gCache = window.SAT_CACHE ? window.SAT_CACHE[gene] : null;
  let q = String(pos);
  if (gCache && gCache[q]) {{
    const obs = gCache[q].s.slice().sort((a,b) => (b[5] - a[5]) || ((b[1]||0) - (a[1]||0)))[0];
    if (obs) q = gCache[q].m[0] + pos + obs[0];
  }}
  document.getElementById('sat-query').value = q;
  inspectSaturationCodon();
}}

function selectMutAA(pos, wt, mutAA) {{
  document.getElementById('sat-query').value = wt + String(pos) + mutAA;
  inspectSaturationCodon();
}}

function focusDomainIn3D(gene, startPos, endPos) {{
  dmsCenterPos = Math.round((startPos + endPos) / 2);
  if (document.getElementById('exp-dms').classList.contains('active')) renderDmsHeatmap();
  if (viewer3d && active3DState.gene === gene) {{
    const span = [];
    for (let p = startPos; p <= endPos; p += 2) span.push(p);
    viewer3d.zoomTo({{model: 0, resi: span}}, 380);
  }}
}}

function syncGeneFromPublicBar(gene) {{
  if (!gene) return;
  renderGeneXrefStrip(gene);
  const lolSel = document.getElementById('lol-gene');
  if (lolSel && lolSel.value !== gene) {{
    lolSel.value = gene;
    document.getElementById('sat-query').value = '';
    renderLollipopAndLoadSat();
  }}
}}

function renderGeneXrefStrip(gene) {{
  const p = ATLAS.panel.find(x => x.geneSymbol === gene);
  if (!p) return;
  const pubSel = document.getElementById('pubbar-gene');
  if (pubSel && pubSel.value !== gene) pubSel.value = gene;

  const loeuf = p.gnomad_LOEUF !== '' ? Number(p.gnomad_LOEUF) : null;
  const lBg = loeuf === null ? '#f3f4f6' : (loeuf < 0.35 ? 'rgba(197,34,31,0.10)' : (loeuf < 0.60 ? 'rgba(176,96,0,0.12)' : 'rgba(0,112,155,0.09)'));
  const lCol = loeuf === null ? '#5f6368' : (loeuf < 0.35 ? '#c5221f' : (loeuf < 0.60 ? '#b06000' : '#00709b'));
  const lTxt = loeuf === null ? 'LOEUF —' : `LOEUF ${{loeuf.toFixed(2)}} · ${{loeuf < 0.35 ? 'Constrained' : (loeuf < 0.60 ? 'Moderate' : 'Tolerant')}}`;
  const pdb1 = (p.pdbComplexIds || '').split(';')[0].trim();
  const smBadge = p.opentargets_sm_tractable === '1'
    ? `<span class="loeuf-pill" data-term="tractability" style="background:rgba(19,115,51,0.10);color:#137333;">SM Ligandable</span>` : '';

  const bioEl = document.getElementById('gene-bio-badges');
  if (bioEl) {{
    bioEl.innerHTML = `<span class="loeuf-pill" data-term="gnomad_loeuf" style="background:${{lBg}};color:${{lCol}};">${{lTxt}}</span>${{smBadge}}`;
  }}

  const pubLinks = document.getElementById('public-datasets-links');
  if (pubLinks) {{
    pubLinks.innerHTML = `
      <a class="xref-pill" href="https://www.uniprot.org/uniprotkb/${{p.uniprotAccession}}" target="_blank" rel="noopener">UniProt ${{p.uniprotAccession}} &#8599;</a>
      <a class="xref-pill" href="https://alphafold.ebi.ac.uk/entry/${{p.uniprotAccession}}" target="_blank" rel="noopener">AlphaFold DB &#8599;</a>
      <a class="xref-pill" href="https://www.ebi.ac.uk/ProtVar/query?search=${{p.uniprotAccession}}" target="_blank" rel="noopener">ProtVar &#8599;</a>
      <a class="xref-pill" href="https://gnomad.broadinstitute.org/gene/${{p.ensemblGeneId}}?dataset=gnomad_r4" target="_blank" rel="noopener">gnomAD v4 &#8599;</a>
      <a class="xref-pill" href="https://platform.opentargets.org/target/${{p.ensemblGeneId}}" target="_blank" rel="noopener">Open Targets &#8599;</a>
      <a class="xref-pill" href="https://www.cbioportal.org/results/mutations?cancer_study_list=pancan_pcawg_2020&amp;gene_list=${{gene}}" target="_blank" rel="noopener">cBioPortal &#8599;</a>
      ${{pdb1 ? `<a class="xref-pill" href="https://www.rcsb.org/structure/${{pdb1}}" target="_blank" rel="noopener">RCSB ${{pdb1}} &#8599;</a>` : ''}}
      <a class="xref-pill" href="https://www.ncbi.nlm.nih.gov/clinvar/?term=${{gene}}%5Bgene%5D" target="_blank" rel="noopener">NCBI ClinVar &#8599;</a>
      <a class="xref-pill" href="https://gtexportal.org/home/gene/${{gene}}" target="_blank" rel="noopener">GTEx v10 &#8599;</a>
      <a class="xref-pill" href="https://doi.org/10.5281/zenodo.24455067" target="_blank" rel="noopener">Zenodo Deposit &#8599;</a>
    `;
  }}
}}

function downloadLollipopSVG() {{
  const svgEl = document.querySelector('#lollipop-svg-wrap svg');
  if (!svgEl) return;
  const src = '<?xml version="1.0" encoding="UTF-8"?>\\n' + svgEl.outerHTML;
  const a = document.createElement('a');
  a.href = URL.createObjectURL(new Blob([src], {{type: 'image/svg+xml;charset=utf-8'}}));
  a.download = `${{active3DState.gene}}_domain_lollipop_track.svg`;
  a.click();
}}

function updateLollipopActiveMarker(pos) {{
  const marker = document.getElementById('lol-active-needle');
  const label = document.getElementById('lol-active-lbl');
  const gene = document.getElementById('lol-gene').value;
  const tr = ATLAS.tracks[gene] || [];
  if (!marker || !tr.length) return;
  const L = tr.length, W = 1420, padL = 94, padR = 26;
  const x = padL + ((Math.min(L, Math.max(1, pos)) - 1) / Math.max(1, L - 1)) * (W - padL - padR);
  marker.setAttribute('x1', x.toFixed(1));
  marker.setAttribute('x2', x.toFixed(1));
  if (label) {{
    label.setAttribute('x', Math.min(W - 46, Math.max(padL + 24, x)).toFixed(1));
    label.textContent = `aa ${{pos}}`;
  }}
}}

function renderLollipopAndLoadSat() {{
  const gene = document.getElementById('lol-gene').value;
  renderGeneXrefStrip(gene);
  const domSel = document.getElementById('dms-domain-sel');
  if (domSel) {{
    const doms = PROTEIN_DOMAINS[gene] || [];
    domSel.innerHTML = '<option value="">Jump to domain...</option>' +
      doms.map(d => `<option value="${{d[0]}}-${{d[1]}}">${{d[2]}} (${{d[0]}}–${{d[1]}})</option>`).join('');
  }}

  const tr = ATLAS.tracks[gene] || [];
  if (!tr.length) return;
  const L = tr.length;
  const W = 1420, H = 206, padL = 94, padR = 26;
  const plotW = W - padL - padR;
  const xOf = pos => padL + ((pos - 1) / Math.max(1, L - 1)) * plotW;

  let maxCnt = 1;
  tr.forEach(c => {{ if (c[7] + c[8] > maxCnt) maxCnt = c[7] + c[8]; }});

  let domCapsules = '';
  (PROTEIN_DOMAINS[gene] || []).forEach(([ds, de, dname, dcol]) => {{
    const x1 = xOf(ds), x2 = xOf(Math.min(L, de));
    const dw = Math.max(6, x2 - x1);
    const showTxt = dw >= 48;
    const shortLbl = dw < 76 && dname.length > 9 ? dname.slice(0, 8) + '…' : dname;
    domCapsules += `<g style="cursor:pointer" onclick="focusDomainIn3D('${{gene}}',${{ds}},${{de}})">
      <rect x="${{x1.toFixed(1)}}" y="112" width="${{dw.toFixed(1)}}" height="20" rx="5" fill="${{dcol}}" fill-opacity="0.88">
        <title>${{gene}} ${{dname}} domain (${{ds}}–${{de}} aa) — Click to focus domain in 3D</title>
      </rect>
      ${{showTxt ? `<text x="${{(x1 + dw/2).toFixed(1)}}" y="125.5" font-size="10" font-weight="600" fill="#ffffff" text-anchor="middle" pointer-events="none">${{shortLbl}}</text>` : ''}}
    </g>`;
  }});

  let pockets = '', amPath = '', stems = '';
  tr.forEach(c => {{
    const [pos, wt, plddt, pkt, am, e2, if1, tcgaN, cvN, pdb] = c;
    const x = xOf(pos);
    if (pkt === 1) {{
      pockets += `<rect x="${{x-1.3}}" y="107" width="2.6" height="5" fill="#c5221f" opacity="0.9" style="cursor:pointer" onclick="selectCodonFromTrack(${{pos}})"><title>${{gene}} ${{wt}}${{pos}} · 3D active pocket &le;6Å (${{pdb}})</title></rect>`;
    }}
    if (am !== null) {{
      const yAm = 182 - am * 28;
      amPath += (amPath ? ' L' : 'M') + x.toFixed(1) + ',' + yAm.toFixed(1);
    }}
    if (tcgaN > 0 || cvN > 0) {{
      const h = 14 + (Math.log2(tcgaN + cvN + 1) / Math.log2(maxCnt + 1)) * 78;
      const yTop = 105 - h;
      const col = pkt === 1 ? '#c5221f' : (cvN > 0 ? '#00709b' : '#5f6368');
      const r = Math.min(6.4, 2.8 + Math.log2(tcgaN + cvN + 1));
      stems += `<g style="cursor:pointer" onclick="selectCodonFromTrack(${{pos}})">
        <line x1="${{x.toFixed(1)}}" y1="105" x2="${{x.toFixed(1)}}" y2="${{yTop.toFixed(1)}}" stroke="${{col}}" stroke-width="1.25" opacity="0.72"/>
        <circle cx="${{x.toFixed(1)}}" cy="${{yTop.toFixed(1)}}" r="${{r.toFixed(1)}}" fill="${{col}}" stroke="#fff" stroke-width="0.8">
          <title>${{gene}} ${{wt}}${{pos}} · TCGA: ${{tcgaN}} · ClinVar P/LP: ${{cvN}}</title>
        </circle>
      </g>`;
    }}
  }});

  let ticks = '';
  const step = L > 1200 ? 200 : (L > 500 ? 100 : 50);
  for (let p = 1; p <= L; p += step) {{
    ticks += `<text x="${{xOf(p).toFixed(1)}}" y="200" font-size="10.5" fill="#80868b" text-anchor="middle">${{p}}</text>`;
  }}
  ticks += `<text x="${{xOf(L).toFixed(1)}}" y="200" font-size="10.5" fill="#80868b" text-anchor="end">${{L}} aa</text>`;

  document.getElementById('lollipop-svg-wrap').innerHTML = `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${{W}} ${{H}}" width="100%" height="206">
    <text x="6" y="24" font-size="10" font-weight="600" fill="#80868b">VARIANTS</text>
    <text x="6" y="126" font-size="10" font-weight="600" fill="#80868b">DOMAINS</text>
    <text x="6" y="174" font-size="10" font-weight="600" fill="#80868b">CONSTRAINT</text>
    <rect x="${{padL}}" y="119" width="${{plotW}}" height="6" rx="3" fill="#dfe3e8"/>
    ${{domCapsules}}${{pockets}}
    <line id="lol-active-needle" x1="${{padL}}" y1="14" x2="${{padL}}" y2="186" stroke="#00709b" stroke-width="1.6" stroke-dasharray="3,2" opacity="0.85"/>
    <text id="lol-active-lbl" x="${{padL + 24}}" y="12" font-size="10" font-family="monospace" font-weight="600" fill="#00709b" text-anchor="middle"></text>
    ${{stems}}
    <path d="${{amPath}}" fill="none" stroke="#00709b" stroke-width="1.15" opacity="0.82"/>
    <line x1="${{padL}}" y1="${{182 - 0.564*28}}" x2="${{W-padR}}" y2="${{182 - 0.564*28}}" stroke="#c5221f" stroke-dasharray="3,3" stroke-width="0.8" opacity="0.65"/>
    ${{ticks}}
  </svg>
  <div class="plot-legend">
    <span class="leg-item"><i class="swatch" style="background:#c5221f;margin:0;"></i>3D Active Pocket (&le;6 &Aring; RNA/Metal)</span>
    <span class="leg-item"><i class="swatch" style="background:#00709b;margin:0;"></i>ClinVar Pathogenic / LP</span>
    <span class="leg-item"><i class="swatch" style="background:#5f6368;margin:0;"></i>TCGA Somatic (stem height = log&#8322; count)</span>
    <span class="leg-item"><i class="leg-box" style="background:#4f46e5;"></i>UniProt Domains (click capsule to zoom 3D)</span>
    <span class="leg-item"><i class="leg-line" style="background:#00709b;"></i>Mean AlphaMissense (<i class="leg-line" style="border-top:2px dashed #c5221f;height:0;width:12px;"></i> 0.564 pathogenic cutoff)</span>
    <span class="leg-item"><span style="color:#00709b;font-weight:700;">&#9478;</span>Active 3D residue</span>
  </div>`;

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
  const kxp = vecCross(axisUnit, p), kdp = vecDot(axisUnit, p);
  return [
    origin[0] + p[0]*c + kxp[0]*s + axisUnit[0]*kdp*(1-c),
    origin[1] + p[1]*c + kxp[1]*s + axisUnit[1]*kdp*(1-c),
    origin[2] + p[2]*c + kxp[2]*s + axisUnit[2]*kdp*(1-c)
  ];
}}

function formatPdbAtom(serial, atomName, res3, chain, resSeq, x, y, z, bfac) {{
  const an = atomName.length < 4 ? (' ' + atomName).padEnd(4, ' ') : atomName.slice(0, 4);
  const el = atomName.replace(/[0-9]/g, '').slice(0, 1);
  return 'ATOM  ' + String(serial).padStart(5, ' ') + ' ' + an + ' ' + res3.padStart(3, ' ') + ' ' + chain +
    String(resSeq).padStart(4, ' ') + '    ' + x.toFixed(3).padStart(8, ' ') + y.toFixed(3).padStart(8, ' ') +
    z.toFixed(3).padStart(8, ' ') + '  1.00' + bfac.toFixed(2).padStart(6, ' ') + '          ' + el.padStart(2, ' ');
}}

function computeMutant3DGeometry(st, pos, wtAA, mutAA) {{
  if (!st || !st._byRes[pos]) return null;
  const resRows = st._byRes[pos], bb = {{}};
  resRows.forEach(r => {{ bb[r[1]] = [r[2], r[3], r[4], r[5]]; }});
  if (!bb.N || !bb.CA || !bb.C) return null;

  const ca = bb.CA, envAtoms = [];
  let caNeighbors8A = 0;
  st._caList.forEach(c => {{
    if (c[0] === pos) return;
    const dCa = Math.hypot(c[1]-ca[0], c[2]-ca[1], c[3]-ca[2]);
    if (dCa <= 8.0) caNeighbors8A++;
    if (dCa <= 12.5) (st._byRes[c[0]] || []).forEach(ar => envAtoms.push([c[0], ar[1], ar[2], ar[3], ar[4]]));
  }});

  const u = vecNorm(vecSub(bb.N, bb.CA));
  const w = vecNorm(vecCross(u, vecSub(bb.C, bb.CA)));
  const v = vecCross(w, u);
  const templates = (ATLAS.rotamers && ATLAS.rotamers[mutAA]) ? ATLAS.rotamers[mutAA] : [[]];
  const spinAngles = (mutAA === 'G' || mutAA === 'A' || mutAA === 'P') ? [0] : [0, -20, 20, -40, 40, 60, -60, 120, -120, 180];

  let bestAtoms = [], bestScore = Infinity, bestMinDist = 99.9, bestClashes = [], bestContacts = [];
  templates.forEach(tpl => {{
    const basePts = tpl.map(([aname, cu, cv, cw]) => [
      aname, ca[0] + cu*u[0] + cv*v[0] + cw*w[0], ca[1] + cu*u[1] + cv*v[1] + cw*w[1], ca[2] + cu*u[2] + cv*v[2] + cw*w[2]
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

      let clashScore = 0.0, minD = 99.9;
      const clashes = [], contacts = [];
      cand.forEach(([aname, x, y, z]) => {{
        envAtoms.forEach(([epos, eaname, ex, ey, ez]) => {{
          if (Math.abs(epos - pos) === 1 && aname === 'CB' && (eaname === 'N' || eaname === 'C')) return;
          const d = Math.hypot(x - ex, y - ey, z - ez);
          if (d < minD) minD = d;
          if (d < 2.40) {{
            clashScore += (2.40 - d) * (2.40 - d) * 25.0;
            clashes.push({{from:[x,y,z], to:[ex,ey,ez], d:d}});
          }} else if (d < 3.0) clashScore += (3.0 - d) * 0.4;
          if (d >= 2.40 && d <= 3.75 && aname !== 'CB') contacts.push({{from:[x,y,z], to:[ex,ey,ez], d:d}});
        }});
      }});
      clashScore += Math.abs(deg) * 0.0005;
      if (clashScore < bestScore) {{
        bestScore = clashScore; bestAtoms = cand; bestMinDist = minD > 50 ? null : minD;
        bestClashes = clashes.sort((a,b)=>a.d-b.d).slice(0, 6);
        bestContacts = contacts.sort((a,b)=>a.d-b.d).slice(0, 5);
      }}
    }});
  }});

  return {{ bb, mutSidechain: bestAtoms, caNeighbors8A, mutMinDist: bestMinDist, clashes: bestClashes, contacts: bestContacts }};
}}

function buildPdbModels(gene, st, pos, wtAA, mutAA, geom, colorMode) {{
  const seq = st.s, trMap = {{}};
  (ATLAS.tracks[gene] || []).forEach(c => {{ trMap[c[0]] = c; }});
  const mode = document.getElementById('mol3d-mode').value;

  const linesWT = [];
  let serial = 1;
  st.a.forEach(([rpos, aname, x, y, z, plddt]) => {{
    if (mode === 'mutant' && rpos === pos && !['N','CA','C','O'].includes(aname)) return;
    const aa1 = seq[rpos - 1] || 'A';
    const res3 = AA_PROPS[aa1] ? AA_PROPS[aa1].three : 'ALA';
    let bval = plddt;
    if (colorMode === 'am') bval = (trMap[rpos] && trMap[rpos][4] !== null) ? trMap[rpos][4] * 100.0 : 20.0;
    else if (colorMode === 'pocket') bval = (trMap[rpos] && trMap[rpos][3] === 1) ? 95.0 : 25.0;
    linesWT.push(formatPdbAtom(serial++, aname, res3, 'A', rpos, x, y, z, bval));
  }});
  linesWT.push('TER', 'END');

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
    geom.mutSidechain.forEach(([aname, mx, my, mz]) => linesMut.push(formatPdbAtom(mser++, aname, mut3, 'M', pos, mx, my, mz, 99.0)));
    linesMut.push('TER', 'END');
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

function update3DLegend(mode, colorMode) {{
  const leg = document.getElementById('mol3d-legend');
  if (!leg) return;
  if (mode === 'rcsb') {{
    leg.innerHTML = `
      <span><strong style="color:#fff;">Experimental PDB:</strong></span>
      <span><i class="swatch" style="background:#64748b;"></i>Protein</span>
      <span><i class="swatch" style="background:#f59e0b;"></i>Bound RNA/DNA</span>
      <span><i class="swatch" style="background:#10b981;"></i>Catalytic Mg&sup2;&#8314;/Zn&sup2;&#8314;</span>
      <span><i class="swatch" style="background:#e11d48;"></i>Selected codon</span>
    `;
    return;
  }}
  let bbKey = '';
  if (colorMode === 'plddt') {{
    bbKey = `<span><strong style="color:#fff;">pLDDT:</strong> <i class="swatch" style="background:#0053d6;"></i>&ge;90 <i class="swatch" style="background:#65cbf3;margin-left:5px;"></i>70–90 <i class="swatch" style="background:#ffdb13;margin-left:5px;"></i>50–70 <i class="swatch" style="background:#ff7d45;margin-left:5px;"></i>&lt;50</span>`;
  }} else if (colorMode === 'domain') {{
    bbKey = `<span><strong style="color:#fff;">Backbone:</strong> UniProt Domains (<i class="swatch" style="background:#475569;margin-left:3px;"></i>Inter-domain)</span>`;
  }} else if (colorMode === 'am') {{
    bbKey = `<span data-term="alphamissense"><strong style="color:#fff;">AM:</strong> <i class="swatch" style="background:#c5221f;"></i>&ge;0.564 Pathogenic <i class="swatch" style="background:#d97706;margin-left:5px;"></i>0.34–0.564 <i class="swatch" style="background:#00709b;margin-left:5px;"></i>&le;0.34 Benign</span>`;
  }} else {{
    bbKey = `<span data-term="dist_na"><strong style="color:#fff;">Pocket:</strong> <i class="swatch" style="background:#e11d48;"></i>&le;6&Aring; Active site <i class="swatch" style="background:#64748b;margin-left:5px;"></i>&gt;6&Aring; Scaffold</span>`;
  }}
  leg.innerHTML = `
    ${{bbKey}}
    <span style="opacity:0.4;">|</span>
    <span><i class="swatch" style="background:#0095c8;"></i>WT</span>
    <span data-term="rotamer_clash"><i class="swatch" style="background:#e11d48;"></i>Mutant rotamer</span>
    <span data-term="burial_8a"><i class="swatch" style="background:#94a3b8;"></i>8&Aring; shell</span>
    <span data-term="dist_na"><i class="swatch" style="background:#f59e0b;"></i>Contact (2.4–3.8&Aring;)</span>
    <span data-term="rotamer_clash"><i class="swatch" style="background:#ff3b30;"></i>Clash (&lt;2.4&Aring;)</span>
  `;
}}

function refresh3DScene(keepCamera) {{
  const v = ensureViewerInitialized();
  if (!v) return;
  const {{ gene, pos, wt, mut }} = active3DState;
  const mode = document.getElementById('mol3d-mode').value;
  const colorMode = document.getElementById('mol3d-color').value;
  update3DLegend(mode, colorMode);

  if (mode === 'rcsb') {{ renderExperimentalRCSB(gene, pos, wt, mut); return; }}

  ensureStruct3DLoaded(gene, st => {{
    if (!st) return;
    const domHit = domainAtPos(gene, pos);
    document.getElementById('mol3d-hud').textContent = `${{gene}} · ${{st.u}} (${{st.s.length}} aa)${{domHit ? ' · ' + domHit.name : ''}}`;
    const geom = computeMutant3DGeometry(st, pos, wt, mut);
    const pdbs = buildPdbModels(gene, st, pos, wt, mut, geom, colorMode);
    active3DState.lastPdbText = pdbs.wtPdb + (mode !== 'wt' ? ('\\n' + pdbs.mutPdb) : '');

    v.removeAllModels(); v.removeAllShapes(); v.removeAllLabels();
    if (typeof v.removeAllSurfaces === 'function') v.removeAllSurfaces();
    v.addModel(pdbs.wtPdb, 'pdb');
    v.setStyle({{model: 0}}, {{cartoon: {{colorfunc: atom => {{
      if (colorMode === 'domain') {{
        const dm = domainAtPos(gene, atom.resi);
        return dm ? dm.color : '#475569';
      }}
      const b = atom.b;
      if (colorMode === 'plddt') return b >= 90 ? '#0053d6' : (b >= 70 ? '#65cbf3' : (b >= 50 ? '#ffdb13' : '#ff7d45'));
      if (colorMode === 'am') return b >= 56.4 ? '#c5221f' : (b >= 34.0 ? '#d97706' : '#00709b');
      return b >= 80 ? '#e11d48' : '#64748b';
    }}, opacity: 0.86}}}});

    let nbrRes = [];
    if (geom && geom.bb.CA) {{
      const [cx, cy, cz] = geom.bb.CA;
      nbrRes = st._caList.filter(c => c[0] !== pos && Math.hypot(c[1]-cx, c[2]-cy, c[3]-cz) <= 8.5).map(c => c[0]);
      if (showNeighbors && nbrRes.length) {{
        v.addStyle({{model: 0, resi: nbrRes}}, {{stick: {{radius: 0.10, colorscheme: 'Jmol', opacity: 0.5}}}});
      }}
      if (showSurface3D && nbrRes.length && typeof $3Dmol !== 'undefined') {{
        v.addSurface($3Dmol.SurfaceType.VDW, {{opacity: 0.32, color: '#cbd5e1'}}, {{model: 0, resi: nbrRes.concat([pos])}});
      }}
    }}
    if (mode === 'overlay' || mode === 'wt') {{
      v.addStyle({{model: 0, resi: pos}}, {{stick: {{radius: 0.21, color: '#0095c8', opacity: mode==='overlay'?0.72:1}}, sphere: {{scale: 0.25, color: '#0095c8', opacity: mode==='overlay'?0.72:1}}}});
    }}
    if ((mode === 'overlay' || mode === 'mutant') && pdbs.mutPdb && wt !== mut) {{
      v.addModel(pdbs.mutPdb, 'pdb');
      v.setStyle({{model: 1}}, {{stick: {{radius: 0.24, color: '#e11d48'}}, sphere: {{scale: 0.28, color: '#e11d48'}}}});
      if (geom) {{
        geom.contacts.forEach(ct => v.addCylinder({{start:{{x:ct.from[0],y:ct.from[1],z:ct.from[2]}}, end:{{x:ct.to[0],y:ct.to[1],z:ct.to[2]}}, radius:0.045, dashed:true, dashLength:0.22, gapLength:0.16, color:'#f59e0b'}}));
        geom.clashes.forEach(cl => v.addCylinder({{start:{{x:cl.from[0],y:cl.from[1],z:cl.from[2]}}, end:{{x:cl.to[0],y:cl.to[1],z:cl.to[2]}}, radius:0.075, dashed:true, dashLength:0.18, gapLength:0.12, color:'#ff3b30'}}));
      }}
    }}
    if (geom && geom.bb.CA) {{
      const [lx, ly, lz] = geom.bb.CA;
      v.addLabel(wt === mut ? `${{gene}} ${{wt}}${{pos}}` : `${{gene}} ${{wt}}${{pos}}${{mut}}`, {{
        position: {{x: lx, y: ly + 1.8, z: lz}}, backgroundColor: '#202124', fontColor: '#fff', fontSize: 12, backgroundOpacity: 0.85
      }});
    }}
    v.setClickable({{model: 0}}, true, atom => {{ if (atom && atom.resi) selectCodonFromTrack(atom.resi); }});
    v.render();
    if (!keepCamera) zoomToMutatedResidue();
    if (showContact2D) renderContactMap2D();
  }});
}}

function renderExperimentalRCSB(gene, pos, wt, mut) {{
  const v = ensureViewerInitialized();
  if (!v) return;
  const gCache = window.SAT_CACHE ? window.SAT_CACHE[gene] : null;
  const entry = gCache ? gCache[String(pos)] : null;
  const pInfo = ATLAS.panel.find(p => p.geneSymbol === gene);
  const pdbRaw = (entry && entry.m[1]) ? entry.m[1].split(':')[0] : ((pInfo && pInfo.pdbComplexIds) ? pInfo.pdbComplexIds.split(';')[0].trim() : '');
  if (!pdbRaw) {{ document.getElementById('mol3d-mode').value = 'overlay'; refresh3DScene(); return; }}

  document.getElementById('mol3d-hud').textContent = `${{gene}} · RCSB PDB ${{pdbRaw}}`;
  const drawCif = cifTxt => {{
    v.removeAllModels(); v.removeAllShapes(); v.removeAllLabels();
    if (typeof v.removeAllSurfaces === 'function') v.removeAllSurfaces();
    v.addModel(cifTxt, 'cif');
    v.setStyle({{}}, {{cartoon: {{color: '#64748b', opacity: 0.82}}}});
    v.addStyle({{resn: ['A','U','G','C','DA','DT','DG','DC']}}, {{cartoon: {{color: '#d97706'}}, stick: {{radius: 0.13, color: '#f59e0b'}}}});
    v.addStyle({{elem: ['MG','ZN','CA','MN']}}, {{sphere: {{radius: 1.3, color: '#10b981'}}}});
    v.addStyle({{resi: pos}}, {{stick: {{radius: 0.26, color: '#e11d48'}}, sphere: {{scale: 0.30, color: '#e11d48'}}}});
    v.addLabel(`PDB ${{pdbRaw}} · ${{gene}} ${{wt}}${{pos}}`, {{resi: pos, backgroundColor: '#202124', fontColor: '#fff', fontSize: 12}});
    v.zoomTo({{resi: pos}}, 400); v.render();
  }};
  if (window.RCSB_CACHE[pdbRaw]) {{ drawCif(window.RCSB_CACHE[pdbRaw]); return; }}
  fetch(`https://files.rcsb.org/download/${{pdbRaw}}.cif`)
    .then(r => {{ if (!r.ok) throw new Error(); return r.text(); }})
    .then(txt => {{ window.RCSB_CACHE[pdbRaw] = txt; drawCif(txt); }})
    .catch(() => {{ document.getElementById('mol3d-mode').value = 'overlay'; refresh3DScene(); }});
}}

function zoomToMutatedResidue() {{ if (viewer3d) {{ const p = active3DState.pos; viewer3d.zoomTo({{model:0, resi:[Math.max(1,p-4),p,p+4]}}, 360); viewer3d.zoom(0.82, 160); }} }}
function resetFullProtein3D() {{ if (viewer3d) viewer3d.zoomTo({{model: 0}}, 400); }}
function toggleNeighbors3D() {{ showNeighbors = !showNeighbors; document.getElementById('btn-nbrs').textContent = `Shell: ${{showNeighbors ? 'On' : 'Off'}}`; refresh3DScene(true); }}
function toggleSurface3D() {{
  showSurface3D = !showSurface3D;
  const btn = document.getElementById('btn-surf');
  if (btn) btn.classList.toggle('active-vp', showSurface3D);
  refresh3DScene(true);
}}
function toggleContactMap2D() {{
  showContact2D = !showContact2D;
  const hud = document.getElementById('contact-hud');
  const btn = document.getElementById('btn-2dmap');
  if (hud) hud.classList.toggle('open', showContact2D);
  if (btn) btn.classList.toggle('active-vp', showContact2D);
  if (showContact2D) renderContactMap2D();
}}
function renderContactMap2D() {{
  const cvs = document.getElementById('contact-canvas');
  if (!cvs) return;
  const {{ gene, pos }} = active3DState;
  const st = window.STRUCT3D_CACHE ? window.STRUCT3D_CACHE[gene] : null;
  if (!st || !st._caList) return;
  const L = st.s.length, N = 192, pad = 6, M = N - pad;
  document.getElementById('contact-max-lbl').textContent = `${{L}} aa`;
  const ctx = cvs.getContext('2d');
  ctx.fillStyle = '#111315';
  ctx.fillRect(0, 0, N, N);

  const caByPos = {{}};
  st._caList.forEach(c => {{ caByPos[c[0]] = c; }});
  const img = ctx.createImageData(M, M);
  const data = img.data;
  for (let py = 0; py < M; py++) {{
    const r1 = 1 + Math.floor((py / M) * L);
    const c1 = caByPos[r1];
    for (let px = 0; px < M; px++) {{
      const r2 = 1 + Math.floor((px / M) * L);
      const c2 = caByPos[r2];
      const off = (py * M + px) * 4;
      if (!c1 || !c2) {{ data[off]=20; data[off+1]=22; data[off+2]=26; data[off+3]=255; continue; }}
      const d = Math.hypot(c1[1]-c2[1], c1[2]-c2[2], c1[3]-c2[3]);
      const conf = Math.min(c1[4], c2[4]) / 100.0;
      if (d <= 24.0) {{
        const t = (1.0 - d / 24.0) * (0.35 + 0.65 * conf);
        data[off]   = Math.round(18 + t * 0);
        data[off+1] = Math.round(26 + t * 165);
        data[off+2] = Math.round(32 + t * 145);
        data[off+3] = 255;
      }} else {{
        data[off]=22; data[off+1]=25; data[off+2]=30; data[off+3]=255;
      }}
    }}
  }}
  ctx.putImageData(img, pad, pad);

  (PROTEIN_DOMAINS[gene] || []).forEach(([ds, de, _, dcol]) => {{
    const a = pad + ((ds - 1) / L) * M, b = pad + (Math.min(L, de) / L) * M;
    ctx.fillStyle = dcol;
    ctx.fillRect(a, 0, Math.max(2, b - a), 4);
    ctx.fillRect(0, a, 4, Math.max(2, b - a));
  }});

  const cp = pad + ((Math.min(L, Math.max(1, pos)) - 1) / L) * M;
  ctx.strokeStyle = '#ff4d6d';
  ctx.lineWidth = 1;
  ctx.beginPath();
  ctx.moveTo(cp, pad); ctx.lineTo(cp, N);
  ctx.moveTo(pad, cp); ctx.lineTo(N, cp);
  ctx.stroke();
}}
function handleContactMapClick(ev) {{
  const cvs = document.getElementById('contact-canvas');
  const st = window.STRUCT3D_CACHE ? window.STRUCT3D_CACHE[active3DState.gene] : null;
  if (!cvs || !st) return;
  const rect = cvs.getBoundingClientRect();
  const x = (ev.clientX - rect.left) * (cvs.width / rect.width);
  const pad = 6, M = cvs.width - pad;
  const targetPos = Math.min(st.s.length, Math.max(1, 1 + Math.round(((x - pad) / M) * (st.s.length - 1))));
  selectCodonFromTrack(targetPos);
}}
function toggleSpin3D() {{ if (viewer3d) {{ spinActive = !spinActive; viewer3d.spin(spinActive ? 'y' : false); document.getElementById('btn-spin').textContent = spinActive ? 'Stop' : 'Spin'; }} }}
function toggleTheatre3D() {{
  theatre3D = !theatre3D;
  const wb = document.getElementById('explorer-workbench');
  const btn = document.getElementById('btn-expand3d');
  if (wb) wb.classList.toggle('theatre', theatre3D);
  if (btn) {{
    btn.textContent = theatre3D ? 'Standard' : 'Expand';
    btn.classList.toggle('active-vp', theatre3D);
  }}
  setTimeout(() => {{
    if (viewer3d) {{ viewer3d.resize(); viewer3d.render(); }}
  }}, 200);
}}
function download3DSnapshotPNG() {{
  if (!viewer3d) return;
  viewer3d.render();
  const uri = viewer3d.pngURI();
  if (!uri) return;
  const {{ gene, wt, pos, mut }} = active3DState;
  const a = document.createElement('a');
  a.href = uri;
  a.download = `${{gene}}_${{wt}}${{pos}}${{mut}}_3D.png`;
  a.click();
}}
function printVariantReport() {{
  if (viewer3d) {{
    viewer3d.render();
    const imgEl = document.getElementById('print-3d-img');
    if (imgEl) imgEl.src = viewer3d.pngURI();
  }}
  const dtEl = document.getElementById('print-report-date');
  if (dtEl) dtEl.textContent = new Date().toISOString().slice(0, 10);
  setTimeout(() => window.print(), 80);
}}
function downloadMutatedPDB() {{
  if (!active3DState.lastPdbText) return;
  const {{ gene, wt, pos, mut }} = active3DState;
  const a = document.createElement('a');
  a.href = URL.createObjectURL(new Blob([active3DState.lastPdbText], {{type: 'chemical/x-pdb'}}));
  a.download = `${{gene}}_${{wt}}${{pos}}${{mut}}.pdb`; a.click();
}}

function buildProtVarRadarSVG(am, e2, if1, gDist, minPocketA, accentCol) {{
  const sAm  = am !== null ? Math.min(1, Math.max(0, am)) : 0;
  const sSeq = e2 !== null ? Math.min(1, Math.max(0, -e2 / 10.0)) : 0;
  const sIf1 = if1 !== null ? Math.min(1, Math.max(0, -if1 / 4.0)) : 0;
  const sChm = Math.min(1, Math.max(0, gDist / 180.0));
  const sPkt = minPocketA !== null ? Math.min(1, Math.max(0, 1.0 - minPocketA / 16.0)) : 0.1;
  const vals = [sAm, sSeq, sIf1, sChm, sPkt];
  const lbls = ['AM', 'Seq', 'Fold', 'Chem', 'Pocket'];
  const cx = 62, cy = 60, R = 36;
  const ptAt = (i, rFrac) => {{
    const ang = -Math.PI / 2 + (i * 2 * Math.PI) / 5;
    return [cx + Math.cos(ang) * R * rFrac, cy + Math.sin(ang) * R * rFrac];
  }};
  let rings = '', axes = '', lblSvg = '';
  [0.5, 1.0].forEach(rf => {{
    const poly = [0,1,2,3,4].map(i => ptAt(i, rf).map(v=>v.toFixed(1)).join(',')).join(' ');
    rings += `<polygon points="${{poly}}" fill="none" stroke="#dadce0" stroke-width="0.8"/>`;
  }});
  lbls.forEach((lb, i) => {{
    const [ax, ay] = ptAt(i, 1.0);
    const [lx, ly] = ptAt(i, 1.32);
    axes += `<line x1="${{cx}}" y1="${{cy}}" x2="${{ax.toFixed(1)}}" y2="${{ay.toFixed(1)}}" stroke="#e2e6ea" stroke-width="0.8"/>`;
    lblSvg += `<text x="${{lx.toFixed(1)}}" y="${{(ly + 3).toFixed(1)}}" font-size="8.5" font-weight="600" fill="#5f6368" text-anchor="middle">${{lb}}</text>`;
  }});
  const dataPoly = vals.map((v, i) => ptAt(i, Math.max(0.08, v)).map(x=>x.toFixed(1)).join(',')).join(' ');
  return `<svg viewBox="0 0 124 116" width="122" height="114" aria-label="5-axis biophysical impact radar">
    ${{rings}}${{axes}}
    <polygon points="${{dataPoly}}" fill="${{accentCol}}" fill-opacity="0.22" stroke="${{accentCol}}" stroke-width="1.6"/>
    ${{lblSvg}}
  </svg>`;
}}

/* Clean Minimal Visual Scorecard */
function renderPropertyDeltaCard(gene, pos, wt, mut, entry, subRow) {{
  const box = document.getElementById('prop-delta-panel');
  const st = window.STRUCT3D_CACHE ? window.STRUCT3D_CACHE[gene] : null;
  const geom = st ? computeMutant3DGeometry(st, pos, wt, mut) : null;
  const wp = AA_PROPS[wt] || AA_PROPS.A, mp = AA_PROPS[mut] || wp;
  const [_, pdb, dna, dmet, dprt, plddt] = entry.m;
  const [alt, am, e2, e1v, if1, tn, cv] = subRow || [mut, null, null, null, null, 0, ''];

  const dq = mp.q - wp.q, dvol = mp.vol - wp.vol, gDist = granthamDist(wt, mut);
  const constrained = (am !== null && am >= 0.564) || (e2 !== null && e2 <= -5.0);
  const regimeTag = wt === mut ? tagDot('dot-emerald', 'Wild-type')
    : (constrained ? ((if1 !== null && if1 > -1.5) ? tagDot('dot-crimson', 'Fold-preserving active site') : tagDot('dot-amber', 'Fold-destabilising'))
                   : tagDot('dot-petrol', 'Neutral'));
  const radarCol = wt === mut ? '#137333' : (constrained ? ((if1 !== null && if1 > -1.5) ? '#c5221f' : '#b06000') : '#00709b');

  const clashCount = geom ? geom.clashes.length : 0;
  const aaCells = entry.s.map(r => {{
    const a = r[0], rAm = r[1], rE2 = r[2], rIf1 = r[4];
    const cns = (rAm !== null && rAm >= 0.564) || (rE2 !== null && rE2 <= -5.0);
    const topCol = cns ? ((rIf1 !== null && rIf1 > -1.5) ? '#c5221f' : '#b06000') : '#00709b';
    return `<button class="aa-cell ${{a===mut?'active':''}}" style="border-top-color:${{topCol}};" onclick="selectMutAA(${{pos}},'${{wt}}','${{a}}')" title="${{wt}}${{pos}}${{a}}">${{a}}</button>`;
  }}).join('');

  const pktDists = [dna, dmet, dprt].filter(x => x !== null);
  const minPktA = pktDists.length ? Math.min(...pktDists) : null;
  const domHit = domainAtPos(gene, pos);
  const domBadge = domHit ? `<span class="tag" style="background:${{domHit.color}}18;color:${{domHit.color}};font-weight:600;">${{domHit.name}}</span>` : '';
  const gInfo = reverseMapVariantToGenomic(gene, pos, wt, mut);
  const genomicBar = gInfo ? `
    <div style="display:flex;justify-content:space-between;align-items:center;gap:6px;flex-wrap:wrap;margin-top:6px;padding:5px 9px;background:var(--canvas);border:1px solid var(--line-soft);border-radius:7px;font-family:var(--mono);font-size:11px;color:var(--ink-2);">
      <span>GRCh38 <strong>${{gInfo.snvVcf}}</strong> &middot; ${{gInfo.cdna}} &middot; Codon ${{gInfo.wtCod}}&rarr;${{gInfo.mutCod}} &middot; Exon ${{gInfo.exon}}/${{gInfo.totalExons}} (${{gInfo.strandSym}})</span>
      <button class="btn-subtle no-print" style="height:21px;padding:0 7px;font-size:10.5px;" onclick="navigator.clipboard.writeText('${{gInfo.snvVcf}}');this.textContent='Copied ✓';setTimeout(()=>this.textContent='Copy VCF',1300)">Copy VCF</button>
    </div>
  ` : '';

  const pocketNote = (dmet !== null && dmet <= 6.0) ? `${{dmet}} Å from catalytic metal`
    : ((dna !== null && dna <= 6.0) ? `${{dna}} Å from bound RNA/DNA` : `${{geom ? geom.caNeighbors8A : '—'}} Cα within 8 Å`);

  const oneLine = wt === mut
    ? `${{wp.name}} at codon ${{pos}} (${{pocketNote}}; pLDDT ${{plddt ?? '—'}}). Pick any amino acid above to model its 3D side chain.`
    : `<strong>${{wp.name}} &rarr; ${{mp.name}}</strong> (${{pocketNote}}${{pdb ? ' in PDB ' + pdb : ''}}) · ${{clashCount === 0 ? '0 steric clashes' : clashCount + ' clash(es) <2.4 Å'}} · TCGA tumours: <strong>${{tn}}</strong>${{cv ? ' · ClinVar: ' + cv : ''}}.`;

  box.innerHTML = `
    <div>
      <div style="display:flex;justify-content:space-between;align-items:center;gap:8px;flex-wrap:wrap;">
        <div style="display:flex;align-items:center;gap:7px;flex-wrap:wrap;">
          <span style="font-size:18px;font-weight:650;letter-spacing:-0.02em;">
            ${{gene}} <code>${{wt}}${{pos}}${{mut}}</code>
          </span>
          <span style="font-size:12.5px;color:var(--ink-3);">p.${{wp.three}}${{pos}}${{mp.three}}</span>
          ${{domBadge}}
        </div>
        <div style="display:flex;align-items:center;gap:6px;flex-wrap:wrap;">
          ${{regimeTag}}
          <button class="btn-subtle no-print" style="height:26px;padding:0 9px;font-size:11.5px;" onclick="shareCurrentVariant(this)" title="Copy direct URL to this 3D variant">Share Link</button>
          <button class="btn-subtle no-print" style="height:26px;padding:0 9px;font-size:11.5px;" onclick="printVariantReport()" title="Print or save 1-page PDF report">Report / PDF</button>
        </div>
      </div>
      ${{genomicBar}}

      <div style="display:flex;justify-content:space-between;align-items:center;flex-wrap:wrap;gap:6px;margin-top:10px;margin-bottom:-6px;font-size:11px;color:var(--ink-3);">
        <span>19-AA Substitution Strip (click any residue to mutate in 3D):</span>
        <span style="display:inline-flex;gap:9px;">
          <span class="leg-item"><i class="leg-line" style="background:#c5221f;"></i>Active-site</span>
          <span class="leg-item"><i class="leg-line" style="background:#b06000;"></i>Unfolding</span>
          <span class="leg-item"><i class="leg-line" style="background:#00709b;"></i>Neutral</span>
        </span>
      </div>
      <div class="aa-strip">${{aaCells}}</div>

      <div class="metric-tiles">
        <div class="m-tile" data-term="alphamissense">
          <span class="m-lbl">AlphaMissense</span>
          <div class="m-val">${{am !== null ? am.toFixed(3) : '—'}}</div>
          <div class="m-sub">${{am !== null ? (am >= 0.564 ? 'Pathogenic' : (am <= 0.34 ? 'Benign' : 'Ambiguous')) : '—'}}</div>
        </div>
        <div class="m-tile" data-term="esm_if1">
          <span class="m-lbl">ESM-IF1 3D &Delta;logP</span>
          <div class="m-val">${{if1 !== null ? if1.toFixed(2) : '—'}}</div>
          <div class="m-sub">${{if1 !== null ? (if1 <= -1.5 ? 'Unfolding' : 'Fold intact') : '—'}}</div>
        </div>
        <div class="m-tile" data-term="esm2_llr">
          <span class="m-lbl">ESM-2 (3B) &Delta;&ell;</span>
          <div class="m-val">${{e2 !== null ? e2.toFixed(2) : '—'}}</div>
          <div class="m-sub">${{e2 !== null ? (e2 <= -5 ? 'Constrained' : 'Tolerated') : '—'}}</div>
        </div>
        <div class="m-tile" data-term="delta_q">
          <span class="m-lbl">Charge &Delta;q</span>
          <div class="m-val">${{dq >= 0 ? '+' : ''}}${{dq.toFixed(1)}} e</div>
          <div class="m-sub">${{wp.cls}} &rarr; ${{mp.cls}}</div>
        </div>
        <div class="m-tile" data-term="delta_vol">
          <span class="m-lbl">Volume &Delta;V</span>
          <div class="m-val">${{dvol >= 0 ? '+' : ''}}${{dvol.toFixed(0)}} &Aring;&sup3;</div>
          <div class="m-sub" data-term="rotamer_clash">${{clashCount === 0 ? '0 clashes (<2.4Å)' : clashCount + ' clash(es)'}}</div>
        </div>
        <div class="m-tile" data-term="grantham">
          <span class="m-lbl">Grantham D</span>
          <div class="m-val">${{gDist}}</div>
          <div class="m-sub">${{gDist <= 50 ? 'Conservative' : (gDist <= 100 ? 'Moderate' : 'Radical')}}</div>
        </div>
      </div>
    </div>

    <div class="score-bottom">
      <div class="summary-pill">${{oneLine}}</div>
      <div class="radar-box" title="ProtVar-style 5-axis normalized fingerprint: AlphaMissense (AM), ESM-2 constraint (Seq), ESM-IF1 unfolding (Fold), Grantham radicality (Chem), and 3D Pocket proximity (Pocket)">
        ${{buildProtVarRadarSVG(am, e2, if1, gDist, minPktA, radarCol)}}
        <span style="font-size:9.5px;color:var(--ink-3);font-family:var(--mono);margin-top:-2px;">0&rarr;1 Impact Radar</span>
      </div>
    </div>
  `;
}}

function inspectSaturationCodon() {{
  const gene = document.getElementById('lol-gene').value;
  const rawQ = normalizeMutToken(document.getElementById('sat-query').value);
  const box = document.getElementById('sat-codon-detail');
  if (!rawQ) return;
  const gCache = window.SAT_CACHE ? window.SAT_CACHE[gene] : null;
  if (!gCache) {{ ensureSatLoaded(gene, inspectSaturationCodon); return; }}

  const mVar = rawQ.match(/^([A-Z])?(\\d+)([A-Z*])?$/);
  if (!mVar || !gCache[mVar[2]]) return;
  const posStr = mVar[2], pos = parseInt(posStr, 10), entry = gCache[posStr];
  const [wt, pdb, dna, dmet, dprt, plddt] = entry.m;
  const targetAlt = (mVar[3] && mVar[3] !== '*') ? mVar[3] : '';

  let activeSubRow = (targetAlt && targetAlt !== wt) ? (entry.s.find(r => r[0] === targetAlt) || null) : null;
  if (!activeSubRow) activeSubRow = entry.s.slice().sort((a,b) => (b[5] - a[5]) || ((b[1]||0) - (a[1]||0)))[0];
  const mut3D = activeSubRow ? activeSubRow[0] : wt;

  active3DState = {{ gene, pos, wt, mut: mut3D, lastPdbText: active3DState.lastPdbText }};
  dmsCenterPos = pos;
  recordRecentVariant(gene, `${{wt}}${{pos}}${{mut3D}}`);
  updateLollipopActiveMarker(pos);
  if (document.getElementById('exp-dms').classList.contains('active')) renderDmsHeatmap();
  if (document.getElementById('tab-variants').classList.contains('active')) {{
    syncUrlHash(`#variant=${{gene}}:${{wt}}${{pos}}${{mut3D}}`);
  }}
  ensureStruct3DLoaded(gene, () => {{
    renderPropertyDeltaCard(gene, pos, wt, mut3D, entry, activeSubRow);
    refresh3DScene(false);
  }});

  let html = `<div class="card-head">
      <h2>All 19 Substitutions at ${{gene}} ${{wt}}${{pos}}</h2>
      <span class="meta-note">PDB ${{pdb || 'AlphaFold v6'}} · RNA ${{dna ?? '—'}} Å · Metal ${{dmet ?? '—'}} Å</span>
    </div>
    <div class="table-guide no-print">
      <span><strong>Thresholds:</strong> <code>AlphaMissense &ge; 0.564</code> Pathogenic &middot; <code>ESM-2 &le; &minus;5.0</code> Constrained &middot; <code>ESM-IF1 &gt; &minus;1.5</code> Fold-preserving (<span class="leg-item"><i class="dot dot-crimson"></i>Active site</span>) vs <code>&le; &minus;1.5</code> (<span class="leg-item"><i class="dot dot-amber"></i>Unfolding</span>) &middot; <code>Grantham &gt; 100</code> Radical.</span>
      <span class="meta-note">Click <strong>View</strong> to model side chain in 3D &middot; Click <code>&#8597;</code> headers to sort</span>
    </div>
    <div class="table-wrap" style="max-height:420px;">
      <table>
        <thead><tr><th>3D</th><th>Variant</th><th data-term="grantham">Grantham</th><th data-term="delta_q">&Delta;q</th><th data-term="delta_vol">&Delta;V (&Aring;&sup3;)</th><th data-term="alphamissense">AlphaMissense</th><th data-term="esm2_llr">ESM-2</th><th data-term="esm1v_llr">ESM-1v</th><th data-term="esm_if1">ESM-IF1</th><th data-term="mech_regime">Mechanism</th><th>TCGA</th><th>ClinVar</th></tr></thead>
        <tbody>`;
  entry.s.forEach(r => {{
    const [alt, am, e2, e1v, if1, tn, cv] = r;
    const wp = AA_PROPS[wt] || AA_PROPS.A, mp = AA_PROPS[alt] || wp;
    const dq = mp.q - wp.q, dv = mp.vol - wp.vol, gd = granthamDist(wt, alt);
    const cns = (am !== null && am >= 0.564) || (e2 !== null && e2 <= -5.0);
    const regime = cns ? ((if1 !== null && if1 > -1.5) ? tagDot('dot-crimson','Active site') : tagDot('dot-amber','Unfolding')) : tagDot('dot-petrol','Neutral');
    html += `<tr>
      <td><button class="btn-subtle" style="height:24px;padding:0 9px;font-size:11.5px;" onclick="selectMutAA(${{pos}},'${{wt}}','${{alt}}')">View</button></td>
      <td><code>${{wt}}${{pos}}${{alt}}</code></td>
      <td>${{gd}}</td><td>${{dq>=0?'+':''}}${{dq.toFixed(1)}}</td><td>${{dv>=0?'+':''}}${{dv.toFixed(0)}}</td>
      <td>${{am !== null ? am.toFixed(3) : '—'}}</td><td>${{e2 !== null ? e2.toFixed(2) : '—'}}</td>
      <td>${{e1v !== null ? e1v.toFixed(2) : '—'}}</td><td>${{if1 !== null ? if1.toFixed(2) : '—'}}</td>
      <td>${{regime}}</td><td>${{tn > 0 ? '<strong>'+tn+'</strong>' : '0'}}</td><td>${{cv || '—'}}</td>
    </tr>`;
  }});
  box.innerHTML = html + '</tbody></table></div>';
  markSortableHeaders(box);
}}

function tierPill(t) {{
  if (t === '3D_Catalytic_or_Substrate_Pocket_Driver') return tagDot('dot-crimson', '3D pocket');
  if (t === 'Pathogenic_Fold_or_Interface_Missense') return tagDot('dot-amber', 'Pathogenic fold');
  if (t === 'Truncating_LoF') return tagDot('dot-slate', 'Truncating');
  if (t === 'Ambiguous_Missense_VUS') return tagDot('dot-slate', 'VUS');
  if (t === 'Likely_Benign_Passenger_Missense') return tagDot('dot-petrol', 'Benign');
  return tagDot('dot-slate', t);
}}
function mechPill(m) {{
  if (m === 'Fold_Preserving_Catalytic_or_Interface') return tagDot('dot-crimson', 'Active site');
  if (m === 'Thermodynamic_Fold_Destabilizing') return tagDot('dot-amber', 'Unfolding');
  if (m === 'Tolerated_Neutral') return tagDot('dot-petrol', 'Neutral');
  return '—';
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
    if (q && !(`${{v.g}} ${{v.p}} ${{v.cohorts}} ${{v.cv}}`.toUpperCase().includes(q))) return false;
    return true;
  }});
  document.getElementById('var-count').textContent = `${{currentFilteredVars.length}} unique variants (${{currentFilteredVars.reduce((a,b)=>a+b.n,0)}} calls)`;
  document.querySelector('#var-table tbody').innerHTML = currentFilteredVars.slice(0, 250).map(v => {{
    const cleanP = v.p.replace(/^p\\./i, '');
    const jumpArg = /^[A-Z]\\d+[A-Z]$/.test(cleanP) ? cleanP : (v.pos > 0 ? String(v.pos) : '');
    return `<tr>
      <td>${{jumpArg ? `<button class="btn-subtle" style="height:23px;padding:0 8px;font-size:11px;" onclick="openVariantIn3D('${{v.g}}','${{jumpArg}}')">View</button>` : '—'}}</td>
      <td><strong>${{v.g}}</strong></td><td><code>${{v.p}}</code></td><td>${{v.n}}</td><td>${{v.cv || '—'}}</td>
      <td>${{tierPill(v.tier)}}</td><td>${{mechPill(v.mech)}}</td><td><code>${{v.pdb || '—'}}</code></td>
      <td>${{v.dna !== null ? v.dna.toFixed(2) : '—'}}</td><td>${{v.dmet !== null ? v.dmet.toFixed(2) : '—'}}</td>
      <td>${{v.am !== null ? v.am.toFixed(3) : '—'}}</td><td>${{v.esm !== null ? v.esm.toFixed(2) : '—'}}</td>
      <td>${{v.esm1v !== null ? v.esm1v.toFixed(2) : '—'}}</td><td>${{v.if1 !== null ? v.if1.toFixed(2) : '—'}}</td>
      <td>${{v.cohorts}}</td>
    </tr>`;
  }}).join('');
}}

function renderGtexBars() {{
  const g = document.getElementById('gtex-gene').value;
  const vals = ATLAS.gtex_tpm[g] || [];
  const pairs = ATLAS.gtex_tissues.map((t, i) => [t, vals[i] || 0]).sort((a,b) => b[1] - a[1]);
  const mx = Math.max(1, ...pairs.map(x => x[1]));
  const gBox = document.getElementById('gtex-bars');
  gBox.innerHTML = '<table><thead><tr><th>Adult Tissue (54 Sites)</th><th data-term="gtex_tpm">Median TPM</th><th>Relative Level</th></tr></thead><tbody>' +
    pairs.map(([t, v]) => `<tr><td>${{t.replace(/_/g, ' ')}}</td><td><strong>${{v.toFixed(2)}}</strong></td><td><div style="width:${{Math.max(2, (v/mx)*240).toFixed(0)}}px;height:7px;background:var(--accent);border-radius:3px;"></div></td></tr>`).join('') + '</tbody></table>';
  markSortableHeaders(gBox);
}}

/* Render Authentic KaTeX Formula Cards with Symbol Keys & Primary Literature Citations */
function renderFormulas() {{
  const cat = document.getElementById('formula-cat').value;
  const q = (document.getElementById('formula-q').value || '').trim().toLowerCase();
  const list = (ATLAS.formulas || []).filter(f => {{
    if (cat !== 'ALL' && f.cat !== cat) return false;
    const varTxt = (f.vars || []).map(v => v.join(' ')).join(' ');
    if (q && !(`${{f.id}} ${{f.title}} ${{f.summary}} ${{f.cat}} ${{f.citation || ''}} ${{varTxt}}`.toLowerCase().includes(q))) return false;
    return true;
  }});
  const wrap = document.getElementById('formula-cards-wrap');
  wrap.innerHTML = list.map((f, idx) => {{
    const refHtml = (f.refs && f.refs.length)
      ? `<div class="eq-cite"><strong>References:</strong> ${{f.refs.map(r => `<a href="${{r.url}}" target="_blank" rel="noopener">${{r.cite}} &#8599;</a>`).join(' &middot; ')}}</div>`
      : (f.citation ? `<div class="eq-cite"><strong>Reference:</strong> <a href="${{f.doi_url || '#'}}" target="_blank" rel="noopener">${{f.citation}} &#8599;</a></div>` : '');
    return `
    <div class="eq-card" id="${{f.id}}">
      <div>
        <div style="display:flex;justify-content:space-between;align-items:center;margin-bottom:6px;">
          <span class="ver-pill" style="color:var(--accent);background:var(--accent-soft);font-weight:600;">${{f.id}}</span>
          <span class="meta-note">${{f.cat}}</span>
        </div>
        <div style="font-size:14px;font-weight:600;margin-bottom:10px;">${{f.title}}</div>
        <div class="eq-math-box" id="katex-eq-${{idx}}"></div>
        <div class="eq-vars" id="katex-eq-vars-${{idx}}"></div>
      </div>
      <div>
        <div style="font-size:12.5px;color:var(--ink-2);">${{f.summary}}</div>
        ${{refHtml}}
      </div>
    </div>
    `;
  }}).join('');
  list.forEach((f, idx) => {{
    renderLatex(f.latex, document.getElementById(`katex-eq-${{idx}}`), true);
    renderFormulaVars(f.vars, document.getElementById(`katex-eq-vars-${{idx}}`), 'Where (Symbols, Alphabets & Operators)');
  }});
}}

function renderDatasetsCatalogue() {{
  const q = (document.getElementById('ds-q').value || '').trim().toLowerCase();
  const list = (ATLAS.datasets || []).filter(d => !q || (`${{d.file}} ${{d.pk}} ${{d.desc}}`.toLowerCase().includes(q)));
  document.querySelector('#datasets-cat-table tbody').innerHTML = list.map(d => `
    <tr>
      <td><code>${{d.file}}</code></td>
      <td><strong>${{Number(d.rows).toLocaleString()}}</strong> &times; ${{d.cols}}</td>
      <td>${{d.size}}</td><td><code>${{d.pk}}</code></td>
      <td style="white-space:normal;min-width:280px;">${{d.desc}}</td>
    </tr>
  `).join('');
}}

function renderGlossaryCards() {{
  const q = (document.getElementById('gloss-q').value || '').trim().toLowerCase();
  const entries = Object.entries(ATLAS.glossary || {{}})
    .map(([k, v]) => ({{key: k, ...v}}))
    .sort((a, b) => a.title.localeCompare(b.title))
    .filter(g => {{
      const varTxt = (g.vars || []).map(v => v.join(' ')).join(' ');
      return !q || (`${{g.title}} ${{g.cat}} ${{g.meaning}} ${{g.thresholds}} ${{varTxt}}`.toLowerCase().includes(q));
    }});
  document.getElementById('gloss-count').textContent = `${{entries.length}} terms`;
  const wrap = document.getElementById('glossary-cards-wrap');
  wrap.innerHTML = entries.map((g, i) => {{
    const fRef = g.anchor ? FORMULA_BY_ID[g.anchor] : null;
    const citeRow = (fRef && fRef.refs && fRef.refs.length)
      ? `<div class="eq-cite"><span>Ref: ${{fRef.refs.slice(0, 2).map(r => `<a href="${{r.url}}" target="_blank" rel="noopener">${{r.cite}} &#8599;</a>`).join(' &middot; ')}}</span><a href="#eq=${{fRef.id}}" onclick="event.preventDefault();openFormulaTab('${{fRef.id}}');" style="margin-left:auto;">Eq. ${{fRef.id}} &rarr;</a></div>`
      : '';
    return `
    <div class="eq-card">
      <div>
        <div style="display:flex;justify-content:space-between;align-items:baseline;margin-bottom:4px;">
          <strong style="font-size:13.5px;">${{g.title}}</strong>
          <span class="meta-note">${{g.cat}}</span>
        </div>
        <div style="font-size:12.5px;color:var(--ink-2);margin-bottom:10px;">${{g.meaning}}</div>
        <div class="eq-math-box" style="min-height:44px;padding:8px 10px;" id="katex-gl-${{i}}"></div>
        <div class="eq-vars" id="katex-gl-vars-${{i}}"></div>
      </div>
      <div>
        <div class="meta-note">${{g.thresholds}}</div>
        ${{citeRow}}
      </div>
    </div>
    `;
  }}).join('');
  entries.forEach((g, i) => {{
    renderLatex(g.latex, document.getElementById(`katex-gl-${{i}}`), false);
    renderFormulaVars(g.vars, document.getElementById(`katex-gl-vars-${{i}}`), 'Symbols in Formula');
  }});
}}

function renderAllTables() {{
  document.querySelector('#cohort-burden-table tbody').innerHTML = (ATLAS.top_tumours || []).map(r => {{
    const n = Number(r.N), c10 = Number(r.patient_prop_core10_strict);
    return `<tr>
      <td>${{r.rank_core10_strict}}</td>
      <td><strong>${{r.tumour_acronym}}${{r.low_coverage_N_lt_100 === '1' ? '*' : ''}}</strong></td>
      <td>${{n.toLocaleString()}}</td>
      <td><strong>${{(c10 * 100).toFixed(1)}}%</strong> (${{Math.round(c10 * n)}})</td>
      <td>${{(Number(r.patient_prop_extended12_strict) * 100).toFixed(1)}}%</td>
      <td>${{(Number(r.patient_prop_panel14_strict) * 100).toFixed(1)}}%</td>
      <td>${{(Number(r.patient_prop_core10_any_alter) * 100).toFixed(1)}}%</td>
      <td>${{(Number(r.ago1_3_4_1p34_concordant_cna_prop) * 100).toFixed(1)}}%</td>
    </tr>`;
  }}).join('');

  document.querySelector('#study-decoder-table tbody').innerHTML = (ATLAS.study_summary || []).map(s => `
    <tr>
      <td><strong>${{s.tumour_acronym}}${{s.low_coverage_flag_N_lt_100 === '1' ? '*' : ''}}</strong></td>
      <td>${{s.studyName.replace(' (TCGA, PanCancer Atlas)', '')}}</td>
      <td><strong>${{Number(s.intersection_sequenced_cna_n).toLocaleString()}}</strong></td>
      <td>${{s.sampleCount_sequenced}}</td><td>${{s.sampleCount_cna}}</td>
      <td>${{s.sampleCount_mrna_matched_n}}</td><td>${{s.sampleCount_mirnaseq_matched_n}}</td>
      <td>${{s.hypermutated_samples_n}}</td><td>${{s.nonhypermutated_samples_n}}</td>
    </tr>
  `).join('');

  document.querySelector('#mir-bench-table tbody').innerHTML = ATLAS.mir_bench.map(b => `
    <tr>
      <td><strong>${{b.benchmark_id}}</strong></td><td><code>${{b.metric}}</code></td>
      <td>${{b.group_A}} (${{b.n_group_A}})</td><td>${{Number(b.median_group_A).toFixed(3)}}</td>
      <td>${{b.group_B}} (${{b.n_group_B}})</td><td>${{Number(b.median_group_B).toFixed(3)}}</td>
      <td><strong>${{(Number(b.median_shift_A_minus_B)>=0?'+':'') + Number(b.median_shift_A_minus_B).toFixed(3)}}</strong></td>
      <td>${{b.p_value_two_sided}}</td>
    </tr>
  `).join('');

  renderDicerMirTable();

  document.querySelector('#top-axes-table tbody').innerHTML = ATLAS.top_axes.map(ax => `
    <tr>
      <td>${{ax.rank}}</td><td><code>${{ax.matureMiRNA}}</code></td><td><strong>${{ax.targetGeneSymbol}}</strong></td>
      <td>${{ax.n_dual_tumours}}</td><td><strong>${{Number(ax.pearson_r_log2).toFixed(4)}}</strong></td><td>${{ax.t_statistic}}</td>
    </tr>
  `).join('');

  document.querySelector('#clinvar-roc-table tbody').innerHTML = ATLAS.clinvar_roc.map(r => `
    <tr>
      <td><code>${{r.benchmark_subset.replace('ClinVar_','')}}</code></td><td><strong>${{r.predictor}}</strong></td>
      <td>${{r.n_pos}}</td><td>${{r.n_neg}}</td>
      <td><strong>${{Number(r.roc_auc).toFixed(4)}}</strong> [${{Number(r.roc_ci95_low).toFixed(3)}}–${{Number(r.roc_ci95_high).toFixed(3)}}]</td>
      <td>${{Number(r.pr_auc).toFixed(4)}} (${{Number(r.pr_baseline).toFixed(3)}})</td><td>${{r.mann_whitney_p}}</td>
    </tr>
  `).join('');

  document.querySelector('#msk-rep-table tbody').innerHTML = ATLAS.msk_rep.map(m => `
    <tr>
      <td><strong>${{m.geneSymbol}}</strong></td><td>${{Number(m.msk_N_profiled).toLocaleString()}}</td>
      <td><strong>${{(Number(m.msk_mut_freq)*100).toFixed(2)}}%</strong> (${{m.msk_mut_samples_n}})</td>
      <td>${{(Number(m.tcga_mut_freq_N10195)*100).toFixed(2)}}%</td>
      <td>${{(Number(m.msk_nonhypermut_mut_freq)*100).toFixed(2)}}%</td>
      <td>${{(Number(m.msk_truncating_freq)*100).toFixed(2)}}%</td>
      <td><strong>${{(Number(m.msk_3dpocket_or_hotspot_freq)*100).toFixed(2)}}%</strong></td>
      <td>${{(Number(m.msk_primary_mut_freq)*100).toFixed(2)}}% / ${{(Number(m.msk_metastasis_mut_freq)*100).toFixed(2)}}%</td>
    </tr>
  `).join('');

  document.querySelector('#msk-dicer-table tbody').innerHTML = ATLAS.msk_dicer.map(d => `
    <tr>
      <td><strong>${{d.msk_cancer_type}}</strong></td><td>${{d.N_samples}}</td>
      <td>${{d.dicer1_mut_n}} (${{(Number(d.dicer1_mut_freq)*100).toFixed(2)}}%)</td>
      <td><strong>${{d.dicer1_rnase3b_hotspot_n}} (${{(Number(d.dicer1_rnase3b_hotspot_freq)*100).toFixed(2)}}%)</strong></td>
      <td>${{d.dicer1_truncating_n}} (${{(Number(d.dicer1_truncating_freq)*100).toFixed(2)}}%)</td>
    </tr>
  `).join('');

  document.querySelector('#gene-ref-table tbody').innerHTML = ATLAS.panel.map(p => `
    <tr>
      <td><strong>${{p.geneSymbol}}</strong></td>
      <td>${{p.gnomad_pLI !== '' ? Number(p.gnomad_pLI).toFixed(2) : '—'}}</td>
      <td><strong>${{p.gnomad_LOEUF !== '' ? Number(p.gnomad_LOEUF).toFixed(3) : '—'}}</strong></td>
      <td>${{p.gnomad_mis_z !== '' ? Number(p.gnomad_mis_z).toFixed(2) : '—'}}</td>
      <td>${{p.hpa_main_subcellular_location || '—'}}</td>
      <td>${{p.gtex_pan54_median_tpm !== '' ? Number(p.gtex_pan54_median_tpm).toFixed(1) : '—'}}</td>
      <td>${{p.gtex_top_tissue}} (${{Number(p.gtex_top_tissue_tpm||0).toFixed(1)}})</td>
      <td>${{p.gtex_tau_tissue_specificity}}</td>
      <td>${{p.opentargets_sm_tractable === '1' ? tagDot('dot-emerald','Ligandable') : '—'}}</td>
      <td>${{p.clinvar_pathogenic_or_likely_pathogenic_n}} / ${{p.clinvar_total_single_gene_n}}</td>
    </tr>
  `).join('');

  const cpMap = {{}}; ATLAS.cptac.forEach(c => {{ cpMap[c.geneSymbol] = c; }});
  const dmMap = {{}}; ATLAS.depmap.forEach(d => {{ dmMap[d.geneSymbol] = d; }});
  document.querySelector('#consortia-table tbody').innerHTML = ATLAS.panel.map(p => {{
    const g = p.geneSymbol, c = cpMap[g] || {{}}, d = dmMap[g] || {{}};
    return `<tr>
      <td><strong>${{g}}</strong></td><td>${{p.panelRole.replace(/_/g,' ')}}</td><td><code>${{p.pdbComplexIds}}</code></td>
      <td>${{c.spearman_CNA_to_mRNA || '—'}}</td><td>${{c.spearman_CNA_to_Protein || '—'}}</td>
      <td>${{c.post_transcriptional_CNA_attenuation || '—'}}</td>
      <td><strong>${{d.median_chronos_gene_effect}}</strong></td>
      <td>${{(Number(d.prop_dependent_le_m05)*100).toFixed(1)}}%</td>
      <td>${{p.depmapIsCommonEssential === '1' ? tagDot('dot-crimson','Essential') : tagDot('dot-petrol','Buffered')}}</td>
    </tr>`;
  }}).join('');

  document.querySelector('#collat-table tbody').innerHTML = (ATLAS.collateral || []).map(c => `
    <tr>
      <td><code>${{c.test_id.split('_')[0]}}</code></td><td><strong>${{c.crispr_target_gene}}</strong></td>
      <td style="white-space:normal;min-width:260px;">${{c.hypothesis}}</td>
      <td>${{c.n_altered_lines}} (${{Number(c.median_chronos_altered).toFixed(3)}})</td>
      <td>${{c.n_control_lines}} (${{Number(c.median_chronos_control).toFixed(3)}})</td>
      <td><strong>${{(Number(c.chronos_shift_altered_minus_control)>=0?'+':'') + Number(c.chronos_shift_altered_minus_control).toFixed(3)}}</strong></td>
      <td>${{c.p_value_two_sided}}</td>
    </tr>
  `).join('');

  document.querySelector('#surv-table tbody').innerHTML = ATLAS.surv_bench.map(s => `
    <tr>
      <td><strong>${{s.test_id.split('_').slice(0,2).join('_')}}</strong></td><td>${{s.cohort}}</td><td><code>${{s.endpoint}}</code></td>
      <td>${{s.events_altered}}/${{s.n_altered}}</td><td>${{s.events_control}}/${{s.n_control}}</td>
      <td><strong>${{Number(s.hazard_ratio_stratified).toFixed(2)}}</strong> [${{Number(s.hr_ci95_low).toFixed(2)}}–${{Number(s.hr_ci95_high).toFixed(2)}}]</td>
      <td>${{s.p_value_two_sided}}</td>
      <td>${{s.multivar_cox_hr ? '<strong>' + Number(s.multivar_cox_hr).toFixed(2) + '</strong> [' + Number(s.multivar_ci95_low).toFixed(2) + '–' + Number(s.multivar_ci95_high).toFixed(2) + ']' : '—'}}</td>
      <td>${{s.multivar_p_value || '—'}}</td>
    </tr>
  `).join('');

  document.querySelector('#imm-table tbody').innerHTML = ATLAS.imm_bench.map(i => `
    <tr>
      <td><strong>${{i.benchmark_id}}</strong></td><td><code>${{i.immune_or_instability_metric}}</code></td>
      <td>${{i.group_A}} (${{i.n_group_A}})</td><td>${{Number(i.median_group_A).toFixed(2)}}</td>
      <td>${{i.group_B}} (${{i.n_group_B}})</td><td>${{Number(i.median_group_B).toFixed(2)}}</td>
      <td><strong>${{(Number(i.median_shift_A_minus_B)>=0?'+':'') + Number(i.median_shift_A_minus_B).toFixed(2)}}</strong></td>
      <td>${{i.p_value_two_sided}}</td>
    </tr>
  `).join('');

  document.querySelector('#excluded-cands-table tbody').innerHTML = (ATLAS.excluded_cands || []).map(e => `
    <tr><td><strong>${{e.candidateGene}}</strong></td><td><code>${{e.cytoband}}</code></td><td>${{e.molecularFunction}}</td><td style="white-space:normal;">${{e.exclusionRationale}}</td></tr>
  `).join('');

  document.querySelector('#consortia-sources-table tbody').innerHTML = (ATLAS.consortia || []).map(c => `
    <tr><td><strong>${{c.consortium}}</strong></td><td>${{c.scale}}</td><td style="white-space:normal;">${{c.layers}}</td><td><code>${{c.join_key}}</code></td></tr>
  `).join('');
}}

/* Interactive Small-RNA DICER1 5p:3p Strand-Collapse SVG Strip Plot */
let selectedDicerSample = '';
function renderDicerMirTable(filterClass) {{
  const rows = [...ATLAS.dicer_pts]
    .filter(d => !filterClass || d.c === filterClass)
    .sort((a,b) => a.d53 - b.d53);
  document.querySelector('#dicer-mir-table tbody').innerHTML = rows.map(d => `
    <tr id="dicer-row-${{d.s}}" class="${{selectedDicerSample===d.s?'row-flash':''}}">
      <td><code>${{d.s}}</code></td><td><strong>${{d.t}}</strong></td>
      <td>${{d.c === 'DICER1_RNaseIIIb_3D_Catalytic_Pocket' ? tagDot('dot-crimson','RNase IIIb pocket') : d.c.replace('DICER1_','')}}</td>
      <td><strong>${{d.d53.toFixed(3)}}</strong></td>
      <td>${{d.f5 !== null ? d.f5.toFixed(3) : '—'}}</td>
      <td>${{d.l7z !== null && d.l7z !== undefined ? (d.l7z>=0?'+':'') + d.l7z.toFixed(2) : '—'}}</td>
    </tr>
  `).join('');
}}

function selectDicerDot(sampleId) {{
  selectedDicerSample = sampleId;
  document.getElementById('dicer-reset-btn').style.display = 'inline-block';
  renderDicerMirTable('');
  const rowEl = document.getElementById('dicer-row-' + sampleId);
  if (rowEl) rowEl.scrollIntoView({{behavior: 'smooth', block: 'center'}});
}}

function filterDicerByClass(cls) {{
  selectedDicerSample = '';
  document.getElementById('dicer-reset-btn').style.display = 'inline-block';
  renderDicerMirTable(cls);
}}

function resetDicerSelection() {{
  selectedDicerSample = '';
  document.getElementById('dicer-reset-btn').style.display = 'none';
  renderDicerMirTable('');
}}

function renderDicerStripPlot() {{
  const wrap = document.getElementById('dicer-svg-plot');
  if (!wrap || !ATLAS.dicer_pts || !ATLAS.dicer_pts.length) return;

  const groups = [
    {{ key: 'DICER1_RNaseIIIb_3D_Catalytic_Pocket', label: 'RNase IIIb 3D Pocket', color: '#c5221f' }},
    {{ key: 'DICER1_Truncating_LoF', label: 'Truncating LoF', color: '#4b5563' }},
    {{ key: 'DICER1_Other_Pathogenic_Missense', label: 'Pathogenic Fold / Interface', color: '#b06000' }},
    {{ key: 'DICER1_Ambiguous_or_Benign_Missense', label: 'Benign / Passenger VUS', color: '#00709b' }}
  ];
  const knownKeys = new Set(groups.map(g => g.key));
  const extraKeys = [...new Set(ATLAS.dicer_pts.map(d => d.c))].filter(k => !knownKeys.has(k));
  extraKeys.forEach(k => groups.push({{ key: k, label: k.replace(/^DICER1_/, '').replace(/_/g, ' '), color: '#00709b' }}));
  const activeGroups = groups.filter(g => ATLAS.dicer_pts.some(d => d.c === g.key));

  const W = 1360, rowH = 44, topPad = 32, botPad = 48, padL = 220, padR = 40;
  const H = topPad + activeGroups.length * rowH + botPad;
  const minX = -2.4, maxX = 1.6;
  const xOf = v => padL + ((Math.min(maxX, Math.max(minX, v)) - minX) / (maxX - minX)) * (W - padL - padR);
  const wtMed = 0.475;

  let svg = `<svg viewBox="0 0 ${{W}} ${{H}}" width="100%" height="${{H}}">`;
  [-2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5].forEach(tk => {{
    const x = xOf(tk).toFixed(1);
    svg += `<line x1="${{x}}" y1="${{topPad - 8}}" x2="${{x}}" y2="${{H - botPad}}" stroke="${{tk === 0.0 ? '#cbd5e1' : '#ebeef2'}}" stroke-width="${{tk === 0.0 ? '1.3' : '1'}}"/>`;
    svg += `<text x="${{x}}" y="${{H - 26}}" font-size="10.5" fill="#80868b" text-anchor="middle">${{tk >= 0 ? '+' + tk.toFixed(1) : tk.toFixed(1)}}</text>`;
  }});
  svg += `<text x="${{xOf(-1.6).toFixed(1)}}" y="${{H - 8}}" font-size="10.5" font-weight="600" fill="#c5221f" text-anchor="middle">&larr; Selective 5p Cleavage Collapse (&lt; 0)</text>`;
  svg += `<text x="${{xOf(0.0).toFixed(1)}}" y="${{H - 8}}" font-size="10.5" fill="#5f6368" text-anchor="middle">Mean Duplex log&#8322;((RPM&#8325;&#8346;+1)/(RPM&#8323;&#8346;+1))</text>`;
  svg += `<text x="${{xOf(1.05).toFixed(1)}}" y="${{H - 8}}" font-size="10.5" font-weight="600" fill="#137333" text-anchor="middle">Canonical 5p Strand Dominance (&gt; 0) &rarr;</text>`;

  const xWt = xOf(wtMed).toFixed(1);
  svg += `<line x1="${{xWt}}" y1="${{topPad - 12}}" x2="${{xWt}}" y2="${{H - botPad}}" stroke="#137333" stroke-width="1.5" stroke-dasharray="4,3"/>`;
  svg += `<text x="${{xWt}}" y="16" font-size="10.5" font-weight="600" fill="#137333" text-anchor="middle">Wild-Type Median (+0.475, n = 9,129)</text>`;

  activeGroups.forEach((grp, idx) => {{
    const yMid = topPad + idx * rowH + rowH * 0.5;
    const pts = ATLAS.dicer_pts.filter(d => d.c === grp.key).sort((a,b) => a.d53 - b.d53);
    const med = pts.length ? pts[Math.floor(pts.length / 2)].d53 : 0;
    svg += `<line x1="${{padL}}" y1="${{yMid}}" x2="${{W - padR}}" y2="${{yMid}}" stroke="#f1f3f5" stroke-width="1"/>`;
    svg += `<text x="${{padL - 14}}" y="${{yMid + 4}}" font-size="12" font-weight="600" fill="#202124" text-anchor="end" style="cursor:pointer;" onclick="filterDicerByClass('${{grp.key}}')">${{grp.label}} (n=${{pts.length}})</text>`;
    const xm = xOf(med).toFixed(1);
    svg += `<line x1="${{xm}}" y1="${{yMid - 12}}" x2="${{xm}}" y2="${{yMid + 12}}" stroke="${{grp.color}}" stroke-width="2.6"/>`;
    pts.forEach((d, j) => {{
      const jitter = ((j % 5) - 2) * 3.4;
      const cx = xOf(d.d53).toFixed(1);
      const cy = (yMid + jitter).toFixed(1);
      svg += `<circle cx="${{cx}}" cy="${{cy}}" r="4.6" fill="${{grp.color}}" fill-opacity="0.78" stroke="#ffffff" stroke-width="0.9" style="cursor:pointer;" onclick="selectDicerDot('${{d.s}}')">
        <title>${{d.s}} (${{d.t}}) · ${{grp.label}} · log2(5p/3p) = ${{d.d53.toFixed(3)}}${{d.l7z!==null?' · let-7 Z = '+d.l7z.toFixed(2):''}}</title>
      </circle>`;
    }});
  }});
  svg += `</svg>
  <div class="plot-legend">
    <span class="leg-item"><i class="swatch" style="background:#c5221f;margin:0;"></i>Circle = 1 TCGA tumour (click dot to jump to sample row)</span>
    <span class="leg-item"><span style="font-weight:800;color:#202124;">&#9475;</span>Thick vertical bar = Class median</span>
    <span class="leg-item"><i class="leg-line" style="border-top:2px dashed #137333;height:0;width:14px;"></i>Pan-cancer wild-type median (+0.475, n = 9,129)</span>
    <span class="leg-item">Click any left-hand class label to filter table</span>
  </div>`;
  wrap.innerHTML = svg;
}}

/* Paste-a-List Batch Variant Triage Across All 213,940 Substitutions */
let currentBatchRows = [];
function openBatchTriage() {{
  switchTab('tab-variants', document.getElementById('btn-tab-variants'));
  switchSub('exp-sub', 'exp-batch', document.getElementById('seg-exp-batch'));
  if (!document.getElementById('batch-input').value.trim()) loadBatchExample();
  document.getElementById('exp-batch').scrollIntoView({{behavior: 'smooth', block: 'start'}});
}}

function loadBatchExample() {{
  document.getElementById('batch-input').value = [
    '##fileformat=VCFv4.2',
    '##reference=GRCh38',
    '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO',
    'chr14\t95094139\t.\tC\tT\t.\tPASS\tGENE=DICER1;HGVSp=p.Glu1705Lys',
    'chr14\t95090846\t.\tC\tA\t.\tPASS\tGENE=DICER1;HGVSp=p.Asp1810Tyr',
    'chr5\t31409348\t.\tC\tT\t.\tPASS\tGENE=DROSHA;HGVSp=p.Glu1147Lys',
    'chr17\t7675088\t.\tC\tT\t.\tPASS\tGENE=TP53;HGVSp=p.Arg175His',
    'chr14\t95108998\t.\tT\tC\t.\tPASS\tGENE=DICER1;HGVSp=p.Lys868Arg'
  ].join('\\n');
  runBatchTriage();
}}

function handleBatchFileSelect(files) {{
  if (!files || !files.length) return;
  const f = files[0];
  if (!/\\.vcf$/i.test(f.name || '')) {{
    document.getElementById('batch-summary').textContent = 'Rejected: Only uncompressed .vcf files are accepted for file upload.';
    return;
  }}
  const reader = new FileReader();
  reader.onload = ev => {{
    const txt = String(ev.target.result || '');
    const lines = txt.split(/\\r?\\n/).filter(l => l.trim() && !l.startsWith('#'));
    const kept = [];
    for (let i = 0; i < lines.length && kept.length < 300; i++) {{
      const u = resolveUniversalVariantQuery(lines[i]);
      if (u && !u.skipVcfNonPanel) kept.push(lines[i]);
    }}
    document.getElementById('batch-input').value = kept.join('\\n');
    if (!kept.length) {{
      document.getElementById('batch-summary').textContent = `Parsed ${{f.name}}: 0 coding SNVs overlapped the 15-gene panel.`;
      return;
    }}
    runBatchTriage();
  }};
  reader.readAsText(f);
}}

function handleBatchFileDrop(ev) {{
  ev.preventDefault();
  if (ev.dataTransfer && ev.dataTransfer.files && ev.dataTransfer.files.length) {{
    handleBatchFileSelect(ev.dataTransfer.files);
  }}
}}

function clearBatchTriage() {{
  document.getElementById('batch-input').value = '';
  document.getElementById('batch-summary').textContent = '';
  document.querySelector('#batch-table tbody').innerHTML = '';
  currentBatchRows = [];
}}

function parseBatchLines(rawText) {{
  const lines = (rawText || '').split(/\\r?\\n/);
  const tokens = [];
  lines.forEach(ln => {{
    const t = ln.trim();
    if (!t || t.startsWith('#')) return;
    if (t.includes('\\t')) tokens.push(t);
    else t.split(/[,;]+/).map(x => x.trim()).filter(Boolean).forEach(x => tokens.push(x));
  }});
  const resolved = [];
  for (let i = 0; i < tokens.length && resolved.length < 250; i++) {{
    const u = resolveUniversalVariantQuery(tokens[i]);
    if (u && !u.skipVcfNonPanel) resolved.push(u);
  }}
  return resolved;
}}

function runBatchTriage() {{
  const parsed = parseBatchLines(document.getElementById('batch-input').value);
  if (!parsed.length) return;
  const neededGenes = [...new Set(parsed.filter(p => p.gene && GENE_SET.has(p.gene)).map(p => p.gene))];
  document.getElementById('batch-summary').textContent = 'Scoring variants...';

  let pending = neededGenes.length;
  const finalize = () => {{
    let nScored = 0, nAct = 0, nUnf = 0, nNeu = 0;
    currentBatchRows = [];
    const tbodyHtml = parsed.map(item => {{
      if (item.err) {{
        currentBatchRows.push({{ input: item.raw, status: item.err }});
        return `<tr><td>—</td><td><code>${{item.raw}}</code></td><td>${{item.gene || '—'}}</td><td colspan="9" style="color:var(--crimson);">${{item.err}}</td></tr>`;
      }}
      const gCache = window.SAT_CACHE ? window.SAT_CACHE[item.gene] : null;
      const entry = gCache ? gCache[String(item.pos)] : null;
      if (!entry) {{
        currentBatchRows.push({{ input: item.raw, gene: item.gene, status: `Codon ${{item.pos}} out of range` }});
        return `<tr><td>—</td><td><code>${{item.raw}}</code></td><td><strong>${{item.gene}}</strong></td><td colspan="9" style="color:var(--crimson);">Codon ${{item.pos}} outside canonical length</td></tr>`;
      }}
      const [trueWt, pdb, dna, dmet] = entry.m;
      if (item.alt === '*') {{
        nScored++;
        currentBatchRows.push({{ input: item.raw, gene: item.gene, variant: `${{trueWt}}${{item.pos}}*`, mechanism: 'Truncating_LoF' }});
        return `<tr>
          <td><button class="btn-subtle" style="height:24px;padding:0 9px;font-size:11.5px;" onclick="openVariantIn3D('${{item.gene}}','${{item.pos}}')">View</button></td>
          <td><code>${{item.raw}}</code></td><td><strong>${{item.gene}}</strong></td><td><code>${{trueWt}}${{item.pos}}*</code></td>
          <td>${{tagDot('dot-slate','Truncating LoF')}}</td><td>—</td><td>—</td><td>—</td><td>—</td>
          <td>${{dmet ?? dna ?? '—'}}</td><td>—</td><td>—</td>
        </tr>`;
      }}
      const sub = entry.s.find(r => r[0] === item.alt);
      if (!sub) {{
        const msg = item.wt !== trueWt ? `Ref AA at ${{item.pos}} is ${{trueWt}} (not ${{item.wt}})` : 'Synonymous / wild-type';
        currentBatchRows.push({{ input: item.raw, gene: item.gene, status: msg }});
        return `<tr>
          <td><button class="btn-subtle" style="height:24px;padding:0 9px;font-size:11.5px;" onclick="openVariantIn3D('${{item.gene}}','${{trueWt}}${{item.pos}}${{trueWt}}')">View</button></td>
          <td><code>${{item.raw}}</code></td><td><strong>${{item.gene}}</strong></td><td><code>${{trueWt}}${{item.pos}}${{item.alt}}</code></td>
          <td colspan="8" style="color:var(--ink-2);">${{msg}}</td>
        </tr>`;
      }}
      const [alt, am, e2, e1v, if1, tn, cv] = sub;
      const gd = granthamDist(trueWt, alt);
      const cns = (am !== null && am >= 0.564) || (e2 !== null && e2 <= -5.0);
      let mechLbl = 'Neutral', mechHtml = tagDot('dot-petrol', 'Neutral');
      if (cns) {{
        if (if1 !== null && if1 > -1.5) {{ mechLbl = 'Fold-preserving active site'; mechHtml = tagDot('dot-crimson', 'Active site'); nAct++; }}
        else {{ mechLbl = 'Fold-destabilising'; mechHtml = tagDot('dot-amber', 'Unfolding'); nUnf++; }}
      }} else {{ nNeu++; }}
      nScored++;
      const minPkt = [dna, dmet].filter(x => x !== null);
      const pktStr = minPkt.length ? Math.min(...minPkt).toFixed(2) : '—';
      const gRev = reverseMapVariantToGenomic(item.gene, item.pos, trueWt, alt);
      const cdnaTag = item.cdna || (gRev ? gRev.cdna : '');
      const vcfTag = item.vcf || (gRev ? gRev.snvVcf : '');
      const wtWarn = item.wt !== trueWt ? ` <span style="color:var(--amber);font-size:11px;" title="Canonical WT residue is ${{trueWt}}">(WT=${{trueWt}})</span>` : '';
      currentBatchRows.push({{
        input: item.raw, grch38_vcf: vcfTag, cdna: cdnaTag, gene: item.gene, variant: `${{trueWt}}${{item.pos}}${{alt}}`,
        mechanism: mechLbl, alphamissense: am, esm_if1_3d_dlogp: if1, esm2_3b_llr: e2,
        grantham_D: gd, min_pocket_A: pktStr, pdb: pdb || 'AFDB_v6', tcga_n: tn, clinvar: cv || ''
      }});
      return `<tr>
        <td><button class="btn-subtle" style="height:24px;padding:0 9px;font-size:11.5px;" onclick="openVariantIn3D('${{item.gene}}','${{trueWt}}${{item.pos}}${{alt}}')">View</button></td>
        <td><code>${{item.raw}}</code>${{item.assembly ? ` <span class="meta-note">[${{item.assembly}}]</span>` : ''}}</td>
        <td><strong>${{item.gene}}</strong></td>
        <td><code>${{trueWt}}${{item.pos}}${{alt}}</code>${{cdnaTag ? ` <span class="meta-note">(${{cdnaTag}})</span>` : ''}}${{wtWarn}}</td>
        <td>${{mechHtml}}</td>
        <td><strong>${{am !== null ? am.toFixed(3) : '—'}}</strong></td>
        <td>${{if1 !== null ? if1.toFixed(2) : '—'}}</td>
        <td>${{e2 !== null ? e2.toFixed(2) : '—'}}</td>
        <td>${{gd}}</td>
        <td>${{pktStr}}${{pdb ? ' ('+pdb.split(':')[0]+')' : ''}}</td>
        <td>${{tn > 0 ? '<strong>'+tn+'</strong>' : '0'}}</td>
        <td>${{cv || '—'}}</td>
      </tr>`;
    }}).join('');

    document.querySelector('#batch-table tbody').innerHTML = tbodyHtml;
    document.getElementById('batch-summary').textContent =
      `${{nScored}} / ${{parsed.length}} scored · ${{nAct}} active-site · ${{nUnf}} unfolding · ${{nNeu}} neutral`;
    markSortableHeaders(document.getElementById('batch-table'));
  }};

  if (!pending) {{ finalize(); return; }}
  neededGenes.forEach(g => ensureSatLoaded(g, () => {{
    pending--;
    if (pending <= 0) finalize();
  }}));
}}

function exportBatchTriageCSV() {{
  if (!currentBatchRows.length) return;
  downloadCSV('miRBioGenAtlas_batch_triage_results.csv', currentBatchRows);
}}

/* Universal One-Click Table Column Sorting */
const NUM_HEAD_RE = /^[+-]?\\d+(?:\\.\\d+)?(?:e[+-]?\\d+)?(?:\\s*(?:%|Å|Å³|e\\b|\\(|\\[|\\/|\\$))?/i;
function parseCellNumeric(raw) {{
  const s = (raw || '').trim().replace(/\\u2212/g, '-').replace(/,/g, '');
  if (!s || s === '—' || s === 'N/A') return null;
  if (!/^[+-]?\\d+(?:\\.\\d+)?(?:e[+-]?\\d+)?(?:\\s*(?:%|Å|Å³|e\\b|\\(|\\[|\\/|$))/i.test(s)) return NaN;
  const m = s.match(/^[+-]?\\d+(?:\\.\\d+)?(?:e[+-]?\\d+)?/i);
  return m ? parseFloat(m[0]) : NaN;
}}

function markSortableHeaders(rootEl) {{
  (rootEl || document).querySelectorAll('table:not(.no-auto-sort) thead th').forEach(th => {{
    const txt = th.textContent.trim();
    if (!txt || txt === '3D' || txt === 'Relative Level') return;
    th.classList.add('sortable');
  }});
}}

function sortTableByTh(th) {{
  const table = th.closest('table');
  const tbody = table ? table.querySelector('tbody') : null;
  if (!tbody) return;
  const ths = Array.from(th.parentElement.children);
  const colIdx = ths.indexOf(th);
  if (colIdx < 0) return;
  const rows = Array.from(tbody.querySelectorAll('tr'));
  if (rows.length < 2 || rows.some(r => r.children.length <= colIdx)) return;

  const curDir = th.classList.contains('sort-asc') ? 'asc' : (th.classList.contains('sort-desc') ? 'desc' : '');
  let numValid = 0, nonBlank = 0;
  const parsed = rows.map((row, idx) => {{
    const txt = (row.children[colIdx] ? row.children[colIdx].textContent : '').trim();
    const nv = parseCellNumeric(txt);
    if (nv !== null) {{
      nonBlank++;
      if (!Number.isNaN(nv)) numValid++;
    }}
    return {{ row, idx, txt, nv }};
  }});
  const isNumeric = nonBlank > 0 && (numValid / nonBlank) >= 0.6;
  const hdrLow = th.textContent.toLowerCase();
  const preferAsc = !isNumeric || /(rank|gene|code|study|variant|sample|subset|rna|metal|pocket|loeuf|\\bp\\b|two-sided p|cox p)/i.test(hdrLow);
  const nextDir = curDir ? (curDir === 'asc' ? 'desc' : 'asc') : (preferAsc ? 'asc' : 'desc');

  ths.forEach(h => h.classList.remove('sort-asc', 'sort-desc'));
  th.classList.add(nextDir === 'asc' ? 'sort-asc' : 'sort-desc');

  const sign = nextDir === 'asc' ? 1 : -1;
  parsed.sort((a, b) => {{
    const aBlank = (a.nv === null || a.txt === '' || a.txt === '—');
    const bBlank = (b.nv === null || b.txt === '' || b.txt === '—');
    if (aBlank && bBlank) return a.idx - b.idx;
    if (aBlank) return 1;
    if (bBlank) return -1;
    if (isNumeric && !Number.isNaN(a.nv) && !Number.isNaN(b.nv)) {{
      if (a.nv !== b.nv) return (a.nv - b.nv) * sign;
      return a.idx - b.idx;
    }}
    const cmp = a.txt.localeCompare(b.txt, undefined, {{ numeric: true, sensitivity: 'base' }});
    return cmp !== 0 ? cmp * sign : (a.idx - b.idx);
  }});
  const frag = document.createDocumentFragment();
  parsed.forEach(p => frag.appendChild(p.row));
  tbody.appendChild(frag);
}}

function initUniversalTableSorting() {{
  markSortableHeaders(document);
  document.addEventListener('click', e => {{
    const th = e.target.closest('th.sortable');
    if (th) sortTableByTh(th);
  }});
}}

/* Deep-Link URL Hash & Query Parameter State Restoration */
function applyUrlState() {{
  suppressHashSync = true;
  try {{
    const params = new URLSearchParams(window.location.search || '');
    const qParam = params.get('q');
    if (qParam) {{
      document.getElementById('global-search').value = qParam;
      runGlobalSearch(qParam);
      return;
    }}
    const rawHash = decodeURIComponent((window.location.hash || '').replace(/^#/, '').trim());
    if (!rawHash) return;
    if (rawHash === 'about-atlas') {{
      const ab = document.getElementById('about-atlas');
      if (ab) ab.scrollIntoView({{behavior: 'smooth'}});
      return;
    }}
    const kv = {{}};
    rawHash.split('&').forEach(pair => {{
      const eq = pair.indexOf('=');
      if (eq > 0) kv[pair.slice(0, eq).toLowerCase()] = pair.slice(eq + 1);
    }});
    if (kv.variant) {{
      const parts = kv.variant.split(/[:\\s]+/);
      if (parts[0] && GENE_SET.has(parts[0].toUpperCase())) {{
        loadPresetMutation(parts[0].toUpperCase(), parts[1] || '');
        return;
      }}
    }}
    if (kv.gene && GENE_SET.has(kv.gene.toUpperCase())) {{
      loadPresetMutation(kv.gene.toUpperCase(), '');
      return;
    }}
    if (kv.eq) {{
      openFormulaTab(kv.eq.toUpperCase());
      return;
    }}
    if (kv.tab) {{
      const tMap = {{
        explorer: ['tab-variants', 'btn-tab-variants'],
        pancancer: ['tab-heatmap', 'btn-tab-heatmap'],
        mirna: ['tab-mirna', 'btn-tab-mirna'],
        clinical: ['tab-clinical', 'btn-tab-clinical'],
        docs: ['tab-docs', 'btn-tab-docs']
      }};
      const hit = tMap[kv.tab.toLowerCase()];
      if (hit) switchTab(hit[0], document.getElementById(hit[1]), true);
      if (kv.tab.toLowerCase() === 'batch') openBatchTriage();
      if (kv.tab.toLowerCase() === 'dms') openDmsLandscape();
      if (kv.tab.toLowerCase() === 'api') openApiDocs();
    }}
  }} finally {{
    setTimeout(() => {{ suppressHashSync = false; }}, 120);
  }}
}}

function openApiDocs() {{
  switchTab('tab-docs', document.getElementById('btn-tab-docs'), true);
  switchSub('doc-sub', 'doc-api', document.getElementById('seg-api'));
  renderApiSnippet();
  syncUrlHash('#tab=api');
}}

function renderApiSnippet() {{
  const sel = document.getElementById('api-gene-sel');
  const g = (sel && sel.value) ? sel.value : 'DICER1';
  const baseUrl = 'https://siddugayam.github.io/miRBioGenAtlas';
  const relPath = `api/v2/genes/${{g}}.json`;
  const pathEl = document.getElementById('api-gene-path');
  const openBtn = document.getElementById('api-open-link');
  const rowLnk = document.getElementById('api-gene-row-link');
  if (pathEl) pathEl.textContent = relPath;
  if (openBtn) openBtn.href = relPath;
  if (rowLnk) rowLnk.href = relPath;
  const codeBox = document.getElementById('api-code-box');
  if (!codeBox) return;
  codeBox.textContent = [
    `# 1. Fetch ${{g}} JSON profile (metadata, GRCh38/37 CDS exons, TCGA variants, ClinVar & 3D track)`,
    `curl -sL ${{baseUrl}}/${{relPath}} | jq '.metadata.gnomad_LOEUF, .tcga_somatic_variants[:2]'`,
    ``,
    `# 2. Python (urllib / requests) — zero dependencies`,
    `import json, urllib.request`,
    `data = json.load(urllib.request.urlopen("${{baseUrl}}/${{relPath}}"))`,
    `print(data["geneSymbol"], len(data["tcga_somatic_variants"]), "observed TCGA variants")`,
    ``,
    `# 3. Local CLI — query by protein, GRCh38/GRCh37 VCF coordinate, or cDNA HGVS`,
    `python3 scripts/query_variant.py DICER1 E1705K`,
    `python3 scripts/query_variant.py chr14:95094139:C>T --json`
  ].join('\\n');
}}

function copyApiCode() {{
  const box = document.getElementById('api-code-box');
  const btn = document.getElementById('btn-copy-api');
  if (!box || !navigator.clipboard) return;
  const orig = btn.textContent;
  navigator.clipboard.writeText(box.textContent).then(() => {{
    btn.textContent = 'Copied ✓';
    setTimeout(() => {{ btn.textContent = orig; }}, 1500);
  }});
}}

/* MaveDB / dms-viz Style 2D Deep Mutational Scanning Heatmap Synced to 3D */
const DMS_AA_ORDER = ['G','A','V','L','I','M','F','Y','W','P','S','T','C','N','Q','H','K','R','D','E'];
function openDmsLandscape() {{
  switchTab('tab-variants', document.getElementById('btn-tab-variants'));
  switchSub('exp-sub', 'exp-dms', document.getElementById('seg-exp-dms'));
  renderDmsHeatmap();
  document.getElementById('exp-dms').scrollIntoView({{behavior: 'smooth', block: 'start'}});
}}

function stepDmsWindow(delta) {{
  const gene = document.getElementById('lol-gene').value;
  const gCache = window.SAT_CACHE ? window.SAT_CACHE[gene] : null;
  if (!gCache) return;
  const L = Object.keys(gCache).length;
  dmsCenterPos = Math.min(L, Math.max(1, (dmsCenterPos || active3DState.pos || 1) + delta));
  renderDmsHeatmap();
}}

function jumpDmsToDomain(val) {{
  if (!val) return;
  const [s, e] = val.split('-').map(Number);
  if (s && e) {{
    dmsCenterPos = Math.round((s + e) / 2);
    renderDmsHeatmap();
  }}
}}

function dmsScoreColor(metric, val) {{
  if (val === null || val === undefined) return '#eef1f4';
  if (metric === 'am') {{
    if (val <= 0.45) {{
      const t = Math.max(0, val / 0.45);
      return `rgb(${{Math.round(35 + t*215)}},${{Math.round(118 + t*132)}},${{Math.round(178 + t*72)}})`;
    }}
    const t = Math.min(1, (val - 0.45) / 0.55);
    return `rgb(${{Math.round(250 - t*53)}},${{Math.round(245 - t*211)}},${{Math.round(245 - t*214)}})`;
  }}
  if (metric === 'if1') {{
    const t = Math.min(1, Math.max(0, -val / 4.2));
    return `rgb(${{Math.round(244 - t*68)}},${{Math.round(247 - t*151)}},${{Math.round(250 - t*240)}})`;
  }}
  const t = Math.min(1, Math.max(0, -val / 11.0));
  return `rgb(${{Math.round(244 - t*47)}},${{Math.round(247 - t*213)}},${{Math.round(250 - t*219)}})`;
}}

function renderDmsHeatmap() {{
  const wrap = document.getElementById('dms-svg-wrap');
  if (!wrap) return;
  const gene = document.getElementById('lol-gene').value;
  const gCache = window.SAT_CACHE ? window.SAT_CACHE[gene] : null;
  if (!gCache) {{ ensureSatLoaded(gene, renderDmsHeatmap); return; }}

  const L = Object.keys(gCache).length;
  const half = parseInt(document.getElementById('dms-span').value || '25', 10);
  const metric = document.getElementById('dms-metric').value || 'am';
  const dmsLeg = document.getElementById('dms-legend-bar');
  if (dmsLeg) {{
    let scaleHtml = '';
    if (metric === 'am') {{
      scaleHtml = `<span class="leg-item"><strong>AlphaMissense:</strong> <code>0.0</code> Benign <span style="display:inline-block;width:96px;height:8px;border-radius:999px;border:1px solid var(--line);background:linear-gradient(90deg,rgb(35,118,178) 0%,rgb(250,250,250) 45%,rgb(197,34,31) 100%);"></span> <code>1.0</code> Pathogenic (&ge;0.564 cutoff)</span>`;
    }} else if (metric === 'if1') {{
      scaleHtml = `<span class="leg-item"><strong>ESM-IF1 3D &Delta;logP:</strong> <code>0.0</code> Fold Intact <span style="display:inline-block;width:96px;height:8px;border-radius:999px;border:1px solid var(--line);background:linear-gradient(90deg,rgb(244,247,250) 0%,rgb(176,96,10) 100%);"></span> <code>&le;&minus;4.2</code> Unfolding (&le;&minus;1.5 cutoff)</span>`;
    }} else {{
      scaleHtml = `<span class="leg-item"><strong>ESM-2 &Delta;&ell;:</strong> <code>0.0</code> Tolerated <span style="display:inline-block;width:96px;height:8px;border-radius:999px;border:1px solid var(--line);background:linear-gradient(90deg,rgb(244,247,250) 0%,rgb(197,34,31) 100%);"></span> <code>&le;&minus;11.0</code> Constrained (&le;&minus;5.0 cutoff)</span>`;
    }}
    dmsLeg.innerHTML = `${{scaleHtml}}<span class="meta-note"><code>&bull;</code> Wild-type AA &middot; <span style="display:inline-block;width:9px;height:9px;border:1.3px solid #202124;border-radius:2px;vertical-align:-1px;"></span> Observed in TCGA / ClinVar P/LP &middot; Column box = Active 3D codon &middot; Top bar = UniProt domain</span>`;
  }}
  const center = Math.min(L, Math.max(1, dmsCenterPos || active3DState.pos || 1));
  let startP = Math.max(1, center - half);
  let endP = Math.min(L, startP + half * 2);
  startP = Math.max(1, endP - half * 2);
  const nCols = endP - startP + 1;

  const cw = half > 30 ? 16 : 24, ch = 18, padL = 42, padT = 32, padB = 26, padR = 16;
  const W = padL + nCols * cw + padR, H = padT + 20 * ch + padB;

  let cells = '', colHdrs = '';
  for (let p = startP; p <= endP; p++) {{
    const entry = gCache[String(p)];
    if (!entry) continue;
    const x = padL + (p - startP) * cw;
    const wt = entry.m[0];
    const dm = domainAtPos(gene, p);
    colHdrs += `<rect x="${{x}}" y="8" width="${{cw - 1}}" height="5" fill="${{dm ? dm.color : '#dfe3e8'}}"><title>${{dm ? dm.name : 'Inter-domain'}} · aa ${{p}}</title></rect>`;
    if (p === active3DState.pos) {{
      colHdrs += `<rect x="${{x - 0.5}}" y="${{padT - 1}}" width="${{cw}}" height="${{20 * ch + 2}}" fill="none" stroke="#202124" stroke-width="1.5"/>`;
    }}
    if (p === startP || p === endP || p % 5 === 0) {{
      colHdrs += `<text x="${{(x + cw/2).toFixed(1)}}" y="25" font-size="9.5" font-family="monospace" fill="${{p===active3DState.pos?'#c5221f':'#5f6368'}}" font-weight="${{p===active3DState.pos?'700':'400'}}" text-anchor="middle">${{p}}</text>`;
    }}
    const subMap = {{}};
    entry.s.forEach(r => {{ subMap[r[0]] = r; }});
    DMS_AA_ORDER.forEach((aa, ri) => {{
      const y = padT + ri * ch;
      if (aa === wt) {{
        cells += `<g style="cursor:pointer" onclick="selectMutAA(${{p}},'${{wt}}','${{wt}}')">
          <rect x="${{x}}" y="${{y}}" width="${{cw - 1}}" height="${{ch - 1}}" fill="#f8f9fa"/>
          <circle cx="${{(x + (cw-1)/2).toFixed(1)}}" cy="${{(y + (ch-1)/2).toFixed(1)}}" r="2.1" fill="#5f6368"/>
          <title>${{gene}} ${{wt}}${{p}} (Wild-type)</title>
        </g>`;
        return;
      }}
      const sr = subMap[aa];
      const val = sr ? (metric === 'am' ? sr[1] : (metric === 'if1' ? sr[4] : sr[2])) : null;
      const bg = dmsScoreColor(metric, val);
      const hasObs = sr && (sr[5] > 0 || (sr[6] && sr[6].toLowerCase().includes('pathogenic')));
      const strokeAttr = (p === active3DState.pos && aa === active3DState.mut)
        ? 'stroke="#111315" stroke-width="1.8"'
        : (hasObs ? 'stroke="#202124" stroke-width="1.1"' : '');
      cells += `<rect x="${{x}}" y="${{y}}" width="${{cw - 1}}" height="${{ch - 1}}" rx="1.5" fill="${{bg}}" ${{strokeAttr}} style="cursor:pointer" onclick="selectMutAA(${{p}},'${{wt}}','${{aa}}')">
        <title>${{gene}} ${{wt}}${{p}}${{aa}} · AM=${{sr&&sr[1]!==null?sr[1].toFixed(3):'—'}} · IF1=${{sr&&sr[4]!==null?sr[4].toFixed(2):'—'}} · ESM2=${{sr&&sr[2]!==null?sr[2].toFixed(2):'—'}}${{sr&&sr[5]>0?' · TCGA n='+sr[5]:''}}${{sr&&sr[6]?' · ClinVar: '+sr[6]:''}}</title>
      </rect>`;
    }});
  }}

  const rowLbls = DMS_AA_ORDER.map((aa, ri) =>
    `<text x="${{padL - 8}}" y="${{padT + ri * ch + 13}}" font-size="10.5" font-family="monospace" font-weight="600" fill="#202124" text-anchor="end">${{aa}}</text>`
  ).join('');

  wrap.innerHTML = `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${{W}} ${{H}}" width="100%" height="${{H}}">
    ${{rowLbls}}${{cells}}${{colHdrs}}
    <text x="${{padL}}" y="${{H - 7}}" font-size="10.5" fill="#80868b">${{gene}} codons ${{startP}}–${{endP}} of ${{L}} aa</text>
  </svg>`;
}}

function downloadCSV(filename, rows) {{
  if (!rows.length) return;
  const keys = Object.keys(rows[0]);
  const csv = [keys.join(','), ...rows.map(r => keys.map(k => JSON.stringify(r[k] ?? '')).join(','))].join('\\n');
  const a = document.createElement('a');
  a.href = URL.createObjectURL(new Blob([csv], {{type: 'text/csv;charset=utf-8;'}}));
  a.download = filename; a.click();
}}
function exportLongCSV() {{ downloadCSV('miRBioGenAtlas_480_gene_study_summary.csv', ATLAS.long); }}
function exportFilteredVariantsCSV() {{ downloadCSV('miRBioGenAtlas_filtered_3D_variants.csv', currentFilteredVars); }}
function exportDicerMirnaCSV() {{ downloadCSV('miRBioGenAtlas_DICER1_mirnaseq_5p3p.csv', ATLAS.dicer_pts); }}

function copyCitationText(fmt) {{
  const plain = 'Gayam, P. K. R., Earny, V. A., Aranjani, J. M. & Mahammad, F. S. miRBioGenAtlas: Pan-cancer alteration landscape of the miRNA biogenesis machinery across 32 TCGA cancer studies (v2.0.2). Zenodo https://doi.org/10.5281/zenodo.24455067 (2026).';
  const bib = `@dataset{{gayam_2026_mirbiogenatlas,\\n  author    = {{Gayam, Prasanna Kumar Reddy and Earny, Venkat Abhiram and Aranjani, Jesil Mathew and Mahammad, Fayaz Shaik}},\\n  title     = {{{{miRBioGenAtlas: Pan-cancer alteration landscape of the miRNA biogenesis machinery across 32 TCGA cancer studies}}}},\\n  year      = {{2026}},\\n  publisher = {{Zenodo}},\\n  version   = {{2.0.2}},\\n  doi       = {{10.5281/zenodo.24455067}},\\n  url       = {{https://doi.org/10.5281/zenodo.24455067}}\\n}}`;
  const txt = fmt === 'bibtex' ? bib : plain;
  const btn = document.getElementById(fmt === 'bibtex' ? 'btn-copy-bib' : 'btn-copy-cite');
  const orig = btn.textContent;
  navigator.clipboard.writeText(txt).then(() => {{
    btn.textContent = 'Copied ✓';
    setTimeout(() => {{ btn.textContent = orig; }}, 1600);
  }});
}}

window.addEventListener('DOMContentLoaded', () => {{
  const selVar = document.getElementById('var-gene');
  const selLol = document.getElementById('lol-gene');
  const selPub = document.getElementById('pubbar-gene');
  const selGtex = document.getElementById('gtex-gene');
  const selApi = document.getElementById('api-gene-sel');
  ATLAS.panel.forEach(p => {{
    const g = p.geneSymbol;
    selVar.insertAdjacentHTML('beforeend', `<option value="${{g}}">${{g}} (${{p.uniprotAccession}})</option>`);
    selLol.insertAdjacentHTML('beforeend', `<option value="${{g}}" ${{g==='DICER1'?'selected':''}}>${{g}} · ${{p.uniprotAccession}} (${{p.proteinLengthAa}} aa)</option>`);
    if (selPub) selPub.insertAdjacentHTML('beforeend', `<option value="${{g}}" ${{g==='DICER1'?'selected':''}}>${{g}} (${{p.uniprotAccession}})</option>`);
    selGtex.insertAdjacentHTML('beforeend', `<option value="${{g}}" ${{g==='LIN28B'?'selected':''}}>${{g}} (&tau; = ${{p.gtex_tau_tissue_specificity}})</option>`);
    if (selApi) selApi.insertAdjacentHTML('beforeend', `<option value="${{g}}" ${{g==='DICER1'?'selected':''}}>${{g}}.json</option>`);
  }});
  document.getElementById('sat-query').value = 'E1705K';
  renderHeatmap();
  renderVariants();
  renderAllTables();
  renderDicerStripPlot();
  renderGtexBars();
  renderFormulas();
  renderDatasetsCatalogue();
  renderApiSnippet();
  renderGlossaryCards();
  initUniversalTableSorting();
  renderLollipopAndLoadSat();
  applyUrlState();
  window.addEventListener('hashchange', applyUrlState);
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
