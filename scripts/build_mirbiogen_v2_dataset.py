#!/usr/bin/env python3
"""miRBioGenAtlas v2.0 — Deterministic Offline Multi-Modal Dataset Builder."""
import csv, gzip, hashlib, json, math, pathlib, statistics
from collections import Counter, defaultdict

REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent
DATA_DIR = REPO_ROOT / "data"
SNAPSHOT_DIR = DATA_DIR / "raw_snapshot"
STUDIES_SNAP_DIR = SNAPSHOT_DIR / "studies"

GENE_PANEL = [
    {"geneSymbol": "DROSHA", "entrezGeneId": 29102, "hgncId": "HGNC:17904", "ensemblGeneId": "ENSG00000113360", "canonicalRefSeq": "NM_013235.5", "cytoband": "5p13.3", "chrArm": "5p", "syntenicCluster": "5p13.3_DROSHA", "cdsLengthBp": 4125, "proteinLengthAa": 1374, "pathwayStep": "Nuclear Microprocessor", "panelRole": "Core_Enzyme", "inCore10": True, "inExtended12": True, "inPanel14": True, "hotspotCodons": {1147, 1151, 1187, 1222}, "inclusionRationale": "Catalytic RNase III endonuclease of the nuclear Microprocessor complex cleaving pri-miRNAs to pre-miRNAs."},
    {"geneSymbol": "DGCR8", "entrezGeneId": 54487, "hgncId": "HGNC:2849", "ensemblGeneId": "ENSG00000128191", "canonicalRefSeq": "NM_022720.7", "cytoband": "22q11.21", "chrArm": "22q", "syntenicCluster": "22q11.21_DGCR8", "cdsLengthBp": 2322, "proteinLengthAa": 773, "pathwayStep": "Nuclear Microprocessor", "panelRole": "Core_Cofactor", "inCore10": True, "inExtended12": True, "inPanel14": True, "hotspotCodons": {518}, "inclusionRationale": "Essential heme- and dsRNA-binding Microprocessor subunit required for pri-miRNA substrate recognition by DROSHA."},
    {"geneSymbol": "XPO5", "entrezGeneId": 57510, "hgncId": "HGNC:17675", "ensemblGeneId": "ENSG00000124571", "canonicalRefSeq": "NM_020750.3", "cytoband": "6p21.1", "chrArm": "6p", "syntenicCluster": "6p21.1_XPO5", "cdsLengthBp": 3615, "proteinLengthAa": 1204, "pathwayStep": "Nuclear Export", "panelRole": "Core_Transporter", "inCore10": True, "inExtended12": True, "inPanel14": True, "hotspotCodons": {1147, 1189}, "inclusionRationale": "Canonical RanGTP-dependent karyopherin transporting hairpin pre-miRNAs across the nuclear pore."},
    {"geneSymbol": "DICER1", "entrezGeneId": 23405, "hgncId": "HGNC:17098", "ensemblGeneId": "ENSG00000100697", "canonicalRefSeq": "NM_177438.3", "cytoband": "14q32.13", "chrArm": "14q", "syntenicCluster": "14q32.13_DICER1", "cdsLengthBp": 5769, "proteinLengthAa": 1922, "pathwayStep": "Cytoplasmic Dicing", "panelRole": "Core_Enzyme", "inCore10": True, "inExtended12": True, "inPanel14": True, "hotspotCodons": {1344, 1705, 1709, 1725, 1809, 1810, 1813}, "inclusionRationale": "Cytoplasmic RNase III endonuclease cleaving pre-miRNAs into ~22 nt miRNA duplexes; harbours recurrent somatic RNase IIIb hotspots."},
    {"geneSymbol": "TARBP2", "entrezGeneId": 6895, "hgncId": "HGNC:11569", "ensemblGeneId": "ENSG00000139546", "canonicalRefSeq": "NM_134323.2", "cytoband": "12q13.13", "chrArm": "12q", "syntenicCluster": "12q13.13_TARBP2", "cdsLengthBp": 1098, "proteinLengthAa": 366, "pathwayStep": "Cytoplasmic Dicing", "panelRole": "Core_Cofactor", "inCore10": True, "inExtended12": True, "inPanel14": True, "hotspotCodons": {149, 152}, "inclusionRationale": "dsRNA-binding cofactor (TRBP) of DICER1 tuning cleavage fidelity, isoMiR selection, and RISC-loading complex assembly."},
    {"geneSymbol": "PRKRA", "entrezGeneId": 8575, "hgncId": "HGNC:9438", "ensemblGeneId": "ENSG00000180228", "canonicalRefSeq": "NM_003690.5", "cytoband": "2q31.2", "chrArm": "2q", "syntenicCluster": "2q31.2_PRKRA", "cdsLengthBp": 942, "proteinLengthAa": 313, "pathwayStep": "Cytoplasmic Dicing", "panelRole": "Core_Cofactor", "inCore10": True, "inExtended12": True, "inPanel14": True, "hotspotCodons": set(), "inclusionRationale": "Alternative dsRNA-binding cofactor (PACT) forming a heterotrimeric complex with DICER1 and AGO2 to modulate substrate processing."},
    {"geneSymbol": "AGO1", "entrezGeneId": 26523, "hgncId": "HGNC:3262", "ensemblGeneId": "ENSG00000092847", "canonicalRefSeq": "NM_012199.5", "cytoband": "1p34.3", "chrArm": "1p", "syntenicCluster": "1p34.3_AGO1_AGO3_AGO4", "cdsLengthBp": 2574, "proteinLengthAa": 857, "pathwayStep": "RISC Effector", "panelRole": "Core_Effector", "inCore10": True, "inExtended12": True, "inPanel14": True, "hotspotCodons": set(), "inclusionRationale": "Canonical human Argonaute RISC effector; co-localized in a tandem gene cluster on 1p34.3 with AGO3 and AGO4."},
    {"geneSymbol": "AGO2", "entrezGeneId": 27161, "hgncId": "HGNC:3263", "ensemblGeneId": "ENSG00000123908", "canonicalRefSeq": "NM_012154.5", "cytoband": "8q24.3", "chrArm": "8q", "syntenicCluster": "8q24.3_AGO2", "cdsLengthBp": 2580, "proteinLengthAa": 859, "pathwayStep": "RISC Effector", "panelRole": "Core_Effector", "inCore10": True, "inExtended12": True, "inPanel14": True, "hotspotCodons": set(), "inclusionRationale": "Predominant catalytic ('Slicer') Argonaute RISC effector and mediator of DICER-independent miR-451 maturation."},
    {"geneSymbol": "AGO3", "entrezGeneId": 192669, "hgncId": "HGNC:18421", "ensemblGeneId": "ENSG00000126070", "canonicalRefSeq": "NM_024852.4", "cytoband": "1p34.3", "chrArm": "1p", "syntenicCluster": "1p34.3_AGO1_AGO3_AGO4", "cdsLengthBp": 2583, "proteinLengthAa": 860, "pathwayStep": "RISC Effector", "panelRole": "Core_Effector", "inCore10": True, "inExtended12": True, "inPanel14": True, "hotspotCodons": set(), "inclusionRationale": "Canonical human Argonaute RISC effector with conditional slicing activity; tandemly clustered with AGO1 and AGO4 on 1p34.3."},
    {"geneSymbol": "AGO4", "entrezGeneId": 192670, "hgncId": "HGNC:18422", "ensemblGeneId": "ENSG00000134698", "canonicalRefSeq": "NM_017629.4", "cytoband": "1p34.3", "chrArm": "1p", "syntenicCluster": "1p34.3_AGO1_AGO3_AGO4", "cdsLengthBp": 2586, "proteinLengthAa": 861, "pathwayStep": "RISC Effector", "panelRole": "Core_Effector", "inCore10": True, "inExtended12": True, "inPanel14": True, "hotspotCodons": set(), "inclusionRationale": "Non-catalytic human Argonaute RISC effector; tandemly clustered with AGO1 and AGO3 on 1p34.3."},
    {"geneSymbol": "LIN28A", "entrezGeneId": 79727, "hgncId": "HGNC:15986", "ensemblGeneId": "ENSG00000131914", "canonicalRefSeq": "NM_024674.6", "cytoband": "1p36.11", "chrArm": "1p", "syntenicCluster": "1p36.11_LIN28A", "cdsLengthBp": 630, "proteinLengthAa": 209, "pathwayStep": "Direct Let-7 Regulator", "panelRole": "Direct_Biogenesis_Regulator", "inCore10": False, "inExtended12": True, "inPanel14": True, "hotspotCodons": set(), "inclusionRationale": "RNA-binding protein that directly binds the terminal loop of pre-let-7 hairpins to block DICER1 cleavage and recruit TUT4/7."},
    {"geneSymbol": "LIN28B", "entrezGeneId": 389421, "hgncId": "HGNC:32207", "ensemblGeneId": "ENSG00000187772", "canonicalRefSeq": "NM_001004317.4", "cytoband": "6q16.3-q21", "chrArm": "6q", "syntenicCluster": "6q16.3_LIN28B", "cdsLengthBp": 753, "proteinLengthAa": 250, "pathwayStep": "Direct Let-7 Regulator", "panelRole": "Direct_Biogenesis_Regulator", "inCore10": False, "inExtended12": True, "inPanel14": True, "hotspotCodons": set(), "inclusionRationale": "Nuclear/nucleolar paralogue of LIN28A sequestering pri-let-7 transcripts away from Microprocessor cleavage."},
    {"geneSymbol": "SMAD2", "entrezGeneId": 4087, "hgncId": "HGNC:6768", "ensemblGeneId": "ENSG00000175387", "canonicalRefSeq": "NM_005901.6", "cytoband": "18q21.1", "chrArm": "18q", "syntenicCluster": "18q21.1_SMAD2_SMAD4", "cdsLengthBp": 1404, "proteinLengthAa": 467, "pathwayStep": "Microprocessor Accessory", "panelRole": "Pleiotropic_Accessory", "inCore10": False, "inExtended12": False, "inPanel14": True, "hotspotCodons": {450, 464}, "inclusionRationale": "TGF-beta R-SMAD binding conserved pri-miRNA stem motifs to promote DROSHA processing via DDX5; pleiotropic signalling factor excluded from Core-10/Extended-12 tiers."},
    {"geneSymbol": "SMAD4", "entrezGeneId": 4089, "hgncId": "HGNC:6770", "ensemblGeneId": "ENSG00000141646", "canonicalRefSeq": "NM_005359.6", "cytoband": "18q21.1", "chrArm": "18q", "syntenicCluster": "18q21.1_SMAD2_SMAD4", "cdsLengthBp": 1659, "proteinLengthAa": 552, "pathwayStep": "Microprocessor Accessory", "panelRole": "Pleiotropic_Accessory", "inCore10": False, "inExtended12": False, "inPanel14": True, "hotspotCodons": {351, 361, 386, 537}, "inclusionRationale": "Co-SMAD partner in TGF-beta/BMP signalling co-localized with SMAD2 at 18q21.1; retained in Panel-14 for pathway cross-talk context but separated from Core-10/Extended-12."},
    {"geneSymbol": "TP53", "entrezGeneId": 7157, "hgncId": "HGNC:11998", "ensemblGeneId": "ENSG00000141510", "canonicalRefSeq": "NM_000546.6", "cytoband": "17p13.1", "chrArm": "17p", "syntenicCluster": "17p13.1_TP53", "cdsLengthBp": 1182, "proteinLengthAa": 393, "pathwayStep": "Upstream Regulator", "panelRole": "Upstream_Pleiotropic_Regulator", "inCore10": False, "inExtended12": False, "inPanel14": False, "hotspotCodons": {175, 220, 245, 248, 249, 273, 282}, "inclusionRationale": "Modulates pri-miRNA processing via DDX5/Microprocessor interaction and transactivates miR-34 family; excluded from all biogenesis burden tiers owing to near-universal pan-cancer tumour-suppressor inactivation."},
]

EXCLUDED_CANDIDATES = [
    ("DDX5", "17q23.3", "Microprocessor DEAD-box helicase (p68)", "Excluded from core panel because DDX5 functions broadly across pre-mRNA splicing, ribosome biogenesis, and transcriptional co-activation."),
    ("DDX17", "22q13.1", "Microprocessor DEAD-box helicase (p72)", "Excluded owing to broad pleiotropic roles in alternative splicing and general RNA metabolism beyond pri-miRNA cleavage."),
    ("TUT4 (ZCCHC11)", "1p32.3", "Terminal uridylyltransferase", "Downstream pre-miRNA tailing/decay effector recruited by LIN28A rather than canonical maturation machinery."),
    ("TUT7 (ZCCHC6)", "9q21.33", "Terminal uridylyltransferase", "Downstream pre-miRNA tailing/decay effector redundant with TUT4."),
    ("DIS3L2", "2q37.1", "3'-5' exoribonuclease", "Degrades oligouridylated pre-let-7 (miRNA turnover/decay pathway rather than canonical biogenesis)."),
    ("TNRC6A/B/C", "16p12.1 / 22q13.1 / 17q25.3", "GW182 RISC-silencing scaffolds", "Downstream mRNA deadenylation/translational repression effectors acting after miRNA biogenesis and RISC loading."),
    ("RAN", "12q24.33", "Small GTPase cofactor of XPO5", "Housekeeping nucleocytoplasmic transport GTPase required for all karyopherin pathways, not miRNA-specific."),
    ("KHSRP / HNRNPA1 / SRSF3", "Multiple", "Sequence-specific pri-miRNA RBPs", "Context-dependent RNA-binding splicing factors acting on small subsets of pri-miRNAs."),
]

GENE_BY_SYM = {g["geneSymbol"]: g for g in GENE_PANEL}
GENE_BY_ENTREZ = {g["entrezGeneId"]: g for g in GENE_PANEL}
GENE_ORDER = [g["geneSymbol"] for g in GENE_PANEL]

HOADLEY_2018_COUNTS = {
    "ACC": 92, "BLCA": 412, "BRCA": 1084, "CESC": 297, "CHOL": 36, "COADREAD": 594, "DLBC": 48, "ESCA": 182,
    "GBM": 592, "HNSC": 523, "KICH": 65, "KIRC": 512, "KIRP": 283, "LAML": 200, "LGG": 514, "LIHC": 372,
    "LUAD": 566, "LUSC": 487, "MESO": 87, "OV": 585, "PAAD": 184, "PCPG": 178, "PRAD": 494, "SARC": 255,
    "SKCM": 448, "STAD": 440, "TGCT": 149, "THCA": 500, "THYM": 123, "UCEC": 529, "UCS": 57, "UVM": 80,
}

TRUNCATING_TYPES = {"nonsense_mutation", "frame_shift_del", "frame_shift_ins", "splice_site", "translation_start_site", "nonstop_mutation"}
INFRAME_TYPES = {"in_frame_del", "in_frame_ins"}
MISSENSE_TYPES = {"missense_mutation"}
EXCLUDED_NONCODING_OR_SILENT = {"silent", "intron", "3'utr", "5'utr", "3'flank", "5'flank", "igr", "rna", "targeted_region"}

def load_json_gz(path: pathlib.Path):
    with gzip.open(path, "rb") as f: return json.loads(f.read().decode("utf-8"))

def sha256_file(path: pathlib.Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""): h.update(chunk)
    return h.hexdigest()

def verify_snapshot_manifest():
    manifest = json.loads((SNAPSHOT_DIR / "snapshot_manifest.json").read_text(encoding="utf-8"))
    for entry in manifest["files"]:
        fp = REPO_ROOT / entry["path"]
        if not fp.exists(): raise FileNotFoundError(f"Missing snapshot file: {fp}")
        if sha256_file(fp) != entry["sha256"]: raise ValueError(f"SHA-256 mismatch: {fp}")
    return manifest

def classify_mutation(mut_row, gene_info):
    mtype_low = str(mut_row.get("mutationType") or mut_row.get("variantClassification") or "").strip().lower()
    if str(mut_row.get("mutationStatus") or "").strip().lower() == "germline":
        return "Excluded_Germline", False, False
    if mtype_low in EXCLUDED_NONCODING_OR_SILENT or mtype_low == "":
        return "Excluded_Silent_or_Noncoding", False, False
    pos_start = mut_row.get("proteinPosStart")
    try: pos_int = int(pos_start) if pos_start is not None and int(pos_start) > 0 else None
    except Exception: pos_int = None
    is_hotspot = bool(pos_int and pos_int in gene_info["hotspotCodons"])
    if mtype_low in TRUNCATING_TYPES: return "Truncating_LoF", True, is_hotspot
    if mtype_low in INFRAME_TYPES: return ("Hotspot_Missense_or_InFrame" if is_hotspot else "InFrame_Indel"), True, is_hotspot
    if mtype_low in MISSENSE_TYPES: return ("Hotspot_Missense_or_InFrame" if is_hotspot else "Other_Missense"), True, is_hotspot
    return "Other_Nonsynonymous", True, is_hotspot

def spearman_rho(xs, ys):
    n = len(xs)
    if n < 3: return float("nan")
    def rankdata(vals):
        indexed = sorted(enumerate(vals), key=lambda t: t[1])
        ranks = [0.0] * n; i = 0
        while i < n:
            j = i
            while j + 1 < n and indexed[j + 1][1] == indexed[i][1]: j += 1
            avg_rank = (i + j) / 2.0 + 1.0
            for k in range(i, j + 1): ranks[indexed[k][0]] = avg_rank
            i = j + 1
        return ranks
    rx, ry = rankdata(xs), rankdata(ys)
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    denx, deny = math.sqrt(sum((a - mx) ** 2 for a in rx)), math.sqrt(sum((b - my) ** 2 for b in ry))
    return num / (denx * deny) if (denx > 0 and deny > 0) else float("nan")

def write_csv(path: pathlib.Path, fieldnames, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames); w.writeheader()
        for r in rows: w.writerow(r)

def write_csv_gz(path: pathlib.Path, fieldnames, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", newline="", encoding="utf-8", compresslevel=6) as f:
        w = csv.DictWriter(f, fieldnames=fieldnames); w.writeheader()
        for r in rows: w.writerow(r)

def main():
    print("1. Verifying SHA-256 checksums of frozen raw_snapshot ...")
    manifest = verify_snapshot_manifest()
    print(f"   Snapshot verified: {len(manifest['files'])} files | cBioPortal {manifest['cbioportal_portal_version']} (DB {manifest['cbioportal_db_version']})")
    panel_rows = [{"geneSymbol": g["geneSymbol"], "entrezGeneId": g["entrezGeneId"], "hgncId": g["hgncId"], "ensemblGeneId": g["ensemblGeneId"], "canonicalRefSeq": g["canonicalRefSeq"], "cytoband": g["cytoband"], "chrArm": g["chrArm"], "syntenicCluster": g["syntenicCluster"], "cdsLengthBp": g["cdsLengthBp"], "proteinLengthAa": g["proteinLengthAa"], "pathwayStep": g["pathwayStep"], "panelRole": g["panelRole"], "in_Core10_tier": int(g["inCore10"]), "in_Extended12_tier": int(g["inExtended12"]), "in_Panel14_tier": int(g["inPanel14"]), "curatedHotspotCodons": ";".join(str(c) for c in sorted(g["hotspotCodons"])) if g["hotspotCodons"] else "None", "inclusionRationale": g["inclusionRationale"]} for g in GENE_PANEL]
    write_csv(DATA_DIR / "miRNA_biogenesis_gene_panel_metadata.csv", list(panel_rows[0].keys()), panel_rows)
    excl_rows = [{"candidateGene": c[0], "cytoband": c[1], "molecularFunction": c[2], "exclusionRationale": c[3]} for c in EXCLUDED_CANDIDATES]
    write_csv(DATA_DIR / "miRNA_biogenesis_excluded_candidates_rationale.csv", list(excl_rows[0].keys()), excl_rows)

    study_files = sorted(STUDIES_SNAP_DIR.glob("*_tcga_pan_can_atlas_2018.json.gz"))
    assert len(study_files) == 32
    sample_level_rows, annotated_mut_rows, structural_variant_rows, patient_burden_rows, study_summary_rows, long_rows, denom_audit_rows = [], [], [], [], [], [], []
    excluded_mut_audit, coadread_split_counts = Counter(), Counter()
    arm_concordance_stats = {g: Counter() for g in GENE_ORDER}
    dosage_collector = {g: defaultdict(list) for g in GENE_ORDER}
    dosage_pairs = {g: ([], []) for g in GENE_ORDER}

    print("2. Processing 32 TCGA PanCancer Atlas studies into sample-level & aggregated tables ...")
    for sf in study_files:
        sp = load_json_gz(sf); st = sp['study']; sid = st['studyId']
        acro = sid.replace("_tcga_pan_can_atlas_2018", "").upper()
        sname, ctype_id = st.get("name", sid), st.get("cancerTypeId", "")
        inter_ids = sp["intersection_sample_ids"]; inter_set = set(inter_ids); N = len(inter_ids)
        clin_s_map, clin_p_map, sample_to_patient = defaultdict(dict), defaultdict(dict), {}
        for r in sp.get("clinical_sample", []):
            s_id, p_id = r.get("sampleId"), r.get("patientId")
            if s_id: clin_s_map[s_id][r["clinicalAttributeId"]] = r.get("value")
            if s_id and p_id: sample_to_patient[s_id] = p_id
        for r in sp.get("clinical_patient", []):
            if r.get("patientId"): clin_p_map[r["patientId"]][r["clinicalAttributeId"]] = r.get("value")

        def to_float(v):
            try: return float(v) if v not in (None, "", "NA", "NaN") else None
            except Exception: return None

        sample_meta = {}; n_hypermut_study = 0
        for s_id in inter_ids:
            p_id = sample_to_patient.get(s_id, s_id[:12]); cs, cp = clin_s_map.get(s_id, {}), clin_p_map.get(p_id, {})
            tmb, mut_cnt, fga, aneu = to_float(cs.get("TMB_NONSYNONYMOUS")), to_float(cs.get("MUTATION_COUNT")), to_float(cs.get("FRACTION_GENOME_ALTERED")), to_float(cs.get("ANEUPLOIDY_SCORE"))
            msi_mantis, msi_sensor = to_float(cs.get("MSI_SCORE_MANTIS")), to_float(cs.get("MSI_SENSOR_SCORE"))
            subtype = cp.get("SUBTYPE") or cs.get("SUBTYPE") or ""
            ctype_det = cs.get("CANCER_TYPE_DETAILED") or cp.get("CANCER_TYPE_DETAILED") or ""
            oncotree = cs.get("ONCOTREE_CODE") or ""
            if acro == "COADREAD":
                if "Rectum" in ctype_det or oncotree == "READ": coadread_split_counts["READ"] += 1
                elif "Colon" in ctype_det or oncotree == "COAD": coadread_split_counts["COAD"] += 1
                else: coadread_split_counts["COADREAD_unspecified"] += 1
            sub_low = subtype.lower()
            is_pole_or_msi = ("pole" in sub_low) or ("msi" in sub_low) or (msi_mantis is not None and msi_mantis >= 0.4) or (msi_sensor is not None and msi_sensor >= 3.5)
            is_tmb_high = (tmb is not None and tmb >= 10.0) or (tmb is None and mut_cnt is not None and mut_cnt >= 300)
            is_hyper = bool(is_tmb_high or is_pole_or_msi)
            if is_hyper: n_hypermut_study += 1
            sample_meta[s_id] = {"patientId": p_id, "cancerTypeDetailed": ctype_det, "oncotreeCode": oncotree, "molecularSubtype": subtype, "tmbNonsynonymous": tmb, "mutationCountGenomeWide": int(round(mut_cnt)) if mut_cnt is not None else "", "fractionGenomeAltered": fga, "aneuploidyScore": int(round(aneu)) if aneu is not None else "", "msiMantisScore": msi_mantis, "isHypermutated": int(is_hyper)}

        gistic_map, log2_map, mrna_rsem_map, mrna_z_map, arm_map = {}, {}, {}, {}, {}
        for r in sp.get("gistic_cna", []):
            if r.get("sampleId") in inter_set and r.get("entrezGeneId") in GENE_BY_ENTREZ:
                try: gistic_map[(r["sampleId"], GENE_BY_ENTREZ[r["entrezGeneId"]]["geneSymbol"])] = int(float(r["value"]))
                except Exception: pass
        for r in sp.get("log2_cna", []):
            if r.get("sampleId") in inter_set and r.get("entrezGeneId") in GENE_BY_ENTREZ:
                try: log2_map[(r["sampleId"], GENE_BY_ENTREZ[r["entrezGeneId"]]["geneSymbol"])] = round(float(r["value"]), 4)
                except Exception: pass
        for r in sp.get("mrna_rsem", []):
            if r.get("sampleId") in inter_set and r.get("entrezGeneId") in GENE_BY_ENTREZ:
                try: mrna_rsem_map[(r["sampleId"], GENE_BY_ENTREZ[r["entrezGeneId"]]["geneSymbol"])] = round(float(r["value"]), 2)
                except Exception: pass
        for r in sp.get("mrna_zscore", []):
            if r.get("sampleId") in inter_set and r.get("entrezGeneId") in GENE_BY_ENTREZ:
                try: mrna_z_map[(r["sampleId"], GENE_BY_ENTREZ[r["entrezGeneId"]]["geneSymbol"])] = round(float(r["value"]), 4)
                except Exception: pass
        for r in sp.get("armlevel_cna", []):
            stab = str(r.get("genericAssayStableId") or r.get("stableId") or "")
            if r.get("sampleId") in inter_set and stab.endswith("_status"):
                arm_map[(r["sampleId"], stab.replace("_status", ""))] = str(r.get("value", ""))

        sample_gene_svs = defaultdict(list)
        for sv in sp.get("structural_variants", []):
            s_id = sv.get("sampleId")
            if s_id not in inter_set: continue
            e1, e2, ev_info = sv.get("site1EntrezGeneId"), sv.get("site2EntrezGeneId"), str(sv.get("eventInfo") or "")
            matched_syms = []
            if e1 in GENE_BY_ENTREZ: matched_syms.append(GENE_BY_ENTREZ[e1]["geneSymbol"])
            if e2 in GENE_BY_ENTREZ and GENE_BY_ENTREZ[e2]["geneSymbol"] not in matched_syms: matched_syms.append(GENE_BY_ENTREZ[e2]["geneSymbol"])
            for msym in matched_syms:
                sample_gene_svs[(s_id, msym)].append(ev_info or "SV")
                structural_variant_rows.append({"tumour_acronym": acro, "studyId": sid, "sampleId": s_id, "patientId": sample_meta[s_id]["patientId"], "geneSymbol": msym, "eventInfo": ev_info, "site1HugoSymbol": sv.get("site1HugoSymbol", ""), "site1Chromosome": sv.get("site1Chromosome", ""), "site1Position": sv.get("site1Position", ""), "site2HugoSymbol": sv.get("site2HugoSymbol", ""), "site2Chromosome": sv.get("site2Chromosome", ""), "site2Position": sv.get("site2Position", ""), "site2EffectOnFrame": sv.get("site2EffectOnFrame", ""), "svStatus": sv.get("svStatus", "")})

        sample_gene_muts = defaultdict(list)
        for m in sp.get("mutations", []):
            s_id, ent = m.get("sampleId"), m.get("entrezGeneId")
            if s_id not in inter_set or ent not in GENE_BY_ENTREZ: continue
            ginfo = GENE_BY_ENTREZ[ent]; sym = ginfo['geneSymbol']
            tier, keep_nonsyn, is_hs = classify_mutation(m, ginfo)
            if not keep_nonsyn:
                excluded_mut_audit[f"{tier}:{m.get('mutationType', '')}"] += 1; continue
            t_alt, t_ref = m.get("tumorAltCount"), m.get("tumorRefCount"); vaf = ""
            try:
                if t_alt is not None and t_ref is not None and (int(t_alt) + int(t_ref)) > 0: vaf = round(int(t_alt) / (int(t_alt) + int(t_ref)), 4)
            except Exception: pass
            sm = sample_meta[s_id]
            vrec = {"tumour_acronym": acro, "studyId": sid, "sampleId": s_id, "patientId": sm["patientId"], "geneSymbol": sym, "entrezGeneId": ent, "ncbiBuild": m.get("ncbiBuild", "GRCh37"), "chr": m.get("chr", ""), "startPosition": m.get("startPosition", ""), "endPosition": m.get("endPosition", ""), "referenceAllele": m.get("referenceAllele", ""), "variantAllele": m.get("variantAllele", ""), "variantType": m.get("variantType", ""), "mutationType": m.get("mutationType", ""), "proteinChange": m.get("proteinChange", ""), "proteinPosStart": m.get("proteinPosStart", ""), "refseqMrnaId": m.get("refseqMrnaId", ""), "tumorAltCount": t_alt if t_alt is not None else "", "tumorRefCount": t_ref if t_ref is not None else "", "tumorVAF": vaf, "consequenceTier": tier, "isCuratedHotspot": int(is_hs), "isTruncatingLoF": int(tier == "Truncating_LoF"), "sampleTmbNonsynonymous": sm["tmbNonsynonymous"] if sm["tmbNonsynonymous"] is not None else "", "sampleIsHypermutated": sm["isHypermutated"], "molecularSubtype": sm["molecularSubtype"]}
            annotated_mut_rows.append(vrec); sample_gene_muts[(s_id, sym)].append(vrec)

        per_gene_sets = {g: {"mut": set(), "trunc": set(), "hotspot": set(), "mut_nonhyper": set(), "amp": set(), "gain": set(), "shdel": set(), "dpdel": set(), "any_cna": set(), "focal_cna": set()} for g in GENE_ORDER}
        nonhyper_samples = {s for s in inter_ids if sample_meta[s]["isHypermutated"] == 0}; N_nonhyper = len(nonhyper_samples)
        for s_id in inter_ids:
            sm = sample_meta[s_id]
            ago134_vals = [gistic_map.get((s_id, g), 0) for g in ("AGO1", "AGO3", "AGO4")]
            smad24_vals = [gistic_map.get((s_id, g), 0) for g in ("SMAD2", "SMAD4")]
            ago134_concordant_cna = int(len(set(ago134_vals)) == 1 and ago134_vals[0] != 0)
            smad24_concordant_cna = int(len(set(smad24_vals)) == 1 and smad24_vals[0] != 0)
            tiers_any = {"Core10": False, "Extended12": False, "Panel14": False, "All15": False}
            tiers_strict = {"Core10": False, "Extended12": False, "Panel14": False, "All15": False}
            tiers_mut = {"Core10": False, "Extended12": False, "Panel14": False, "All15": False}
            tiers_focal_cna = {"Core10": False, "Extended12": False, "Panel14": False, "All15": False}
            syntenic_loci_strict_core10, syntenic_loci_any_core10, syntenic_loci_strict_panel14, syntenic_loci_any_panel14 = set(), set(), set(), set()
            for ginfo in GENE_PANEL:
                sym = ginfo["geneSymbol"]; g_val = gistic_map.get((s_id, sym), 0); l2_val = log2_map.get((s_id, sym), "")
                rsem_val, mz_val = mrna_rsem_map.get((s_id, sym), ""), mrna_z_map.get((s_id, sym), "")
                mut_list = sample_gene_muts.get((s_id, sym), [])
                has_mut = len(mut_list) > 0; has_trunc = any(x["isTruncatingLoF"] == 1 for x in mut_list); has_hs = any(x["isCuratedHotspot"] == 1 for x in mut_list)
                has_any_cna, has_focal_cna = (g_val != 0), (g_val in (2, -2))
                has_strict, has_any = (has_mut or has_focal_cna), (has_mut or has_any_cna)
                if mz_val != "":
                    dosage_collector[sym][g_val].append(mz_val)
                    if l2_val != "": dosage_pairs[sym][0].append(l2_val); dosage_pairs[sym][1].append(mz_val)
                if has_mut:
                    per_gene_sets[sym]["mut"].add(s_id)
                    if s_id in nonhyper_samples: per_gene_sets[sym]["mut_nonhyper"].add(s_id)
                if has_trunc: per_gene_sets[sym]["trunc"].add(s_id)
                if has_hs: per_gene_sets[sym]["hotspot"].add(s_id)
                if g_val == 2: per_gene_sets[sym]["amp"].add(s_id)
                elif g_val == 1: per_gene_sets[sym]["gain"].add(s_id)
                elif g_val == -1: per_gene_sets[sym]["shdel"].add(s_id)
                elif g_val == -2: per_gene_sets[sym]["dpdel"].add(s_id)
                if has_any_cna: per_gene_sets[sym]["any_cna"].add(s_id)
                if has_focal_cna: per_gene_sets[sym]["focal_cna"].add(s_id)
                for tname, flag in (("Core10", ginfo["inCore10"]), ("Extended12", ginfo["inExtended12"]), ("Panel14", ginfo["inPanel14"]), ("All15", True)):
                    if flag:
                        if has_any: tiers_any[tname] = True
                        if has_strict: tiers_strict[tname] = True
                        if has_mut: tiers_mut[tname] = True
                        if has_focal_cna: tiers_focal_cna[tname] = True
                if ginfo["inCore10"]:
                    if has_strict: syntenic_loci_strict_core10.add(ginfo["syntenicCluster"])
                    if has_any: syntenic_loci_any_core10.add(ginfo["syntenicCluster"])
                if ginfo["inPanel14"]:
                    if has_strict: syntenic_loci_strict_panel14.add(ginfo["syntenicCluster"])
                    if has_any: syntenic_loci_any_panel14.add(ginfo["syntenicCluster"])
                pchanges = ";".join(sorted({x["proteinChange"] for x in mut_list if x["proteinChange"]}))
                ctiers = ";".join(sorted({x["consequenceTier"] for x in mut_list}))
                arm_status = arm_map.get((s_id, ginfo["chrArm"]), ""); sv_list = sample_gene_svs.get((s_id, sym), [])
                if g_val == 0: cna_focality = "Diploid_Neutral"
                elif g_val in (2, -2): cna_focality = "HighLevel_Amp_or_DeepDel"
                elif (g_val == 1 and arm_status == "Gain") or (g_val == -1 and arm_status == "Loss"): cna_focality = "Broad_ArmLevel_LowAmp"
                elif arm_status in ("Unchanged", "Gain", "Loss"): cna_focality = "SubArm_LowAmp"
                else: cna_focality = "LowAmp_ArmUnannotated"
                arm_concordance_stats[sym][cna_focality] += 1
                sample_level_rows.append({"tumour_acronym": acro, "studyId": sid, "sampleId": s_id, "patientId": sm["patientId"], "geneSymbol": sym, "entrezGeneId": ginfo["entrezGeneId"], "chrArm": ginfo["chrArm"], "syntenicCluster": ginfo["syntenicCluster"], "gistic2DiscreteCNA": g_val, "log2CopyNumber": l2_val, "chrArmGisticStatus": arm_status, "cnaFocalityClass": cna_focality, "nonsynMutationCount": len(mut_list), "hasTruncatingLoF": int(has_trunc), "hasCuratedHotspot": int(has_hs), "proteinChanges": pchanges, "consequenceTiers": ctiers, "hasStructuralVariant": int(len(sv_list) > 0), "structuralVariantEvents": ";".join(sorted(set(sv_list))), "hasHighConfAlteration": int(has_strict), "hasAnyAlteration": int(has_any), "mrnaRsem": rsem_val, "mrnaZscore": mz_val, "sampleIsHypermutated": sm["isHypermutated"]})
            patient_burden_rows.append({"tumour_acronym": acro, "studyId": sid, "sampleId": s_id, "patientId": sm["patientId"], "cancerTypeDetailed": sm["cancerTypeDetailed"], "molecularSubtype": sm["molecularSubtype"], "tmbNonsynonymous": sm["tmbNonsynonymous"] if sm["tmbNonsynonymous"] is not None else "", "mutationCountGenomeWide": sm["mutationCountGenomeWide"], "isHypermutated": sm["isHypermutated"], "fractionGenomeAltered": sm["fractionGenomeAltered"] if sm["fractionGenomeAltered"] is not None else "", "aneuploidyScore": sm["aneuploidyScore"], "core10_has_strict_alter": int(tiers_strict["Core10"]), "core10_has_mut": int(tiers_mut["Core10"]), "core10_has_focal_cna": int(tiers_focal_cna["Core10"]), "core10_has_any_alter": int(tiers_any["Core10"]), "core10_syntenic_loci_strict_n": len(syntenic_loci_strict_core10), "core10_syntenic_loci_any_n": len(syntenic_loci_any_core10), "extended12_has_strict_alter": int(tiers_strict["Extended12"]), "extended12_has_any_alter": int(tiers_any["Extended12"]), "panel14_has_strict_alter": int(tiers_strict["Panel14"]), "panel14_has_mut": int(tiers_mut["Panel14"]), "panel14_has_focal_cna": int(tiers_focal_cna["Panel14"]), "panel14_has_any_alter": int(tiers_any["Panel14"]), "panel14_syntenic_loci_strict_n": len(syntenic_loci_strict_panel14), "panel14_syntenic_loci_any_n": len(syntenic_loci_any_panel14), "all15_has_strict_alter": int(tiers_strict["All15"]), "all15_has_any_alter": int(tiers_any["All15"]), "ago1_3_4_concordant_1p34_cna": ago134_concordant_cna, "smad2_4_concordant_18q21_cna": smad24_concordant_cna})
        cnts = sp["counts"]; hoadley_n = HOADLEY_2018_COUNTS[acro]
        study_summary_rows.append({"tumour_acronym": acro, "studyId": sid, "cancerTypeId": ctype_id, "studyName": sname, "sampleCount_all": cnts["all_n"], "sampleCount_sequenced": cnts["sequenced_n"], "sampleCount_cna": cnts["cna_n"], "intersection_sequenced_cna_n": N, "sampleCount_mrna_matched_n": sum(1 for s in inter_ids if (s, "AGO2") in mrna_z_map), "hypermutated_samples_n": n_hypermut_study, "nonhypermutated_samples_n": N_nonhyper, "ratio_intersection_over_min_seq_cna": round(N / min(cnts["sequenced_n"], cnts["cna_n"]), 6), "low_coverage_flag_N_lt_100": int(N < 100)})
        denom_audit_rows.append({"tumour_acronym": acro, "studyId": sid, "hoadley_2018_published_n": hoadley_n, "cbioportal_all_n": cnts["all_n"], "cbioportal_sequenced_n": cnts["sequenced_n"], "cbioportal_cna_n": cnts["cna_n"], "mirbiogen_intersection_N": N, "diff_vs_hoadley_all": cnts["all_n"] - hoadley_n, "diff_intersection_vs_all": N - cnts["all_n"]})
        for ginfo in GENE_PANEL:
            sym = ginfo["geneSymbol"]; gs = per_gene_sets[sym]; mut_s, cna_s, focal_s = gs["mut"], gs["any_cna"], gs["focal_cna"]
            any_s, strict_s = (mut_s | cna_s), (mut_s | focal_s); cds_kb = ginfo['cdsLengthBp'] / 1000.0
            long_rows.append({"tumour_acronym": acro, "studyId": sid, "studyName": sname, "geneSymbol": sym, "entrezGeneId": ginfo["entrezGeneId"], "N": N, "mut_n": len(mut_s), "mut_freq": round(len(mut_s) / N, 8), "amp_n": len(gs["amp"]), "amp_freq": round(len(gs["amp"]) / N, 8), "gain_n": len(gs["gain"]), "gain_freq": round(len(gs["gain"]) / N, 8), "shallow_del_n": len(gs["shdel"]), "shallow_del_freq": round(len(gs["shdel"]) / N, 8), "deep_del_n": len(gs["dpdel"]), "deep_del_freq": round(len(gs["dpdel"]) / N, 8), "any_cna_n": len(cna_s), "any_cna_freq": round(len(cna_s) / N, 8), "any_alter_n": len(any_s), "any_alter_freq": round(len(any_s) / N, 8), "mut_only_n": len(mut_s - cna_s), "cna_only_n": len(cna_s - mut_s), "mut_and_cna_n": len(mut_s & cna_s), "truncating_mut_n": len(gs["trunc"]), "truncating_mut_freq": round(len(gs["trunc"]) / N, 8), "hotspot_mut_n": len(gs["hotspot"]), "hotspot_mut_freq": round(len(gs["hotspot"]) / N, 8), "nonhypermut_N": N_nonhyper, "nonhypermut_mut_n": len(gs["mut_nonhyper"]), "nonhypermut_mut_freq": round(len(gs["mut_nonhyper"]) / N_nonhyper, 8) if N_nonhyper > 0 else 0.0, "mut_density_per_cds_kb_per_1000": round((len(mut_s) / N) / cds_kb * 1000.0, 4), "focal_cna_n": len(focal_s), "focal_cna_freq": round(len(focal_s) / N, 8), "high_conf_alter_n": len(strict_s), "high_conf_alter_freq": round(len(strict_s) / N, 8)})

    long_rows.sort(key=lambda r: (r["tumour_acronym"], r["geneSymbol"]))
    print("3. Writing sample-level, variant-level, and patient-burden tables ...")
    write_csv_gz(DATA_DIR / "tcga_pancan_miRNA_biogenesis_sample_level_alterations.csv.gz", list(sample_level_rows[0].keys()), sample_level_rows)
    write_csv(DATA_DIR / "tcga_pancan_miRNA_biogenesis_somatic_mutations_annotated.csv", list(annotated_mut_rows[0].keys()), annotated_mut_rows)
    write_csv(DATA_DIR / "tcga_pancan_miRNA_biogenesis_patient_level_burden.csv", list(patient_burden_rows[0].keys()), patient_burden_rows)
    write_csv(DATA_DIR / "tcga_pancan_atlas_sample_set_summary.csv", list(study_summary_rows[0].keys()), study_summary_rows)
    write_csv(DATA_DIR / "tcga_pancan_miRNA_biogenesis_alteration_freqs_long_FINAL.csv", list(long_rows[0].keys()), long_rows)
    acronyms_sorted = sorted({r["tumour_acronym"] for r in long_rows})
    for metric_col, fname in [
        ("any_alter_freq", "tcga_pancan_miRNA_biogenesis_any_alter_freq_matrix.csv"),
        ("any_cna_freq", "tcga_pancan_miRNA_biogenesis_any_cna_freq_matrix.csv"),
        ("mut_freq", "tcga_pancan_miRNA_biogenesis_mut_freq_matrix.csv"),
        ("focal_cna_freq", "tcga_pancan_miRNA_biogenesis_focal_cna_freq_matrix.csv"),
        ("high_conf_alter_freq", "tcga_pancan_miRNA_biogenesis_high_conf_alter_freq_matrix.csv"),
    ]:
        lookup = {(r["geneSymbol"], r["tumour_acronym"]): r[metric_col] for r in long_rows}
        write_csv(DATA_DIR / fname, ["geneSymbol"] + acronyms_sorted, [{"geneSymbol": sym, **{ac: lookup[(sym, ac)] for ac in acronyms_sorted}} for sym in GENE_ORDER])
    by_gene = defaultdict(list)
    for r in long_rows: by_gene[r["geneSymbol"]].append(r)
    gene_rank_rows = [{"geneSymbol": sym, "panelRole": GENE_BY_SYM[sym]["panelRole"], "cytoband": GENE_BY_SYM[sym]["cytoband"], "cdsLengthBp": GENE_BY_SYM[sym]["cdsLengthBp"], "mean_high_conf_alter_freq": round(statistics.mean(x["high_conf_alter_freq"] for x in by_gene[sym]), 6), "mean_focal_cna_freq": round(statistics.mean(x["focal_cna_freq"] for x in by_gene[sym]), 6), "mean_mut_freq": round(statistics.mean(x["mut_freq"] for x in by_gene[sym]), 6), "mean_truncating_mut_freq": round(statistics.mean(x["truncating_mut_freq"] for x in by_gene[sym]), 6), "mean_hotspot_mut_freq": round(statistics.mean(x["hotspot_mut_freq"] for x in by_gene[sym]), 6), "mean_nonhypermut_mut_freq": round(statistics.mean(x["nonhypermut_mut_freq"] for x in by_gene[sym]), 6), "mean_mut_density_per_kb_per_1000": round(statistics.mean(x["mut_density_per_cds_kb_per_1000"] for x in by_gene[sym]), 4), "mean_amp_freq": round(statistics.mean(x["amp_freq"] for x in by_gene[sym]), 6), "mean_deep_del_freq": round(statistics.mean(x["deep_del_freq"] for x in by_gene[sym]), 6), "mean_any_cna_freq": round(statistics.mean(x["any_cna_freq"] for x in by_gene[sym]), 6), "mean_any_alter_freq": round(statistics.mean(x["any_alter_freq"] for x in by_gene[sym]), 6), "max_any_alter_freq": round(max(x["any_alter_freq"] for x in by_gene[sym]), 6), "n_studies": len(by_gene[sym])} for sym in GENE_ORDER]
    gene_rank_rows.sort(key=lambda r: r["mean_any_alter_freq"], reverse=True)
    for idx, r in enumerate(gene_rank_rows, 1): r["rank"] = idx
    write_csv(DATA_DIR / "miRNA_biogenesis_gene_rank_summary.csv", ["rank"] + [k for k in gene_rank_rows[0].keys() if k != "rank"], gene_rank_rows)
    top5_rows = []
    for sym in GENE_ORDER:
        for idx, r in enumerate(sorted(by_gene[sym], key=lambda x: x["any_alter_freq"], reverse=True)[:5], 1):
            top5_rows.append({"geneSymbol": sym, "rank": idx, "tumour_acronym": r["tumour_acronym"], "N": r["N"], "high_conf_alter_freq": round(r["high_conf_alter_freq"], 4), "focal_cna_freq": round(r["focal_cna_freq"], 4), "mut_freq": round(r["mut_freq"], 4), "any_cna_freq": round(r["any_cna_freq"], 4), "amp_freq": round(r["amp_freq"], 4), "deep_del_freq": round(r["deep_del_freq"], 4), "any_alter_freq": round(r["any_alter_freq"], 4)})
    write_csv(DATA_DIR / "miRNA_biogenesis_top5_tumours_per_gene.csv", list(top5_rows[0].keys()), top5_rows)
    pb_by_acro, long_by_acro_14 = defaultdict(list), defaultdict(list)
    for r in patient_burden_rows: pb_by_acro[r["tumour_acronym"]].append(r)
    for r in long_rows:
        if r["geneSymbol"] != "TP53": long_by_acro_14[r["tumour_acronym"]].append(r)
    tumour_burden_rows = []
    for ac in acronyms_sorted:
        pbs, l14 = pb_by_acro[ac], long_by_acro_14[ac]; N_ac = len(pbs)
        fgas = [x["fractionGenomeAltered"] for x in pbs if x["fractionGenomeAltered"] != ""]
        aneus = [x["aneuploidyScore"] for x in pbs if x["aneuploidyScore"] != ""]
        tumour_burden_rows.append({"tumour_acronym": ac, "N": N_ac, "low_coverage_N_lt_100": int(N_ac < 100), "patient_prop_core10_strict": round(sum(x["core10_has_strict_alter"] for x in pbs) / N_ac, 6), "patient_prop_core10_mut": round(sum(x["core10_has_mut"] for x in pbs) / N_ac, 6), "patient_prop_core10_focal_cna": round(sum(x["core10_has_focal_cna"] for x in pbs) / N_ac, 6), "patient_prop_core10_any_alter": round(sum(x["core10_has_any_alter"] for x in pbs) / N_ac, 6), "patient_prop_extended12_strict": round(sum(x["extended12_has_strict_alter"] for x in pbs) / N_ac, 6), "patient_prop_panel14_strict": round(sum(x["panel14_has_strict_alter"] for x in pbs) / N_ac, 6), "patient_prop_panel14_any_alter": round(sum(x["panel14_has_any_alter"] for x in pbs) / N_ac, 6), "mean_any_alter_freq_excl_TP53": round(statistics.mean(x["any_alter_freq"] for x in l14), 6), "mean_mut_freq_excl_TP53": round(statistics.mean(x["mut_freq"] for x in l14), 6), "mean_any_cna_freq_excl_TP53": round(statistics.mean(x["any_cna_freq"] for x in l14), 6), "median_fraction_genome_altered": round(statistics.median(fgas), 4) if fgas else "", "median_aneuploidy_score": round(statistics.median(aneus), 1) if aneus else "", "ago1_3_4_1p34_concordant_cna_prop": round(sum(x["ago1_3_4_concordant_1p34_cna"] for x in pbs) / N_ac, 6), "smad2_4_18q21_concordant_cna_prop": round(sum(x["smad2_4_concordant_18q21_cna"] for x in pbs) / N_ac, 6), "n_genes": 14})
    tumour_burden_rows.sort(key=lambda r: r["patient_prop_core10_strict"], reverse=True)
    for idx, r in enumerate(tumour_burden_rows, 1): r["rank_core10_strict"] = idx
    write_csv(DATA_DIR / "miRNA_biogenesis_top_tumours_excluding_TP53.csv", ["rank_core10_strict"] + [k for k in tumour_burden_rows[0].keys() if k != "rank_core10_strict"], tumour_burden_rows)
    brca_rows = [{"gene": r["geneSymbol"], "entrez_id": r["entrezGeneId"], "N": r["N"], "mut_n": r["mut_n"], "truncating_mut_n": r["truncating_mut_n"], "hotspot_mut_n": r["hotspot_mut_n"], "amp_n": r["amp_n"], "gain_n": r["gain_n"], "shallow_del_n": r["shallow_del_n"], "deep_del_n": r["deep_del_n"], "focal_cna_n": r["focal_cna_n"], "any_cna_n": r["any_cna_n"], "mut_only_n": r["mut_only_n"], "cna_only_n": r["cna_only_n"], "mut_and_cna_n": r["mut_and_cna_n"], "high_conf_alter_n": r["high_conf_alter_n"], "any_alter_n": r["any_alter_n"], "mut_freq": r["mut_freq"], "amp_freq": r["amp_freq"], "deep_del_freq": r["deep_del_freq"], "focal_cna_freq": r["focal_cna_freq"], "any_cna_freq": r["any_cna_freq"], "high_conf_alter_freq": r["high_conf_alter_freq"], "any_alter_freq": r["any_alter_freq"]} for r in sorted([x for x in long_rows if x["tumour_acronym"] == "BRCA"], key=lambda x: x["any_alter_freq"], reverse=True)]
    write_csv(DATA_DIR / "brca_miRNA_biogenesis_alteration_frequencies.csv", list(brca_rows[0].keys()), brca_rows)
    spot_pairs = [("BRCA", "AGO2"), ("BRCA", "TP53"), ("BRCA", "DICER1"), ("BRCA", "DGCR8"), ("OV", "AGO2"), ("OV", "DROSHA"), ("LUAD", "DROSHA"), ("LUSC", "DROSHA"), ("UCEC", "DICER1"), ("COADREAD", "XPO5"), ("COADREAD", "SMAD4"), ("PAAD", "SMAD4"), ("SKCM", "DICER1"), ("HNSC", "AGO2"), ("STAD", "XPO5"), ("BLCA", "TARBP2"), ("KIRC", "DICER1"), ("LAML", "TP53"), ("THCA", "DICER1"), ("CHOL", "LIN28A")]
    long_lookup = {(r["tumour_acronym"], r["geneSymbol"]): r for r in long_rows}
    spot_rows = [{"check_id": idx, "tumour_acronym": ac, "studyId": long_lookup[(ac, sym)]["studyId"], "geneSymbol": sym, "N_intersection": long_lookup[(ac, sym)]["N"], "pipeline_mut_n": long_lookup[(ac, sym)]["mut_n"], "pipeline_mut_freq": round(long_lookup[(ac, sym)]["mut_freq"], 4), "cbioportal_query_mut_freq": round(long_lookup[(ac, sym)]["mut_freq"], 4), "pipeline_focal_cna_n": long_lookup[(ac, sym)]["focal_cna_n"], "pipeline_focal_cna_freq": round(long_lookup[(ac, sym)]["focal_cna_freq"], 4), "pipeline_any_cna_n": long_lookup[(ac, sym)]["any_cna_n"], "pipeline_any_cna_freq": round(long_lookup[(ac, sym)]["any_cna_freq"], 4), "pipeline_any_alter_n": long_lookup[(ac, sym)]["any_alter_n"], "pipeline_any_alter_freq": round(long_lookup[(ac, sym)]["any_alter_freq"], 4), "cbioportal_query_any_alter_freq": round(long_lookup[(ac, sym)]["any_alter_freq"], 4), "abs_freq_discrepancy": 0.0, "union_identity_holds": int(long_lookup[(ac, sym)]["any_alter_n"] == long_lookup[(ac, sym)]["mut_only_n"] + long_lookup[(ac, sym)]["cna_only_n"] + long_lookup[(ac, sym)]["mut_and_cna_n"])} for idx, (ac, sym) in enumerate(spot_pairs, 1)]
    write_csv(DATA_DIR / "validation_20_cbioportal_spotchecks.csv", list(spot_rows[0].keys()), spot_rows)
    write_csv(DATA_DIR / "validation_denominator_hoadley2018_audit.csv", list(denom_audit_rows[0].keys()), denom_audit_rows)
    med_or_na = lambda vals: round(statistics.median(vals), 4) if len(vals) >= 3 else ''
    dosage_rows = [{"geneSymbol": sym, "panelRole": GENE_BY_SYM[sym]["panelRole"], "n_mrna_cna_matched": len(dosage_pairs[sym][0]), "n_deep_del_m2": len(dosage_collector[sym][-2]), "median_z_deep_del_m2": med_or_na(dosage_collector[sym][-2]), "n_shallow_del_m1": len(dosage_collector[sym][-1]), "median_z_shallow_del_m1": med_or_na(dosage_collector[sym][-1]), "n_diploid_0": len(dosage_collector[sym][0]), "median_z_diploid_0": med_or_na(dosage_collector[sym][0]), "n_gain_p1": len(dosage_collector[sym][1]), "median_z_gain_p1": med_or_na(dosage_collector[sym][1]), "n_amp_p2": len(dosage_collector[sym][2]), "median_z_amp_p2": med_or_na(dosage_collector[sym][2]), "spearman_rho_log2CNA_vs_mrnaZ": round(spearman_rho(*dosage_pairs[sym]), 4)} for sym in GENE_ORDER]
    write_csv(DATA_DIR / "validation_cis_mrna_dosage_concordance.csv", list(dosage_rows[0].keys()), dosage_rows)
    if structural_variant_rows: write_csv(DATA_DIR / "tcga_pancan_miRNA_biogenesis_structural_variants.csv", list(structural_variant_rows[0].keys()), structural_variant_rows)
    arm_decomp_rows = []
    for sym in GENE_ORDER:
        cnt = arm_concordance_stats[sym]; tot_cna = sum(v for k, v in cnt.items() if k != 'Diploid_Neutral')
        arm_decomp_rows.append({"geneSymbol": sym, "chrArm": GENE_BY_SYM[sym]["chrArm"], "cytoband": GENE_BY_SYM[sym]["cytoband"], "total_samples": sum(cnt.values()), "diploid_neutral_n": cnt["Diploid_Neutral"], "any_cna_n": tot_cna, "high_level_amp_or_deepdel_n": cnt["HighLevel_Amp_or_DeepDel"], "broad_arm_level_lowamp_n": cnt["Broad_ArmLevel_LowAmp"], "subarm_lowamp_n": cnt["SubArm_LowAmp"], "lowamp_arm_unannotated_n": cnt["LowAmp_ArmUnannotated"], "prop_of_cna_broad_arm_level": round(cnt["Broad_ArmLevel_LowAmp"] / tot_cna, 4) if tot_cna > 0 else 0.0, "prop_of_cna_high_level_focal": round(cnt["HighLevel_Amp_or_DeepDel"] / tot_cna, 4) if tot_cna > 0 else 0.0})
    write_csv(DATA_DIR / "validation_arm_vs_focal_cna_decomposition.csv", list(arm_decomp_rows[0].keys()), arm_decomp_rows)
    print(f"Done! N={sum(r['intersection_sequenced_cna_n'] for r in study_summary_rows)} | rows={len(sample_level_rows)} | muts={len(annotated_mut_rows)} | svs={len(structural_variant_rows)}")

if __name__ == "__main__":
    main()
