#!/usr/bin/env python3
"""miRBioGenAtlas v2.0 — Frozen Multi-Modal cBioPortal Snapshot Retriever."""
import concurrent.futures, datetime, gzip, hashlib, json, pathlib, time, urllib.error, urllib.parse, urllib.request

BASE_URL = "https://www.cbioportal.org/api"
GENE_MAP = {
    "DROSHA": 29102, "DGCR8": 54487, "XPO5": 57510, "DICER1": 23405, "TARBP2": 6895,
    "PRKRA": 8575, "AGO1": 26523, "AGO2": 27161, "AGO3": 192669, "AGO4": 192670,
    "LIN28A": 79727, "LIN28B": 389421, "SMAD2": 4087, "SMAD4": 4089, "TP53": 7157,
}
ENTREZ_IDS = list(GENE_MAP.values())
KEEP_CLINICAL_ATTRS = {
    "MUTATION_COUNT", "TMB_NONSYNONYMOUS", "ANEUPLOIDY_SCORE", "FRACTION_GENOME_ALTERED",
    "MSI_SCORE_MANTIS", "MSI_SENSOR_SCORE", "BUFFA_HYPOXIA_SCORE", "SAMPLE_TYPE",
    "ONCOTREE_CODE", "CANCER_TYPE", "CANCER_TYPE_DETAILED", "SUBTYPE", "SEX", "AGE",
    "AJCC_PATHOLOGIC_TUMOR_STAGE", "OS_MONTHS", "OS_STATUS", "DSS_MONTHS", "DSS_STATUS", "PFS_MONTHS", "PFS_STATUS",
}
REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent
SNAPSHOT_DIR = REPO_ROOT / "data" / "raw_snapshot"
STUDIES_SNAP_DIR = SNAPSHOT_DIR / "studies"

def http_json(method, path_or_url, params=None, body=None, timeout=90, retries=4):
    url = path_or_url if path_or_url.startswith("http") else f"{BASE_URL}{path_or_url}"
    if params:
        url = f"{url}?{urllib.parse.urlencode(params)}"
    headers = {"Accept": "application/json", "User-Agent": "miRBioGenAtlas-Snapshot/2.0 (Python-urllib)"}
    data_bytes = None
    if body is not None:
        data_bytes = json.dumps(body).encode("utf-8")
        headers["Content-Type"] = "application/json"
    for attempt in range(retries):
        req = urllib.request.Request(url, data=data_bytes, headers=headers, method=method)
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                raw = resp.read().decode("utf-8")
                return json.loads(raw) if raw else []
        except urllib.error.HTTPError as e:
            if e.code in (429, 500, 502, 503, 504) and attempt < retries - 1:
                time.sleep(2 ** (attempt + 1)); continue
            raise
        except Exception:
            if attempt < retries - 1:
                time.sleep(2 ** (attempt + 1)); continue
            raise
    return []

def save_json_gz(path: pathlib.Path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    raw = json.dumps(obj, sort_keys=True, separators=(",", ":")).encode("utf-8")
    with gzip.open(tmp, "wb", compresslevel=6) as f: f.write(raw)
    tmp.replace(path)

def sha256_file(path: pathlib.Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""): h.update(chunk)
    return h.hexdigest()

def fetch_armlevel_cna(arm_profile_id, sample_ids):
    if not arm_profile_id: return []
    for method, path, params, body in [
        ("GET", f"/generic_assay_data/{arm_profile_id}", {"projection": "SUMMARY"}, None),
        ("POST", f"/generic_assay_data/{arm_profile_id}/fetch", {"projection": "SUMMARY"}, {"sampleIds": sample_ids}),
        ("GET", f"/molecular-profiles/{arm_profile_id}/generic-assay-data", {"projection": "SUMMARY"}, None),
        ("POST", "/generic-assay-data/fetch", {"projection": "SUMMARY"}, {"molecularProfileIds": [arm_profile_id]}),
    ]:
        try:
            res = http_json(method, path, params=params, body=body, timeout=60, retries=2)
            if isinstance(res, list) and len(res) > 0:
                sset = set(sample_ids)
                return [r for r in res if r.get("sampleId") in sset]
        except Exception: continue
    return []

def fetch_structural_variants(sv_profile_id, sample_ids):
    if not sv_profile_id: return []
    sids = [{"molecularProfileId": sv_profile_id, "sampleId": s} for s in sample_ids]
    sset = set(sample_ids)
    for proj, body in [
        ("DETAILED", {"entrezGeneIds": ENTREZ_IDS, "molecularProfileIds": [sv_profile_id]}),
        ("SUMMARY", {"entrezGeneIds": ENTREZ_IDS, "sampleMolecularIdentifiers": sids}),
    ]:
        try:
            res = http_json("POST", "/structural-variant/fetch", params={"projection": proj}, body=body, timeout=60, retries=2)
            if isinstance(res, list):
                return [r for r in res if r.get("sampleId") in sset]
        except Exception: continue
    return []

def fetch_one_study(study_row):
    sid = study_row["studyId"]
    out_file = STUDIES_SNAP_DIR / f"{sid}.json.gz"
    if out_file.exists() and out_file.stat().st_size > 500:
        return sid, "cached", out_file.stat().st_size
    mol_profiles = http_json("GET", f"/studies/{sid}/molecular-profiles", params={"projection": "SUMMARY"})
    mp_by_id = {m["molecularProfileId"]: m for m in mol_profiles}
    def pick_profile(suffix, alt_type=None, datatype=None):
        cand = f"{sid}_{suffix}"
        if cand in mp_by_id: return cand
        for m in mol_profiles:
            if alt_type and m.get("molecularAlterationType") != alt_type: continue
            if datatype and m.get("datatype") != datatype: continue
            return m["molecularProfileId"]
        return None
    mut_prof = pick_profile("mutations", "MUTATION_EXTENDED", "MAF")
    gistic_prof = pick_profile("gistic", "COPY_NUMBER_ALTERATION", "DISCRETE")
    log2_prof = pick_profile("log2CNA", "COPY_NUMBER_ALTERATION", "LOG2-VALUE")
    arm_prof = f"{sid}_armlevel_cna" if f"{sid}_armlevel_cna" in mp_by_id else None
    sv_prof = pick_profile("structural_variants", "STRUCTURAL_VARIANT", None)
    mrna_rsem_prof = f"{sid}_rna_seq_v2_mrna" if f"{sid}_rna_seq_v2_mrna" in mp_by_id else None
    mrna_z_prof = f"{sid}_rna_seq_v2_mrna_median_all_sample_Zscores" if f"{sid}_rna_seq_v2_mrna_median_all_sample_Zscores" in mp_by_id else (f"{sid}_rna_seq_v2_mrna_median_Zscores" if f"{sid}_rna_seq_v2_mrna_median_Zscores" in mp_by_id else None)
    sample_lists = http_json("GET", f"/studies/{sid}/sample-lists", params={"projection": "SUMMARY"})
    sl_ids_by_suffix = {}
    for sl in sample_lists:
        slid = sl["sampleListId"]
        for suf in ("_all", "_sequenced", "_cna", "_cnaseq", "_3way_complete", "_rna_seq_v2_mrna", "_sv"):
            if slid == f"{sid}{suf}": sl_ids_by_suffix[suf] = slid
    sample_sets = {}
    for suf, slid in sl_ids_by_suffix.items():
        sids = http_json("GET", f"/sample-lists/{slid}/sample-ids")
        sample_sets[suf] = sorted(set(str(x).strip() for x in sids if str(x).strip()))
    seq_set, cna_set = set(sample_sets.get("_sequenced", [])), set(sample_sets.get("_cna", []))
    intersection_ids = sorted(seq_set & cna_set)
    mutations = http_json("POST", f"/molecular-profiles/{mut_prof}/mutations/fetch", params={"projection": "DETAILED"}, body={"entrezGeneIds": ENTREZ_IDS, "sampleIds": intersection_ids}, timeout=120) if (mut_prof and intersection_ids) else []
    gistic_cna = http_json("POST", f"/molecular-profiles/{gistic_prof}/molecular-data/fetch", params={"projection": "SUMMARY"}, body={"entrezGeneIds": ENTREZ_IDS, "sampleIds": intersection_ids}, timeout=120) if (gistic_prof and intersection_ids) else []
    log2_cna = http_json("POST", f"/molecular-profiles/{log2_prof}/molecular-data/fetch", params={"projection": "SUMMARY"}, body={"entrezGeneIds": ENTREZ_IDS, "sampleIds": intersection_ids}, timeout=120) if (log2_prof and intersection_ids) else []
    mrna_rsem, mrna_zscore = [], []
    if mrna_rsem_prof and intersection_ids:
        try: mrna_rsem = http_json("POST", f"/molecular-profiles/{mrna_rsem_prof}/molecular-data/fetch", params={"projection": "SUMMARY"}, body={"entrezGeneIds": ENTREZ_IDS, "sampleIds": intersection_ids}, timeout=120)
        except Exception: mrna_rsem = []
    if mrna_z_prof and intersection_ids:
        try: mrna_zscore = http_json("POST", f"/molecular-profiles/{mrna_z_prof}/molecular-data/fetch", params={"projection": "SUMMARY"}, body={"entrezGeneIds": ENTREZ_IDS, "sampleIds": intersection_ids}, timeout=120)
        except Exception: mrna_zscore = []
    armlevel_cna = fetch_armlevel_cna(arm_prof, intersection_ids)
    structural_variants = fetch_structural_variants(sv_prof, intersection_ids)
    clinical_sample, clinical_patient = [], []
    try:
        raw_clin_s = http_json("GET", f"/studies/{sid}/clinical-data", params={"clinicalDataType": "SAMPLE", "projection": "SUMMARY", "pageSize": 100000}, timeout=90)
        inter_set = set(intersection_ids)
        clinical_sample = [{"sampleId": r.get("sampleId"), "patientId": r.get("patientId"), "clinicalAttributeId": r.get("clinicalAttributeId"), "value": r.get("value")} for r in raw_clin_s if r.get("sampleId") in inter_set and r.get("clinicalAttributeId") in KEEP_CLINICAL_ATTRS]
    except Exception: pass
    try:
        raw_clin_p = http_json("GET", f"/studies/{sid}/clinical-data", params={"clinicalDataType": "PATIENT", "projection": "SUMMARY", "pageSize": 100000}, timeout=90)
        patient_set = {s[:12] for s in intersection_ids} | {r.get("patientId") for r in clinical_sample if r.get("patientId")}
        clinical_patient = [{"patientId": r.get("patientId"), "clinicalAttributeId": r.get("clinicalAttributeId"), "value": r.get("value")} for r in raw_clin_p if r.get("patientId") in patient_set and r.get("clinicalAttributeId") in KEEP_CLINICAL_ATTRS]
    except Exception: pass
    payload = {
        "study": study_row,
        "profiles_used": {"mutations": mut_prof, "gistic": gistic_prof, "log2CNA": log2_prof, "armlevel_cna": arm_prof, "structural_variants": sv_prof, "mrna_rsem": mrna_rsem_prof, "mrna_zscore": mrna_z_prof},
        "sample_lists_metadata": sample_lists, "sample_sets": sample_sets, "intersection_sample_ids": intersection_ids,
        "counts": {"all_n": len(sample_sets.get("_all", [])), "sequenced_n": len(seq_set), "cna_n": len(cna_set), "intersection_n": len(intersection_ids), "cnaseq_n": len(sample_sets.get("_cnaseq", [])), "threeway_n": len(sample_sets.get("_3way_complete", [])), "mrna_n": len(sample_sets.get("_rna_seq_v2_mrna", [])), "mutations_raw_rows": len(mutations), "gistic_rows": len(gistic_cna), "log2_cna_rows": len(log2_cna), "mrna_rsem_rows": len(mrna_rsem), "mrna_zscore_rows": len(mrna_zscore), "armlevel_cna_rows": len(armlevel_cna), "sv_rows": len(structural_variants), "clinical_sample_rows": len(clinical_sample), "clinical_patient_rows": len(clinical_patient)},
        "mutations": mutations, "gistic_cna": gistic_cna, "log2_cna": log2_cna, "mrna_rsem": mrna_rsem, "mrna_zscore": mrna_zscore, "armlevel_cna": armlevel_cna, "structural_variants": structural_variants, "clinical_sample": clinical_sample, "clinical_patient": clinical_patient,
    }
    save_json_gz(out_file, payload)
    return sid, payload["counts"], out_file.stat().st_size

def main():
    SNAPSHOT_DIR.mkdir(parents=True, exist_ok=True); STUDIES_SNAP_DIR.mkdir(parents=True, exist_ok=True)
    info = http_json("GET", "/info")
    genes_meta = http_json("POST", "/genes/fetch", params={"geneIdType": "ENTREZ_GENE_ID", "projection": "DETAILED"}, body=ENTREZ_IDS)
    all_studies, page = [], 0
    while True:
        chunk = http_json("GET", "/studies", params={"projection": "DETAILED", "pageSize": 1000, "pageNumber": page, "sortBy": "studyId", "direction": "ASC"})
        if not chunk: break
        all_studies.extend(chunk)
        if len(chunk) < 1000: break
        page += 1
    pancan_studies = sorted([s for s in all_studies if s.get("studyId", "").endswith("_tcga_pan_can_atlas_2018")], key=lambda x: x["studyId"])
    save_json_gz(SNAPSHOT_DIR / "portal_and_genes_metadata.json.gz", {"retrieved_utc": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"), "base_url": BASE_URL, "api_info": info, "target_gene_map": GENE_MAP, "genes_detailed": genes_meta, "pancan_studies": pancan_studies})
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        for fut in concurrent.futures.as_completed({pool.submit(fetch_one_study, st): st["studyId"] for st in pancan_studies}):
            sid_out, _, nbytes = fut.result()
            print(f"   {sid_out}: {nbytes/1024:.1f} KB")
    manifest_files = [{"path": p.relative_to(REPO_ROOT).as_posix(), "size_bytes": p.stat().st_size, "sha256": sha256_file(p)} for p in sorted(SNAPSHOT_DIR.rglob("*.gz"))]
    manifest = {"dataset": "miRBioGenAtlas", "snapshot_version": "2.0.0", "retrieved_utc": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"), "cbioportal_portal_version": info.get("portalVersion"), "cbioportal_db_version": info.get("dbVersion"), "n_studies": len(pancan_studies), "n_genes": len(ENTREZ_IDS), "files": manifest_files}
    (SNAPSHOT_DIR / "snapshot_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")

if __name__ == "__main__":
    main()
