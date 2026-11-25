#!/usr/bin/env python3
"""Command-line lookup for variant and codon annotations across the 213,940-variant atlas."""

import argparse
import csv
import gzip
import json
import pathlib
import re
import sys

AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}
HGVS3_RE = re.compile(r"^(?:p\.)?([A-Za-z]{3})(\d+)([A-Za-z]{3})$")
HGVS1_RE = re.compile(r"^(?:p\.)?([A-Z])(\d+)([A-Z])$")
VCF_RE = re.compile(r"^(?:(hg19|grch37|hg38|grch38)[:\s-]+)?(?:chr)?([0-9XY]+)[:\s-]+(\d+)[:\s-]+([ACGT])(?:>|:|-|\s+)([ACGT])$", re.I)
CDNA_RE = re.compile(r"^c\.(\d+)([ACGT])>([ACGT])$", re.I)
DNA_COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
CODON_TO_AA = {
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L","CTT":"L","CTC":"L","CTA":"L","CTG":"L",
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M","GTT":"V","GTC":"V","GTA":"V","GTG":"V",
    "TCT":"S","TCC":"S","TCA":"S","TCG":"S","CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "ACT":"T","ACC":"T","ACA":"T","ACG":"T","GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*","CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K","GAT":"D","GAC":"D","GAA":"E","GAG":"E",
    "TGT":"C","TGC":"C","TGA":"*","TGG":"W","CGT":"R","CGC":"R","CGA":"R","CGG":"R",
    "AGT":"S","AGC":"S","AGA":"R","AGG":"R","GGT":"G","GGC":"G","GGA":"G","GGG":"G",
}

REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent
DATA_DIR = REPO_ROOT / "data"


def resolve_genomic_or_cdna(arg1: str, arg2: str = None):
    cds_file = DATA_DIR / "canonical_cds_genomic_15genes.json"
    if not cds_file.exists():
        return None
    cds_db = json.loads(cds_file.read_text(encoding="utf-8")).get("genes", {})

    if arg2 and CDNA_RE.match(arg2.strip()):
        gene = arg1.strip().upper()
        gd = cds_db.get(gene)
        if not gd:
            return None
        m = CDNA_RE.match(arg2.strip())
        c_idx = int(m.group(1)) - 1
        c_ref, c_alt = m.group(2).upper(), m.group(3).upper()
        seq = gd["cds_seq"]
        if 0 <= c_idx < len(seq) and seq[c_idx] == c_ref:
            cs = (c_idx // 3) * 3
            wt_c = seq[cs:cs+3]
            mut_c = wt_c[:c_idx % 3] + c_alt + wt_c[c_idx % 3 + 1:]
            return gene, f"{CODON_TO_AA[wt_c]}{c_idx // 3 + 1}{CODON_TO_AA[mut_c]}"

    token = f"{arg1} {arg2}".strip() if (arg2 and not HGVS1_RE.match(arg2) and not HGVS3_RE.match(arg2)) else arg1.strip()
    mv = VCF_RE.match(token)
    if mv:
        pref = "GRCh37" if (mv.group(1) and ("19" in mv.group(1) or "37" in mv.group(1))) else "GRCh38"
        chrom = mv.group(2).upper()
        pos = int(mv.group(3))
        ref, alt = mv.group(4).upper(), mv.group(5).upper()
        builds = ["GRCh37", "GRCh38"] if pref == "GRCh37" else ["GRCh38", "GRCh37"]
        for bld in builds:
            ek = "cds_exons_grch38" if bld == "GRCh38" else "cds_exons_grch37"
            sk = "strand" if bld == "GRCh38" else "strand_grch37"
            for gene, gd in cds_db.items():
                if str(gd["chr"]).upper() != chrom:
                    continue
                cum = 0
                strand = gd.get(sk, gd["strand"])
                for _, cs, ce in gd.get(ek, []):
                    if cs <= pos <= ce:
                        off = (pos - cs) if strand == 1 else (ce - pos)
                        c_idx = cum + off
                        c_ref = ref if strand == 1 else DNA_COMP[ref]
                        c_alt = alt if strand == 1 else DNA_COMP[alt]
                        seq = gd["cds_seq"]
                        if 0 <= c_idx < len(seq) and seq[c_idx] == c_ref:
                            cod_s = (c_idx // 3) * 3
                            wt_c = seq[cod_s:cod_s+3]
                            mut_c = wt_c[:c_idx % 3] + c_alt + wt_c[c_idx % 3 + 1:]
                            return gene, f"{CODON_TO_AA[wt_c]}{c_idx // 3 + 1}{CODON_TO_AA[mut_c]}"
                    cum += (ce - cs + 1)
    return None


def normalize_variant(raw: str) -> str:
    s = raw.strip()
    m1 = HGVS1_RE.match(s.upper())
    if m1:
        return f"{m1.group(1)}{m1.group(2)}{m1.group(3)}"
    m3 = HGVS3_RE.match(s)
    if m3:
        wt = AA3_TO_1.get(m3.group(1).upper())
        mt = AA3_TO_1.get(m3.group(3).upper())
        if wt and mt:
            return f"{wt}{m3.group(2)}{mt}"
    raise ValueError(f"Unrecognised missense notation '{raw}' (expected e.g. E1705K or p.Glu1705Lys)")


def classify_biophysical_regime(am_val, e2_val, if1_val) -> str:
    constrained = (am_val is not None and am_val >= 0.564) or (e2_val is not None and e2_val <= -5.0)
    if not constrained:
        return "Tolerated_Neutral"
    if if1_val is not None and if1_val > -1.5:
        return "Fold_Preserving_Catalytic_or_Interface"
    return "Thermodynamic_Fold_Destabilizing"


def load_gene_metadata(gene: str) -> dict:
    with open(DATA_DIR / "miRNA_biogenesis_gene_panel_metadata.csv", "r", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            if row["geneSymbol"].upper() == gene.upper():
                return row
    raise KeyError(f"Gene '{gene}' is not in the 15-gene miRBioGenAtlas panel.")


def load_clinvar_hits(gene: str, variant_or_codon) -> list:
    p = DATA_DIR / "clinvar_clinical_variants_15genes_annotated.csv"
    if not p.exists():
        return []
    hits = []
    with open(p, "r", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            if r["geneSymbol"] != gene:
                continue
            if isinstance(variant_or_codon, int):
                if int(r["codonPos"]) == variant_or_codon:
                    hits.append(r)
            elif r["proteinChange"] == variant_or_codon:
                hits.append(r)
    return hits


def load_msk_counts(gene: str) -> dict:
    p = DATA_DIR / "raw_snapshot" / "external" / "msk_impact_10945_replication.json.gz"
    if not p.exists():
        return {}
    with gzip.open(p, "rt", encoding="utf-8") as f:
        snap = json.load(f)
    counts = {}
    for m in snap.get("mutations", []):
        if m.get("g") == gene:
            pch = (m.get("p") or "").lstrip("p.")
            counts[pch] = counts.get(pch, 0) + 1
    return counts


def query_single_variant(gene: str, raw_var: str, as_json: bool = False):
    gmeta = load_gene_metadata(gene)
    sym = gmeta["geneSymbol"]
    var = normalize_variant(raw_var)

    target = None
    with gzip.open(DATA_DIR / "miRNA_biogenesis_saturation_mutagenesis_213940_atlas.csv.gz", "rt", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            if r["geneSymbol"] == sym and r["variant"] == var:
                target = r
                break

    if target is None:
        sys.exit(
            f"Variant {sym} {var} not found in 213,940-variant saturation atlas "
            f"(check canonical UniProt {gmeta['uniprotAccession']} wild-type residue)."
        )

    cv_hits = load_clinvar_hits(sym, var)
    msk_counts = load_msk_counts(sym)
    msk_n = msk_counts.get(var, 0) if sym in ("DICER1", "SMAD2", "SMAD4", "TP53") else "not_on_IM3_IM5_panel"

    am_f = float(target["alphamissenseScore"]) if target["alphamissenseScore"] != "" else None
    e2_f = float(target["esm2_3b_llr"])
    e1_f = float(target["esm1v_650m_llr"])
    if1_f = float(target["esm_if1_3d_llr"])
    regime = classify_biophysical_regime(am_f, e2_f, if1_f)

    payload = {
        "geneSymbol": sym,
        "uniprotAccession": gmeta["uniprotAccession"],
        "pathwayStep": gmeta["pathwayStep"],
        "variant": f"p.{var}",
        "codonPos": int(target["codonPos"]),
        "pdbComplexId": target["pdbComplexId"],
        "afdbPlddt": float(target["afdbPlddt"]) if target["afdbPlddt"] != "" else None,
        "distSubstrateA": float(target["distSubstrateA"]) if target["distSubstrateA"] != "" else None,
        "distMetalA": float(target["distMetalA"]) if target["distMetalA"] != "" else None,
        "distPartnerA": float(target["distPartnerA"]) if target["distPartnerA"] != "" else None,
        "alphamissenseScore": am_f,
        "alphamissenseClass": target["alphamissenseClass"],
        "esm2_3b_llr": e2_f,
        "esm1v_650m_llr": e1_f,
        "esm_if1_3d_llr": if1_f,
        "biophysicalRegime": regime,
        "tcgaPanCanSomaticCount_N10195": int(target["observedTcgaSomaticCount"]),
        "mskImpactSomaticCount_N10945": msk_n,
        "gnomad_pLI": float(gmeta["gnomad_pLI"]) if gmeta.get("gnomad_pLI", "") != "" else None,
        "gnomad_LOEUF": float(gmeta["gnomad_LOEUF"]) if gmeta.get("gnomad_LOEUF", "") != "" else None,
        "gnomad_mis_z": float(gmeta["gnomad_mis_z"]) if gmeta.get("gnomad_mis_z", "") != "" else None,
        "clinvarRecords": [
            {
                "clinvarVariationId": c["clinvarVariationId"],
                "clinicalSignificance": c["clinicalSignificance"],
                "reviewStars": int(c["reviewStars"]),
                "phenotypes": c["phenotypes"],
                "rsid": c["rsid"],
                "grch38Vcf": c["grch38Vcf"],
            }
            for c in cv_hits
        ],
    }

    if as_json:
        print(json.dumps(payload, indent=2))
        return

    print("=" * 78)
    print(f"miRBioGenAtlas — Structural & Clinical Variant Triage: {sym} p.{var}")
    print("=" * 78)
    print(f"  Gene / UniProt / Step : {sym} ({payload['uniprotAccession']}) | {payload['pathwayStep']}")
    print(f"  3D Structure & Conf   : Experimental PDB: {payload['pdbComplexId'] or 'AlphaFold v6'} | AFDB pLDDT={payload['afdbPlddt']}")
    print(f"  3D Atomic Distances   : RNA/DNA substrate={payload['distSubstrateA']} A | Catalytic Mg2+/Zn2+={payload['distMetalA']} A | Partner={payload['distPartnerA']} A")
    print(f"  Deep Mutational Scores: AlphaMissense={payload['alphamissenseScore']} ({payload['alphamissenseClass']}) | ESM-2 (3B) LLR={payload['esm2_3b_llr']:+.3f} | ESM-1v LLR={payload['esm1v_650m_llr']:+.3f}")
    print(f"  3D Inverse Folding    : ESM-IF1 3D LLR={payload['esm_if1_3d_llr']:+.3f} => Biophysical regime: {payload['biophysicalRegime']}")
    print(f"  Germline Constraint   : gnomAD v4 pLI={payload['gnomad_pLI']} | LOEUF={payload['gnomad_LOEUF']} | Missense Z={payload['gnomad_mis_z']}")
    print(f"  Somatic Tumour Counts : TCGA PanCanAtlas (N=10,195)={payload['tcgaPanCanSomaticCount_N10195']} | MSK-IMPACT (N=10,945)={payload['mskImpactSomaticCount_N10945']}")
    if payload["clinvarRecords"]:
        for c in payload["clinvarRecords"]:
            vcf_str = f" | GRCh38 {c['grch38Vcf']}" if c["grch38Vcf"] else ""
            print(f"  NCBI ClinVar          : VCV/{c['clinvarVariationId']} | {c['clinicalSignificance']} ({c['reviewStars']}*) | {c['phenotypes']}{vcf_str}")
    else:
        print("  NCBI ClinVar          : Not yet reported in ClinVar single-gene submissions")
    print("=" * 78)


def query_codon_profile(gene: str, codon: int, as_json: bool = False):
    gmeta = load_gene_metadata(gene)
    sym = gmeta["geneSymbol"]
    rows = []
    with gzip.open(DATA_DIR / "miRNA_biogenesis_saturation_mutagenesis_213940_atlas.csv.gz", "rt", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            if r["geneSymbol"] == sym and int(r["codonPos"]) == codon:
                am_f = float(r["alphamissenseScore"]) if r["alphamissenseScore"] != "" else None
                r["biophysicalRegime"] = classify_biophysical_regime(am_f, float(r["esm2_3b_llr"]), float(r["esm_if1_3d_llr"]))
                rows.append(r)
    if not rows:
        sys.exit(f"Codon {codon} out of range for {sym} (length {gmeta['proteinLengthAa']} aa).")
    rows.sort(key=lambda r: -(float(r["alphamissenseScore"]) if r["alphamissenseScore"] != "" else -1.0))
    if as_json:
        print(json.dumps({"geneSymbol": sym, "codonPos": codon, "substitutions": rows}, indent=2))
        return
    r0 = rows[0]
    print("=" * 88)
    print(f"miRBioGenAtlas Codon Saturation Scan — {sym} Codon {codon} (WT {r0['wtAA']}, PDB {r0['pdbComplexId'] or 'AlphaFold v6'})")
    print(f"3D Distances: Substrate={r0['distSubstrateA'] or 'NA'} A | Metal={r0['distMetalA'] or 'NA'} A | Partner={r0['distPartnerA'] or 'NA'} A | pLDDT={r0['afdbPlddt']}")
    print("-" * 88)
    print(f"{'Variant':<10} {'AM_Score':>9} {'AM_Class':<20} {'ESM2_3B':>9} {'ESM1v':>9} {'ESM_IF1':>9} {'TCGA':>5}  {'Biophysical_Regime'}")
    print("-" * 88)
    for r in rows:
        print(
            f"p.{r['variant']:<8} {r['alphamissenseScore'] or 'NA':>9} {r['alphamissenseClass'] or 'NA':<20} "
            f"{float(r['esm2_3b_llr']):+9.3f} {float(r['esm1v_650m_llr']):+9.3f} {float(r['esm_if1_3d_llr']):+9.3f} "
            f"{int(r['observedTcgaSomaticCount']):>5}  {r['biophysicalRegime']}"
        )
    print("=" * 88)


def main():
    ap = argparse.ArgumentParser(description="Query miRBioGenAtlas v2.0.2 structural, PLM, clinical, and somatic profiles.")
    ap.add_argument("gene", nargs="?", help="Gene symbol (e.g. DICER1, DROSHA, AGO2)")
    ap.add_argument("variant", nargs="?", help="Missense substitution (e.g. E1705K or p.Glu1705Lys)")
    ap.add_argument("--codon", nargs=2, metavar=("GENE", "POS"), help="List all 19 single-AA substitutions at a codon (e.g. --codon DICER1 1705)")
    ap.add_argument("--json", action="store_true", help="Emit structured JSON output")
    args = ap.parse_args()

    if args.codon:
        query_codon_profile(args.codon[0], int(args.codon[1]), as_json=args.json)
    elif args.gene:
        resolved = resolve_genomic_or_cdna(args.gene, args.variant)
        if resolved:
            query_single_variant(resolved[0], resolved[1], as_json=args.json)
        elif args.variant:
            query_single_variant(args.gene, args.variant, as_json=args.json)
        else:
            ap.print_help()
            sys.exit(1)
    else:
        ap.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()
