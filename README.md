# miRBioGenAtlas (v2.0.0)

**Sample-level and aggregated multi-modal somatic alteration dataset of the miRNA biogenesis machinery across 32 TCGA cancer studies (`N = 10,195` tumours)**

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24337287.svg)](https://doi.org/10.5281/zenodo.24337287)
[![Code License: MIT](https://img.shields.io/badge/Code_License-MIT-blue.svg)](LICENSE)
[![Data License: CC BY 4.0](https://img.shields.io/badge/Data_License-CC_BY_4.0-green.svg)](LICENSE-data.md)
[![cBioPortal Snapshot](https://img.shields.io/badge/cBioPortal-v7.2.3_(DB_3.1.0)-informational)](https://www.cbioportal.org)

---

## Overview

**miRBioGenAtlas (v2.0.0)** provides a frozen, SHA-256-verified, sample-level and study-level multi-omic dataset covering **15 canonical microRNA (miRNA) biogenesis and regulatory genes** across **32 TCGA PanCancer Atlas (2018) studies** (spanning 33 cancer types; `COAD = 336` and `READ = 196` combined in `COADREAD = 532`; total matched `sequenced ∩ CNA` denominator **`N = 10,195` tumours**).

### What is included in v2.0.0
1. **Frozen input archive (`data/raw_snapshot/`, 9.07 MB):** Versioned, SHA-256-checksummed cBioPortal REST API (`v7.2.3`, DB `3.1.0`) responses across 6 modalities (`mutations`, discrete `gistic`, continuous `log2CNA`, chromosome-arm `armlevel_cna`, `structural_variants`, matched `rna_seq_v2_mrna` RSEM/Z-scores, and sample/patient clinical/TMB/aneuploidy attributes).
2. **Sample-level multi-modal table (`152,925` rows = `10,195` samples × `15` genes):** Discrete GISTIC2 (`-2..+2`), continuous `log2` copy-number ratio, matched chromosome-arm CNA status (`Broad_ArmLevel_LowAmp` vs `HighLevel_Amp_or_DeepDel`), functional variant counts/flags, structural variants, matched mRNA expression (`n = 9,620` tumours), and sample hypermutator status.
3. **Annotated somatic variants (`6,292` coding calls) & structural variants (`72` calls):** Stratified into `Truncating_LoF`, `Hotspot_Missense_or_InFrame` (e.g. `DICER1` RNase IIIa `S1344` & RNase IIIb `E1705/D1709/G1809/D1810/E1813`, `SMAD4` MH2 `R361/D351/D537`), `InFrame_Indel`, and unannotated `Other_Missense` (VUS), paired with tumour VAF and sample `TMB`/MSI/`POLE` hypermutator flags.
4. **Synteny-aware patient-level pathway burden (`10,195` samples):** Proportion of patients carrying $\ge 1$ high-confidence (`mutation ∪ GISTIC2 ±2`) or broad alteration across three nested gene-set tiers (**Core-10**, **Extended-12**, **Panel-14**), collapsing tandemly clustered loci (`1p34.3 AGO1/AGO3/AGO4`; `18q21.1 SMAD2/SMAD4`) to prevent double-counting.
5. **Aggregated frequency tables & 4 technical validation audit tables:** Full backward compatibility with v1.0 gene × study matrices plus cBioPortal 20-pair verification, Hoadley et al. (2018) denominator audit, *cis*-CNA-to-mRNA dosage concordance, and arm-vs-focal CNA decomposition.

---

## Repository Structure

```text
miRBioGenAtlas/
├── data/
│   ├── raw_snapshot/                                                  # Frozen SHA-256-verified API inputs (9.07 MB)
│   ├── tcga_pancan_miRNA_biogenesis_sample_level_alterations.csv.gz   # 152,925 sample x gene multi-modal rows
│   ├── tcga_pancan_miRNA_biogenesis_somatic_mutations_annotated.csv   # 6,292 annotated somatic variant calls
│   ├── tcga_pancan_miRNA_biogenesis_structural_variants.csv           # 72 somatic structural variants / fusions
│   ├── tcga_pancan_miRNA_biogenesis_patient_level_burden.csv          # 10,195 patient-level burden & TMB/FGA records
│   ├── miRNA_biogenesis_gene_panel_metadata.csv                       # 15-gene coordinates, CDS length, tiers & rationale
│   ├── miRNA_biogenesis_excluded_candidates_rationale.csv             # Exclusion rationale for non-core candidates
│   ├── tcga_pancan_miRNA_biogenesis_alteration_freqs_long_FINAL.csv   # 480-row gene x study master summary table
│   ├── tcga_pancan_miRNA_biogenesis_*_matrix.csv                      # 5 wide (15 x 32) frequency matrices
│   ├── miRNA_biogenesis_gene_rank_summary.csv                         # Gene-level summary & CDS-normalized densities
│   ├── miRNA_biogenesis_top_tumours_excluding_TP53.csv                # Study-level patient-burden & FGA summary
│   ├── miRNA_biogenesis_top5_tumours_per_gene.csv                     # Top 5 studies per gene
│   ├── brca_miRNA_biogenesis_alteration_frequencies.csv               # Breast carcinoma cohort summary (N = 1,052)
│   ├── tcga_pancan_atlas_sample_set_summary.csv                       # Study denominators & hypermutator counts
│   ├── tcga_low_coverage_tumours_for_miRNA_biogenesis.txt             # Flagged low-N cohorts (N < 100)
│   └── validation_*.csv                                               # 4 technical validation audit tables
├── scripts/
│   ├── fetch_cbioportal_snapshot.py                                   # Live API snapshot retriever (Python stdlib)
│   └── build_mirbiogen_v2_dataset.py                                  # Deterministic offline dataset builder (<10 s)
├── notebooks/
│   └── notebook.ipynb                                                 # Portable standalone Jupyter pipeline
├── figures/
│   └── Figure_1..Figure_4 (.pdf, .tiff, .png)                         # Vector PDF + 450-DPI TIFF + 300-DPI PNG figures
├── miRBioGenAtlas_figures.R                                           # R >= 4.2 figure-generation script
├── DATA_DICTIONARY.md                                                 # Complete column-by-column schema
├── CITATION.cff                                                       # Citation metadata
├── .zenodo.json                                                       # Zenodo deposit metadata
├── requirements.txt                                                   # Python environment specification
├── LICENSE                                                            # MIT License (code)
└── LICENSE-data.md                                                    # CC BY 4.0 License (data)
```

---

## Reproducing the Dataset & Figures

### 1. Deterministic offline rebuild from frozen snapshot (~8 seconds, no internet required)
```bash
python3 scripts/build_mirbiogen_v2_dataset.py
```
Or open and run `notebooks/notebook.ipynb` (`EXECUTION_MODE = "FROZEN_SNAPSHOT"`). Both scripts use only the Python $\ge 3.9$ standard library with relative `pathlib` paths.

### 2. Optional live re-query against cBioPortal REST API (~60–90 seconds)
```bash
python3 scripts/fetch_cbioportal_snapshot.py
python3 scripts/build_mirbiogen_v2_dataset.py
```

### 3. Regenerating all publication figures (Vector PDF + 450-DPI LZW TIFF)
Requires R $\ge 4.2$ with `tidyverse`, `cowplot`, `patchwork`, `scales`, and `ragg`:
```bash
Rscript miRBioGenAtlas_figures.R
```

---

## Licensing

- **Data files (`data/`):** Creative Commons Attribution 4.0 International ([CC BY 4.0](LICENSE-data.md))
- **Source code (`scripts/`, `notebooks/`, `miRBioGenAtlas_figures.R`):** [MIT License](LICENSE)
