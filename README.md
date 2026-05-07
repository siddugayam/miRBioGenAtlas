# miRBioGenAtlas

**Pan-cancer alteration landscape of the miRNA biogenesis machinery**

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Data: TCGA](https://img.shields.io/badge/Data-TCGA%20PanCancer%20Atlas-green)](https://www.cbioportal.org)
[![R version](https://img.shields.io/badge/R-%3E%3D4.2-blue)](https://www.r-project.org)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)

---

## Overview

**miRBioGenAtlas** is a reproducible data resource and analysis pipeline that systematically quantifies copy number alteration (CNA) and somatic mutation frequencies for **15 miRNA biogenesis pathway genes** across **32 TCGA PanCancer Atlas tumour types** (n = 10,195 samples).

This repository accompanies the manuscript:

> Gayam PK. *Pan-cancer alteration landscape of the miRNA biogenesis machinery: a TCGA analysis across 32 tumour types.* [Journal name, 2026] [DOI pending]

---

## Key findings

- **AGO2** is the most frequently altered biogenesis gene across cancers (mean **48.2%**), driven predominantly by copy number gain at 8q24
- **Copy number alteration dominates** over somatic mutation across all 15 biogenesis genes
- **Ovarian serous cystadenocarcinoma (OV)** shows the highest pan-biogenesis alteration burden (**65.8%**)
- **Thyroid carcinoma (THCA)** and **acute myeloid leukaemia (LAML)** show the lowest burden (3.1% and 3.0%), suggesting tissue-specific pathway preservation
- In breast cancer (BRCA, n=996), **AGO2 amplification** was detected in **11.3%** of samples

---

## Repository structure

```
miRBioGenAtlas/
│
├── data/                          # All CSV data files
│   ├── tcga_pancan_miRNA_biogenesis_alteration_freqs_long_FINAL.csv
│   ├── tcga_pancan_miRNA_biogenesis_any_alter_freq_matrix.csv
│   ├── tcga_pancan_miRNA_biogenesis_any_cna_freq_matrix.csv
│   ├── tcga_pancan_miRNA_biogenesis_mut_freq_matrix.csv
│   ├── miRNA_biogenesis_gene_rank_summary.csv
│   ├── miRNA_biogenesis_top5_tumours_per_gene.csv
│   ├── miRNA_biogenesis_top_tumours_excluding_TP53.csv
│   ├── brca_miRNA_biogenesis_alteration_frequencies.csv
│   ├── tcga_pancan_atlas_sample_set_summary.csv
│   └── tcga_low_coverage_tumours_for_miRNA_biogenesis.txt
│
├── notebooks/
│   └── notebook.ipynb             # Python data retrieval pipeline (cBioPortal API)
│
├── figures/                       # Generated TIFF figures (300 DPI)
│   ├── Fig1_pancan_heatmap.tiff
│   ├── Fig2_gene_ranking.tiff
│   ├── Fig3_tumour_burden_ranking.tiff
│   ├── Fig4_BRCA_deepdive.tiff
│   └── Fig_combined_panel.tiff
│
├── miRBioGenAtlas_figures.R       # R figure generation script
├── README.md
└── LICENSE
```

---

## Genes analysed

| Category | Genes |
|---|---|
| Microprocessor complex | DROSHA, DGCR8, SMAD2, SMAD4 |
| Nuclear export | XPO5 |
| Cytoplasmic processing | DICER1, TARBP2 (TRBP), PRKRA (PACT) |
| RISC / Argonaute | AGO1, AGO2, AGO3, AGO4 |
| LIN28 regulators | LIN28A, LIN28B |
| Upstream regulator | TP53 |

---

## Reproducing the figures

### Requirements

- R ≥ 4.2
- Packages: `tidyverse`, `pheatmap`, `RColorBrewer`, `ggpubr`, `scales`, `cowplot`, `viridis`, `ggtext`, `patchwork`

### Steps

```r
# 1. Clone this repository
# git clone https://github.com/siddugayam/miRBioGenAtlas.git
# cd miRBioGenAtlas

# 2. Open R and run
source("miRBioGenAtlas_figures.R")

# Output: 5 TIFF files in the figures/ directory
```

All packages are automatically installed if not already present.

---

## Reproducing the data retrieval

The Python notebook (`notebooks/notebook.ipynb`) retrieves data directly from the **cBioPortal public API** and is fully reproducible.

### Requirements

```bash
pip install requests pandas jupyter
```

### Steps

```bash
cd notebooks
jupyter notebook notebook.ipynb
# Run all cells — data is downloaded from cBioPortal API
```

> No login or API key required. cBioPortal data is publicly available.

---

## Data description

| File | Description | Rows | Columns |
|---|---|---|---|
| `*_long_FINAL.csv` | Master long-format table | ~480 | gene, tumour, N, frequencies |
| `*_any_alter_freq_matrix.csv` | Gene × tumour matrix (any alteration) | 15 | 32 |
| `*_any_cna_freq_matrix.csv` | Gene × tumour matrix (CNA only) | 15 | 32 |
| `*_mut_freq_matrix.csv` | Gene × tumour matrix (mutation only) | 15 | 32 |
| `gene_rank_summary.csv` | Gene ranking across 32 studies | 15 | 8 |
| `top5_tumours_per_gene.csv` | Top 5 tumours per gene | 75 | 9 |
| `top_tumours_excluding_TP53.csv` | Tumour burden ranking (14 genes) | 32 | 7 |
| `brca_*.csv` | BRCA deep-dive (n=996) | 15 | 14 |
| `sample_set_summary.csv` | Study metadata and sample counts | 32 | 14 |

---

## Tumour type abbreviations

| Code | Cancer type | N |
|---|---|---|
| OV | Ovarian serous cystadenocarcinoma | 511 |
| BRCA | Breast invasive carcinoma | 1,052 |
| LUSC | Lung squamous cell carcinoma | 484 |
| LUAD | Lung adenocarcinoma | 511 |
| COADREAD | Colorectal adenocarcinoma | 532 |
| STAD | Stomach adenocarcinoma | 434 |
| ESCA | Oesophageal adenocarcinoma | 182 |
| HNSC | Head and neck squamous cell carcinoma | 509 |
| UCEC | Uterine corpus endometrial carcinoma | 511 |
| GBM | Glioblastoma multiforme | 380 |
| *...and 22 more* | | |

Full list in `data/tcga_pancan_atlas_sample_set_summary.csv`

---

## Citation

If you use this data or code, please cite:

```bibtex
@article{gayam2026mirbiogenatlas,
  title   = {Pan-cancer alteration landscape of the {miRNA} biogenesis machinery:
             a {TCGA} analysis across 32 tumour types},
  author  = {Gayam, Prasanna Kumar Reddy},
  journal = {[Journal name]},
  year    = {2026},
  doi     = {[DOI pending]}
}
```

Please also cite the cBioPortal:

> Cerami E, et al. The cBio Cancer Genomics Portal: An open platform for exploring multidimensional cancer genomics data. *Cancer Discov.* 2012;2(5):401–404.

---

## Data availability

All underlying data were obtained from publicly available TCGA PanCancer Atlas studies via the cBioPortal API. No individual patient data are included. The processed dataset is also archived at Zenodo: [DOI pending].

---

## License

This project is licensed under the **MIT License** — see [LICENSE](LICENSE) for details.

---

## Contact

**Gayam Prasanna Kumar Reddy**
Manipal Academy of Higher Education, India
gayam.prasanna@learner.manipal.edu
siddugayam1994@gmail.com

*Issues and pull requests welcome.*
