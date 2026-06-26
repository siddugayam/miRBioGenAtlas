# Data Dictionary — miRBioGenAtlas

**Dataset:** Pan-cancer alteration landscape of the miRNA biogenesis machinery (TCGA PanCancer Atlas, 32 tumour types, n = 10,195).
**Version:** 1.0.0 · **Source:** cBioPortal public REST API (`https://www.cbioportal.org/api`).

## Conventions (apply to all files)

- **Encoding / delimiter:** UTF-8, comma-delimited, one header row.
- **Frequencies (`*_freq`)** are **fractions in [0, 1]**, *not* percentages. Multiply by 100 for a percentage.
- **Counts (`*_n`)** are non-negative integers; the denominator is always **`N`** for the corresponding study/cohort.
- **`N` (denominator)** = number of samples in the **intersection of the `_sequenced` and CNA sample lists** for that study (every counted sample contributes both a mutation and a CNA call).
- **Copy-number classes** follow cBioPortal discrete (GISTIC2) calls:
  `amp` = high-level amplification (+2), `gain` = low-level gain (+1), `shallow_del` = shallow/heterozygous deletion (−1), `deep_del` = deep/homozygous deletion (−2).
- **`any_cna`** = sample has **any** of {gain, shallow_del, amp, deep_del} (CN ≠ 0).
- **`any_alter`** = **union** of mutation and any CNA (a sample with both is counted **once**, not twice).
- **Genes (15):** DROSHA, DGCR8 (Microprocessor); XPO5 (nuclear export); DICER1, TARBP2, PRKRA (cytoplasmic processing); AGO1–AGO4 (RISC/Argonaute); LIN28A, LIN28B (let-7 regulators); SMAD2, SMAD4 (Microprocessor **accessory** factors); TP53 (**upstream regulator — accessory, not a core biogenesis enzyme**; excluded from biogenesis-burden rankings).

---

## 1. `tcga_pancan_miRNA_biogenesis_alteration_freqs_long_FINAL.csv`
**Master long-format table — 480 rows (15 genes × 32 studies) × 23 columns.** Primary file for downstream analysis.

| Column | Type | Units / range | Definition |
|---|---|---|---|
| `tumour_acronym` | string | — | TCGA cancer-type code (e.g. `BRCA`). |
| `studyId` | string | — | cBioPortal study identifier (e.g. `brca_tcga_pan_can_atlas_2018`). |
| `studyName` | string | — | Full study name. |
| `geneSymbol` | string | HGNC | Gene symbol. |
| `entrezGeneId` | integer | NCBI | Entrez Gene ID. |
| `N` | integer | samples | Intersection denominator (sequenced ∩ CNA). |
| `mut_n` | integer | count | Samples with ≥1 non-synonymous somatic mutation. |
| `mut_freq` | float | [0,1] | `mut_n / N`. |
| `amp_n` | integer | count | Samples with high-level amplification (+2). |
| `amp_freq` | float | [0,1] | `amp_n / N`. |
| `gain_n` | integer | count | Samples with low-level gain (+1). |
| `gain_freq` | float | [0,1] | `gain_n / N`. |
| `shallow_del_n` | integer | count | Samples with shallow deletion (−1). |
| `shallow_del_freq` | float | [0,1] | `shallow_del_n / N`. |
| `deep_del_n` | integer | count | Samples with deep deletion (−2). |
| `deep_del_freq` | float | [0,1] | `deep_del_n / N`. |
| `any_cna_n` | integer | count | Samples with any CNA (gain ∪ shallow_del ∪ amp ∪ deep_del). |
| `any_cna_freq` | float | [0,1] | `any_cna_n / N`. |
| `any_alter_n` | integer | count | Samples with any alteration (mutation ∪ any CNA). |
| `any_alter_freq` | float | [0,1] | `any_alter_n / N`. |
| `mut_only_n` | integer | count | Samples mutated but no CNA. |
| `cna_only_n` | integer | count | Samples with CNA but no mutation. |
| `mut_and_cna_n` | integer | count | Samples with both a mutation and a CNA. |

**Internal check (holds for all 480 rows):** `any_alter_n == mut_only_n + cna_only_n + mut_and_cna_n`.

---

## 2–4. Matrix files (gene × tumour)
`tcga_pancan_miRNA_biogenesis_any_alter_freq_matrix.csv` · `..._any_cna_freq_matrix.csv` · `..._mut_freq_matrix.csv`
**15 data rows (genes) × 33 physical columns** = `geneSymbol` index column **+ 32 tumour-type columns**.

| Column | Type | Units / range | Definition |
|---|---|---|---|
| `geneSymbol` | string | HGNC | Gene symbol (row index). |
| `ACC`, `BLCA`, … `UVM` (32 columns) | float | [0,1] | Frequency for that gene in that tumour type. Metric depends on the file: **any-alteration** / **CNA-only** / **mutation-only** respectively. |

*Tumour columns (alphabetical):* ACC, BLCA, BRCA, CESC, CHOL, COADREAD, DLBC, ESCA, GBM, HNSC, KICH, KIRC, KIRP, LAML, LGG, LIHC, LUAD, LUSC, MESO, OV, PAAD, PCPG, PRAD, SARC, SKCM, STAD, TGCT, THCA, THYM, UCEC, UCS, UVM.

---

## 5. `miRNA_biogenesis_gene_rank_summary.csv`
**Gene ranking — 15 rows × 9 columns.** Means/maxima taken across the 32 studies.

| Column | Type | Range | Definition |
|---|---|---|---|
| `rank` | integer | 1–15 | Rank by `mean_any_alter_freq` (descending). |
| `geneSymbol` | string | HGNC | Gene symbol. |
| `mean_any_alter_freq` | float | [0,1] | Mean any-alteration frequency across 32 studies. |
| `mean_mut_freq` | float | [0,1] | Mean mutation frequency. |
| `mean_any_cna_freq` | float | [0,1] | Mean any-CNA frequency. |
| `mean_amp_freq` | float | [0,1] | Mean amplification frequency. |
| `mean_deep_del_freq` | float | [0,1] | Mean deep-deletion frequency. |
| `max_any_alter_freq` | float | [0,1] | Maximum any-alteration frequency across studies. |
| `n_studies` | integer | 32 | Number of studies averaged. |

> **Note:** TP53 ranks #1 here (mean 0.5507). Among **biogenesis** genes (TP53 excluded), **AGO2** is highest (mean 0.4816 ≈ 48.2%). Wording in the manuscript should always read "most-altered biogenesis gene (TP53 excluded)."

---

## 6. `miRNA_biogenesis_top5_tumours_per_gene.csv`
**Top-5 tumour types per gene — 75 rows (15 genes × 5) × 9 columns.**

| Column | Type | Range | Definition |
|---|---|---|---|
| `geneSymbol` | string | HGNC | Gene symbol. |
| `rank` | integer | 1–5 | Rank of the tumour type for this gene (by `any_alter_freq`). |
| `tumour_acronym` | string | — | Tumour-type code. |
| `N` | integer | samples | Intersection denominator for that study. |
| `any_alter_freq` | float | [0,1] | Any-alteration frequency. |
| `mut_freq` | float | [0,1] | Mutation frequency. |
| `any_cna_freq` | float | [0,1] | Any-CNA frequency. |
| `amp_freq` | float | [0,1] | Amplification frequency. |
| `deep_del_freq` | float | [0,1] | Deep-deletion frequency. |

---

## 7. `miRNA_biogenesis_top_tumours_excluding_TP53.csv`
**Tumour-type burden ranking (TP53 excluded) — 32 rows × 7 columns.**

| Column | Type | Range | Definition |
|---|---|---|---|
| `rank` | integer | 1–32 | Rank by `mean_any_alter_freq_excl_TP53` (descending; OV = 1, LAML = 32). |
| `tumour_acronym` | string | — | Tumour-type code. |
| `mean_any_alter_freq_excl_TP53` | float | [0,1] | Mean any-alteration frequency across the **14** biogenesis genes (TP53 excluded). |
| `mean_mut_freq_excl_TP53` | float | [0,1] | Mean mutation frequency across 14 genes. |
| `mean_any_cna_freq_excl_TP53` | float | [0,1] | Mean any-CNA frequency across 14 genes. |
| `N` | integer | samples | Intersection denominator for that study. |
| `n_genes` | integer | 14 | Genes averaged (TP53 excluded). |

---

## 8. `brca_miRNA_biogenesis_alteration_frequencies.csv`
**Breast invasive carcinoma deep-dive (15 rows × 15 columns). `N = 1,052`.**

> **Note:** this deep-dive cohort uses the same BRCA denominator (N = 1,052) as the pan-cancer master table and sample-set summary. Column names in this file are `gene`/`entrez_id` (the master table uses `geneSymbol`/`entrezGeneId`).

| Column | Type | Range | Definition |
|---|---|---|---|
| `gene` | string | HGNC | Gene symbol. |
| `entrez_id` | integer | NCBI | Entrez Gene ID. |
| `mut_n` | integer | count | Samples with a somatic mutation. |
| `amp_n` | integer | count | Samples with amplification (+2). |
| `gain_n` | integer | count | Samples with gain (+1). |
| `shallow_del_n` | integer | count | Samples with shallow deletion (−1). |
| `deep_del_n` | integer | count | Samples with deep deletion (−2). |
| `any_cna_n` | integer | count | Samples with any CNA. |
| `any_alter_n` | integer | count | Samples with any alteration (mutation ∪ CNA). |
| `N` | integer | 1052 | Cohort denominator. |
| `mut_freq` | float | [0,1] | `mut_n / N`. |
| `amp_freq` | float | [0,1] | `amp_n / N`. |
| `deep_del_freq` | float | [0,1] | `deep_del_n / N`. |
| `any_cna_freq` | float | [0,1] | `any_cna_n / N`. |
| `any_alter_freq` | float | [0,1] | `any_alter_n / N`. |

---

## 9. `tcga_pancan_atlas_sample_set_summary.csv` (CORRECTED — 32 rows × 8 columns)
Use the regenerated file `tcga_pancan_atlas_sample_set_summary_CORRECTED.csv`. **Six columns were removed** from the original because they contained broadcast placeholder values copied from unrelated studies (`cnaseq_sampleListId`, `cnaseq_n`, `threeway_sampleListId`, `threeway_n`) and two ratios derived from them (`ratio_cnaseq_over_intersection`, `ratio_threeway_over_intersection`). The retained columns are all verified.

| Column | Type | Range | Definition |
|---|---|---|---|
| `studyId` | string | — | cBioPortal study identifier. |
| `cancerTypeId` | string | OncoTree | OncoTree cancer-type code (note: GBM and LGG both map to `difg` = *diffuse glioma*). |
| `studyName` | string | — | Full study name. |
| `sampleCount_all` | integer | samples | All samples in the study. |
| `sampleCount_sequenced` | integer | samples | Samples in the `_sequenced` list. |
| `sampleCount_cna` | integer | samples | Samples with CNA data. |
| `intersection_sequenced_cna_n` | integer | samples | **`N`** used throughout (sequenced ∩ CNA). Sum across 32 studies = **10,195**. |
| `ratio_intersection_over_min_seq_cna` | float | (0,1] | `intersection_sequenced_cna_n / min(sampleCount_sequenced, sampleCount_cna)`. |

---

## 10. `tcga_low_coverage_tumours_for_miRNA_biogenesis.txt`
Plain-text note listing the **seven** tumour types with `N < 100` (CHOL 36, DLBC 41, UCS 56, KICH 65, UVM 80, MESO 86, ACC 89) and recommending COSMIC supplementation. Not tabular.
