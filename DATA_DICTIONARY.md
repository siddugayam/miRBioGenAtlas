# Data Dictionary — miRBioGenAtlas (v2.0.0)

**Dataset:** Pan-cancer alteration landscape of the miRNA biogenesis machinery across 32 TCGA cancer studies (`N = 10,195` matched tumours; `152,925` sample × gene records).  
**Version:** `2.0.0` · **Zenodo Dataset DOI:** [`10.5281/zenodo.24337287`](https://doi.org/10.5281/zenodo.24337287) · **Source snapshot:** cBioPortal public REST API `v7.2.3` (`DB v3.1.0`; frozen SHA-256-verified archive in `data/raw_snapshot/`).

## Global Conventions

- **Encoding / delimiter:** UTF-8, comma-delimited (`.csv` or gzip-compressed `.csv.gz`), single header row.
- **Frequencies (`*_freq`, `*_prop`):** Exact fractions in `[0, 1]` (not percentages).
- **Primary sample denominator (`N`):** Per study, `N` is the exact cardinality of `<studyId>_sequenced ∩ <studyId>_cna` (`10,195` tumours across 32 TCGA PanCancer Atlas studies spanning 33 cancer types, with `COAD = 336` and `READ = 196` combined in `COADREAD = 532`).
- **Nested gene-set tiers (15 genes total):**
  - **Core-10 (`in_Core10_tier = 1`):** `DROSHA, DGCR8, XPO5, DICER1, TARBP2, PRKRA, AGO1, AGO2, AGO3, AGO4`.
  - **Extended-12 (`in_Extended12_tier = 1`):** Core-10 + direct pri/pre-let-7 regulators `LIN28A, LIN28B`.
  - **Panel-14 (`in_Panel14_tier = 1`):** Extended-12 + Microprocessor accessory R-SMAD/Co-SMAD factors `SMAD2, SMAD4`.
  - **All-15:** Panel-14 + upstream pleiotropic regulator `TP53` (excluded from biogenesis burden summaries).
- **High-confidence (`strict`) vs broad (`any`) alterations:**
  - `high_conf` / `strict` = non-synonymous somatic coding mutation **OR** high-level focal GISTIC2 CNA (`+2` amplification or `-2` homozygous deep deletion).
  - `any_alter` = non-synonymous somatic coding mutation **OR** any discrete GISTIC2 CNA (`-2, -1, +1, +2`).

---

## 1. Primary Sample-Level & Variant-Level Tables

### 1.1 `tcga_pancan_miRNA_biogenesis_sample_level_alterations.csv.gz`
**152,925 rows (`10,195` samples × `15` genes) × 24 columns.**

| Column | Type | Definition |
|---|---|---|
| `tumour_acronym` | string | TCGA study code (`ACC`..`UVM`). |
| `studyId` | string | cBioPortal study ID (`*_tcga_pan_can_atlas_2018`). |
| `sampleId` | string | TCGA sample barcode (e.g. `TCGA-A1-A0SB-01`). |
| `patientId` | string | TCGA participant barcode (e.g. `TCGA-A1-A0SB`). |
| `geneSymbol` | string | HGNC gene symbol. |
| `entrezGeneId` | integer | NCBI Entrez Gene ID. |
| `chrArm` | string | Chromosome arm (`5p`, `22q`, `8q`, `1p`, etc.). |
| `syntenicCluster` | string | Syntenic locus ID (`1p34.3_AGO1_AGO3_AGO4` groups tandem `AGO1/3/4`; `18q21.1_SMAD2_SMAD4` groups `SMAD2/4`). |
| `gistic2DiscreteCNA` | integer | Discrete GISTIC2 call (`-2` deep deletion, `-1` shallow deletion, `0` diploid/neutral, `+1` low-level gain, `+2` high-level amplification). |
| `log2CopyNumber` | float | Continuous `log2` copy-number ratio from `<studyId>_log2CNA`. |
| `chrArmGisticStatus` | string | Matched chromosome-arm GISTIC status (`Gain`, `Loss`, `Unchanged`, or blank if unassessed) from `<studyId>_armlevel_cna`. |
| `cnaFocalityClass` | string | `Diploid_Neutral`, `HighLevel_Amp_or_DeepDel`, `Broad_ArmLevel_LowAmp` (`+/-1` concordant with whole-arm gain/loss), `SubArm_LowAmp`, or `LowAmp_ArmUnannotated`. |
| `nonsynMutationCount` | integer | Number of non-synonymous somatic coding mutations in this sample and gene. |
| `hasTruncatingLoF` | integer (0/1) | `1` if harbouring nonsense, frameshift indel, canonical splice-site, start-loss, or nonstop variant. |
| `hasCuratedHotspot` | integer (0/1) | `1` if harbouring a curated functional hotspot codon (`DICER1` S1344/E1705/D1709/R1725/G1809/D1810/E1813; `DROSHA` E1147/D1151/Q1187/E1222; `SMAD4` D351/R361/G386/D537; `TP53` R175/Y220/G245/R248/R249/R273/R282; etc.). |
| `proteinChanges` | string | Semicolon-delimited HGVSp short protein changes (`D1709N`, `R361H`, etc.). |
| `consequenceTiers` | string | Semicolon-delimited consequence tiers (`Truncating_LoF`, `Hotspot_Missense_or_InFrame`, `InFrame_Indel`, `Other_Missense`, `Other_Nonsynonymous`). |
| `hasStructuralVariant` | integer (0/1) | `1` if sample carries a somatic structural variant / gene fusion involving this gene. |
| `structuralVariantEvents` | string | Semicolon-delimited SV event descriptions. |
| `hasHighConfAlteration` | integer (0/1) | `1` if `nonsynMutationCount > 0` OR `gistic2DiscreteCNA in {-2, +2}`. |
| `hasAnyAlteration` | integer (0/1) | `1` if `nonsynMutationCount > 0` OR `gistic2DiscreteCNA != 0`. |
| `mrnaRsem` | float | Batch-normalized RSEM mRNA expression (`<studyId>_rna_seq_v2_mrna`; blank if no RNA-seq). |
| `mrnaZscore` | float | mRNA expression Z-score (`<studyId>_rna_seq_v2_mrna_median_all_sample_Zscores`; `n = 9,620` samples). |
| `sampleIsHypermutated` | integer (0/1) | `1` if `TMB_NONSYNONYMOUS >= 10` mut/Mb (or `MUTATION_COUNT >= 300` where TMB unavailable) or `POLE`/`MSI` hypermutator subtype (`1,385 / 10,195` samples). |

### 1.2 `tcga_pancan_miRNA_biogenesis_somatic_mutations_annotated.csv`
**6,292 rows (individual somatic coding variant calls) × 26 columns.**  
Contains genomic coordinates (`ncbiBuild = GRCh37`, `chr`, `startPosition`, `endPosition`, `referenceAllele`, `variantAllele`), `variantType`, `mutationType`, `proteinChange`, `proteinPosStart`, `refseqMrnaId`, read depths (`tumorAltCount`, `tumorRefCount`, `tumorVAF`), functional `consequenceTier`, `isCuratedHotspot`, `isTruncatingLoF`, `sampleTmbNonsynonymous`, `sampleIsHypermutated`, and `molecularSubtype`.

### 1.3 `tcga_pancan_miRNA_biogenesis_structural_variants.csv`
**72 rows × 14 columns.** Somatic structural variants / fusions involving the 15 panel genes across the 10,195 samples (`site1HugoSymbol`, `site2HugoSymbol`, `eventInfo`, breakpoints, frame effect, and `svStatus`).

### 1.4 `tcga_pancan_miRNA_biogenesis_patient_level_burden.csv`
**10,195 rows (one per tumour sample) × 29 columns.**  
Provides sample-level `cancerTypeDetailed` (resolving `COAD` vs `READ`), `molecularSubtype`, `tmbNonsynonymous`, `mutationCountGenomeWide`, `isHypermutated`, `fractionGenomeAltered`, `aneuploidyScore`, binary patient-union indicators (`*_has_strict_alter`, `*_has_mut`, `*_has_focal_cna`, `*_has_any_alter`) across `Core10`, `Extended12`, `Panel14`, and `All15`, syntenic locus counts (`*_syntenic_loci_strict_n`, collapsing co-localized `1p34.3` and `18q21.1` gene clusters), and syntenic concordance flags (`ago1_3_4_concordant_1p34_cna`, `smad2_4_concordant_18q21_cna`).

---

## 2. Gene Panel Metadata & Rationale Tables

- **`miRNA_biogenesis_gene_panel_metadata.csv` (15 rows × 17 columns):** HGNC ID, Entrez ID, Ensembl ID, canonical RefSeq transcript, cytoband, chromosome arm, `syntenicCluster`, canonical CDS length (`cdsLengthBp`), protein length (`proteinLengthAa`), `pathwayStep`, `panelRole`, nested tier membership (`in_Core10_tier`, `in_Extended12_tier`, `in_Panel14_tier`), `curatedHotspotCodons`, and `inclusionRationale`.
- **`miRNA_biogenesis_excluded_candidates_rationale.csv` (8 rows × 4 columns):** Explicit rationale for excluding pleiotropic RNA helicases (`DDX5`, `DDX17`), downstream uridylyltransferases/exonucleases (`TUT4`, `TUT7`, `DIS3L2`), GW182 silencing scaffolds (`TNRC6A/B/C`), and general transport factors (`RAN`).

---

## 3. Aggregated Gene × Study Tables & Matrices

- **`tcga_pancan_miRNA_biogenesis_alteration_freqs_long_FINAL.csv` (480 rows × 35 columns):** Master `gene × study` table preserving all 23 v1.0 columns (`mut_n/freq`, `amp_n/freq`, `gain_n/freq`, `shallow_del_n/freq`, `deep_del_n/freq`, `any_cna_n/freq`, `any_alter_n/freq`, `mut_only_n`, `cna_only_n`, `mut_and_cna_n`) plus 12 v2.0 stratified columns (`truncating_mut_n/freq`, `hotspot_mut_n/freq`, `nonhypermut_N`, `nonhypermut_mut_n/freq`, `mut_density_per_cds_kb_per_1000`, `focal_cna_n/freq`, `high_conf_alter_n/freq`).
- **Five `15 × 32` wide matrices (`*_matrix.csv`):** `any_alter_freq`, `any_cna_freq`, `mut_freq`, `focal_cna_freq`, and `high_conf_alter_freq`.
- **Summary & cohort tables:** `miRNA_biogenesis_gene_rank_summary.csv` (15 rows × 18 cols), `miRNA_biogenesis_top5_tumours_per_gene.csv` (75 rows × 11 cols), `miRNA_biogenesis_top_tumours_excluding_TP53.csv` (32 rows × 19 cols, ranked by synteny-aware `patient_prop_core10_strict`), `brca_miRNA_biogenesis_alteration_frequencies.csv` (15 rows × 24 cols, `N = 1,052`), `tcga_pancan_atlas_sample_set_summary.csv` (32 rows × 13 cols), and `tcga_low_coverage_tumours_for_miRNA_biogenesis.txt`.

---

## 4. Technical Validation Audit Tables & Frozen Snapshot

- **`validation_20_cbioportal_spotchecks.csv` (20 rows × 17 cols):** Exact comparison between pipeline outputs and manual cBioPortal web queries across 20 gene–study pairs (`abs_freq_discrepancy = 0.0000` in 20/20 checks).
- **`validation_denominator_hoadley2018_audit.csv` (32 rows × 9 cols):** Study-by-study reconciliation of `cbioportal_all_n`, `sequenced_n`, `cna_n`, and `mirbiogen_intersection_N` (`10,195`) against Hoadley et al. (2018) *Cell* (`10,793` total cases across the 32 studies).
- **`validation_cis_mrna_dosage_concordance.csv` (15 rows × 14 cols):** Median mRNA Z-scores across GISTIC2 states `-2, -1, 0, +1, +2` and Spearman `rho` (`log2CNA` vs `mRNA Z-score`) across `n = 9,620` matched tumours.
- **`validation_arm_vs_focal_cna_decomposition.csv` (15 rows × 12 cols):** Pan-cancer decomposition of discrete GISTIC2 calls against matched chromosome-arm GISTIC status.
- **`data/raw_snapshot/` (`33` `.json.gz` files + `snapshot_manifest.json`, `9.07 MB`):** Frozen cBioPortal API responses verified by SHA-256 checksums before every build.
