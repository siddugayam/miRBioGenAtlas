# Data Dictionary — miRBioGenAtlas (v2.0.2)

**Dataset:** Multi-consortium, structurally and functionally calibrated alteration landscape of the miRNA biogenesis machinery across 32 TCGA cancer studies (`N = 10,195` matched tumours; `152,925` sample × gene records; `9,578` matched mature `miRNA-seq` profiles; `9,353` dual-profiled `miRNA × mRNA` tumours across `13,369,752` pairs; `213,940` proteome-wide saturation mutagenesis substitutions across `11,260` codons; `6,292` 3D Cryo-EM/PLM/AlphaMissense/ClinVar-annotated TCGA somatic variants; `16,569` NCBI ClinVar clinical variant records; `10,945` MSK-IMPACT clinical sequencing tumours; `810` GTEx v10 54-tissue healthy adult baselines; `1,002` CPTAC mass-spectrometry proteomes; `1,208` Broad DepMap CRISPR-Cas9 cell lines; TCGA-CDR curated survival and Thorsson Pan-Cancer Immune/HRD profiles).  
**Version:** `2.0.2` · **Zenodo Dataset DOI:** [`10.5281/zenodo.24455067`](https://doi.org/10.5281/zenodo.24455067) · **Interactive Explorer:** `docs/index.html` + `docs/sat/*.js` (`https://siddugayam.github.io/miRBioGenAtlas/`).

## Global Conventions

- **Encoding / delimiter:** UTF-8, comma-delimited (`.csv` or gzip-compressed `.csv.gz`), single header row.
- **Frequencies (`*_freq`, `*_prop`):** Exact fractions in `[0, 1]` (not percentages).
- **Primary sample denominator (`N`):** Per TCGA study, `N` is the exact cardinality of `<studyId>_sequenced ∩ <studyId>_cna` (`10,195` tumours across 32 TCGA PanCancer Atlas studies spanning 33 cancer types, with `COAD = 336` and `READ = 196` combined in `COADREAD = 532`).
- **Nested gene-set tiers (15 genes total):**
  - **Core-10 (`in_Core10_tier = 1`):** `DROSHA, DGCR8, XPO5, DICER1, TARBP2, PRKRA, AGO1, AGO2, AGO3, AGO4`.
  - **Extended-12 (`in_Extended12_tier = 1`):** Core-10 + direct pri/pre-let-7 regulators `LIN28A, LIN28B`.
  - **Panel-14 (`in_Panel14_tier = 1`):** Extended-12 + Microprocessor accessory R-SMAD/Co-SMAD factors `SMAD2, SMAD4`.
  - **All-15:** Panel-14 + upstream pleiotropic regulator `TP53` (excluded from biogenesis burden summaries).
- **High-confidence (`strict`) vs broad (`any`) alterations:**
  - `high_conf` / `strict` = non-synonymous somatic coding mutation **OR** high-level focal GISTIC2 CNA (`+2` amplification or `-2` homozygous deep deletion).
  - `any_alter` = non-synonymous somatic coding mutation **OR** any discrete GISTIC2 CNA (`-2, -1, +1, +2`).

---

## 1. Primary Sample-Level, Variant-Level & Mature miRNA-seq Tables

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
| `chrArmGisticStatus` | string | Matched chromosome-arm GISTIC status (`Gain`, `Loss`, `Unchanged`, or blank) from `<studyId>_armlevel_cna`. |
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
**6,292 rows (individual somatic coding variant calls) × 47 columns.**  
Contains `GRCh37` genomic coordinates, `variantType`, `mutationType`, `proteinChange`, `proteinPosStart`, `refseqMrnaId`, read depths (`tumorAltCount`, `tumorRefCount`, `tumorVAF`), `consequenceTier`, `isCuratedHotspot`, `isTruncatingLoF`, `sampleTmbNonsynonymous`, `sampleIsHypermutated`, `molecularSubtype`, plus **21 structural, inverse-folding, protein-language-model & NCBI ClinVar columns**:
- `uniprotAccession`: Canonical reviewed human UniProtKB accession (`Q9NRR4`..`P04637`).
- `afdbResiduePlddt`: Per-residue structural confidence (`0–100`) from AlphaFold DB v6 (`AF-<UniProt>-F1`).
- `pdbComplexId`: Experimental human RCSB PDB Cryo-EM/X-ray structure resolving the mutated codon (`9ASM`, `6V5B`, `3A6P`, `7XW2`, `5ZAK`, `5N8L`, `8ZU6`, `4KRF`, `4W5N`, `5VM9`, `6OON`, `5UDZ`, `4A4I`, `1U7V`, `5MEY`, `1TUP`).
- `distToSubstrateA`: Minimum heavy-atom 3D Euclidean distance (`Å`) from the mutated residue to bound pri-/pre-/mature RNA or DNA substrate (`blank` if unresolved or apo).
- `distToCatalyticMetalA`: Minimum heavy-atom 3D Euclidean distance (`Å`) to bound active-site catalytic/structural divalent metal ions (`Mg2+`, `Ca2+`, `Zn2+`).
- `distToPartnerInterfaceA`: Minimum heavy-atom 3D Euclidean distance (`Å`) to co-complex protein partner chains (`DROSHA–DGCR8`, `DICER1–TARBP2`, `XPO5–RanGTP`, `SMAD2–SMAD4`).
- `is3dPocketOrInterface_le6A`: `1` if within `<= 6.0 Å` of bound nucleic acid substrate, metal ion, or obligate partner protein.
- `alphamissenseScore` & `alphamissenseClass`: DeepMind AlphaMissense pathogenicity score (`0.0000–1.0000`) and tier (`likely_pathogenic`, `ambiguous`, `likely_benign`, or `non_missense`).
- `alphamissenseResidueMean`: Mean AlphaMissense pathogenicity across all 19 single-AA substitutions at that codon position (`0.0000–1.0000`).
- `functionalVariantTierV2`: Integrated 3D + AlphaMissense variant classification (`Truncating_LoF` [`n = 1,716`], `3D_Catalytic_or_Substrate_Pocket_Driver` [`n = 1,251`], `Pathogenic_Fold_or_Interface_Missense` [`n = 2,370`], `Ambiguous_Missense_VUS` [`n = 177`], `Likely_Benign_Passenger_Missense` [`n = 590`], `InFrame_Indel` [`n = 99`], `Other_Nonsynonymous` [`n = 89`]).
- `esm2_3b_llr` & `esm1v_650m_llr`: Zero-shot masked-marginal evolutionary log-likelihood ratios from `esm2_t36_3B_UR50D` (3B parameters) and `esm1v_t33_650M_UR90S_1` (650M parameters).
- `esm_if1_3d_dlogp` & `esm_if1_3d_marginal_llr`: Full-chain autoregressive 3D inverse-folding $\Delta\log P$ and marginal 3D LLR from `esm_if1_gvp4_t16_142M_UR50` conditioned on AlphaFold v6 backbone coordinates (`N, CA, C`).
- `biophysicalMechanismClass`: Mechanistic decomposition into `Fold_Preserving_Catalytic_or_Interface`, `Thermodynamic_Fold_Destabilizing`, or `Tolerated_Neutral`.
- `clinvarVariationId`, `clinvarClinicalSignificance`, `clinvarReviewStars`, `clinvarPhenotypes`: Cross-annotated NCBI ClinVar Variation ID, clinical assertion, review-status star rating (`0–4`), and associated human Mendelian/hereditary cancer phenotypes where the exact protein change has been reported in ClinVar.

### 1.3 `tcga_pancan_miRNA_biogenesis_mirnaseq_functional_readouts.csv`
**9,578 rows (all `miRNA-seq`-matched tumours in the `10,195` cohort; `93.95%` overlap) × 21 columns.**  
Derived from NCI GDC PanCanAtlas batch-corrected mature miRNA-seq (`pancanMiRs`, GDC UUID `1c6174d9-8ffb-466e-b5ee-07b204c15cf8`) and `EBPlusPlusAdjustPANCAN` 20,531-gene mRNA-seq (`N_dual = 9,353` dual-profiled tumours; `13,369,752` transcriptome-wide `miRNA × mRNA` pairs):
- `dicer1FunctionalClass`: Sample *DICER1* somatic status (`DICER1_RNaseIIIb_3D_Catalytic_Pocket`, `DICER1_Truncating_LoF`, `DICER1_Other_Domain_Pathogenic_Missense`, `DICER1_Likely_Benign_Passenger_VUS`, `DICER1_Ambiguous_VUS`, `DICER1_WildType`).
- `mirna_5p_arm_fraction`: Linear RPM ratio $\sum \text{RPM}_{5p} / \sum (\text{RPM}_{5p} + \text{RPM}_{3p})$ across 192 batch-corrected canonical `-5p`/`-3p` duplex pairs.
- `mirna_duplex_5p_3p_log2_ratio`: Mean paired $\log_2((\text{RPM}_{5p}+1)/(\text{RPM}_{3p}+1))$ across expressed canonical duplex pairs ($\text{RPM}_{5p}+\text{RPM}_{3p} \ge 10$; `n_expressed_canonical_duplexes`).
- `let7_family_mean_log2_rpm`, `mir200_family_mean_log2_rpm`, `mir34_family_mean_log2_rpm`: Mean $\log_2(\text{RPM}+1)$ across canonical `let-7` (8 members), `miR-200` (5 members), and `miR-34` (3 members) mature families.
- `mir451a_log2_rpm` & `mir451a_minus_mir144_5p_log2`: Endogenous **DICER1-independent / AGO2-Slicer-dependent** bicistronic processing control (`hsa-miR-451a` vs `hsa-miR-144-5p`).
- `global_canonical_median_log2_rpm`, `n_detected_mature_mirnas_rpm_ge1`, `mirnaome_shannon_entropy_bits`, `mirnaome_pc1_score`: Global mature miRNA output, Shannon diversity (`bits`), and singular-value decomposition principal component 1.
- `let7_target_oncofetal_mrna_z`, `mir200_target_emt_mrna_z`, `empirical_5p_target_derepression_z`, `global_risc_target_derepression_z`: Downstream RISC target mRNA expression Z-scores across `n = 9,353` dual-profiled tumours (`HMGA2, IGF2BP1–3, LIN28B, MYCN, NRAS, CCND2, CDC25A, CDK6` for `let-7`; `ZEB1, ZEB2, FN1, VIM, SNAI2, TWIST1` for `miR-200`).

### 1.4 `miRNA_biogenesis_saturation_mutagenesis_213940_atlas.csv.gz` & `miRNA_biogenesis_residue_biophysical_11260_profile.csv.gz`
- **`miRNA_biogenesis_saturation_mutagenesis_213940_atlas.csv.gz` (213,940 rows × 17 columns; `4.31 MB`):** Complete proteome-scale single-amino-acid saturation mutagenesis lookup across all `11,260` codons × `19` substitutions of the 15 human miRNA biogenesis proteins, combining experimental 3D PDB distances (`distSubstrateA`, `distMetalA`, `distPartnerA`), AlphaFold v6 `afdbPlddt`, `alphamissenseScore`/`Class`, `esm2_3b_llr`, `esm1v_650m_llr`, `esm_if1_3d_llr`, and `observedTcgaSomaticCount`.
- **`miRNA_biogenesis_residue_biophysical_11260_profile.csv.gz` (11,260 rows × 14 columns; `262.9 KB`):** Codon-level biophysical constraint profile (`mean_alphamissense`, `mean_esm2_3b_llr`, `mean_esm1v_650m_llr`, `mean_esm_if1_3d_llr`, 3D PDB distances, and `observedTcgaCodonMutCount`).
- **`tcga_somatic_variant_plm_if1_scores.csv` (1,960 rows × 8 columns)**, **`tcga_pancan_sample_risc_target_derepression.csv` (9,353 rows × 6 columns)**, **`tcga_pancan_sample_mirna_pca_scores.csv` (9,578 rows × 4 columns)**, **`plm_inference_provenance.json`**, and **`tcga_pancan_mirna_mrna_top250_repressive_axes.csv` (250 rows × 6 columns)**.

### 1.5 `tcga_pancan_miRNA_biogenesis_structural_variants.csv`
**72 rows × 14 columns.** Somatic structural variants / fusions involving the 15 panel genes across the 10,195 samples.

### 1.6 `tcga_pancan_miRNA_biogenesis_patient_level_burden.csv`
**10,195 rows (one per tumour sample) × 62 columns.**  
Integrates `cancerTypeDetailed` (resolving `COAD` vs `READ`), `molecularSubtype`, `tmbNonsynonymous`, `mutationCountGenomeWide`, `isHypermutated`, `fractionGenomeAltered`, `aneuploidyScore`, synteny-aware patient union indicators across `Core10`, `Extended12`, `Panel14`, and `All15`, `core10_has_3d_or_am_pathogenic_mut`, syntenic concordance flags (`ago1_3_4_concordant_1p34_cna`, `smad2_4_concordant_18q21_cna`), curated TCGA-CDR clinical survival & staging attributes (`ageAtDiagnosis`, `sex`, `ajccPathologicStage`, `osEvent`, `osMonths`, `dssEvent`, `dssMonths`, `pfsEvent`, `pfsMonths`, `buffaHypoxiaScore`), matched mature `miRNA-seq` summary columns, downstream target mRNA Z-scores (`let7_target_oncofetal_mrna_z`, `mir200_target_emt_mrna_z`, `empirical_5p_target_derepression_z`, `global_risc_target_derepression_z`), and official NCI GDC PanCanAtlas `Thorsson et al. (2018)` immune/instability columns (`immuneSubtypeC1C6`, `leukocyteFraction`, `ifngResponseScore`, `tgfbResponseScore`, `lymphocyteInfiltrationScore`, `macrophageCsf1Score`, `woundHealingScore`, `hrdScore`).

---

## 2. Multi-Consortium Clinical, Germline, Healthy-Tissue, Proteogenomic & CRISPR Tables

- **`clinvar_clinical_variants_15genes_annotated.csv` (7,669 rows × 21 columns):** Unique coding missense and loss-of-function variants across the 15 genes extracted from `16,569` single-gene NCBI ClinVar records, annotated with `clinicalSignificance`, `significanceCategory`, `reviewStars` (`0–4`), `phenotypes`, dbSNP `rsid` and GRCh38 VCF coordinates (`grch38Vcf` `<chrom>:<pos>:<ref>><alt>` where resolved), `tcgaSomaticCount`, 3D PDB distances, `alphamissenseScore`, `esm2_3b_llr`, `esm1v_650m_llr`, and `esm_if1_3d_llr`.
- **`msk_impact_10945_clinical_cohort_replication.csv` (4 rows × 19 columns) & `msk_impact_dicer1_cancer_type_enrichment.csv` (40 rows × 8 columns):** Independent clinical sequencing cohort replication across **`N = 10,945` prospective clinical patients** (`msk_impact_2017`; `DICER1`, `SMAD2`, `SMAD4`, `TP53` profiled across 100% of `IM3` + `IM5` tumours), stratified by non-hypermutator status, truncating LoF, 3D pocket/hotspot status, focal CNA (`+2`/`-2`), primary vs metastatic specimen type, and clinical cancer type.
- **`gtex_v10_normal_54tissues_miRNA_biogenesis_tpm.csv` (810 rows × 5 columns):** Healthy human adult median mRNA expression (`medianTPM`, `gencodeId`, `ontologyId`, `tissueSiteDetailId`) across all `54` GTEx v10 (`GENCODE v39`) non-diseased tissue sites for all 15 panel genes.
- **`cptac_pancan_miRNA_biogenesis_proteogenomic_validation.csv` (15 rows × 8 columns):** Within-cohort standardised Spearman correlations (`spearman_CNA_to_mRNA`, `spearman_CNA_to_Protein`, `spearman_mRNA_to_Protein`, and `post_transcriptional_CNA_attenuation = rho_CNA_mRNA - rho_CNA_Protein`) pooled across **8 CPTAC Pan-Cancer mass-spectrometry cohorts (`N = 1,002` primary tumours)**.
- **`cptac_cohort_complex_stoichiometry.csv` (8 rows × 7 columns):** Per-cohort CPTAC sample counts and heterodimeric co-complex protein correlations (`DROSHA` vs `DGCR8` Microprocessor; `DICER1` vs `TARBP2` RLC).
- **`depmap_crispr_miRNA_biogenesis_essentiality.csv` (15 rows × 12 columns):** Genome-wide CRISPR-Cas9 Chronos `geneEffect` essentiality across **`N = 1,208` human cancer cell lines (`29` tissues)** via Open Targets Platform v4 / Broad DepMap (`median_chronos_gene_effect`, `iqr_chronos_gene_effect`, `% <= -0.5`, `% <= -1.0`, `spearman_expression_vs_chronos`, and `opentargetsIsCommonEssential`).

---

## 3. Gene Panel Metadata, Aggregated Matrices & Validation Audit Tables

- **`miRNA_biogenesis_gene_panel_metadata.csv` (15 rows × 42 columns):** HGNC/Entrez/Ensembl/RefSeq/UniProt IDs, cytobands, syntenic clusters, CDS/protein lengths, tiers, curated hotspots, AlphaFold v6 global pLDDT & domain count, experimental PDB complex IDs & residue coverage, DepMap essentiality flag, gnomAD v4 germline constraint (`gnomad_pLI`, `gnomad_LOEUF`, `gnomad_oe_mis`, `gnomad_mis_z`), Human Protein Atlas subcellular localisation & IHC tissue counts, GTEx v10 pan-54-tissue median TPM & tissue-specificity index (`gtex_tau_tissue_specificity`), Open Targets small-molecule/PROTAC tractability & clinical candidate counts, and ClinVar pathogenic/benign summary counts.
- **`miRNA_biogenesis_excluded_candidates_rationale.csv` (8 rows × 4 columns):** Rationale for excluding pleiotropic/downstream factors (`DDX5`, `DDX17`, `TUT4`, `TUT7`, `DIS3L2`, `TNRC6A/B/C`, `RAN`).
- **`tcga_pancan_miRNA_biogenesis_alteration_freqs_long_FINAL.csv` (480 rows × 35 columns)** & **5 wide (`15 × 33`) matrices (`*_matrix.csv`)** + `miRNA_biogenesis_gene_rank_summary.csv` (15 × 18), `miRNA_biogenesis_top_tumours_excluding_TP53.csv` (32 × 19), `miRNA_biogenesis_top5_tumours_per_gene.csv` (75 × 11), `brca_miRNA_biogenesis_alteration_frequencies.csv` (15 × 24), and `tcga_pancan_atlas_sample_set_summary.csv` (32 × 14).
- **Ten validation benchmark tables (`validation_*.csv`):**
  1. `validation_20_cbioportal_spotchecks.csv` (20 × 17; `0.0000` discrepancy across 20/20 pairs)
  2. `validation_denominator_hoadley2018_audit.csv` (32 × 9)
  3. `validation_cis_mrna_dosage_concordance.csv` (15 × 14; `n = 9,620` TCGA RNA-seq tumours)
  4. `validation_arm_vs_focal_cna_decomposition.csv` (15 × 12)
  5. `validation_mirnaseq_enzymatic_benchmarks.csv` (12 × 14; `B1–B12` across `n = 9,578` miRNA-seq and `n = 9,353` dual `miRNA × mRNA` tumours)
  6. `validation_depmap_collateral_vulnerability.csv` (10 × 15; `n = 1,208` DepMap / `667` CCLE CN-matched cell lines)
  7. `validation_plm_alphamissense_concordance.csv` (16 × 11; `ESM-2 3B`, `ESM-1v`, `ESM-IF1 3D` vs `AlphaMissense` & 3D PDB pocket calibration)
  8. `validation_clinvar_roc_auc_benchmarks.csv` (20 × 10; non-circular human clinical `ClinVar` Pathogenic vs Benign missense ROC-AUC & PR-AUC across `n = 612` assertions; `AlphaMissense` ROC-AUC = `0.9784`)
  9. `validation_clinical_survival_benchmarks.csv` (12 × 23; cancer-type-stratified Mantel–Haenszel log-rank hazard ratios AND exact Breslow-tied multivariable Cox proportional hazards (`multivar_cox_hr`, `multivar_ci95_low`, `multivar_ci95_high`, `multivar_wald_z`, `multivar_p_value`, `multivar_covariates`) adjusting simultaneously for age at diagnosis, biological sex, AJCC pathologic stage, and hypermutator status across 32 TCGA-CDR studies and metastatic/TMB/panel covariates across 43 MSK-IMPACT cancer types)
  10. `validation_panimmune_and_hrd_benchmarks.csv` (8 × 12; `Thorsson` PanImmune `HRD`, `IFN-gamma`, `Leukocyte Fraction`, `TGF-beta`, and `Wound Healing` benchmarks)
- **CLI & Automated Test Suite (`scripts/query_variant.py` & `tests/test_atlas_integrity.py`):** Zero-dependency command-line lookup across all `213,940` saturation variants and `unittest` integrity verification across all `45` SHA-256 archives and exported tables.
- **`data/raw_snapshot/` (`43` `.json.gz` files + `snapshot_manifest.json` tracking `45` SHA-256-verified archives):** Frozen TCGA cBioPortal, 3D PDB/AFDB, AlphaMissense, GDC PanCanAtlas miRNA-seq, PLM/IF1/target-mRNA snapshot, CPTAC, DepMap/CCLE, gnomAD, GTEx v10, HPA, Open Targets, ClinVar, MSK-IMPACT, and Thorsson PanImmune snapshots.
