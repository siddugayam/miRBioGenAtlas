# Publication and technical validation figure generator for miRBioGenAtlas.
# Exports vector PDF, 450-DPI LZW TIFF, and 300-DPI PNG files into figures/.

suppressPackageStartupMessages({
  library(tidyverse)
  library(cowplot)
  library(patchwork)
  library(scales)
  library(ragg)
})

theme_set(theme_cowplot(font_size = 10, font_family = "sans"))

data_dir   <- "data"
output_dir <- "figures"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

save_pub_figure <- function(plot_obj, stem, width_in, height_in) {
  ggsave(
    filename = file.path(output_dir, paste0(stem, ".pdf")),
    plot = plot_obj, width = width_in, height = height_in,
    units = "in", device = pdf, useDingbats = FALSE
  )
  agg_tiff(
    filename = file.path(output_dir, paste0(stem, ".tiff")),
    width = width_in, height = height_in, units = "in",
    res = 450, compression = "lzw", background = "white"
  )
  print(plot_obj)
  invisible(dev.off())
  agg_png(
    filename = file.path(output_dir, paste0(stem, ".png")),
    width = width_in, height = height_in, units = "in",
    res = 300, background = "white"
  )
  print(plot_obj)
  invisible(dev.off())
  cat(sprintf("  Saved %s (.pdf, .tiff [450 DPI], .png [300 DPI])\n", stem))
}


long_df      <- read_csv(file.path(data_dir, "tcga_pancan_miRNA_biogenesis_alteration_freqs_long_FINAL.csv"), show_col_types = FALSE)
gene_meta    <- read_csv(file.path(data_dir, "miRNA_biogenesis_gene_panel_metadata.csv"), show_col_types = FALSE)
gene_rank    <- read_csv(file.path(data_dir, "miRNA_biogenesis_gene_rank_summary.csv"), show_col_types = FALSE)
tumour_rank  <- read_csv(file.path(data_dir, "miRNA_biogenesis_top_tumours_excluding_TP53.csv"), show_col_types = FALSE)
arm_decomp   <- read_csv(file.path(data_dir, "validation_arm_vs_focal_cna_decomposition.csv"), show_col_types = FALSE)
dosage_df    <- read_csv(file.path(data_dir, "validation_cis_mrna_dosage_concordance.csv"), show_col_types = FALSE)
mut_annot    <- read_csv(file.path(data_dir, "tcga_pancan_miRNA_biogenesis_somatic_mutations_annotated.csv"), show_col_types = FALSE)

LOW_N_TUMOURS <- c("CHOL", "DLBC", "UCS", "KICH", "UVM", "MESO", "ACC")
GENE_ORDER <- c("DROSHA", "DGCR8", "XPO5", "DICER1", "TARBP2", "PRKRA",
                "AGO1", "AGO2", "AGO3", "AGO4",
                "LIN28A", "LIN28B", "SMAD2", "SMAD4", "TP53")


cat("Generating Figure 1 (Methods workflow)...\n")

boxes <- tibble(
  xmin = c(0.2, 3.6, 7.0,   0.2, 5.3),
  xmax = c(3.2, 6.6, 10.0,  4.9, 10.0),
  ymin = c(3.3, 3.3, 3.3,   0.3, 0.3),
  ymax = c(5.8, 5.8, 5.8,   2.6, 2.6),
  hdr_fill = c("#0F5499", "#0F5499", "#0F5499", "#1B7340", "#7A3E00"),
  box_fill = c("#F4F8FD", "#F4F8FD", "#F4F8FD", "#F3FAF6", "#FDF8F3"),
  title = c(
    "1. Panel & Cohort Scope",
    "2. Multi-Consortium Frozen Snapshot",
    "3. 3D, Functional & Sample Harmonisation",
    "4. Deposited Data Records & Web Portal (v2.0.1)",
    "5. Multi-Modal Technical Validation Suite"
  ),
  body = c(
    paste(
      "- 15 selected genes in 3 nested tiers:",
      "  Core-10 enzymes/cofactors/AGO1-4,",
      "  Extended-12 (+LIN28A/B),",
      "  Panel-14 (+SMAD2/4), +TP53 context",
      "- 32 TCGA PanCancer Atlas studies",
      "  (33 cancer types; COAD+READ=532)",
      "- N = 10,195 TCGA + 10,945 MSK-IMPACT",
      "  (21,140 combined human tumours)",
      sep = "\n"
    ),
    paste(
      "- SHA-256 raw_snapshot/ (45 files, verified):",
      "  1) TCGA cBioPortal v7.2.3 + TCGA-CDR survival",
      "  2) GDC miRNA-seq (9,578), mRNA (9,353) &",
      "     Thorsson PanImmune / HRD (9,899)",
      "  3) 16 PDB complexes, AFDB v6, AlphaMissense,",
      "     ESM-2 3B, ESM-1v & ESM-IF1 3D",
      "  4) ClinVar (16,569), gnomAD, GTEx (54), HPA,",
      "     OpenTargets, MSK-IMPACT, CPTAC & DepMap",
      sep = "\n"
    ),
    paste(
      "- 152,925 sample x gene multi-omic records",
      "- 213,940-variant saturation atlas (ESM-2 3B,",
      "  ESM-1v, ESM-IF1 3D inverse folding, AM)",
      "- 13.37M mature miRNA x mRNA target pairs",
      "- Mature 5p:3p duplexes (192) & let-7/EMT Z",
      "- Synteny-aware patient union burden +",
      "  survival, immune, germline & tissue context",
      sep = "\n"
    ),
    paste(
      "- Sample-level genomic/mRNA (.csv.gz, 152,925) & miRNA/target Z (9,578)",
      "- 213,940-variant saturation (.csv.gz), 6,292 TCGA & ClinVar/MSK tables",
      "- Patient burden + TCGA-CDR survival & Thorsson immune/HRD (10,195)",
      "- GTEx 54-tissue TPM (810), CPTAC (1,002) & DepMap CRISPR (1,208)",
      "- 480-row master table, 5 matrices, 10 validation_*.csv audit tables",
      "- Offline stdlib builder (~10 s) & docs/index.html + sat/*.js explorer",
      sep = "\n"
    ),
    paste(
      "- 20-pair cBioPortal audit (0.0000 diff) & MSK-IMPACT (N = 10,945) check",
      "- ClinVar human clinical missense ROC-AUC = 0.978 (n = 612 assertions)",
      "- Mature miRNA-seq 5p:3p shift in DICER1 3D pocket (p = 1.5e-15)",
      "- 5p let-7 target mRNA concordance (p = 6.5e-5; LIN28 p = 1.7e-31)",
      "- ESM-IF1 3D fold-preserving pocket vs destabilising core (p = 0.012)",
      "- CPTAC complex buffering, GTEx specificity & DepMap CRISPR (p = 5.4e-4)",
      sep = "\n"
    )
  )
)

arrows_df <- tibble(
  x    = c(3.22, 6.62, 8.50, 4.92),
  xend = c(3.58, 6.98, 8.50, 5.28),
  y    = c(4.55, 4.55, 3.26, 1.45),
  yend = c(4.55, 4.55, 2.66, 1.45)
)

fig1 <- ggplot() +
  geom_rect(data = boxes, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = box_fill),
            colour = "#B8C4D0", linewidth = 0.5) +
  geom_rect(data = boxes, aes(xmin = xmin, xmax = xmax, ymin = ymax - 0.48, ymax = ymax, fill = hdr_fill),
            colour = NA) +
  scale_fill_identity() +
  geom_text(data = boxes, aes(x = xmin + 0.14, y = ymax - 0.24, label = title),
            hjust = 0, vjust = 0.5, colour = "white", fontface = "bold", size = 3.7, family = "sans") +
  geom_text(data = boxes, aes(x = xmin + 0.14, y = ymax - 0.62, label = body),
            hjust = 0, vjust = 1, colour = "#1A2430", size = 3.05, lineheight = 1.18, family = "sans") +
  geom_segment(data = arrows_df, aes(x = x, xend = xend, y = y, yend = yend),
               arrow = arrow(length = unit(0.20, "cm"), type = "closed"),
               colour = "#334E68", linewidth = 0.7) +
  coord_cartesian(xlim = c(0.1, 10.1), ylim = c(0.15, 5.95), expand = FALSE) +
  theme_void() +
  theme(plot.margin = margin(6, 6, 6, 6))

save_pub_figure(fig1, "Figure_1", width_in = 11.5, height_in = 6.2)



cat("Generating Figure 2 (Data Overview)...\n")

# Prepare tumour labels with '*' for N < 100
study_labels <- long_df %>%
  distinct(tumour_acronym, N) %>%
  arrange(tumour_acronym) %>%
  mutate(tumour_label = if_else(tumour_acronym %in% LOW_N_TUMOURS,
                                paste0(tumour_acronym, "*"),
                                tumour_acronym))

label_map <- setNames(study_labels$tumour_label, study_labels$tumour_acronym)

heat_df <- long_df %>%
  mutate(
    tumour_label = factor(label_map[tumour_acronym], levels = study_labels$tumour_label),
    geneSymbol   = factor(geneSymbol, levels = rev(GENE_ORDER)),
    txt_col      = if_else(any_alter_freq >= 0.52, "white", "#102A43")
  )

tier_strip_df <- gene_meta %>%
  mutate(
    geneSymbol = factor(geneSymbol, levels = rev(GENE_ORDER)),
    tier_group = case_when(
      in_Core10_tier == 1     ~ "Core-10 machinery",
      in_Extended12_tier == 1 ~ "Let-7 regulators",
      in_Panel14_tier == 1    ~ "SMAD accessory",
      TRUE                    ~ "Upstream TP53"
    ),
    tier_group = factor(tier_group, levels = c("Core-10 machinery", "Let-7 regulators", "SMAD accessory", "Upstream TP53"))
  )

p2a_strip <- ggplot(tier_strip_df, aes(x = 1, y = geneSymbol, fill = tier_group)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  scale_fill_manual(
    values = c(
      "Core-10 machinery" = "#185FA5",
      "Let-7 regulators"  = "#D4537E",
      "SMAD accessory"    = "#D97706",
      "Upstream TP53"     = "#6B7280"
    ),
    name = "Panel tier"
  ) +
  scale_y_discrete(labels = function(x) parse(text = paste0("italic('", x, "')"))) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 9.5, base_family = "sans") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 9, colour = "#102A43"),
    legend.position = "bottom",
    legend.title = element_text(size = 8.5, face = "bold"),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.35, "cm"),
    plot.margin = margin(4, 0, 4, 4)
  )

p2a_heat <- ggplot(heat_df, aes(x = tumour_label, y = geneSymbol, fill = any_alter_freq)) +
  geom_tile(colour = "grey92", linewidth = 0.25) +
  geom_text(aes(label = sprintf("%.2f", any_alter_freq), colour = txt_col), size = 2.25, family = "sans") +
  scale_colour_identity() +
  scale_fill_gradientn(
    colours = c("#FFFFFF", "#E6F1FB", "#85B7EB", "#185FA5", "#042C53"),
    limits = c(0, 1),
    breaks = c(0, 0.25, 0.50, 0.75, 1.00),
    labels = c("0.00", "0.25", "0.50", "0.75", "1.00"),
    name = "Alteration frequency (0-1)"
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 9.5, base_family = "sans") +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8.2, colour = "#102A43"),
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(size = 8.5, face = "bold"),
    legend.text = element_text(size = 8),
    legend.key.width = unit(1.0, "cm"),
    legend.key.height = unit(0.32, "cm"),
    plot.margin = margin(4, 6, 4, 1)
  )

# Panel b & c: Side-by-side split for the 14 biogenesis/accessory genes (TP53 reported in text)
# Left (b): Mean CNA frequency decomposed into Focal (+/-2) vs Broad low-level (+/-1) on [0.00, 0.50] scale
# Right (c): Mean Somatic Mutation frequency decomposed into Truncating LoF, Hotspot, Other coding on [0.000, 0.028] zoomed scale

gene_order_14 <- gene_rank %>%
  filter(geneSymbol != "TP53") %>%
  arrange(mean_any_cna_freq) %>%
  pull(geneSymbol)

cna_split_df <- gene_rank %>%
  filter(geneSymbol != "TP53") %>%
  mutate(
    geneSymbol = factor(geneSymbol, levels = gene_order_14),
    low_level_cna = pmax(0, mean_any_cna_freq - mean_amp_freq - mean_deep_del_freq)
  ) %>%
  select(geneSymbol, mean_amp_freq, mean_deep_del_freq, low_level_cna) %>%
  pivot_longer(-geneSymbol, names_to = "cna_class", values_to = "freq") %>%
  mutate(cna_class = factor(
    cna_class,
    levels = c("low_level_cna", "mean_deep_del_freq", "mean_amp_freq"),
    labels = c("Low-level gain/loss (GISTIC +/-1)", "Deep deletion (GISTIC -2)", "Focal amplification (GISTIC +2)")
  ))

p2b_cna <- ggplot(cna_split_df, aes(x = freq, y = geneSymbol, fill = cna_class)) +
  geom_col(width = 0.70, position = position_stack(reverse = TRUE)) +
  scale_fill_manual(
    values = c(
      "Focal amplification (GISTIC +2)"  = "#2E7D32",
      "Deep deletion (GISTIC -2)"        = "#D9383A",
      "Low-level gain/loss (GISTIC +/-1)"= "#90CAF9"
    ),
    name = "CNA class"
  ) +
  scale_x_continuous(
    limits = c(0, 0.52),
    breaks = seq(0, 0.50, by = 0.10),
    labels = number_format(accuracy = 0.01),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_y_discrete(labels = function(x) parse(text = paste0("italic('", x, "')"))) +
  labs(
    title = "Copy-number alteration layer (0.00-0.50 scale)",
    x = "Mean CNA frequency across 32 studies (fraction)",
    y = NULL
  ) +
  theme_cowplot(font_size = 9.5) +
  theme(
    plot.title = element_text(size = 9.5, face = "bold"),
    axis.text.y = element_text(size = 9),
    axis.text.x = element_text(size = 8.5),
    axis.title.x = element_text(size = 9),
    panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.35),
    legend.position = "bottom",
    legend.title = element_text(size = 8, face = "bold"),
    legend.text = element_text(size = 7.8),
    legend.key.size = unit(0.32, "cm")
  ) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

mut_split_df <- gene_rank %>%
  filter(geneSymbol != "TP53") %>%
  mutate(
    geneSymbol = factor(geneSymbol, levels = gene_order_14),
    other_coding_mut = pmax(0, mean_mut_freq - mean_truncating_mut_freq - mean_hotspot_mut_freq)
  ) %>%
  select(geneSymbol, mean_truncating_mut_freq, mean_hotspot_mut_freq, other_coding_mut) %>%
  pivot_longer(-geneSymbol, names_to = "mut_class", values_to = "freq") %>%
  mutate(mut_class = factor(
    mut_class,
    levels = c("other_coding_mut", "mean_hotspot_mut_freq", "mean_truncating_mut_freq"),
    labels = c("Other missense / in-frame (VUS)", "Curated domain hotspot", "Truncating LoF (nonsense/frameshift/splice)")
  ))

p2c_mut <- ggplot(mut_split_df, aes(x = freq, y = geneSymbol, fill = mut_class)) +
  geom_col(width = 0.70, position = position_stack(reverse = TRUE)) +
  scale_fill_manual(
    values = c(
      "Truncating LoF (nonsense/frameshift/splice)" = "#8B1E1E",
      "Curated domain hotspot"                      = "#E67E22",
      "Other missense / in-frame (VUS)"             = "#F2B880"
    ),
    name = "Somatic mutation class"
  ) +
  scale_x_continuous(
    limits = c(0, 0.026),
    breaks = seq(0, 0.025, by = 0.005),
    labels = number_format(accuracy = 0.001),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    title = "Somatic mutation layer (zoomed 0.000-0.025 scale)",
    x = "Mean somatic mutation frequency across 32 studies (fraction)",
    y = NULL
  ) +
  theme_cowplot(font_size = 9.5) +
  theme(
    plot.title = element_text(size = 9.5, face = "bold"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 8.5),
    axis.title.x = element_text(size = 9),
    panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.35),
    legend.position = "bottom",
    legend.title = element_text(size = 8, face = "bold"),
    legend.text = element_text(size = 7.8),
    legend.key.size = unit(0.32, "cm")
  ) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

top_row_fig2 <- (p2a_strip + p2a_heat + plot_layout(widths = c(0.045, 0.955), guides = "collect")) &
  theme(legend.position = "bottom", legend.box = "horizontal", legend.margin = margin(2, 10, 2, 10))
bot_row_fig2 <- (p2b_cna | p2c_mut)

fig2 <- plot_grid(
  top_row_fig2,
  bot_row_fig2,
  ncol = 1,
  rel_heights = c(1.12, 0.98),
  labels = c("a", "b"),
  label_size = 12,
  label_fontface = "bold",
  label_fontfamily = "sans"
)

save_pub_figure(fig2, "Figure_2", width_in = 12.0, height_in = 9.8)



cat("Generating Figure 3 (Technical Validation I: CNA & cis-mRNA dosage)...\n")

# Panel a: Chromosome-arm vs focal CNA decomposition across all 10,195 samples
arm_plot_df <- arm_decomp %>%
  mutate(
    gene_arm = paste0(geneSymbol, " (", chrArm, ")"),
    gene_arm = factor(gene_arm, levels = rev(paste0(GENE_ORDER, " (", gene_meta$chrArm[match(GENE_ORDER, gene_meta$geneSymbol)], ")"))),
    Broad_arm = broad_arm_level_lowamp_n / total_samples,
    Subarm_low = (subarm_lowamp_n + lowamp_arm_unannotated_n) / total_samples,
    Focal_high = high_level_amp_or_deepdel_n / total_samples
  ) %>%
  select(gene_arm, Focal_high, Subarm_low, Broad_arm) %>%
  pivot_longer(-gene_arm, names_to = "focality", values_to = "prop") %>%
  mutate(focality = factor(
    focality,
    levels = c("Broad_arm", "Subarm_low", "Focal_high"),
    labels = c("Arm-concordant +/-1 (broad)", "Sub-arm / other +/-1", "High-level focal +/-2 (amp/deep del)")
  ))

p3a <- ggplot(arm_plot_df, aes(x = prop, y = gene_arm, fill = focality)) +
  geom_col(width = 0.72, position = position_stack(reverse = TRUE)) +
  scale_fill_manual(
    values = c(
      "High-level focal +/-2 (amp/deep del)" = "#1B5E20",
      "Sub-arm / other +/-1"                 = "#64B5F6",
      "Arm-concordant +/-1 (broad)"          = "#BBDEFB"
    ),
    name = "CNA focality"
  ) +
  scale_x_continuous(limits = c(0, 0.55), labels = number_format(accuracy = 0.05), expand = expansion(mult = c(0, 0.02))) +
  labs(
    title = "Arm-level vs high-level focal CNA",
    x = "Proportion of pan-cancer cohort (n = 10,195)",
    y = NULL
  ) +
  theme_cowplot(font_size = 9.5) +
  theme(
    plot.title = element_text(size = 10, face = "bold"),
    axis.text.y = element_text(size = 8.5),
    panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.35),
    legend.position = "bottom",
    legend.title = element_text(size = 8, face = "bold"),
    legend.text = element_text(size = 7.6),
    legend.key.size = unit(0.30, "cm")
  ) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

# Panel b: Orthogonal cis-CNA-to-mRNA dosage concordance across n = 9,620 tumours
dosage_long <- dosage_df %>%
  mutate(geneSymbol = factor(geneSymbol, levels = rev(GENE_ORDER))) %>%
  select(geneSymbol, spearman_rho_log2CNA_vs_mrnaZ,
         `-2` = median_z_deep_del_m2,
         `-1` = median_z_shallow_del_m1,
         `0`  = median_z_diploid_0,
         `+1` = median_z_gain_p1,
         `+2` = median_z_amp_p2) %>%
  pivot_longer(c(`-2`, `-1`, `0`, `+1`, `+2`), names_to = "gistic_state", values_to = "median_z") %>%
  mutate(
    gistic_state = factor(gistic_state, levels = c("-2", "-1", "0", "+1", "+2")),
    median_z = as.numeric(median_z)
  )

rho_labels <- dosage_df %>%
  mutate(
    geneSymbol = factor(geneSymbol, levels = rev(GENE_ORDER)),
    rho_lbl = sprintf("rho = %+.2f", spearman_rho_log2CNA_vs_mrnaZ)
  )

p3b <- ggplot(dosage_long, aes(x = gistic_state, y = geneSymbol, fill = median_z)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = ifelse(is.na(median_z), "n<3", sprintf("%+.2f", median_z)),
                colour = ifelse(!is.na(median_z) & abs(median_z) > 1.15, "white", "#102A43")),
            size = 2.6, family = "sans") +
  geom_text(data = rho_labels, aes(x = 5.65, y = geneSymbol, label = rho_lbl),
            inherit.aes = FALSE, hjust = 0, size = 2.75, fontface = "bold", colour = "#102A43", family = "sans") +
  scale_colour_identity() +
  scale_fill_gradient2(
    low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
    midpoint = 0, limits = c(-2.5, 2.5), oob = squish,
    name = "Median mRNA Z-score"
  ) +
  scale_x_discrete(expand = expansion(add = c(0.5, 1.55))) +
  scale_y_discrete(labels = function(x) parse(text = paste0("italic('", x, "')"))) +
  labs(
    title = "Orthogonal cis-CNA vs mRNA dosage (n = 9,620)",
    x = "Discrete GISTIC2 copy-number state",
    y = NULL
  ) +
  theme_cowplot(font_size = 9.5) +
  theme(
    plot.title = element_text(size = 10, face = "bold"),
    axis.text.y = element_text(size = 8.8),
    legend.position = "bottom",
    legend.title = element_text(size = 8, face = "bold"),
    legend.text = element_text(size = 7.8),
    legend.key.width = unit(0.85, "cm"),
    legend.key.height = unit(0.30, "cm")
  )

# Panel c: Broad CNA burden vs genome-wide Fraction Genome Altered (FGA), contrasted with Core-10 strict burden
fga_df <- tumour_rank %>%
  filter(!is.na(median_fraction_genome_altered)) %>%
  mutate(
    disp = if_else(low_coverage_N_lt_100 == 1, paste0(tumour_acronym, "*"), tumour_acronym),
    highlight = tumour_acronym %in% c("TGCT", "OV", "UCS", "LAML", "UCEC", "SKCM", "BRCA", "LUSC")
  )

tgct_row <- filter(fga_df, tumour_acronym == "TGCT")

p3c <- ggplot(fga_df, aes(x = median_fraction_genome_altered)) +
  geom_smooth(aes(y = patient_prop_core10_any_alter), method = "lm", formula = y ~ x, se = FALSE, colour = "#90CAF9", linetype = "dashed", linewidth = 0.7) +
  geom_smooth(aes(y = patient_prop_core10_strict), method = "lm", formula = y ~ x, se = FALSE, colour = "#1B5E20", linetype = "solid", linewidth = 0.7) +
  geom_segment(aes(xend = median_fraction_genome_altered, y = patient_prop_core10_strict, yend = patient_prop_core10_any_alter),
               colour = "grey80", linewidth = 0.45) +
  geom_point(aes(y = patient_prop_core10_any_alter, colour = "Broad Core-10 union (mut OR GISTIC != 0)"), size = 2.3) +
  geom_point(aes(y = patient_prop_core10_strict, colour = "Strict Core-10 union (mut OR GISTIC +/-2)"), size = 2.5) +
  geom_text(
    data = filter(fga_df, highlight),
    aes(y = patient_prop_core10_strict, label = disp),
    nudge_y = -0.048, size = 2.65, fontface = "bold", colour = "#102A43"
  ) +
  annotate("segment", x = 0.455, xend = tgct_row$median_fraction_genome_altered + 0.005,
           y = 0.80, yend = tgct_row$patient_prop_core10_any_alter - 0.015,
           colour = "#0F5499", linewidth = 0.45) +
  annotate("text", x = 0.46, y = 0.79, label = "TGCT: 0.99 broad vs 0.04 strict",
           hjust = 0, size = 2.75, colour = "#0F5499", fontface = "italic") +
  scale_colour_manual(
    values = c(
      "Broad Core-10 union (mut OR GISTIC != 0)" = "#42A5F5",
      "Strict Core-10 union (mut OR GISTIC +/-2)" = "#1B5E20"
    ),
    name = "Patient-level metric"
  ) +
  scale_x_continuous(breaks = seq(0, 0.6, 0.2), labels = number_format(accuracy = 0.1)) +
  scale_y_continuous(breaks = seq(0, 1, 0.25), labels = number_format(accuracy = 0.01)) +
  coord_cartesian(xlim = c(-0.02, 0.63), ylim = c(-0.07, 1.04)) +
  labs(
    title = "Genome-wide FGA vs patient-level Core-10 burden",
    x = "Study median Fraction Genome Altered (FGA)",
    y = "Proportion of patients altered (0-1)"
  ) +
  theme_cowplot(font_size = 9.5) +
  theme(
    plot.title = element_text(size = 10, face = "bold"),
    panel.grid.major = element_line(colour = "grey92", linewidth = 0.35),
    legend.position = "bottom",
    legend.title = element_text(size = 8, face = "bold"),
    legend.text = element_text(size = 7.8)
  ) +
  guides(colour = guide_legend(nrow = 1))

fig3 <- plot_grid(
  plot_grid(p3a, p3b, ncol = 2, rel_widths = c(0.95, 1.05), labels = c("a", "b"), label_size = 12),
  p3c,
  ncol = 1,
  rel_heights = c(1.15, 0.88),
  labels = c("", "c"),
  label_size = 12
)

save_pub_figure(fig3, "Figure_3", width_in = 12.0, height_in = 9.5)



cat("Generating Figure 4 (Technical Validation II: Mutations & nested tiers)...\n")

# Panel a: Full cohort vs Non-hypermutated mutation frequency + CDS-length-normalized density (14 genes excl TP53)
mut_cal_order <- gene_rank %>%
  filter(geneSymbol != "TP53") %>%
  arrange(mean_mut_freq) %>%
  pull(geneSymbol)

mut_cal_df <- gene_rank %>%
  filter(geneSymbol != "TP53") %>%
  mutate(
    gene_cds_lbl = sprintf("%s (%.1f kb)", geneSymbol, cdsLengthBp / 1000),
    geneSymbol = factor(geneSymbol, levels = mut_cal_order)
  )

lvl_cds <- mut_cal_df$gene_cds_lbl[match(mut_cal_order, mut_cal_df$geneSymbol)]
mut_cal_df$gene_cds_lbl <- factor(mut_cal_df$gene_cds_lbl, levels = lvl_cds)

mut_cal_long <- mut_cal_df %>%
  select(gene_cds_lbl,
         `All samples (n = 10,195)` = mean_mut_freq,
         `Non-hypermutated only (n = 8,810)` = mean_nonhypermut_mut_freq) %>%
  pivot_longer(-gene_cds_lbl, names_to = "cohort_subset", values_to = "freq")

p4a_left <- ggplot(mut_cal_long, aes(x = freq, y = gene_cds_lbl, fill = cohort_subset)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.65) +
  scale_fill_manual(
    values = c(
      "All samples (n = 10,195)"          = "#C0392B",
      "Non-hypermutated only (n = 8,810)" = "#2C3E50"
    ),
    name = "Sample stratification"
  ) +
  scale_x_continuous(limits = c(0, 0.026), labels = number_format(accuracy = 0.005), expand = expansion(mult = c(0, 0.03))) +
  labs(
    title = "Hypermutator stratification",
    x = "Mean mutation frequency (0-1)",
    y = NULL
  ) +
  theme_cowplot(font_size = 9) +
  theme(
    plot.title = element_text(size = 9.5, face = "bold"),
    axis.text.y = element_text(size = 8.3),
    panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.35),
    legend.position = "bottom",
    legend.title = element_text(size = 7.8, face = "bold"),
    legend.text = element_text(size = 7.5),
    legend.key.size = unit(0.30, "cm")
  ) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

p4a_right <- ggplot(mut_cal_df, aes(x = mean_mut_density_per_kb_per_1000, y = gene_cds_lbl, fill = panelRole)) +
  geom_vline(xintercept = median(mut_cal_df$mean_mut_density_per_kb_per_1000[mut_cal_df$panelRole != "Pleiotropic_Accessory"]),
             linetype = "dashed", colour = "grey45", linewidth = 0.45) +
  geom_col(width = 0.65) +
  scale_fill_manual(
    values = c(
      "Core_Enzyme"                 = "#185FA5",
      "Core_Cofactor"               = "#42A5F5",
      "Core_Transporter"            = "#26A69A",
      "Core_Effector"               = "#5C6BC0",
      "Direct_Biogenesis_Regulator" = "#D4537E",
      "Pleiotropic_Accessory"       = "#D97706"
    ),
    labels = c(
      "Core_Enzyme"                 = "Core enzyme",
      "Core_Cofactor"               = "Core cofactor",
      "Core_Transporter"            = "Nuclear export",
      "Core_Effector"               = "RISC effector",
      "Direct_Biogenesis_Regulator" = "Let-7 regulator",
      "Pleiotropic_Accessory"       = "SMAD accessory"
    ),
    name = "Role"
  ) +
  scale_x_continuous(limits = c(0, 16), expand = expansion(mult = c(0, 0.03))) +
  labs(
    title = "CDS-length normalized density",
    x = "Mutations / kb CDS / 1,000 samples",
    y = NULL
  ) +
  theme_cowplot(font_size = 9) +
  theme(
    plot.title = element_text(size = 9.5, face = "bold"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.35),
    legend.position = "bottom",
    legend.title = element_text(size = 7.8, face = "bold"),
    legend.text = element_text(size = 7.3),
    legend.key.size = unit(0.28, "cm"),
    plot.margin = margin(6, 12, 6, 2)
  ) +
  guides(fill = guide_legend(ncol = 2, byrow = TRUE))

# Panel b: Study-level patient-level high-confidence alteration burden across Core-10 vs Extended-12 vs Panel-14
tb_plot_df <- tumour_rank %>%
  arrange(patient_prop_core10_strict) %>%
  mutate(
    tumour_lbl = if_else(low_coverage_N_lt_100 == 1, paste0(tumour_acronym, "*"), tumour_acronym),
    tumour_lbl = factor(tumour_lbl, levels = tumour_lbl)
  )

tb_long <- tb_plot_df %>%
  select(
    tumour_lbl,
    `Core-10`       = patient_prop_core10_strict,
    `Extended-12`   = patient_prop_extended12_strict,
    `Panel-14`      = patient_prop_panel14_strict
  ) %>%
  pivot_longer(-tumour_lbl, names_to = "tier", values_to = "prop") %>%
  mutate(tier = factor(tier, levels = c("Core-10", "Extended-12", "Panel-14")))

p4b <- ggplot(tb_long, aes(x = prop, y = tumour_lbl, colour = tier, shape = tier)) +
  geom_line(aes(group = tumour_lbl), colour = "grey80", linewidth = 0.5) +
  geom_point(size = 2.35, alpha = 0.92) +
  scale_colour_manual(
    values = c(
      "Core-10"     = "#0F5499",
      "Extended-12" = "#D4537E",
      "Panel-14"    = "#D97706"
    ),
    name = "Nested tier"
  ) +
  scale_shape_manual(
    values = c(16, 17, 15),
    name = "Nested tier"
  ) +
  scale_x_continuous(limits = c(0, 0.56), breaks = seq(0, 0.50, 0.10), labels = number_format(accuracy = 0.01), expand = expansion(mult = c(0, 0.03))) +
  labs(
    title = "Synteny-aware patient burden (* N < 100)",
    x = "Patients with >=1 mutation or GISTIC2 +/-2 (0-1)",
    y = NULL
  ) +
  theme_cowplot(font_size = 9) +
  theme(
    plot.title = element_text(size = 9.5, face = "bold"),
    axis.text.y = element_text(size = 8.0),
    panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.35),
    panel.grid.major.y = element_line(colour = "grey96", linewidth = 0.25),
    legend.position = "bottom",
    legend.title = element_text(size = 8, face = "bold"),
    legend.text = element_text(size = 7.8),
    plot.margin = margin(6, 10, 6, 10)
  )

left_col_fig4 <- (p4a_left | p4a_right) + plot_layout(widths = c(1.08, 0.92))

fig4 <- plot_grid(
  left_col_fig4,
  p4b,
  ncol = 2,
  rel_widths = c(1.08, 0.92),
  labels = c("a", "b"),
  label_size = 12
)

save_pub_figure(fig4, "Figure_4", width_in = 12.6, height_in = 7.4)



cat("Generating Figure 5 (Technical Validation III: miRNA-seq, 3D/AlphaMissense, CPTAC & DepMap)...\n")

mir_df   <- read_csv(file.path(data_dir, "tcga_pancan_miRNA_biogenesis_mirnaseq_functional_readouts.csv"), show_col_types = FALSE)
mut_df   <- read_csv(file.path(data_dir, "tcga_pancan_miRNA_biogenesis_somatic_mutations_annotated.csv"), show_col_types = FALSE)
cptac_df <- read_csv(file.path(data_dir, "cptac_pancan_miRNA_biogenesis_proteogenomic_validation.csv"), show_col_types = FALSE)
depmap_df <- read_csv(file.path(data_dir, "depmap_crispr_miRNA_biogenesis_essentiality.csv"), show_col_types = FALSE)

# Panel a: Matched TCGA mature miRNA-seq 5p:3p duplex ratio across DICER1 3D/AlphaMissense classes
mir_dicer <- mir_df %>%
  filter(!is.na(mirna_duplex_5p_3p_log2_ratio)) %>%
  filter(dicer1FunctionalClass %in% c(
    "DICER1_RNaseIIIb_3D_Catalytic_Pocket",
    "DICER1_Truncating_LoF",
    "DICER1_Likely_Benign_Passenger_VUS",
    "DICER1_WildType"
  )) %>%
  mutate(
    cls_lbl = case_when(
      dicer1FunctionalClass == "DICER1_RNaseIIIb_3D_Catalytic_Pocket" ~ "RNase IIIb 3D pocket\n(<=8 A Mg2+/RNA, n = 25)",
      dicer1FunctionalClass == "DICER1_Truncating_LoF"                ~ "Truncating LoF\n(monoallelic, n = 50)",
      dicer1FunctionalClass == "DICER1_Likely_Benign_Passenger_VUS"   ~ "Likely benign VUS\n(AlphaMissense, n = 59)",
      TRUE                                                            ~ "DICER1 wild-type\n(n = 9,332)"
    )
  )

cnt_tbl <- mir_dicer %>% count(dicer1FunctionalClass)
n_r3b <- cnt_tbl$n[cnt_tbl$dicer1FunctionalClass == "DICER1_RNaseIIIb_3D_Catalytic_Pocket"]
n_trc <- cnt_tbl$n[cnt_tbl$dicer1FunctionalClass == "DICER1_Truncating_LoF"]
n_ben <- cnt_tbl$n[cnt_tbl$dicer1FunctionalClass == "DICER1_Likely_Benign_Passenger_VUS"]
n_wt  <- cnt_tbl$n[cnt_tbl$dicer1FunctionalClass == "DICER1_WildType"]

lvl_dicer <- c(
  sprintf("3D dicing pocket\n(<=8 A, n = %d)", n_r3b),
  sprintf("Truncating LoF\n(n = %d)", n_trc),
  sprintf("Likely benign VUS\n(n = %d)", n_ben),
  sprintf("DICER1 wild-type\n(n = %s)", comma(n_wt))
)

mir_dicer <- mir_dicer %>%
  mutate(
    cls_lbl = case_when(
      dicer1FunctionalClass == "DICER1_RNaseIIIb_3D_Catalytic_Pocket" ~ lvl_dicer[1],
      dicer1FunctionalClass == "DICER1_Truncating_LoF"                ~ lvl_dicer[2],
      dicer1FunctionalClass == "DICER1_Likely_Benign_Passenger_VUS"   ~ lvl_dicer[3],
      TRUE                                                            ~ lvl_dicer[4]
    ),
    cls_lbl = factor(cls_lbl, levels = rev(lvl_dicer))
  )

p5a <- ggplot(mir_dicer, aes(x = mirna_duplex_5p_3p_log2_ratio, y = cls_lbl, fill = cls_lbl)) +
  geom_vline(xintercept = median(mir_dicer$mirna_duplex_5p_3p_log2_ratio[mir_dicer$dicer1FunctionalClass == "DICER1_WildType"]),
             linetype = "dashed", colour = "grey45", linewidth = 0.45) +
  geom_boxplot(width = 0.52, outlier.size = 0.55, outlier.alpha = 0.25, alpha = 0.88, colour = "grey20", linewidth = 0.4) +
  geom_jitter(
    data = mir_dicer %>% filter(dicer1FunctionalClass != "DICER1_WildType"),
    height = 0.16, size = 1.35, alpha = 0.78, colour = "#111827"
  ) +
  scale_fill_manual(
    values = setNames(c("#DC2626", "#D97706", "#64748B", "#CBD5E1"), lvl_dicer),
    guide = "none"
  ) +
  scale_x_continuous(breaks = seq(-1.5, 1.0, 0.5)) +
  coord_cartesian(xlim = c(-1.55, 1.25)) +
  labs(
    title = "TCGA mature miRNA-seq 5p:3p strand balance (n = 9,578)",
    subtitle = "3D pocket vs WT: p = 1.5e-15 | vs Benign VUS: p = 2.5e-11 | Benign VUS vs WT: p = 0.35",
    x = "Mean paired duplex log2((5p + 1) / (3p + 1)) across 192 canonical pairs",
    y = NULL
  ) +
  theme_cowplot(font_size = 9) +
  theme(
    plot.title = element_text(size = 9.5, face = "bold"),
    plot.subtitle = element_text(size = 7.6, colour = "#991B1B", face = "bold"),
    axis.text.y = element_text(size = 8.0),
    panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.35),
    plot.margin = margin(6, 10, 6, 6)
  )

# Panel b: 3D Cryo-EM + AlphaMissense resolution of somatic variants across the 14 non-TP53 genes
v2_levels <- c(
  "Truncating_LoF",
  "3D_Catalytic_or_Substrate_Pocket_Driver",
  "Pathogenic_Fold_or_Interface_Missense",
  "Ambiguous_Missense_VUS",
  "Likely_Benign_Passenger_Missense",
  "Other"
)

mut_non_tp53 <- mut_df %>%
  filter(geneSymbol != "TP53") %>%
  mutate(
    tier_plot = if_else(functionalVariantTierV2 %in% v2_levels[1:5], functionalVariantTierV2, "Other"),
    tier_plot = factor(tier_plot, levels = rev(v2_levels))
  )

gene_mut_totals <- mut_non_tp53 %>% count(geneSymbol, name = "tot") %>% arrange(tot)

mut_prop_df <- mut_non_tp53 %>%
  count(geneSymbol, tier_plot) %>%
  left_join(gene_mut_totals, by = "geneSymbol") %>%
  mutate(
    prop = n / tot,
    gene_lbl = factor(sprintf("%s (n = %d)", geneSymbol, tot),
                      levels = sprintf("%s (n = %d)", gene_mut_totals$geneSymbol, gene_mut_totals$tot))
  )

p5b <- ggplot(mut_prop_df, aes(x = prop, y = gene_lbl, fill = tier_plot)) +
  geom_col(width = 0.68, colour = "white", linewidth = 0.2) +
  scale_fill_manual(
    values = c(
      "Truncating_LoF"                          = "#1E293B",
      "3D_Catalytic_or_Substrate_Pocket_Driver" = "#DC2626",
      "Pathogenic_Fold_or_Interface_Missense"   = "#EA580C",
      "Ambiguous_Missense_VUS"                  = "#FBBF24",
      "Likely_Benign_Passenger_Missense"        = "#60A5FA",
      "Other"                                   = "#CBD5E1"
    ),
    labels = c(
      "Truncating_LoF"                          = "Truncating LoF",
      "3D_Catalytic_or_Substrate_Pocket_Driver" = "3D pocket <=6 A (PDB)",
      "Pathogenic_Fold_or_Interface_Missense"   = "Pathogenic fold (AM)",
      "Ambiguous_Missense_VUS"                  = "Ambiguous VUS",
      "Likely_Benign_Passenger_Missense"        = "Likely benign VUS",
      "Other"                                   = "In-frame / splice"
    ),
    breaks = v2_levels,
    name = NULL
  ) +
  scale_x_continuous(limits = c(0, 1.001), labels = number_format(accuracy = 0.25), expand = expansion(mult = c(0, 0.02))) +
  labs(
    title = "3D Cryo-EM / PDB & AlphaMissense variant calibration",
    x = "Proportion of non-synonymous somatic mutations (0-1)",
    y = NULL
  ) +
  theme_cowplot(font_size = 9) +
  theme(
    plot.title = element_text(size = 9.5, face = "bold"),
    axis.text.y = element_text(size = 7.9),
    legend.position = "bottom",
    legend.text = element_text(size = 7.2),
    legend.key.size = unit(0.28, "cm"),
    plot.margin = margin(6, 10, 6, 6)
  ) +
  guides(fill = guide_legend(ncol = 3, byrow = TRUE))

# Panel c: Sequence conservation (ESM-2 3B) vs 3D fold stability (ESM-IF1) across 1,960 unique TCGA missense variants
plm_uniq <- mut_df %>%
  filter(!is.na(esm2_3b_llr), !is.na(esm_if1_3d_dlogp), !is.na(biophysicalMechanismClass)) %>%
  distinct(geneSymbol, proteinChange, .keep_all = TRUE) %>%
  mutate(
    mech_lbl = case_when(
      biophysicalMechanismClass == "Fold_Preserving_Catalytic_or_Interface" ~ "Fold-preserving active/interface (n = 471)",
      biophysicalMechanismClass == "Thermodynamic_Fold_Destabilizing"       ~ "Predicted fold-destabilising (n = 1,041)",
      TRUE                                                                  ~ "Tolerated neutral (n = 448)"
    )
  )

exemplars <- tibble(
  geneSymbol = c("DICER1", "DICER1", "DICER1", "TP53", "SMAD4", "TP53", "TP53"),
  proteinChange = c("E1705K", "D1709N", "E1813G", "R248Q", "R361H", "R175H", "V157F"),
  lbl = c("DICER1 E1705K", "DICER1 D1709N", "DICER1 E1813G", "TP53 R248Q", "SMAD4 R361H", "TP53 R175H", "TP53 V157F"),
  lx = c(-14.1, -13.4, -14.2, -11.4, -4.3, -10.2, -11.2),
  ly = c(  0.9,   2.5,  -2.3,  -4.2,  2.1,  -9.5, -11.9)
) %>%
  inner_join(plm_uniq %>% select(geneSymbol, proteinChange, esm2_3b_llr, esm_if1_3d_dlogp), by = c("geneSymbol", "proteinChange"))

p5c <- ggplot(plm_uniq, aes(x = esm2_3b_llr, y = esm_if1_3d_dlogp, colour = mech_lbl)) +
  geom_hline(yintercept = -1.5, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  geom_vline(xintercept = -5.0, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  geom_point(size = 1.15, alpha = 0.62) +
  geom_segment(
    data = exemplars,
    aes(x = lx, y = ly, xend = esm2_3b_llr, yend = esm_if1_3d_dlogp),
    inherit.aes = FALSE, colour = "#111827", linewidth = 0.35
  ) +
  geom_point(
    data = exemplars,
    aes(x = esm2_3b_llr, y = esm_if1_3d_dlogp),
    inherit.aes = FALSE, shape = 21, fill = "#FACC15", colour = "#111827", size = 2.3, stroke = 0.55
  ) +
  geom_label(
    data = exemplars,
    aes(x = lx, y = ly, label = lbl),
    inherit.aes = FALSE, size = 2.25, fontface = "bold", colour = "#0F172A",
    fill = alpha("white", 0.88), label.padding = unit(0.10, "lines"), linewidth = 0.2
  ) +
  scale_colour_manual(
    values = c(
      "Fold-preserving active/interface (n = 471)" = "#DC2626",
      "Predicted fold-destabilising (n = 1,041)"   = "#D97706",
      "Tolerated neutral (n = 448)"                = "#3B82F6"
    ),
    name = NULL
  ) +
  coord_cartesian(xlim = c(-16.2, 2.5), ylim = c(-13.8, 3.2)) +
  labs(
    title = "Sequence conservation (ESM-2 3B) vs 3D stability (ESM-IF1)",
    subtitle = "1,960 unique TCGA missense variants | 3D pocket vs core fold: p = 0.012",
    x = "Evolutionary constraint: ESM-2 (3B) masked marginal LLR",
    y = "3D backbone stability: ESM-IF1 delta logP"
  ) +
  theme_cowplot(font_size = 9) +
  theme(
    plot.title = element_text(size = 9.5, face = "bold"),
    plot.subtitle = element_text(size = 7.5, colour = "#991B1B", face = "bold"),
    panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
    legend.position = "bottom",
    legend.text = element_text(size = 7.0),
    plot.margin = margin(6, 10, 6, 6)
  ) +
  guides(colour = guide_legend(ncol = 2, byrow = TRUE, override.aes = list(size = 2.4, alpha = 1)))

# Panel d: Downstream miRNA target mRNA concordance across n = 9,353 dual-profiled TCGA tumours
sl_lin28 <- read_csv(file.path(data_dir, "tcga_pancan_miRNA_biogenesis_sample_level_alterations.csv.gz"), show_col_types = FALSE) %>%
  filter(geneSymbol %in% c("LIN28A", "LIN28B"), !is.na(mrnaZscore)) %>%
  group_by(sampleId) %>%
  summarise(max_lin28_z = max(mrnaZscore, na.rm = TRUE), .groups = "drop")

derep_base <- mir_df %>%
  filter(!is.na(let7_target_oncofetal_mrna_z), !is.na(mir200_target_emt_mrna_z)) %>%
  left_join(sl_lin28, by = "sampleId")

q25_m200 <- quantile(derep_base$mir200_family_mean_log2_rpm, 0.25, na.rm = TRUE)
q75_m200 <- quantile(derep_base$mir200_family_mean_log2_rpm, 0.75, na.rm = TRUE)

derep_contrasts <- bind_rows(
  derep_base %>%
    filter(dicer1FunctionalClass %in% c("DICER1_RNaseIIIb_3D_Catalytic_Pocket", "DICER1_WildType")) %>%
    transmute(
      axis = "DICER1 3D pocket ->\nlet-7 (5p) targets (p = 6.5e-5)",
      state = if_else(dicer1FunctionalClass == "DICER1_RNaseIIIb_3D_Catalytic_Pocket", "Altered / low miRNA", "Wild-type / control"),
      z = let7_target_oncofetal_mrna_z
    ),
  derep_base %>%
    filter(dicer1FunctionalClass %in% c("DICER1_RNaseIIIb_3D_Catalytic_Pocket", "DICER1_WildType")) %>%
    transmute(
      axis = "DICER1 3D pocket ->\nmiR-200 (3p) control (p = 0.34)",
      state = if_else(dicer1FunctionalClass == "DICER1_RNaseIIIb_3D_Catalytic_Pocket", "Altered / low miRNA", "Wild-type / control"),
      z = mir200_target_emt_mrna_z
    ),
  derep_base %>%
    filter(!is.na(max_lin28_z), max_lin28_z >= 2.0 | max_lin28_z <= 0.0) %>%
    transmute(
      axis = "LIN28A/B Z >= 2 ->\nlet-7 oncofetal targets (p = 1.7e-31)",
      state = if_else(max_lin28_z >= 2.0, "Altered / low miRNA", "Wild-type / control"),
      z = let7_target_oncofetal_mrna_z
    ),
  derep_base %>%
    filter(!is.na(mir200_family_mean_log2_rpm), mir200_family_mean_log2_rpm <= q25_m200 | mir200_family_mean_log2_rpm >= q75_m200) %>%
    transmute(
      axis = "miR-200 Q1 vs Q4 ->\nEMT target mRNAs (p = 1.6e-176)",
      state = if_else(mir200_family_mean_log2_rpm <= q25_m200, "Altered / low miRNA", "Wild-type / control"),
      z = mir200_target_emt_mrna_z
    )
) %>%
  mutate(
    axis = factor(axis, levels = rev(c(
      "DICER1 3D pocket ->\nlet-7 (5p) targets (p = 6.5e-5)",
      "DICER1 3D pocket ->\nmiR-200 (3p) control (p = 0.34)",
      "LIN28A/B Z >= 2 ->\nlet-7 oncofetal targets (p = 1.7e-31)",
      "miR-200 Q1 vs Q4 ->\nEMT target mRNAs (p = 1.6e-176)"
    ))),
    state = factor(state, levels = c("Wild-type / control", "Altered / low miRNA"))
  )

p5d <- ggplot(derep_contrasts, aes(x = z, y = axis, fill = state)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey45", linewidth = 0.45) +
  geom_boxplot(width = 0.62, position = position_dodge(width = 0.72), outlier.shape = NA, alpha = 0.90, colour = "grey20", linewidth = 0.38) +
  scale_fill_manual(
    values = c("Altered / low miRNA" = "#DC2626", "Wild-type / control" = "#CBD5E1"),
    breaks = c("Altered / low miRNA", "Wild-type / control"),
    name = NULL
  ) +
  coord_cartesian(xlim = c(-1.55, 2.05)) +
  labs(
    title = "Downstream miRNA target mRNA concordance (n = 9,353)",
    subtitle = "5p let-7 targets vs 3p miR-200 control across 13.37M miRNA-mRNA pairs",
    x = "Mean target mRNA expression Z-score",
    y = NULL
  ) +
  theme_cowplot(font_size = 9) +
  theme(
    plot.title = element_text(size = 9.5, face = "bold"),
    plot.subtitle = element_text(size = 7.5, colour = "#1E3A8A", face = "italic"),
    axis.text.y = element_text(size = 7.8),
    panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.35),
    legend.position = "bottom",
    legend.text = element_text(size = 7.4),
    legend.key.size = unit(0.30, "cm"),
    plot.margin = margin(6, 10, 6, 6)
  )

# Panel e: CPTAC Pan-Cancer mass-spec proteomics CNA -> mRNA vs CNA -> Protein dosage buffering (n = 1,002)
cptac_plot <- cptac_df %>%
  filter(n_cptac_mrna_protein_paired >= 100) %>%
  arrange(spearman_CNA_to_Protein) %>%
  mutate(geneSymbol = factor(geneSymbol, levels = geneSymbol))

cptac_long <- cptac_plot %>%
  select(geneSymbol,
         `CNA -> mRNA`    = spearman_CNA_to_mRNA,
         `CNA -> Protein` = spearman_CNA_to_Protein) %>%
  pivot_longer(-geneSymbol, names_to = "layer", values_to = "rho")

p5e <- ggplot(cptac_long, aes(x = rho, y = geneSymbol, colour = layer, shape = layer)) +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey55", linewidth = 0.4) +
  geom_line(aes(group = geneSymbol), colour = "grey72", linewidth = 0.75) +
  geom_point(size = 2.6, alpha = 0.95) +
  scale_colour_manual(
    values = c("CNA -> mRNA" = "#0284C7", "CNA -> Protein" = "#BE185D"),
    name = "CPTAC dosage coupling (8 cohorts)"
  ) +
  scale_shape_manual(
    values = c("CNA -> mRNA" = 16, "CNA -> Protein" = 17),
    name = "CPTAC dosage coupling (8 cohorts)"
  ) +
  scale_x_continuous(limits = c(-0.14, 0.52), breaks = seq(-0.1, 0.5, 0.1), labels = number_format(accuracy = 0.1)) +
  labs(
    title = "CPTAC proteogenomic dosage attenuation (n = 1,002)",
    subtitle = "Co-complex protein rho: DICER1-TARBP2 +0.47 (n = 655) | DROSHA-DGCR8 +0.37 (n = 598)",
    x = "Within-cohort Spearman rho",
    y = NULL
  ) +
  theme_cowplot(font_size = 9) +
  theme(
    plot.title = element_text(size = 9.5, face = "bold"),
    plot.subtitle = element_text(size = 7.6, colour = "#1E3A8A", face = "italic"),
    axis.text.y = element_text(size = 8.0),
    panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.35),
    panel.grid.major.y = element_line(colour = "grey96", linewidth = 0.25),
    legend.position = "bottom",
    legend.title = element_text(size = 7.6, face = "bold"),
    legend.text = element_text(size = 7.4),
    plot.margin = margin(6, 10, 6, 6)
  )

# Panel f: Broad DepMap CRISPR-Cas9 Chronos essentiality across 1,208 human cancer cell lines
depmap_plot <- depmap_df %>%
  arrange(desc(median_chronos_gene_effect)) %>%
  mutate(
    gene_dep_lbl = sprintf("%s (%.0f%%)", geneSymbol, prop_dependent_le_m05 * 100),
    gene_dep_lbl = factor(gene_dep_lbl, levels = gene_dep_lbl),
    ess_group = case_when(
      geneSymbol == "TP53"                  ~ "Growth-suppressive TSG",
      opentargetsIsCommonEssential == 1     ~ "Pan-essential trunk enzyme/cofactor",
      TRUE                                  ~ "Paralogue-buffered / context-specific"
    ),
    q_low  = median_chronos_gene_effect - iqr_chronos_gene_effect / 2,
    q_high = median_chronos_gene_effect + iqr_chronos_gene_effect / 2
  )

p5f <- ggplot(depmap_plot, aes(x = median_chronos_gene_effect, y = gene_dep_lbl, colour = ess_group)) +
  geom_vline(xintercept = -0.5, linetype = "dashed", colour = "#DC2626", linewidth = 0.45) +
  geom_vline(xintercept = 0.0, linetype = "dotted", colour = "grey50", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = q_low, xmax = q_high), height = 0.28, linewidth = 0.65, alpha = 0.85) +
  geom_point(size = 2.5) +
  scale_colour_manual(
    values = c(
      "Pan-essential trunk enzyme/cofactor"   = "#DC2626",
      "Paralogue-buffered / context-specific" = "#2563EB",
      "Growth-suppressive TSG"                = "#059669"
    ),
    name = NULL
  ) +
  scale_x_continuous(limits = c(-1.05, 0.60), breaks = seq(-1.0, 0.5, 0.25)) +
  labs(
    title = "DepMap CRISPR-Cas9 essentiality (n = 1,208 cell lines)",
    x = "Chronos gene effect (median +/- IQR/2; label = % lines <= -0.5)",
    y = NULL
  ) +
  theme_cowplot(font_size = 9) +
  theme(
    plot.title = element_text(size = 9.5, face = "bold"),
    axis.text.y = element_text(size = 7.9),
    panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.35),
    legend.position = "bottom",
    legend.text = element_text(size = 7.1),
    plot.margin = margin(6, 10, 6, 6)
  ) +
  guides(colour = guide_legend(ncol = 2, byrow = TRUE))

fig5 <- plot_grid(
  p5a, p5b, p5c, p5d, p5e, p5f,
  ncol = 2, nrow = 3,
  labels = c("a", "b", "c", "d", "e", "f"),
  label_size = 12,
  rel_heights = c(0.95, 1.02, 1.03)
)

save_pub_figure(fig5, "Figure_5", width_in = 13.2, height_in = 13.8)



cat("Generating Figure 6 (Clinical, Germline, MSK-IMPACT, GTEx, Immune & Survival)...\n")

cv_roc_df <- read_csv(file.path(data_dir, "validation_clinvar_roc_auc_benchmarks.csv"), show_col_types = FALSE)
msk_df    <- read_csv(file.path(data_dir, "msk_impact_10945_clinical_cohort_replication.csv"), show_col_types = FALSE)
surv_df   <- read_csv(file.path(data_dir, "validation_clinical_survival_benchmarks.csv"), show_col_types = FALSE)
pb_df     <- read_csv(file.path(data_dir, "tcga_pancan_miRNA_biogenesis_patient_level_burden.csv"), show_col_types = FALSE)
sl_ago2   <- read_csv(file.path(data_dir, "tcga_pancan_miRNA_biogenesis_sample_level_alterations.csv.gz"), show_col_types = FALSE) %>%
  filter(geneSymbol == "AGO2") %>%
  select(sampleId, ago2_gistic = gistic2DiscreteCNA)

# Panel a: gnomAD v4 Germline LoF Intolerance (LOEUF) vs Missense Z-Score across 15 genes
gm_plot <- gene_meta %>%
  mutate(
    grp = case_when(
      geneSymbol %in% c("AGO1", "AGO2", "AGO3", "AGO4") ~ "Argonaute RISC effectors",
      in_Core10_tier == 1                               ~ "Core biogenesis enzymes/cofactors",
      in_Extended12_tier == 1                           ~ "Oncofetal let-7 regulators",
      TRUE                                              ~ "SMAD / TP53 regulators"
    )
  )

p6a <- ggplot(gm_plot, aes(x = gnomad_LOEUF, y = gnomad_mis_z, colour = grp)) +
  geom_vline(xintercept = 0.35, linetype = "dashed", colour = "#DC2626", linewidth = 0.4) +
  geom_hline(yintercept = 3.09, linetype = "dotted", colour = "grey45", linewidth = 0.4) +
  geom_point(size = 2.9, alpha = 0.92) +
  geom_text(aes(label = geneSymbol), nudge_y = 0.34, size = 2.55, fontface = "bold", show.legend = FALSE) +
  scale_colour_manual(
    values = c(
      "Argonaute RISC effectors"          = "#7C3AED",
      "Core biogenesis enzymes/cofactors" = "#0F5499",
      "Oncofetal let-7 regulators"        = "#DB2777",
      "SMAD / TP53 regulators"            = "#D97706"
    ),
    name = NULL
  ) +
  scale_x_continuous(limits = c(0.02, 1.06), breaks = seq(0.1, 1.0, 0.2)) +
  scale_y_continuous(limits = c(0.6, 8.8), breaks = seq(1, 8, 2)) +
  labs(
    title = "Human germline constraint (gnomAD v4)",
    subtitle = "Dashed: LOEUF < 0.35 LoF-intolerant | Dotted: missense Z > 3.09",
    x = "gnomAD LOEUF (lower = stronger LoF selection)",
    y = "gnomAD missense constraint Z-score"
  ) +
  theme_cowplot(font_size = 9) +
  theme(
    plot.title = element_text(size = 9.5, face = "bold"),
    plot.subtitle = element_text(size = 7.4, colour = "#334155"),
    panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
    legend.position = "bottom",
    legend.text = element_text(size = 7.0),
    plot.margin = margin(6, 10, 6, 6)
  ) +
  guides(colour = guide_legend(ncol = 2, byrow = TRUE))

# Panel b: ClinVar Pathogenic vs Benign Missense ROC-AUC Benchmark (n = 612 assertions)
cv_sub <- cv_roc_df %>%
  filter(predictor != "Joint_Evolutionary_3D_Composite") %>%
  mutate(
    subset_lbl = case_when(
      benchmark_subset == "ClinVar_All15_Missense_AllStars"         ~ "All-15 ClinVar\n(n = 612)",
      benchmark_subset == "ClinVar_All15_Missense_HighConf_ge2Stars" ~ ">=2-star expert\n(n = 384)",
      benchmark_subset == "ClinVar_NonTP53_14Genes_Missense_AllStars" ~ "14 non-TP53\n(n = 222)",
      TRUE                                                          ~ "Core-10 only\n(n = 170)"
    ),
    subset_lbl = factor(subset_lbl, levels = rev(c(
      "All-15 ClinVar\n(n = 612)",
      ">=2-star expert\n(n = 384)",
      "14 non-TP53\n(n = 222)",
      "Core-10 only\n(n = 170)"
    ))),
    pred_lbl = case_when(
      predictor == "AlphaMissense"      ~ "AlphaMissense",
      predictor == "neg_ESM1v_650M_LLR" ~ "ESM-1v (650M)",
      predictor == "neg_ESM2_3B_LLR"    ~ "ESM-2 (3B)",
      TRUE                              ~ "ESM-IF1 (3D fold)"
    ),
    pred_lbl = factor(pred_lbl, levels = c("AlphaMissense", "ESM-1v (650M)", "ESM-2 (3B)", "ESM-IF1 (3D fold)"))
  )

p6b <- ggplot(cv_sub, aes(x = roc_auc, y = subset_lbl, colour = pred_lbl, shape = pred_lbl)) +
  geom_vline(xintercept = 0.5, linetype = "dotted", colour = "grey50", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = roc_ci95_low, xmax = roc_ci95_high), position = position_dodge(width = 0.62), height = 0.25, linewidth = 0.55) +
  geom_point(position = position_dodge(width = 0.62), size = 2.4) +
  scale_colour_manual(
    values = c(
      "AlphaMissense"     = "#DC2626",
      "ESM-1v (650M)"     = "#0F5499",
      "ESM-2 (3B)"        = "#0D9488",
      "ESM-IF1 (3D fold)" = "#D97706"
    ),
    name = NULL
  ) +
  scale_shape_manual(values = c(16, 17, 15, 18), name = NULL) +
  scale_x_continuous(limits = c(0.48, 1.01), breaks = seq(0.5, 1.0, 0.1)) +
  labs(
    title = "Human clinical ClinVar missense benchmark",
    subtitle = "AlphaMissense AUC = 0.978 (all) / 0.980 (>=2-star) / 0.968 (Core-10)",
    x = "ClinVar Pathogenic vs Benign ROC-AUC (95% CI)",
    y = NULL
  ) +
  theme_cowplot(font_size = 9) +
  theme(
    plot.title = element_text(size = 9.5, face = "bold"),
    plot.subtitle = element_text(size = 7.4, colour = "#991B1B", face = "bold"),
    panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.35),
    legend.position = "bottom",
    legend.text = element_text(size = 7.1),
    plot.margin = margin(6, 10, 6, 6)
  ) +
  guides(colour = guide_legend(ncol = 2, byrow = TRUE))

# Panel c: Independent Clinical Replication in MSK-IMPACT (N = 10,945) vs TCGA (N = 10,195)
msk_non_tp53 <- msk_df %>%
  filter(geneSymbol != "TP53") %>%
  select(
    geneSymbol,
    `TCGA all mut (N=10,195)`       = tcga_mut_freq_N10195,
    `MSK-IMPACT all mut (N=10,945)` = msk_mut_freq,
    `TCGA 3D pocket/hotspot`        = tcga_3dpocket_or_hotspot_freq,
    `MSK 3D pocket/hotspot`         = msk_3dpocket_or_hotspot_freq
  ) %>%
  pivot_longer(-geneSymbol, names_to = "metric", values_to = "freq") %>%
  mutate(
    geneSymbol = factor(geneSymbol, levels = c("SMAD2", "DICER1", "SMAD4")),
    metric = factor(metric, levels = c(
      "TCGA all mut (N=10,195)", "MSK-IMPACT all mut (N=10,945)",
      "TCGA 3D pocket/hotspot", "MSK 3D pocket/hotspot"
    ))
  )

p6c <- ggplot(msk_non_tp53, aes(x = freq * 100, y = geneSymbol, fill = metric)) +
  geom_col(position = position_dodge(width = 0.76), width = 0.68, colour = "white", linewidth = 0.2) +
  scale_fill_manual(
    values = c(
      "TCGA all mut (N=10,195)"       = "#93C5FD",
      "MSK-IMPACT all mut (N=10,945)" = "#1D4ED8",
      "TCGA 3D pocket/hotspot"        = "#FCA5A5",
      "MSK 3D pocket/hotspot"         = "#DC2626"
    ),
    name = NULL
  ) +
  scale_x_continuous(limits = c(0, 4.8), breaks = seq(0, 4.5, 1.0)) +
  labs(
    title = "MSK-IMPACT (N = 10,945) vs TCGA clinical replication",
    subtitle = "TP53 also concordant: MSK 41.5% (14.0% 3D/hs) vs TCGA 37.0% (12.7%)",
    x = "Somatic alteration prevalence (% of profiled patients)",
    y = NULL
  ) +
  theme_cowplot(font_size = 9) +
  theme(
    plot.title = element_text(size = 9.5, face = "bold"),
    plot.subtitle = element_text(size = 7.4, colour = "#1E3A8A", face = "italic"),
    panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.35),
    legend.position = "bottom",
    legend.text = element_text(size = 7.0),
    legend.key.size = unit(0.28, "cm"),
    plot.margin = margin(6, 10, 6, 6)
  ) +
  guides(fill = guide_legend(ncol = 2, byrow = TRUE))

# Panel d: GTEx v10 54-Tissue Healthy Baseline Specificity (tau vs median TPM)
gm_6d <- gm_plot %>%
  mutate(
    lbl_x = case_when(
      geneSymbol == "DROSHA" ~ gtex_tau_tissue_specificity - 0.038,
      geneSymbol == "SMAD4"  ~ gtex_tau_tissue_specificity - 0.036,
      geneSymbol == "DICER1" ~ gtex_tau_tissue_specificity + 0.042,
      geneSymbol == "LIN28B" ~ gtex_tau_tissue_specificity - 0.055,
      geneSymbol == "LIN28A" ~ gtex_tau_tissue_specificity + 0.005,
      TRUE                   ~ gtex_tau_tissue_specificity
    ),
    lbl_y = case_when(
      geneSymbol == "DROSHA" ~ gtex_pan54_median_tpm - 1.35,
      geneSymbol == "SMAD4"  ~ gtex_pan54_median_tpm + 1.35,
      geneSymbol == "DICER1" ~ gtex_pan54_median_tpm + 1.25,
      geneSymbol == "LIN28B" ~ gtex_pan54_median_tpm + 2.55,
      geneSymbol == "LIN28A" ~ gtex_pan54_median_tpm + 1.15,
      TRUE                   ~ gtex_pan54_median_tpm + 1.25
    )
  )

p6d <- ggplot(gm_6d, aes(x = gtex_tau_tissue_specificity, y = gtex_pan54_median_tpm, colour = grp)) +
  geom_vline(xintercept = 0.80, linetype = "dashed", colour = "#DB2777", linewidth = 0.45) +
  geom_point(size = 2.9, alpha = 0.92) +
  geom_text(aes(x = lbl_x, y = lbl_y, label = geneSymbol), size = 2.55, fontface = "bold", show.legend = FALSE) +
  scale_colour_manual(
    values = c(
      "Argonaute RISC effectors"          = "#7C3AED",
      "Core biogenesis enzymes/cofactors" = "#0F5499",
      "Oncofetal let-7 regulators"        = "#DB2777",
      "SMAD / TP53 regulators"            = "#D97706"
    ),
    guide = "none"
  ) +
  scale_x_continuous(limits = c(0.16, 1.05), breaks = seq(0.2, 1.0, 0.2)) +
  scale_y_continuous(limits = c(-0.8, 28.0), breaks = seq(0, 25, 5)) +
  labs(
    title = "Healthy adult expression & specificity (GTEx v10, 54 tissues)",
    subtitle = "LIN28A (tau = 0.997) & LIN28B (tau = 0.978) restricted to adult testis",
    x = "Tissue-specificity index tau (0 = ubiquitous, 1 = tissue-restricted)",
    y = "Pan-54-tissue median expression (TPM)"
  ) +
  theme_cowplot(font_size = 9) +
  theme(
    plot.title = element_text(size = 9.5, face = "bold"),
    plot.subtitle = element_text(size = 7.4, colour = "#9D174D", face = "bold"),
    panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
    plot.margin = margin(6, 10, 6, 6)
  )

# Panel e: Thorsson Pan-Cancer HRD & IFN-gamma Immune Landscape across AGO2 CNA states
pb_ago2 <- pb_df %>%
  inner_join(sl_ago2, by = "sampleId") %>%
  filter(!is.na(hrdScore), ago2_gistic %in% c(-1, 0, 1, 2)) %>%
  mutate(
    cna_lbl = factor(
      case_when(
        ago2_gistic == -1 ~ "Shallow loss (-1)",
        ago2_gistic ==  0 ~ "Diploid (0)",
        ago2_gistic ==  1 ~ "Low gain (+1)",
        ago2_gistic ==  2 ~ "Focal amp (+2)"
      ),
      levels = c("Shallow loss (-1)", "Diploid (0)", "Low gain (+1)", "Focal amp (+2)")
    )
  )

p6e <- ggplot(pb_ago2, aes(x = hrdScore, y = cna_lbl, fill = cna_lbl)) +
  geom_boxplot(width = 0.56, outlier.shape = NA, alpha = 0.90, colour = "grey20", linewidth = 0.4) +
  scale_fill_manual(
    values = c(
      "Shallow loss (-1)" = "#93C5FD",
      "Diploid (0)"       = "#CBD5E1",
      "Low gain (+1)"     = "#FDBA74",
      "Focal amp (+2)"    = "#DC2626"
    ),
    guide = "none"
  ) +
  coord_cartesian(xlim = c(0, 72)) +
  labs(
    title = "AGO2 8q24.3 copy state vs HRD & immune state (n = 9,869)",
    subtitle = "Focal +2 vs diploid: HRD 33 vs 7 (p = 1.4e-183) | IFN-gamma +0.34 vs -0.19 (p = 2.7e-26)",
    x = "Thorsson Homologous Recombination Deficiency (HRD) score",
    y = NULL
  ) +
  theme_cowplot(font_size = 9) +
  theme(
    plot.title = element_text(size = 9.5, face = "bold"),
    plot.subtitle = element_text(size = 7.3, colour = "#991B1B", face = "bold"),
    panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.35),
    plot.margin = margin(6, 10, 6, 6)
  )

# Panel f: Cancer-Type-Stratified Survival Forest Plot (TCGA-CDR & MSK-IMPACT)
surv_plot <- surv_df %>%
  filter(test_id %in% c(
    "S_MSK_IMPACT_TP53_Pathogenic_Stratified_OS",
    "S6_LIN28AB_Oncofetal_Reactivation_Stratified_OS",
    "S9_TP53_PositiveControl_Stratified_OS",
    "S4_AGO2_FocalAmp_p2_Stratified_OS",
    "S7_1p34_AGO1_3_4_CoDeletion_Stratified_OS",
    "S1_Core10_StrictAlter_Stratified_OS"
  )) %>%
  mutate(
    lbl = case_when(
      test_id == "S_MSK_IMPACT_TP53_Pathogenic_Stratified_OS"      ~ "MSK-IMPACT TP53 control (p = 2.9e-23)",
      test_id == "S6_LIN28AB_Oncofetal_Reactivation_Stratified_OS" ~ "LIN28A/B reactivation Z >= 2 (p = 8.0e-4)",
      test_id == "S9_TP53_PositiveControl_Stratified_OS"           ~ "TCGA TP53 control (p = 1.4e-3)",
      test_id == "S4_AGO2_FocalAmp_p2_Stratified_OS"               ~ "AGO2 focal amp +2 (p = 0.18)",
      test_id == "S7_1p34_AGO1_3_4_CoDeletion_Stratified_OS"       ~ "1p34.3 AGO1/3/4 co-deletion (p = 0.032)",
      TRUE                                                         ~ "Pooled Core-10 any strict (p = 0.66)"
    ),
    lbl = factor(lbl, levels = rev(c(
      "LIN28A/B reactivation Z >= 2 (p = 8.0e-4)",
      "1p34.3 AGO1/3/4 co-deletion (p = 0.032)",
      "AGO2 focal amp +2 (p = 0.18)",
      "Pooled Core-10 any strict (p = 0.66)",
      "TCGA TP53 control (p = 1.4e-3)",
      "MSK-IMPACT TP53 control (p = 2.9e-23)"
    ))),
    sig_col = if_else(as.numeric(p_value_two_sided) < 0.05, "Significant (p < 0.05)", "Null / non-significant")
  )

p6f <- ggplot(surv_plot, aes(x = hazard_ratio_stratified, y = lbl, colour = sig_col)) +
  geom_vline(xintercept = 1.0, linetype = "dashed", colour = "grey40", linewidth = 0.45) +
  geom_errorbarh(aes(xmin = hr_ci95_low, xmax = hr_ci95_high), height = 0.25, linewidth = 0.65) +
  geom_point(size = 2.6) +
  scale_colour_manual(
    values = c("Significant (p < 0.05)" = "#DC2626", "Null / non-significant" = "#64748B"),
    name = NULL
  ) +
  scale_x_continuous(limits = c(0.80, 1.92), breaks = seq(0.8, 1.8, 0.2)) +
  labs(
    title = "Cancer-type-stratified overall survival (Mantel-Haenszel)",
    subtitle = "Stratified across 32 TCGA studies (TCGA-CDR) & 43 MSK-IMPACT types",
    x = "Stratified Hazard Ratio for Death (95% CI)",
    y = NULL
  ) +
  theme_cowplot(font_size = 9) +
  theme(
    plot.title = element_text(size = 9.5, face = "bold"),
    plot.subtitle = element_text(size = 7.4, colour = "#334155"),
    panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.35),
    legend.position = "bottom",
    legend.text = element_text(size = 7.3),
    plot.margin = margin(6, 10, 6, 6)
  )

fig6 <- plot_grid(
  p6a, p6b, p6c, p6d, p6e, p6f,
  ncol = 2, nrow = 3,
  labels = c("a", "b", "c", "d", "e", "f"),
  label_size = 12,
  rel_heights = c(1.0, 1.0, 1.0)
)

save_pub_figure(fig6, "Figure_6", width_in = 13.2, height_in = 13.8)

cat("\nAll 6 publication figures generated cleanly in:", output_dir, "\n")
