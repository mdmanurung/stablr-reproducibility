# Paper Figure 4 / S4 — Tem GZMK+/GZMB+ (TemK/TemB) ratio biology
# (Terekhova et al. 2026).
#
# Ported from the authors' `Figures_pipeline/figures_pipeline.R`
# (repo teresho4/Pan_atlas_immune_aging_TemK_TemB_ratio, section
# "Figure 4, S4", lines 962-1290). Ratio of Tem GZMK+ to Tem GZMB+ cell
# counts per donor, across sex, ancestry, BMI category, age bands, and
# CMV status; the CMV-vs-CMV+ population-change screen (SoundLife + KCL,
# Young/Old Wilcoxon tests); and the CMV-colored CD8 PCA panel.
#
# NOT ported (data not shipped in this compendium): the Aquino et al.
# external-validation panel (authors' `Aquino_immune_percentages.csv`).
# See the README in this directory.
#
# Run: ./run_fig4.sh

figdir <- "fig4_tem_ratio"
.args <- commandArgs(trailingOnly = FALSE); .fa <- grep("^--file=", .args, value = TRUE)
.root <- if (length(.fa)) dirname(dirname(normalizePath(sub("^--file=", "", .fa[1])))) else normalizePath(getwd())
source(file.path(.root, "fig_support.R"))
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(rstatix); library(forcats); library(tibble)
})

cells <- read_paper_cells()
minor_pct <- as.data.frame(paper_minor_percentages(cells))
# PORT NOTE: the authors' pipeline carries age/sex/ethnicity/CMV/BMI on every
# per-cell row upstream; the compendium's helper returns counts only, so merge
# donor-level covariates here (one row per donor).
minor_pct <- merge(
  minor_pct,
  cells %>% select(donor_id, age, sex, ethnicity, CMV, BMI) %>%
    # authors' BMI parse: raw file mixes "22" and "22.0" spellings
    mutate(BMI = gsub(".0", "", BMI, fixed = TRUE)) %>% distinct(),
  by = "donor_id"
)
out <- fig_output_dir(figdir)

cd8_ratio <- minor_pct %>%
  filter(Major == "CD8 T", Annotation %in% c("Tem GZMK+", "Tem GZMB+")) %>%
  select(donor_id, Annotation, n_cells, percent) %>%
  pivot_wider(names_from = Annotation, values_from = c(n_cells, percent)) %>%
  # PORT NOTE: values_fill is unreliable across tidyr builds; zero-fill after pivoting
  mutate(across(where(is.numeric), ~ tidyr::replace_na(.x, 0))) %>%
  left_join(
    cells %>% select(donor_id, sex, BMI, CMV, age, ethnicity, Dataset) %>%
      mutate(BMI = gsub(".0", "", BMI, fixed = TRUE)) %>% distinct(),
    by = "donor_id"
  ) %>%
  mutate(
    Ratio = `n_cells_Tem GZMK+` / `n_cells_Tem GZMB+`,
    Ratio = ifelse(is.infinite(Ratio),
                   max(Ratio[is.finite(Ratio)], na.rm = TRUE),
                   Ratio)
  )

ratio_lab <- ylab("log1p(Tem GZMK+/Tem GZMB+)")

# ---- Ratio by sex (verbatim) -----------------------------------------------------
p_sex <- ggplot(cd8_ratio, aes(x = sex, y = log1p(Ratio), fill = sex, color = sex)) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
  geom_boxplot(outlier.shape = NA, width = 0.5, color = "black") +
  scale_fill_manual(values = c("#7293A0", "#A37986")) +
  scale_color_manual(values = c("#7293A0", "#A37986")) +
  xlab("") + ratio_lab + theme_minimal() +
  theme(axis.text.x = element_text(color = "black"), legend.position = "none")
ggsave(file.path(out, "fig4_ratio_by_sex.png"), p_sex, width = 5, height = 5, dpi = 200, bg = "white")

# ---- Ratio by ancestry (verbatim) -------------------------------------------------
ratio_ethnicity <- cd8_ratio %>% filter(ethnicity %in% c("Asian", "European"))
p_eth <- ggplot(ratio_ethnicity, aes(x = ethnicity, y = log1p(Ratio), fill = ethnicity, color = ethnicity)) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
  geom_boxplot(outlier.shape = NA, width = 0.5, color = "black") +
  scale_fill_manual(values = c("#96af7e", "#bfbd7d")) +
  scale_color_manual(values = c("#96af7e", "#bfbd7d")) +
  xlab("") + ratio_lab + theme_minimal() +
  theme(axis.text.x = element_text(color = "black"), legend.position = "none")
ggsave(file.path(out, "fig4_ratio_by_ancestry.png"), p_eth, width = 5, height = 5, dpi = 200, bg = "white")

# ---- Ratio by BMI category (verbatim) ----------------------------------------------
ratio_bmi <- cd8_ratio %>%
  filter(BMI != "No_info") %>%
  mutate(BMI = as.numeric(gsub(".0", "", BMI, fixed = TRUE)),
         BMI_category = case_when(
           BMI < 18.5 ~ "Underweight",
           BMI < 25 ~ "Normal",
           BMI < 30 ~ "Overweight",
           BMI >= 30 ~ "Obese"))
ratio_bmi$BMI_category <- factor(ratio_bmi$BMI_category,
                                 levels = c("Underweight", "Normal", "Overweight", "Obese"))
p_bmi <- ggplot(ratio_bmi, aes(x = BMI_category, y = log1p(Ratio), fill = BMI_category, color = BMI_category)) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
  geom_boxplot(outlier.shape = NA, width = 0.5, color = "black") +
  scale_fill_manual(values = c("#b0c4b1", "#88a096", "#d9b382", "#a65e3f")) +
  scale_color_manual(values = c("#b0c4b1", "#88a096", "#d9b382", "#a65e3f")) +
  xlab("") + ratio_lab + theme_minimal() +
  theme(axis.text.x = element_text(color = "black"), legend.position = "none")
ggsave(file.path(out, "fig4_ratio_by_bmi.png"), p_bmi, width = 6, height = 5, dpi = 200, bg = "white")

# ---- Ratio by age band (verbatim) ---------------------------------------------------
ratio_age <- cd8_ratio %>%
  mutate(Age_group = case_when(
    age < 35 ~ "A", age < 45 ~ "B", age < 55 ~ "C",
    age < 65 ~ "D", age < 75 ~ "E", age >= 75 ~ "F"))
ratio_age$Age_group <- factor(ratio_age$Age_group, levels = c("A", "B", "C", "D", "E", "F"))
p_age <- ggplot(ratio_age, aes(x = Age_group, y = log1p(Ratio), fill = Age_group, color = Age_group)) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
  geom_boxplot(outlier.shape = NA, width = 0.5, color = "black") +
  scale_fill_manual(values = c("#d4b483", "#a29dae", "#7d8ca3", "#b2836d", "#8e8c68", "#9c7a97")) +
  scale_color_manual(values = c("#d4b483", "#a29dae", "#7d8ca3", "#b2836d", "#8e8c68", "#9c7a97")) +
  xlab("") + ratio_lab + theme_minimal() +
  theme(axis.text.x = element_text(color = "black"), legend.position = "none")
ggsave(file.path(out, "fig4_ratio_by_age_band.png"), p_age, width = 7, height = 5, dpi = 200, bg = "white")

# ---- Ratio by CMV, faceted by dataset (verbatim) -------------------------------------
ratio_cmv <- cd8_ratio %>% filter(CMV != "No_info")
p_cmv <- ggplot(ratio_cmv, aes(x = CMV, y = log1p(Ratio), fill = CMV)) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  geom_jitter(width = 0, alpha = 0.7, size = 1.5, shape = 21, color = "#706e6e") +
  scale_fill_manual(values = c("#8fb996", "#c26a6a")) +
  xlab("") + ratio_lab + theme_minimal() +
  theme(axis.text.x = element_text(color = "black"), legend.position = "none") +
  facet_wrap(~Dataset)
ggsave(file.path(out, "fig4_ratio_by_cmv_faceted.png"), p_cmv, width = 12, height = 6, dpi = 200, bg = "white")
cat("wrote: outputs/fig4_ratio_by_*.png\n")

# ---- CMV vs CMV+ population changes (verbatim; SoundLife + KCL) -----------------------
cmv_df <- minor_pct %>%
  filter(Dataset %in% c("SoundLife", "KCL")) %>%
  # minor_pct already carries donor-level CMV (merged above); the authors'
  # version joins it here from their per-cell table
  filter(CMV %in% c("Negative", "Positive")) %>%
  mutate(Age_group = case_when(age < 40 ~ "Young", age > 50 ~ "Old", TRUE ~ NA_character_)) %>%
  filter(!is.na(Age_group)) %>% filter(age < 75 & CMV != "No_info")
cmv_df$Major_Annotation <- paste(cmv_df$Major, cmv_df$Annotation)

result <- cmv_df %>%
  group_by(Dataset, Age_group, Major_Annotation) %>%
  rstatix::wilcox_test(percent ~ CMV) %>%
  mutate(p.adj = p * 55)
result_left <- result %>%
  group_by(Age_group, Major_Annotation) %>%
  summarise(n_sig = sum(p.adj < 0.05, na.rm = TRUE)) %>%
  filter(n_sig >= 1)
result <- result %>% filter(Major_Annotation %in% result_left$Major_Annotation)
result[result$p.adj > 1, "p.adj"] <- 1

plot_df <- result %>%
  mutate(neglog_padj = -log10(p.adj),
         neglog_padj = ifelse(is.infinite(neglog_padj), NA, neglog_padj))
order_df <- plot_df %>%
  group_by(Major_Annotation) %>%
  summarise(max_p = max(neglog_padj, na.rm = TRUE), .groups = "drop")
plot_df <- plot_df %>%
  left_join(order_df, by = "Major_Annotation") %>%
  mutate(Major_Annotation = fct_reorder(Major_Annotation, max_p),
         Age_group = factor(Age_group, levels = c("Young", "Old")))
bar_df <- plot_df %>%
  group_by(Major_Annotation, Age_group) %>%
  summarise(mean_neglog_padj = mean(neglog_padj, na.rm = TRUE),
            sem = sd(neglog_padj, na.rm = TRUE) / sqrt(sum(!is.na(neglog_padj))),
            upper = mean_neglog_padj + sem, .groups = "drop")

p_cmv_change <- ggplot(bar_df, aes(y = Major_Annotation, group = Age_group)) +
  geom_col(aes(x = mean_neglog_padj, fill = Age_group), position = position_dodge(width = 0.75),
           width = 0.65, color = "black", alpha = 0.8) +
  geom_errorbarh(aes(xmin = mean_neglog_padj, xmax = upper), position = position_dodge(width = 0.75),
                 height = 0.18, color = "black") +
  geom_point(data = plot_df, aes(x = neglog_padj, y = Major_Annotation, group = Age_group, fill = Age_group),
             position = position_dodge(width = 0.75), size = 2, alpha = 0.8,
             inherit.aes = FALSE, shape = 21) +
  scale_fill_manual(values = c("Young" = "#D9A75F", "Old" = "#8FA4BF")) +
  labs(title = "CMV- vs CMV+ Population Changes", x = expression(-log[10]("p.adj")),
       y = NULL, fill = "Age group") +
  theme_bw()
ggsave(file.path(out, "fig4_cmv_population_changes.png"), p_cmv_change, width = 9, height = 8, dpi = 200, bg = "white")
write.csv(as.data.frame(result), file.path(out, "fig4_cmv_wilcoxon_tests.csv"), row.names = FALSE)
cat("wrote: outputs/fig4_cmv_*.png/csv\n")

# ---- CMV-colored CD8 PCA panel (verbatim; recomputes the Fig 3 PCA) --------------------
cd8 <- c("RTE", "Naive", "Naive-IFN", "Tcm CCR4-", "Tcm CCR4+", "Trm", "Tem GZMK+",
         "Tem GZMB+", "Temra", "NKT-like", "Tmem KLRC2+", "HLA-DR+", "Proliferative")
cd8_mat <- minor_pct %>%
  filter(Major == "CD8 T", Annotation %in% cd8) %>%
  select(donor_id, Annotation, percent) %>%
  pivot_wider(names_from = Annotation, values_from = percent) %>%
  mutate(across(where(is.numeric), ~ tidyr::replace_na(.x, 0)))
pca <- prcomp(column_to_rownames(cd8_mat, "donor_id"), scale. = FALSE, center = TRUE)
pca_df <- as.data.frame(pca$x) %>%
  rownames_to_column("donor_id") %>%
  left_join(cells %>% select(donor_id, age, Dataset, CMV) %>% distinct(), by = "donor_id")
global_hull <- pca_df %>% slice(chull(PC1, PC2))
p_cmv_pca <- ggplot() +
  geom_point(data = pca_df %>% filter(CMV == "No_info"), aes(x = PC1, y = PC2), color = "#d4d4d4", size = 1.5) +
  geom_point(data = pca_df %>% filter(CMV != "No_info"), aes(x = PC1, y = PC2, color = CMV), size = 1.5) +
  theme_minimal() +
  labs(x = paste0("PC1 (", round(summary(pca)$importance[2, 1] * 100, 1), "% variance)"),
       y = paste0("PC2 (", round(summary(pca)$importance[2, 2] * 100, 1), "% variance)")) +
  theme(legend.position = "bottom") +
  scale_color_manual(values = c("#8fb996", "#c26a6a")) +
  geom_polygon(data = global_hull, aes(x = PC1, y = PC2), fill = NA, color = "#696969", linewidth = 0.6)
ggsave(file.path(out, "fig4_pca_cmv.png"), p_cmv_pca, width = 8, height = 6, dpi = 200, bg = "white")
cat("wrote: outputs/fig4_pca_cmv.png\n")

cat("fig4 done. (Aquino et al. external-validation panel not ported: data not shipped; see README.)\n")
