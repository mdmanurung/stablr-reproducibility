# Paper Figure 2 / S2 / S3 — Cell-type age associations (Terekhova et al. 2026).
#
# Ported from the authors' `Figures_pipeline/figures_pipeline.R`
# (repo teresho4/Pan_atlas_immune_aging_TemK_TemB_ratio, section
# "Figure 2, S2, S3", lines 183-842). Model and plotting code follows the
# authors' version: per (dataset, population) Pearson r of percent vs age,
# lm(percent ~ Age + Sex + Ethnicity) with BH adjustment, the association
# heatmap + boxplot, major-population boxplots and alluvial, CD8 minor
# scatters, the oldest-old age-band comparison (Dunn tests), the selected
# sex/ancestry panels, and the variance partitioning of minor populations.
#
# Run: ./run_fig2.sh

figdir <- "fig2_age_associations"
.args <- commandArgs(trailingOnly = FALSE); .fa <- grep("^--file=", .args, value = TRUE)
.root <- if (length(.fa)) dirname(dirname(normalizePath(sub("^--file=", "", .fa[1])))) else normalizePath(getwd())
source(file.path(.root, "fig_support.R"))
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(patchwork)
  library(ggalluvial); library(rstatix); library(broom); library(stringr); library(tibble)
  library(variancePartition)
})

cells <- read_paper_cells()
cat(sprintf("cells: %d rows\n", nrow(cells)))

major_pct <- as.data.frame(paper_major_percentages(cells))
minor_pct <- as.data.frame(paper_minor_percentages(cells))
cat(sprintf("major rows: %d; minor rows: %d\n", nrow(major_pct), nrow(minor_pct)))

# ---- Authors' age-association model (verbatim) ---------------------------------
prepare_age_association_df <- function(data, major_percentages) {
  df_summary <- major_percentages
  make_ratio <- function(df, labels, annotation_name, denominator_cols) {
    df_wide <- df %>%
      filter(Annotation %in% names(labels)) %>%
      mutate(Type = recode(Annotation, !!!labels)) %>%
      group_by(Dataset, donor_id, Type) %>%
      summarise(n = sum(n_cells), .groups = "drop") %>%
      pivot_wider(names_from = Type, values_from = n, values_fill = 0)
    df_wide %>%
      mutate(denominator = rowSums(across(all_of(denominator_cols))),
             percent = ifelse(denominator == 0, NA, CD4 / denominator),
             Annotation = annotation_name) %>%
      select(Dataset, donor_id, Annotation, percent)
  }
  ratio_rows <- make_ratio(df_summary,
                           labels = c("CD4 T" = "CD4", "CD8 T" = "CD8", "MAIT" = "MAIT"),
                           annotation_name = "CD4/CD8 T cells",
                           denominator_cols = c("CD8", "MAIT"))
  ratio_rows_2 <- make_ratio(df_summary,
                             labels = c("CD4 T" = "CD4", "CD8 T" = "CD8"),
                             annotation_name = "CD4/TRAV1-2- CD8 T cells",
                             denominator_cols = c("CD8"))
  meta <- data %>% select(donor_id, sex, age, ethnicity) %>% distinct()
  bind_rows(df_summary, ratio_rows, ratio_rows_2) %>%
    left_join(meta, by = "donor_id") %>%
    mutate(Dataset = factor(Dataset), Sex = factor(sex),
           Ethnicity = factor(ethnicity), Age = age) %>%
    as.data.frame()
}

calculate_age_associations <- function(df) {
  model_results <- df %>%
    group_by(Dataset, Annotation) %>%
    do({
      subdf <- .
      r <- if (length(unique(subdf$Age)) > 1 && length(unique(subdf$percent)) > 1) {
        cor(subdf$Age, subdf$percent, method = "pearson", use = "complete.obs")
      } else NA
      if (n_distinct(subdf$Ethnicity) > 1) {
        model <- lm(percent ~ Age + Sex + Ethnicity, data = subdf)
      } else {
        model <- lm(percent ~ Age + Sex, data = subdf)
      }
      out <- tidy(model); out$r <- r; out
    }) %>%
    filter(term == "Age") %>%
    ungroup() %>%
    group_by(Dataset) %>%
    mutate(p.adj = p.adjust(p.value, method = "BH")) %>%
    ungroup() %>%
    mutate(sig = case_when(
      is.na(p.adj) ~ " ",
      p.adj < 0.001 & estimate > 0 ~ "+++",
      p.adj < 0.001 & estimate < 0 ~ "---",
      p.adj < 0.01  & estimate > 0 ~ "++",
      p.adj < 0.01  & estimate < 0 ~ "--",
      p.adj < 0.05  & estimate > 0 ~ "+",
      p.adj < 0.05  & estimate < 0 ~ "-",
      TRUE ~ "ns")) %>%
    select(Dataset, Annotation, estimate, std.error, statistic, p.value, p.adj, r, sig)

  model_results$Annotation[model_results$Annotation == "CD8 T"] <- "TRAV1-2- CD8 T cells"
  model_results$Dataset <- factor(model_results$Dataset, levels = PAPER_DATASETS)
  model_results$Annotation <- ifelse(
    !grepl("cells", model_results$Annotation) & model_results$Annotation != "ILCs",
    paste(model_results$Annotation, "cells"), model_results$Annotation)
  model_results$Annotation <- factor(model_results$Annotation, levels = rev(c(
    "MAIT cells", "TRAV1-2- CD8 T cells", "CD4 T cells", "gd T cells", "DN T cells",
    "ILCs", "NK cells", "Myeloid cells", "B cells", "Progenitor cells",
    "CD4/CD8 T cells", "CD4/TRAV1-2- CD8 T cells")))
  model_results
}

plot_age_association_summary <- function(model_results, celltype_colors_major) {
  p_box <- ggplot(model_results, aes(y = Annotation, x = r, fill = Annotation)) +
    geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.6) +
    geom_point(position = position_jitterdodge(jitter.width = 0), size = 1, alpha = 0.8) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.4) +
    scale_x_continuous(name = "Correlation", limits = c(-0.5, 0.5)) +
    scale_fill_manual(values = celltype_colors_major) +
    theme_bw(base_size = 10) +
    theme(legend.position = "none", axis.title.y = element_blank(),
          axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
    coord_fixed(ratio = 0.18)
  heatmap_data <- model_results %>% mutate(r_plot = ifelse(sig == "ns", 0, r))
  p_heatmap <- ggplot(heatmap_data, aes(x = Dataset, y = Annotation, fill = r_plot)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sig), size = 3) +
    scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027", midpoint = 0,
                         name = "Pearson cor", limits = c(-0.5, 0.5), na.value = "white") +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
          axis.text.y = element_text(color = "black"), axis.title = element_blank(),
          panel.grid = element_blank(), legend.position = "top",
          legend.key.size = unit(0.5, "cm"), legend.title = element_text(size = 8),
          legend.text = element_text(size = 7)) +
    coord_fixed(ratio = 0.9)
  p_heatmap + p_box
}

df_age <- prepare_age_association_df(data = as.data.frame(cells), major_percentages = major_pct)
model_results <- calculate_age_associations(df_age)
out <- fig_output_dir(figdir)
ggsave(file.path(out, "fig2_age_association_summary.png"),
       plot_age_association_summary(model_results, PAPER_CELLTYPE_COLORS_MAJOR),
       width = 10, height = 6, dpi = 200, bg = "white")
cat("wrote: outputs/fig2_age_association_summary.png\n")
write.csv(as.data.frame(model_results), file.path(out, "fig2_age_association_models.csv"), row.names = FALSE)

# ---- Major-population boxplots + alluvial (verbatim) ---------------------------
plot_major_percent_boxplot <- function(major_percentages,
                                       annotation_levels = c("CD4 T", "CD8 T", "MAIT", "gd T",
                                                             "DN T", "ILCs", "NK", "Myeloid",
                                                             "B", "Progenitor")) {
  df <- major_percentages %>%
    mutate(Dataset = factor(Dataset, levels = PAPER_DATASETS),
           Annotation = factor(Annotation, levels = annotation_levels))
  donor_n_text <- df %>% distinct(Dataset, donor_id) %>% count(Dataset) %>%
    mutate(label = paste0("n = ", n))
  donor_caption <- paste(donor_n_text$Dataset, donor_n_text$label, sep = ": ", collapse = "   ")
  ggplot(df, aes(x = Dataset, y = percent, fill = Dataset, color = Dataset)) +
    geom_jitter(width = 0.15, size = 0.2, alpha = 0.8) +
    stat_boxplot(geom = "errorbar", width = 0.4) +
    geom_boxplot(width = 0.5, outlier.shape = NA, color = "black", alpha = 0.6) +
    scale_fill_manual(values = PAPER_DATASET_COLORS) +
    scale_color_manual(values = PAPER_DATASET_COLORS) +
    labs(y = "% of CD45+ cells", x = NULL, caption = donor_caption) +
    facet_wrap(~Annotation, scales = "free", ncol = 5) +
    theme_classic(base_size = 12) +
    theme(legend.position = "top", axis.ticks.x = element_blank(),
          panel.grid.major.y = element_line(color = "grey40", linewidth = 0.5),
          axis.text.x = element_text(angle = 45, hjust = 1, color = "black"))
}

plot_major_percent_alluvial <- function(major_percentages,
                                        annotation_levels = c("CD4 T", "CD8 T", "MAIT", "gd T",
                                                              "DN T", "ILCs", "NK", "Myeloid",
                                                              "B", "Progenitor")) {
  df_plot <- major_percentages %>%
    mutate(Dataset = factor(Dataset, levels = PAPER_DATASETS),
           Annotation = factor(Annotation, levels = annotation_levels)) %>%
    group_by(Dataset, Annotation) %>%
    summarise(mean_percent = mean(percent), .groups = "drop")
  ggplot(df_plot, aes(x = Dataset, stratum = Annotation, alluvium = Annotation,
                      y = mean_percent, fill = Annotation, label = Annotation)) +
    geom_flow(stat = "alluvium", lode.guidance = "forward", color = "black", alpha = 0.3, width = 0.7) +
    geom_stratum(width = 0.7, color = "black", alpha = 0.6) +
    scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
    scale_fill_manual(values = PAPER_CELLTYPE_COLORS_MAJOR) +
    xlab("Dataset") + ylab("") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "right", axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
          panel.grid.major.y = element_line(color = "gray70"))
}

ggsave(file.path(out, "fig2_major_percent_boxplots.png"),
       plot_major_percent_boxplot(major_pct), width = 16, height = 8, dpi = 200, bg = "white")
ggsave(file.path(out, "fig2_major_percent_alluvial.png"),
       plot_major_percent_alluvial(major_pct), width = 10, height = 6, dpi = 200, bg = "white")
cat("wrote: outputs/fig2_major_percent_*.png\n")

# ---- Minor-population age associations + CD8 scatters (verbatim) ---------------
calculate_age_associations_minor <- function(minor_percentages) {
  minor_percentages %>%
    group_by(Major, Dataset, Annotation) %>%
    do({
      subdf <- .
      r <- if (length(unique(subdf$age)) > 1 && length(unique(subdf$percent)) > 1) {
        cor(subdf$age, subdf$percent, method = "pearson", use = "complete.obs")
      } else NA
      model <- if (n_distinct(subdf$Ethnicity) > 1) {
        lm(percent ~ age + sex + Ethnicity, data = subdf)
      } else {
        lm(percent ~ age + sex, data = subdf)
      }
      out <- broom::tidy(model); out$r <- r; out
    }) %>%
    filter(term == "age") %>%
    ungroup() %>%
    group_by(Major, Dataset) %>%
    mutate(p.adj = p.adjust(p.value, method = "BH")) %>%
    ungroup()
}

donor_meta_df <- cells %>% select(donor_id, age, sex, ethnicity) %>% distinct()
minor_pct <- merge(minor_pct, donor_meta_df, by = "donor_id")
model_minor_results <- calculate_age_associations_minor(minor_pct)
write.csv(as.data.frame(model_minor_results), file.path(out, "fig2_minor_age_models.csv"), row.names = FALSE)

plot_minor_age_scatter <- function(minor_percentages, model_results, major = "CD8 T",
                                   annotations = c("RTE", "Naive", "Tcm CCR4+", "Tem GZMK+",
                                                   "Tem GZMB+", "Tmem KLRC2+", "HLA-DR+")) {
  plot_list <- list()
  for (subset in annotations) {
    df_plot <- minor_percentages %>%
      filter(Major == major, Annotation == subset) %>%
      left_join(model_results %>% filter(Major == major) %>%
                  select(Annotation, Dataset, p.adj), by = c("Annotation", "Dataset")) %>%
      mutate(Dataset = factor(Dataset, levels = PAPER_DATASETS),
             p_label = ifelse(is.na(p.adj), "adj. p = NA",
                              paste0("adj. p = ", signif(p.adj, 2))))
    labels_df <- df_plot %>% group_by(Dataset) %>% slice(1) %>% ungroup()
    plot_list[[subset]] <- ggplot(df_plot, aes(x = age, y = percent, fill = Dataset)) +
      geom_point(alpha = 0.5, shape = 21) +
      geom_smooth(method = "lm", se = FALSE, color = "black") +
      facet_wrap(~Dataset, scales = "free_y", ncol = 4) +
      geom_text(data = labels_df, aes(x = -Inf, y = Inf, label = p_label),
                hjust = -0.1, vjust = 1.1, inherit.aes = FALSE, size = 3.5) +
      scale_fill_manual(values = PAPER_DATASET_COLORS) +
      theme_bw() +
      labs(x = "Age", y = paste0("% from ", major)) +
      theme(axis.text.x = element_text(color = "black"),
            axis.text.y = element_text(color = "black"),
            legend.position = "none", plot.title = element_text(hjust = 0.5)) +
      ggtitle(subset)
  }
  plot_list
}

plot_list <- plot_minor_age_scatter(minor_percentages = minor_pct, model_results = model_minor_results)
ggsave(file.path(out, "fig2_cd8_tem_gzmk_scatter.png"), plot_list[["Tem GZMK+"]],
       width = 12, height = 6, dpi = 200, bg = "white")
ggsave(file.path(out, "fig2_cd8_tem_gzmb_scatter.png"), plot_list[["Tem GZMB+"]],
       width = 12, height = 6, dpi = 200, bg = "white")
cat("wrote: outputs/fig2_cd8_*_scatter.png\n")

# ---- Oldest-old age bands (verbatim; Dunn tests) --------------------------------
meta_keep <- cells %>% select(donor_id, sex, age, ethnicity, Dataset) %>% distinct()
donors_keep <- meta_keep %>%
  mutate(Age_group = case_when(
    age < 35 ~ "< 35",
    age >= 60 & age <= 70 & Dataset %in% c("OneK1K", "SPAC", "KCL") ~ "60-70",
    age >= 85 ~ "85-97",
    TRUE ~ NA_character_)) %>%
  filter(!is.na(Age_group))
percentages_age <- minor_pct %>%
  filter(donor_id %in% donors_keep$donor_id) %>%
  left_join(donors_keep %>% select(donor_id, Age_group), by = "donor_id") %>%
  mutate(Age_group = factor(Age_group, levels = c("< 35", "60-70", "85-97")))
dunn_results <- percentages_age %>%
  group_by(Annotation) %>%
  rstatix::dunn_test(percent ~ Age_group, p.adjust.method = "bonferroni") %>%
  filter(group1 %in% c("< 35", "60-70") & group2 %in% c("60-70", "85-97")) %>%
  mutate(y.position = 12.5)
p_old <- ggplot(percentages_age, aes(x = Age_group, y = percent, fill = Age_group, color = Age_group)) +
  geom_jitter(width = 0.2, alpha = 0.8, size = 1) +
  geom_boxplot(outlier.shape = NA, width = 0.5, color = "black") +
  facet_wrap(~Annotation, scales = "free") +
  theme_minimal() +
  scale_fill_manual(values = c("< 35" = "#d9b782", "60-70" = "#8a9167", "85-97" = "#a96f98")) +
  scale_color_manual(values = c("< 35" = "#d9b782", "60-70" = "#8a9167", "85-97" = "#a96f98"))
ggsave(file.path(out, "fig2_oldest_old_bands.png"), p_old, width = 16, height = 10, dpi = 200, bg = "white")
write.csv(as.data.frame(dunn_results), file.path(out, "fig2_oldest_old_dunn.csv"), row.names = FALSE)
cat("wrote: outputs/fig2_oldest_old_*.png/csv\n")

# ---- Variance partition of minor populations (verbatim, Asian+European) ---------
percentages_vp <- minor_pct %>% filter(ethnicity %in% c("Asian", "European"))
meta_vp <- cells %>% filter(ethnicity %in% c("Asian", "European")) %>%
  select(donor_id, sex, age, ethnicity, Dataset) %>% distinct() %>%
  rename(Sex = sex, Age = age, Ethnicity = ethnicity)
expr_df <- percentages_vp %>%
  mutate(Annotation = paste(Major, Annotation, sep = "_"),
         Annotation = str_replace_all(Annotation, " T cells_", "_"),
         Annotation = str_replace_all(Annotation, " cells_", "_")) %>%
  select(Annotation, donor_id, percent) %>%
  pivot_wider(names_from = donor_id, values_from = percent, values_fill = 0)
expr <- as.data.frame(expr_df)
rownames(expr) <- expr$Annotation
expr$Annotation <- NULL
rownames(meta_vp) <- meta_vp$donor_id
expr <- expr[, rownames(meta_vp)]
expr[is.na(expr)] <- 0
form <- ~ (1 | Sex) + (1 | Ethnicity) + (1 | Dataset) + Age
varPart <- fitExtractVarPartModel(expr, form, meta_vp)
vp_long <- as.data.frame(varPart) %>%
  rownames_to_column("population") %>%
  pivot_longer(cols = c("Age", "Sex", "Ethnicity"), names_to = "Variable", values_to = "Variance") %>%
  mutate(Variance = Variance * 100,
         Cluster = case_when(
           str_detect(population, "^CD4 T_") ~ "CD4",
           str_detect(population, "^CD8 T_") ~ "CD8",
           str_detect(population, "^gd T_") ~ "gd",
           str_detect(population, "^NK_") ~ "NK",
           str_detect(population, "^Myeloid_") ~ "Myeloid",
           str_detect(population, "^B_") ~ "B",
           TRUE ~ "Other"),
         CleanName = population %>%
           str_remove("^CD4 T_") %>% str_remove("^CD8 T_") %>% str_remove("^gd T_") %>%
           str_remove("^NK_") %>% str_remove("^Myeloid_") %>% str_remove("^B_"),
         CleanName = ifelse(population == "B_CD5+ B cells", "CD5+ B cells", paste(Cluster, CleanName)))
final_order <- vp_long %>%
  filter(Variable == "Age") %>%
  group_by(CleanName) %>%
  summarise(mean_var = mean(Variance, na.rm = TRUE), .groups = "drop") %>%
  arrange(mean_var) %>%
  pull(CleanName)
vp_bar <- vp_long %>%
  filter(Variance >= 1) %>%
  mutate(CleanName = factor(CleanName, levels = final_order), Variance_label = round(Variance, 1))
p_vp <- ggplot(vp_bar, aes(x = Variance, y = CleanName, fill = Cluster)) +
  geom_col(width = 0.85, color = "black", alpha = 0.7) +
  geom_text(aes(label = Variance_label), hjust = -0.2, size = 3) +
  facet_grid(. ~ Variable, scales = "free_x", space = "free") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.25))) +
  scale_fill_manual(values = c("B" = "#ff8899", "CD4" = "#eb883c", "CD8" = "#298283",
                               "gd" = "#6f7dc6", "Myeloid" = "#94755e", "NK" = "#9d6db8")) +
  theme_bw() +
  theme(strip.text = element_text(face = "bold"), axis.text.y = element_text(size = 8),
        axis.title = element_blank(), panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(fill = "Cluster")
ggsave(file.path(out, "fig2_variance_partition.png"), p_vp, width = 12, height = 10, dpi = 200, bg = "white")
write.csv(as.data.frame(vp_long), file.path(out, "fig2_variance_partition.csv"), row.names = FALSE)
cat("wrote: outputs/fig2_variance_partition_*.png/csv\n")

cat("fig2 done.\n")
