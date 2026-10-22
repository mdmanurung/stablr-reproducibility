# Paper Figure 3 / S3 — CD8 T-compartment composition PCA (Terekhova et al. 2026).
#
# Ported from the authors' `Figures_pipeline/figures_pipeline.R`
# (repo teresho4/Pan_atlas_immune_aging_TemK_TemB_ratio, section
# "Figure 3, S3", lines 843-961, plus the PCA panels at lines 880-960).
# CD8 minor-population percentages per donor -> PCA; scree plot, biplot,
# PC2 contribution, age-gradient coloring, Tem GZMK+/GZMB+ ratio
# coloring per dataset, and the per-population age-association summary.
#
# Run: ./run_fig3.sh

figdir <- "fig3_cd8_pca"
.args <- commandArgs(trailingOnly = FALSE); .fa <- grep("^--file=", .args, value = TRUE)
.root <- if (length(.fa)) dirname(dirname(normalizePath(sub("^--file=", "", .fa[1])))) else normalizePath(getwd())
source(file.path(.root, "fig_support.R"))
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(factoextra); library(ggrepel); library(tibble)
})

cells <- read_paper_cells()
minor_pct <- as.data.frame(paper_minor_percentages(cells))
donor_meta_df <- cells %>% select(donor_id, age, sex, ethnicity, CMV, Dataset) %>% distinct()
out <- fig_output_dir(figdir)

cd8 <- c("RTE", "Naive", "Naive-IFN", "Tcm CCR4-", "Tcm CCR4+", "Trm", "Tem GZMK+",
         "Tem GZMB+", "Temra", "NKT-like", "Tmem KLRC2+", "HLA-DR+", "Proliferative")
cd8_mat <- minor_pct %>%
  filter(Major == "CD8 T", Annotation %in% cd8) %>%
  select(donor_id, Annotation, percent) %>%
  pivot_wider(names_from = Annotation, values_from = percent, values_fill = 0)
t_percentages <- cd8_mat %>% column_to_rownames("donor_id")
pca <- prcomp(t_percentages, scale. = FALSE, center = TRUE)
pca_df <- as.data.frame(pca$x) %>%
  rownames_to_column("donor_id") %>%
  left_join(
    cells %>% select(donor_id, sex, age, ethnicity, Dataset, CMV) %>% distinct(),
    by = "donor_id"
  ) %>%
  left_join(cd8_mat, by = "donor_id")

pc_lab <- function(pca, i) {
  paste0("PC", i, " (", round(summary(pca)$importance[2, i] * 100, 1), "% variance)")
}

# ---- Scree / biplot / contribution (factoextra, authors' calls) -----------------
ggsave(file.path(out, "fig3_pca_scree.png"),
       fviz_eig(pca) + ggtitle("% of CD8 T cell populations for each individual"),
       width = 8, height = 5, dpi = 200, bg = "white")
ggsave(file.path(out, "fig3_pca_biplot.png"),
       fviz_pca_biplot(pca, geom = "point", col.var = "#1F78B4", col.ind = "#bdbbbb") +
         ggtitle("Loadings of CD8 T cell subpopulations"),
       width = 8, height = 6, dpi = 200, bg = "white")
ggsave(file.path(out, "fig3_pca_pc2_contribution.png"),
       fviz_contrib(pca, choice = "var", axes = 2, top = 5) +
         ggtitle("Contribution of variables to PC2"),
       width = 8, height = 5, dpi = 200, bg = "white")
cat("wrote: outputs/fig3_pca_*.png\n")

# ---- Age-gradient PCA (verbatim) -------------------------------------------------
p_age <- ggplot(pca_df, aes(PC1, PC2, color = Naive)) +
  geom_point(size = 1.5, alpha = 0.8) +
  theme_minimal() +
  labs(x = pc_lab(pca, 1), y = pc_lab(pca, 2)) +
  theme(legend.position = "bottom") +
  scale_color_gradient2(low = "#569cb4", mid = "#b88084", high = "#ff001e", midpoint = 40)
ggsave(file.path(out, "fig3_pca_age_gradient.png"), p_age, width = 8, height = 6, dpi = 200, bg = "white")

# ---- Tem GZMK+/GZMB+ ratio coloring, faceted by dataset (verbatim) ---------------
global_hull <- pca_df %>% slice(chull(PC1, PC2))
global_hull_expanded <- bind_rows(lapply(unique(pca_df$Dataset), function(ds) {
  global_hull %>% mutate(Dataset = ds)
}))
pca_df <- pca_df %>% mutate(emK_emB = `Tem GZMK+` / `Tem GZMB+`)
pca_df[is.infinite(pca_df$emK_emB), "emK_emB"] <-
  max(pca_df[!is.infinite(pca_df$emK_emB), "emK_emB"], na.rm = TRUE)
p_ratio <- ggplot(pca_df, aes(x = PC1 * -1, y = PC2)) +
  geom_polygon(data = global_hull_expanded, aes(group = Dataset), fill = NA,
               color = "#696969", linewidth = 0.6) +
  geom_point(aes(color = log1p(emK_emB)), size = 1.5, alpha = 0.8) +
  scale_color_gradient2(low = "#dadfe0", mid = "#f27466", high = "#c71602", midpoint = 3) +
  facet_wrap(~Dataset, ncol = 4) +
  labs(x = pc_lab(pca, 1), y = pc_lab(pca, 2)) +
  theme_minimal()
ggsave(file.path(out, "fig3_pca_ratio_faceted.png"), p_ratio, width = 14, height = 7, dpi = 200, bg = "white")
cat("wrote: outputs/fig3_pca_age_gradient.png, fig3_pca_ratio_faceted.png\n")

# ---- Per-population age-association summary (verbatim; needs minor models) -------
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
minor_with_meta <- merge(minor_pct, donor_meta_df %>% select(donor_id, age, sex, ethnicity), by = "donor_id")
model_minor_results <- calculate_age_associations_minor(minor_with_meta)

model_cd8_summary <- model_minor_results %>%
  filter(Major == "CD8 T") %>%
  group_by(Annotation) %>%
  summarise(mean_r = mean(r, na.rm = TRUE),
            mean_logp = mean(-log10(p.adj), na.rm = TRUE),
            n_sig = sum(p.adj < 0.05, na.rm = TRUE), .groups = "drop")
p_sum <- ggplot(model_cd8_summary, aes(x = mean_r, y = mean_logp)) +
  geom_point(aes(size = n_sig, color = mean_r), alpha = 0.8) +
  geom_text_repel(data = model_cd8_summary %>% filter(n_sig > 0), aes(label = Annotation), size = 3) +
  scale_color_gradient2(low = "blue", mid = "grey85", high = "red", midpoint = 0) +
  scale_size_continuous(name = "N significant datasets", range = c(1.5, 6)) +
  theme_bw() +
  labs(x = "Mean Pearson r across datasets", y = "Mean -log10(adj. p)")
ggsave(file.path(out, "fig3_cd8_association_summary.png"), p_sum, width = 8, height = 6, dpi = 200, bg = "white")
write.csv(as.data.frame(model_cd8_summary), file.path(out, "fig3_cd8_association_summary.csv"), row.names = FALSE)
write.csv(as.data.frame(pca_df[, c("donor_id", "PC1", "PC2", "age", "Dataset", "CMV", "emK_emB")]),
          file.path(out, "fig3_pca_coordinates.csv"), row.names = FALSE)
cat("wrote: outputs/fig3_cd8_association_summary.*, fig3_pca_coordinates.csv\n")

cat("fig3 done.\n")
