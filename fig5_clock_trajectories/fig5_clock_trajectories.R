# Paper Figure 5 / S5 — Aging-clock and health-trajectory components
# (Terekhova et al. 2026).
#
# Two of the paper's Figure-5 components are reproduced here from the
# compendium's committed donor tables:
#
#   1. The age-predicting clock (authors' `Age_clock_model/
#      clean_notebook_for_aging_models.ipynb`): ridge regression on the
#      20-population reduced panel (`selected_by_marina`, notebook cell 32)
#      with one-hot sex/ethnicity covariates, evaluated leave-one-cohort-out.
#      The paper reports MAE 8.2 y (full panel) / 8.9 y (reduced panel);
#      this script recomputes both under the compendium's LOCO engine
#      (cv.glmnet ridge, alpha = 0, seed 20261013).
#   2. The all-feature vs reduced-panel comparison behind the clock figure.
#
# NOT ported (data not shipped in this compendium): the Olink proteomics
# panel (`Olink_proteins_all_cohorts.csv`) and the UK Biobank protein-model
# survival analysis (authors' `UKBB_protein_model/`). See the README here.
#
# Run: ./run_fig5.sh

figdir <- "fig5_clock_trajectories"
.args <- commandArgs(trailingOnly = FALSE); .fa <- grep("^--file=", .args, value = TRUE)
.root <- if (length(.fa)) dirname(dirname(normalizePath(sub("^--file=", "", .fa[1])))) else normalizePath(getwd())
source(file.path(.root, "fig_support.R"))
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(glmnet)
})

set.seed(20261013)
out <- fig_output_dir(figdir)

d <- paper_harmonized_matrix()
X_all <- d$X
meta <- d$meta
y_age <- meta$age
env <- meta$dataset
cat(sprintf("donors: %d; features: %d; cohorts: %s\n",
            nrow(X_all), ncol(X_all), paste(d$cohorts, collapse = ", ")))

# ---- The paper's reduced 20-population panel (notebook cell 32) ----------------
# Authors' names -> compendium harmonized names (verified 1:1 by name
# inspection; the compendium's 59 features are the `annotation` levels).
P20_MAP <- list(
  "CD4 RTE"                = "CD4+ RTE",
  "CD4 Th2"                = "CD4+ Th2",
  "CD4 HLA-DR+ memory"     = "CD4+ T HLA-DR+ memory",
  "CD4 Treg memory"        = "CD4+ Treg memory",   # authors' " CD4 Treg memory"
  "CD8 RTE"                = "CD8+ RTE",
  "CD8 Naive"              = "CD8+ T Naive",
  "CD8 Tcm CCR4+"          = "CD8+ Tcm CCR4+",
  "CD8 Tem GZMK+"          = "CD8+ Tem GZMK+",
  "CD8 Tem GZMB+"          = "CD8+ Tem GZMB+",
  "CD8 HLA-DR+"            = "CD8+ HLA-DR+",
  "CD8 Tmem KLRC2+"        = "CD8+ Tmem KLRC2+",
  "gd gd naive"            = "gd T naive",
  "gd Vd1 GZMK+"           = "gd T Vd1 GZMK+",
  "gd Vd1 GZMB+"           = "gd T Vd1 GZMB+",
  "gd Vd2 GZMK+"           = "gd T Vd2 GZMK+",
  "NK CD56bright"          = "NK CD56bright",
  "Non-classical monocytes" = "Myeloid Non-classical monocytes",
  "CD8 Temra"              = "CD8+ Temra",
  "CD4 Temra"              = "CD4+ Temra",
  "CD4 Terminal effector"  = "CD4+ Terminal effector"
)
p20 <- unlist(P20_MAP)
missing <- p20[!p20 %in% colnames(X_all)]
if (length(missing)) stop("P20 populations not found in harmonized features: ",
                          paste(missing, collapse = ", "))
X_p20 <- X_all[, p20, drop = FALSE]

# One-hot sex/ethnicity covariates (authors' get_dummies, drop_first)
onehot <- function(v) {
  lev <- sort(unique(na.omit(v)))[-1]
  sapply(lev, function(l) as.numeric(v == l))
}
cov_p20 <- cbind(onehot(meta$sex), onehot(meta$ethnicity))
colnames(cov_p20) <- paste0("COV:", colnames(cov_p20))
X_p20_cov <- cbind(X_p20, cov_p20)

# ---- LOCO ridge clock (compendium engine: cv.glmnet alpha = 0) -----------------
loco_clock <- function(Xt, yt, ev) {
  cohs <- sort(unique(ev))
  do.call(rbind, lapply(cohs, function(h) {
    tr <- ev != h; te <- ev == h
    cvfit <- cv.glmnet(Xt[tr, , drop = FALSE], yt[tr], alpha = 0, nfolds = 5)
    pr <- as.numeric(predict(cvfit, newx = Xt[te, , drop = FALSE], s = "lambda.min"))
    data.frame(cohort = h, n = sum(te), pred = pr, obs = yt[te],
               mae = mean(abs(pr - yt[te])))
  }))
}

cat("== LOCO ridge clock: all 59 features ==\n")
loco_all <- loco_clock(X_all, y_age, env)
cat("== LOCO ridge clock: P20 reduced panel ==\n")
loco_p20 <- loco_clock(X_p20, y_age, env)
cat("== LOCO ridge clock: P20 + sex/ethnicity (authors' spec) ==\n")
loco_p20_cov <- loco_clock(X_p20_cov, y_age, env)

mae_all <- mean(loco_all$mae); mae_p20 <- mean(loco_p20$mae); mae_p20c <- mean(loco_p20_cov$mae)
cat(sprintf("LOCO MAE: all-59 %.2f y | P20 %.2f y | P20+demog %.2f y (paper: 8.2 full / 8.9 reduced)\n",
            mae_all, mae_p20, mae_p20c))

# ---- Figures --------------------------------------------------------------------
plot_obs_pred <- function(res, title) {
  res$cohort <- factor(res$cohort, levels = sort(unique(res$cohort)))
  ggplot(res, aes(x = obs, y = pred, color = cohort)) +
    geom_point(size = 1.2, alpha = 0.6) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
    facet_wrap(~cohort) +
    theme_minimal(base_size = 11) +
    coord_fixed(xlim = c(18, 100), ylim = c(18, 100)) +
    labs(title = title, x = "Chronological age (y)", y = "Predicted age (y)", color = "Cohort") +
    theme(legend.position = "none", strip.text = element_text(face = "bold"))
}
ggsave(file.path(out, "fig5_clock_all59_loco.png"),
       plot_obs_pred(loco_all, sprintf("All-59 ridge clock — LOCO MAE %.1f y (paper: 8.2 y)", mae_all)),
       width = 12, height = 9, dpi = 200, bg = "white")
ggsave(file.path(out, "fig5_clock_p20_loco.png"),
       plot_obs_pred(loco_p20, sprintf("P20 reduced-panel ridge clock — LOCO MAE %.1f y (paper: 8.9 y)", mae_p20)),
       width = 12, height = 9, dpi = 200, bg = "white")
ggsave(file.path(out, "fig5_clock_p20_demog_loco.png"),
       plot_obs_pred(loco_p20_cov, sprintf("P20 + sex/ethnicity ridge clock — LOCO MAE %.1f y", mae_p20c)),
       width = 12, height = 9, dpi = 200, bg = "white")

mae_tab <- bind_rows(
  data.frame(model = "all-59", loco_all, paper_mae_ref = 8.2),
  data.frame(model = "P20", loco_p20, paper_mae_ref = 8.9),
  data.frame(model = "P20+demog", loco_p20_cov, paper_mae_ref = NA)
)
write.csv(mae_tab, file.path(out, "fig5_loco_mae.csv"), row.names = FALSE)

mae_sum <- mae_tab %>%
  group_by(model) %>%
  summarise(mean_mae = mean(mae), .groups = "drop") %>%
  mutate(paper_ref = c(8.2, 8.9, NA)[match(model, c("all-59", "P20", "P20+demog"))])
p_mae <- ggplot(mae_sum, aes(x = model, y = mean_mae, fill = model)) +
  geom_col(color = "black", width = 0.6) +
  geom_hline(yintercept = 8.2, linetype = "dashed", color = "grey30") +
  geom_hline(yintercept = 8.9, linetype = "dotted", color = "grey50") +
  geom_text(aes(label = sprintf("%.1f", mean_mae)), vjust = -0.4, size = 4) +
  scale_fill_manual(values = c("#4477AA", "#EE6677", "#228833")) +
  theme_minimal(base_size = 12) + theme(legend.position = "none") +
  labs(x = NULL, y = "LOCO MAE (years)",
       title = "Age-clock accuracy vs the paper's reported clock",
       subtitle = "dashed: paper full-panel 8.2 y; dotted: paper reduced-panel 8.9 y")
ggsave(file.path(out, "fig5_clock_mae_comparison.png"), p_mae, width = 7, height = 5, dpi = 200, bg = "white")
cat("wrote: outputs/fig5_clock_*.png, fig5_loco_mae.csv\n")

cat("fig5 done. (Olink proteomics and UK Biobank survival panels not ported: data not shipped; see README.)\n")
