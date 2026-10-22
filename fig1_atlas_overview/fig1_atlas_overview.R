# Paper Figure 1 / S1 — Atlas overview (Terekhova et al. 2026, Immunity).
#
# Ported from the authors' `Figures_pipeline/figures_pipeline.R`
# (repo teresho4/Pan_atlas_immune_aging_TemK_TemB_ratio, section
# "Figure 1, S1", lines 32-182). Panel functions are the authors' code
# essentially verbatim; the data source is this compendium's MD5-pinned
# raw per-cell table (identical schema to the authors'
# ALL_DATASETS_ANNOTATION.csv). Aggregation for the large per-cell
# counts uses data.table (fig_support.R) for speed; plotted semantics
# are unchanged.
#
# Panels: donor sex by dataset; total cells by dataset; age histograms
# by dataset; ancestry distribution; Bad-quality cell fraction;
# annotation alluvial (harmonized vs author cell-type labels).
#
# Run: ./run_fig1.sh   (or Rscript fig1_atlas_overview/fig1_atlas_overview.R)

figdir <- "fig1_atlas_overview"
.args <- commandArgs(trailingOnly = FALSE); .fa <- grep("^--file=", .args, value = TRUE)
.root <- if (length(.fa)) dirname(dirname(normalizePath(sub("^--file=", "", .fa[1])))) else normalizePath(getwd())
source(file.path(.root, "fig_support.R"))
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(patchwork)
  library(ggbreak); library(ggalluvial)
})

cells <- read_paper_cells()
cat(sprintf("cells: %d rows; donors: %s\n", nrow(cells),
            format(uniqueN(cells$donor_id), big.mark = ",")))

# ---- Authors' panel functions (verbatim) --------------------------------------
# PORT NOTE: the authors' donor-sex panel uses ggpattern::geom_col_pattern
# (bars pattern-filled by sex). ggpattern's dependency chain (sf/units) is
# not installable here, so the panel is rendered as sex-dodged bars — same
# counts, same axis, different fill encoding.
plot_donor_sex_by_dataset <- function(data) {
  donor_n_text <- data %>%
    distinct(Dataset, donor_id, sex) %>%
    count(Dataset, sex) %>%
    mutate(label = paste0("n = ", n))
  donor_n_text$Dataset <- factor(donor_n_text$Dataset, levels = PAPER_DATASETS)
  ggplot(donor_n_text, aes(x = Dataset, y = n, fill = sex)) +
    geom_col(position = position_dodge(width = 0.9), color = "black") +
    geom_text(aes(label = label, group = sex), position = position_dodge(width = 0.9),
              vjust = -0.3, size = 2.4) +
    theme_minimal() + ylab("") + xlab("") +
    theme(axis.text.x = element_text(colour = "black"),
          axis.text.y = element_text(colour = "black")) +
    scale_fill_manual(values = c(Female = "#BF4294", Male = "#6BA6A8"))
}

plot_total_cells_by_dataset <- function(data) {
  n_cells <- as.data.frame(table(data$Dataset))
  colnames(n_cells) <- c("Dataset", "N_cells")
  n_cells$Dataset <- factor(n_cells$Dataset, levels = PAPER_DATASETS)
  ggplot(n_cells, aes(x = Dataset, y = N_cells, fill = Dataset)) +
    geom_col(color = "black") +
    scale_y_break(c(2000000, 5000000), scales = 1) +
    scale_fill_manual(values = PAPER_DATASET_COLORS) +
    theme_minimal(base_size = 11) +
    labs(x = NULL, y = NULL) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
          axis.text.y = element_text(color = "black"),
          legend.key.size = unit(0.7, "cm"))
}

plot_age_histograms_by_dataset <- function(data) {
  df <- data %>% distinct(Dataset, donor_id, age)
  x_limits <- c(18, 99)
  make_hist <- function(dataset_name) {
    ggplot(df %>% filter(Dataset == dataset_name), aes(x = age)) +
      geom_histogram(binwidth = 3, fill = unname(PAPER_DATASET_COLORS[[dataset_name]]),
                     color = "black", boundary = 0, closed = "left") +
      geom_density(aes(y = after_stat(count * 3)),
                   color = unname(PAPER_DENSITY_COLORS[[dataset_name]]),
                   linewidth = 1.1) +
      scale_x_continuous(limits = x_limits, breaks = seq(18, 98, 5)) +
      theme_minimal(base_size = 11) +
      labs(x = NULL, y = NULL, title = dataset_name)
  }
  plots <- lapply(names(PAPER_DATASET_COLORS), make_hist)
  patchwork::wrap_plots(plots, ncol = 1)
}

plot_ethnicity_distribution <- function(data) {
  df <- data %>% distinct(donor_id, ethnicity) %>% count(ethnicity, name = "Count") %>%
    rename(Ethnicity = ethnicity)
  ggplot(df, aes(x = Ethnicity, y = Count, fill = Ethnicity)) +
    geom_col(color = "black") +
    scale_y_break(c(50, 600), scales = 1) +
    theme_minimal(base_size = 11) +
    labs(title = "Ancestry Distribution", x = NULL, y = "Count") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
          legend.position = "none")
}

plot_bad_quality <- function(data) {
  df_clean <- data %>%
    mutate(quality = ifelse(major_celltypes %in% c("Bad quality"), "Bad quality cells", "Rest"))
  df_percent <- df_clean %>%
    group_by(Dataset, quality) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(Dataset) %>%
    mutate(percentage = n / sum(n) * 100)
  ggplot(df_percent, aes(x = Dataset, y = percentage, fill = quality)) +
    geom_col(position = "stack", color = "black") +
    ylab("Percentage") +
    scale_fill_manual(values = c("Rest" = "#00bf7dff", "Bad quality cells" = "#ffe680ff")) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(colour = "black", angle = 45, hjust = 1))
}

plot_annotation_alluvial <- function(data, celltype_colors) {
  df <- data %>% filter(major_celltypes != "Bad quality")
  df <- df %>% filter(Dataset %in% c("OneK1K", "AIDA", "CIMA", "UCSF", "SoundLife", "SPAC"))
  df_count <- df %>% count(author_celltype, minor_celltypes)
  mapping <- df_count %>%
    group_by(author_celltype) %>%
    slice_max(order_by = n, n = 1, with_ties = FALSE) %>%
    ungroup()
  ordered_author <- mapping %>% arrange(minor_celltypes) %>% pull(author_celltype)
  df$author_celltype <- factor(df$author_celltype, levels = ordered_author)
  df_alluv <- df %>% count(Dataset, minor_celltypes, author_celltype) %>% filter(n > 10)
  ggplot(df_alluv, aes(axis1 = minor_celltypes, axis2 = author_celltype, y = n)) +
    geom_alluvium(aes(fill = minor_celltypes), width = 1 / 12) +
    geom_stratum(width = 1 / 12, fill = "gray80", color = "black") +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 2.5) +
    scale_x_discrete(limits = c("New", "Author cell type"), expand = c(0.05, 0.05)) +
    facet_wrap(~Dataset, scales = "free") +
    theme_minimal(base_size = 11) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          panel.grid = element_blank(), strip.text = element_text(size = 10),
          axis.text.x = element_text(size = 10, color = "black")) +
    labs(y = "Cell count", x = NULL)
}

# ---- Render panels + composed figure ------------------------------------------
out <- fig_output_dir(figdir)

p_sex    <- plot_donor_sex_by_dataset(cells)
p_cells  <- plot_total_cells_by_dataset(cells)
p_age    <- plot_age_histograms_by_dataset(cells)
p_eth    <- plot_ethnicity_distribution(cells)
p_bad    <- plot_bad_quality(cells)
p_alluv  <- plot_annotation_alluvial(cells)

ggsave(file.path(out, "fig1_panel_donor_sex.png"), p_sex, width = 8, height = 4, dpi = 200, bg = "white")
ggsave(file.path(out, "fig1_panel_total_cells.png"), p_cells, width = 8, height = 4.5, dpi = 200, bg = "white")
ggsave(file.path(out, "fig1_panel_age_histograms.png"), p_age, width = 6, height = 12, dpi = 200, bg = "white")
ggsave(file.path(out, "fig1_panel_ancestry.png"), p_eth, width = 7, height = 5, dpi = 200, bg = "white")
ggsave(file.path(out, "fig1_panel_bad_quality.png"), p_bad, width = 8, height = 4.5, dpi = 200, bg = "white")
ggsave(file.path(out, "fig1_panel_annotation_alluvial.png"), p_alluv, width = 14, height = 6, dpi = 200, bg = "white")
cat("wrote 6 panel PNGs to", out, "\n")

# Cohort summary stats behind the panels
stats <- cells %>% distinct(Dataset, donor_id, sex, age, ethnicity) %>%
  group_by(Dataset) %>%
  summarise(donors = n(), females = sum(sex == "Female"), males = sum(sex == "Male"),
            age_min = min(age, na.rm = TRUE), age_max = max(age, na.rm = TRUE),
            median_age = median(age, na.rm = TRUE), .groups = "drop")
cells_by_ds <- cells[, .(cells = .N), by = Dataset]
stats <- stats %>% left_join(as.data.frame(cells_by_ds), by = "Dataset")
write.csv(stats, file.path(out, "fig1_cohort_summary.csv"), row.names = FALSE)
cat("wrote: outputs/fig1_cohort_summary.csv\n")

cat("fig1 done.\n")
