# Shared data layer for the paper-figure reproduction scripts (fig1-fig5).
#
# Source data: Terekhova et al. (2026), Immunity, DOI 10.1016/j.immuni.2026.09.012
# — raw per-cell annotation table `data/raw/all_cells_metadata.csv.gz`
# (13.1M cells, MD5 50651968855bf6ca3efcdf199e2a9fd5), fetched by
# ./00_download_data.sh.
#
# Percentage semantics ported EXACTLY from the authors'
# `Figures_pipeline/figures_pipeline.R` (repo teresho4/Pan_atlas_immune_aging_TemK_TemB_ratio):
#   - major populations: % of all typed (non-"Bad quality") cells per donor;
#   - minor populations: % of the corresponding major population per donor,
#     excluding majors c("Bad quality", "ILCs", "MAIT", "DN T", "Progenitor").
#
# Aggregation uses data.table for speed on the 13.1M-row table; the
# donor-level outputs are identical to the authors' dplyr version.

suppressPackageStartupMessages(library(data.table))

# ---- Authors' configuration (figures_pipeline.R lines 19-31) ------------------
PAPER_DATASETS <- c("OneK1K", "AIDA", "CIMA", "KCL", "ABF300", "UCSF", "SoundLife", "SPAC")
PAPER_DATASET_COLORS <- c("OneK1K" = "#6BA6A8", "AIDA" = "#E69A45", "SoundLife" = "#C4B454",
                          "SPAC" = "#A78CC1", "UCSF" = "#C15B58", "ABF300" = "#8FB98B",
                          "CIMA" = "#A68A74", "KCL" = "#BF4294")
PAPER_DENSITY_COLORS <- c("OneK1K" = "#487173", "AIDA" = "#c97e28", "SoundLife" = "#968833",
                          "SPAC" = "#5a3082", "UCSF" = "#9c322f", "ABF300" = "#31732a",
                          "KCL" = "#8c2a6a", "CIMA" = "#75573f")
PAPER_CELLTYPE_COLORS_MAJOR <- rev(c("#00A087", "#4DBBD5", "#4DBBD5", "#7CAE00", "#c7345d",
                                     "#8491B4", "#F4CAE4", "#E64B35", "#3C5488", "#a39b83",
                                     "#EEC900", "#B15928"))

# ---- Paths --------------------------------------------------------------------
repo_root <- function() {
  # Resolve from the invoking script path (Rscript figX/figX.R) or fall back
  # to the working directory; walk up until the compendium root is found.
  args <- commandArgs(trailingOnly = FALSE)
  fa <- grep("^--file=", args, value = TRUE)
  starts <- if (length(fa)) dirname(normalizePath(sub("^--file=", "", fa[1]))) else normalizePath(getwd())
  d <- starts
  while (!file.exists(file.path(d, "00_download_data.sh"))) {
    p <- dirname(d)
    if (p == d) stop("Compendium root (00_download_data.sh) not found from ", starts)
    d <- p
  }
  d
}

raw_atlas_path <- function(root = repo_root()) {
  p <- file.path(root, "data", "raw", "all_cells_metadata.csv.gz")
  if (!file.exists(p)) {
    stop("Raw atlas not found at ", p, " — run ./00_download_data.sh first.")
  }
  p
}

# ---- Raw per-cell table -------------------------------------------------------
read_paper_cells <- function(root = repo_root()) {
  cols <- c("donor_id", "age", "sex", "ethnicity", "CMV", "Dataset",
            "major_celltypes", "minor_celltypes", "author_celltype", "BMI")
  dt <- fread(raw_atlas_path(root), select = cols, showProgress = FALSE)
  dt[, age := as.numeric(age)]
  dt
}

# ---- One metadata row per donor -----------------------------------------------
paper_donor_meta <- function(cells) {
  meta <- unique(cells[, .(donor_id, age, sex, ethnicity, CMV, Dataset)])
  meta
}

# ---- Major-population percentages (authors' calculate_major_percentages) ------
paper_major_percentages <- function(cells) {
  df <- cells[major_celltypes != "Bad quality"]
  counts <- df[, .(n_cells = .N), by = .(Dataset, donor_id, major_celltypes)]
  totals <- df[, .(total_cells = .N), by = .(Dataset, donor_id)]
  grid <- unique(df[, .(Dataset, donor_id)])
  majors <- unique(df$major_celltypes)   # authors: global level set x every donor
  grid <- grid[, .(major_celltypes = majors), by = .(Dataset, donor_id)]
  out <- merge(grid, counts, by = c("Dataset", "donor_id", "major_celltypes"), all.x = TRUE)
  out <- merge(out, totals, by = c("Dataset", "donor_id"), all.x = TRUE)
  out[is.na(n_cells), n_cells := 0]
  out[, percent := 100 * n_cells / total_cells]
  setnames(out, "major_celltypes", "Annotation")[, .(Dataset, donor_id, Annotation, n_cells, total_cells, percent)]
}

# ---- Minor-population percentages (authors' calculate_minor_percentages) ------
paper_minor_percentages <- function(cells) {
  `%ni%` <- Negate(`%in%`)
  df <- cells[major_celltypes %ni% c("Bad quality", "ILCs", "MAIT", "DN T", "Progenitor") &
                minor_celltypes != "Bad quality"]
  setnames(df, "major_celltypes", "Major"); setnames(df, "minor_celltypes", "Annotation")
  counts <- df[, .(n_cells = .N), by = .(Dataset, donor_id, Major, Annotation)]
  totals <- df[, .(total_major_cells = .N), by = .(Dataset, donor_id, Major)]
  # full grid: every (donor, Major) crossed with every Annotation within that Major
  grid <- unique(df[, .(Dataset, donor_id, Major)])
  ann_by_major <- unique(df[, .(Major, Annotation)])
  grid <- grid[ann_by_major, on = "Major", allow.cartesian = TRUE]
  out <- merge(grid, counts, by = c("Dataset", "donor_id", "Major", "Annotation"), all.x = TRUE)
  out <- merge(out, totals, by = c("Dataset", "donor_id", "Major"), all.x = TRUE)
  out[is.na(n_cells), n_cells := 0]
  out[, percent := 100 * n_cells / total_major_cells]
  out[, .(Dataset, donor_id, Major, Annotation, n_cells, total_major_cells, percent)]
}

# ---- Donor x cell-type matrix in the compendium's harmonized space ------------
# Mirrors stablr/benchmarks/ws8_v2/data_prepare.R: Bad quality excluded from
# numerator and denominator, add-half smoothing, donor-level logit transform.
paper_harmonized_matrix <- function(root = repo_root()) {
  counts <- fread(file.path(root, "data", "derived", "donor_celltype_counts.tsv"),
                  header = FALSE, sep = "\t", col.names = c("key", "n"))
  counts[, donor := sub("\\|.*$", "", key)][, celltype := sub("^[^|]*\\|", "", key)]
  meta_raw <- fread(file.path(root, "data", "derived", "donor_meta.tsv"),
                    header = FALSE, sep = "\t", col.names = c("type", "donor", "info"))
  meta_raw <- meta_raw[type == "META"]
  parts <- tstrsplit(meta_raw$info, "|", fixed = TRUE)
  meta <- data.frame(donor = meta_raw$donor, dataset = parts[[1]], age = as.numeric(parts[[2]]),
                     sex = parts[[3]], ethnicity = parts[[4]], cmv = parts[[5]])
  ct_levels <- sort(unique(counts$celltype))
  donor_levels <- sort(intersect(meta$donor, unique(counts$donor)))
  counts <- counts[donor %in% donor_levels & celltype %in% ct_levels]
  ct_mat <- matrix(0L, nrow = length(donor_levels), ncol = length(ct_levels),
                   dimnames = list(donor_levels, ct_levels))
  idx <- cbind(match(counts$donor, donor_levels), match(counts$celltype, ct_levels))
  ct_mat[idx] <- counts$n
  total <- rowSums(ct_mat)
  prop <- (ct_mat + 0.5) / (total + 0.5)
  logit <- function(p) log(p / (1 - p))
  X <- logit(prop)
  meta <- meta[match(rownames(X), meta$donor), , drop = FALSE]
  list(X = X, meta = meta, cohorts = sort(unique(meta$dataset)))
}

# ---- Output helper ------------------------------------------------------------
fig_output_dir <- function(figdir) {
  d <- file.path(figdir, "outputs")
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  d
}

save_png <- function(p, filename, width = 10, height = 7, dpi = 200) {
  ggplot2::ggsave(file.path(fig_output_dir(figdir), filename), p,
                  width = width, height = height, dpi = dpi, bg = "white")
  cat("wrote:", file.path("outputs", filename), "\n")
}
