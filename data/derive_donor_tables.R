# derive_donor_tables.R — Rebuild the v1 donor-level tables from the raw
# portal file.
#
# Reconstructs `data/derived/donor_celltype_counts.tsv` and
# `data/derived/donor_meta.tsv` from
# `data/raw/all_cells_metadata.csv.gz` (Terekhova et al., Immunity 2026;
# MD5-pinned, see `00_download_data.sh`).  The original v1 derivation
# script was never archived; this script is the reconstructed, verified
# replacement.  It applies the same `Bad quality` annotation exclusion
# that the committed tables embody (verified: zero `Bad quality` rows in
# the committed counts) and writes the exact v1 file formats consumed by
# `stablr/benchmarks/ws8_aging_showcase.R`:
#
#   donor_celltype_counts.tsv:  "<donor>|<celltype>\t<n>"   (no header)
#   donor_meta.tsv:             "META\t<donor>\t<Dataset|age|sex|ethnicity|CMV>"
#                               (no header)
#
# Rows are written in canonical sorted order (donor, then cell type).
# The original traversal order was not archived and is not reconstructible;
# the committed files differ from this script's output in ROW ORDER ONLY.
# Content is identical — verified by `--verify`, which order-normalizes
# both sides and diffs them (0 mismatches, 0 unmatched rows, on the
# MD5-pinned raw file).  All consumers (`ws8_aging_showcase.R` and the v2
# pipeline) read these tables order-agnostically via keyed matching.
#
# Usage:
#   Rscript data/derive_donor_tables.R            # write + verify
#   Rscript data/derive_donor_tables.R --verify   # verify only
#
# Run from the compendium root.  Requires data.table.

suppressPackageStartupMessages(library(data.table))

verify_only <- "--verify" %in% commandArgs(trailingOnly = TRUE)

raw_file <- file.path("data", "raw", "all_cells_metadata.csv.gz")
counts_out <- file.path("data", "derived", "donor_celltype_counts.tsv")
meta_out   <- file.path("data", "derived", "donor_meta.tsv")
stopifnot(file.exists(raw_file))

expected_md5 <- "50651968855bf6ca3efcdf199e2a9fd5"
md5 <- unname(tools::md5sum(raw_file))
if (md5 != expected_md5) {
  stop("MD5 mismatch for ", raw_file, ": expected ", expected_md5,
       ", got ", md5, call. = FALSE)
}
cat(sprintf("input: %s\nmd5:   %s (OK)\n", raw_file, md5))

keep_cols <- c("donor_id", "age", "sex", "ethnicity", "CMV", "Dataset",
               "annotation")
dt <- fread(raw_file, select = keep_cols, showProgress = FALSE)
cat(sprintf("cell rows read: %d\n", nrow(dt)))

# ---- Donor x cell-type counts (Bad quality excluded) -------------------------
counts <- dt[!grepl("Bad quality", annotation, fixed = TRUE),
             .(n = .N), by = .(donor_id, annotation)]
setorderv(counts, c("donor_id", "annotation"))
cat(sprintf("donor x cell-type pairs: %d (Bad-quality rows excluded)\n",
            nrow(counts)))

# ---- One metadata row per donor ----------------------------------------------
n_combo <- dt[, .(n_combo = uniqueN(paste(age, sex, ethnicity, CMV,
                                          Dataset))),
              by = donor_id]
n_conflicts <- sum(n_combo$n_combo > 1L)
if (n_conflicts > 0L) {
  stop(n_conflicts, " donors have inconsistent metadata rows; ",
       "refusing to pick arbitrarily.", call. = FALSE)
}
meta <- unique(dt[, .(donor_id, Dataset, age, sex, ethnicity, CMV)],
               by = "donor_id")
setorderv(meta, "donor_id")
cat(sprintf("donors: %d; metadata conflicts: %d\n", nrow(meta),
            n_conflicts))
rm(dt); gc(verbose = FALSE)

# ---- Writers (exact v1 formats) ----------------------------------------------
write_counts <- function(ct, path) {
  fwrite(data.frame(key = paste(ct$donor_id, ct$annotation, sep = "|"),
                    n = ct$n),
         path, sep = "\t", col.names = FALSE)
}
write_meta <- function(m, path) {
  fwrite(data.frame(tag = "META", donor = m$donor_id,
                    field = paste(m$Dataset, m$age, m$sex, m$ethnicity,
                                  m$CMV, sep = "|")),
         path, sep = "\t", col.names = FALSE)
}

# ---- Verification (order-normalized content diff) -----------------------------
read_counts <- function(path) {
  x <- fread(path, header = FALSE, sep = "\t", colClasses = c("character",
                                                              "integer"))
  setnames(x, c("key", "n"))
  setkey(x, key)
  x
}
read_meta <- function(path) {
  x <- fread(path, header = FALSE, sep = "\t", colClasses = "character")
  setnames(x, c("tag", "donor", "field"))
  setkey(x, donor)
  x
}

verify <- function() {
  tmp <- tempfile()
  on.exit(unlink(tmp), add = TRUE)
  write_counts(counts, tmp)
  same_counts <- identical(read_counts(tmp), read_counts(counts_out))
  write_meta(meta, tmp)
  same_meta <- identical(read_meta(tmp), read_meta(meta_out))
  cat(sprintf("verify: counts content %s; meta content %s\n",
              if (same_counts) "IDENTICAL" else "DIFFERS",
              if (same_meta) "IDENTICAL" else "DIFFERS"))
  if (!same_counts || !same_meta) {
    stop("Regenerated tables differ in content from the committed ones; ",
         "inspect before overwriting.", call. = FALSE)
  }
  invisible(TRUE)
}

if (verify_only) {
  verify()
} else {
  if (file.exists(counts_out) && file.exists(meta_out)) {
    verify()   # prove content identity before touching the committed files
  }
  write_counts(counts, counts_out)
  write_meta(meta, meta_out)
  cat(sprintf("wrote %s and %s (canonical sorted row order; content\n",
              counts_out, meta_out))
  cat("identical to the previously committed files — see this script's\n")
  cat("header and data/README.md).\n")
}
