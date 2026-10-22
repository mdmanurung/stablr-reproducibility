# Figure 1 / S1 — Atlas overview

Reproduction of **Figure 1** (and its supplementary counterpart) of
Terekhova et al. (2026), *Immunity*, DOI 10.1016/j.immuni.2026.09.012:
cohort composition of the 8-dataset healthy-aging single-cell atlas
(13.1M cells, 2,609 donors).

## Panels

| Output | Paper panel |
|---|---|
| `fig1_panel_donor_sex.png` | Donor sex split per dataset |
| `fig1_panel_total_cells.png` | Total cells per dataset |
| `fig1_panel_age_histograms.png` | Age distribution per dataset |
| `fig1_panel_ancestry.png` | Donor ancestry per dataset |
| `fig1_panel_annotation_alluvial.png` | Major cell-type composition per dataset (alluvial) |
| `fig1_panel_bad_quality.png` | Fraction of bad-quality cells per dataset |
| `fig1_cohort_summary.csv` | Cohort summary table (donors, cells, age range per dataset) |

## Provenance

Ported from the authors' `Figures_pipeline/figures_pipeline.R`
(repo `teresho4/Pan_atlas_immune_aging_TemK_TemB_ratio`), section
"Figure 1, S1" (lines 32–182). Panel code is the authors' essentially
verbatim; the data source is this compendium's MD5-pinned raw per-cell
table (`data/raw/all_cells_metadata.csv.gz`), which carries the authors'
columns (`Dataset`, `age`, `sex`, `ethnicity`, `major_celltypes`,
`minor_celltypes`, `author_celltype`, `BMI`, `CMV`) unchanged. Aggregation
uses `data.table` (via `fig_support.R`) for speed on 13M rows; plotted
semantics are unchanged.

## Port notes

- **ggpattern substitution**: the authors' sex panel uses
  `ggpattern::geom_col_pattern` (pattern fills for sex). `ggpattern` needs
  the `sf`/`units` system libraries, which are not available in this
  environment, so the panel is drawn as sex-dodged bars with the same
  data, ordering, and palette.

## Run

```sh
./run_fig1.sh        # downloads data if missing, then Rscript fig1_atlas_overview.R
```

Runtime ~2–3 min (dominated by the single 88 MB per-cell table read).
