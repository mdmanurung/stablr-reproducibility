# Figure 2 / S2 / S3 — Cell-type age associations

Reproduction of **Figure 2** (with supplementary panels) of
Terekhova et al. (2026), *Immunity*, DOI 10.1016/j.immuni.2026.09.012:
per-population associations of cell-type composition with age across the
8 atlas datasets.

## What is reproduced

- Per (dataset, population) Pearson correlation of cell-type percentage
  with age, and `lm(percent ~ Age + Sex + Ethnicity)` (Ethnicity dropped
  where it has a single level) with Benjamini–Hochberg adjustment per
  dataset — for major populations and minor populations.
- The association heatmap + effect-size boxplot panel.
- Major-population percentage boxplots by age and the major-composition
  alluvial.
- CD8 minor-population percentage-vs-age scatters (Tem GZMK+, Tem GZMB+).
- The oldest-old comparison: age bands A–F (<35, <45, <55, <65, <75, ≥75)
  with Dunn tests.
- Variance partitioning of minor populations (variancePartition,
  `fitExtractVarPartModel` with donor random effects).

## Outputs

| Output | Content |
|---|---|
| `fig2_age_association_summary.png` | Heatmap + boxplot of age associations (major populations) |
| `fig2_age_association_models.csv` | Full model table: r, lm coefficients, BH-adjusted p per dataset × major population |
| `fig2_minor_age_models.csv` | Same model table for minor populations |
| `fig2_major_percent_boxplots.png` | Major-population percentages by age |
| `fig2_major_percent_alluvial.png` | Major-composition alluvial by dataset |
| `fig2_cd8_tem_gzmk_scatter.png` | CD8 Tem GZMK+ % vs age, per dataset |
| `fig2_cd8_tem_gzmb_scatter.png` | CD8 Tem GZMB+ % vs age, per dataset |
| `fig2_oldest_old_bands.png` | Oldest-old age-band comparison |
| `fig2_oldest_old_dunn.csv` | Dunn test results for the oldest-old bands |
| `fig2_variance_partition.png` | variancePartition barplot of minor populations (Age / Ethnicity / Sex shares) |
| `fig2_variance_partition.csv` | Variance fractions per population and covariate |

The variance-partitioning panel reproduces the paper's headline structure:
naive populations carry the largest age share (CD8 Naive ~37%, CD8/CD4 RTE
~18–19%, Tem GZMK+/GZMB+ ~16%), with cohort/ethnicity and sex contributing
far less for most populations.

## Provenance

Ported from the authors' `Figures_pipeline/figures_pipeline.R`
(repo `teresho4/Pan_atlas_immune_aging_TemK_TemB_ratio`), section
"Figure 2, S2, S3" (lines 183–842). Model formulas, significance codes
(`+++`/`++`/`+`/`−`/`−−`/`−−−` at .001/.01/.05), dataset ordering, and
palettes follow the authors' version; the data source is this
compendium's MD5-pinned raw per-cell table via `fig_support.R`
(`paper_major_percentages` / `paper_minor_percentages` reproduce the
authors' percentage semantics: % of typed cells for majors, % of the
parent major for minors, full donor × population grid with zeros).

## Run

```sh
./run_fig2.sh
```

Runtime ~4–6 min (13M-row read + per-dataset model fits + variance
partitioning).
