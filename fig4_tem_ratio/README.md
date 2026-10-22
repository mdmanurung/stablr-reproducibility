# Figure 4 / S4 — Tem GZMK+/GZMB+ (TemK/TemB) ratio biology

Reproduction of **Figure 4** of Terekhova et al. (2026), *Immunity*,
DOI 10.1016/j.immuni.2026.09.012: the paper's headline Tem GZMK+ /
Tem GZMB+ ratio analyses across sex, ancestry, BMI, age bands, and CMV
status.

## What is reproduced

- Per-donor ratio `n_cells(Tem GZMK+) / n_cells(Tem GZMB+)` within the
  CD8 T major compartment (infinite ratios clipped to the max finite
  value, plotted as `log1p`), across:
  - sex, ancestry (Asian/European), BMI category
    (<18.5 / <25 / <30 / ≥30, parsed as in the authors' code), and
    age bands A–F (<35, <45, <55, <65, <75, ≥75);
  - CMV status, faceted by dataset.
- The CMV-vs-population-change screen: within SoundLife + KCL
  (Young <40 / Old >50), Wilcoxon tests of each minor population's
  percentage by CMV status within CD4 T / CD8 T / B / NK / Myeloid
  majors, with `p.adj = p × 55` (authors' correction).
- The CMV-colored CD8 PCA panel.

## Not ported (data not shipped)

- **Aquino et al. external-validation panel** — the authors'
  `Aquino_immune_percentages.csv` is not distributed with this
  compendium and is not downloadable from a pinned archive; the panel is
  documented as a stub.

## Outputs

| Output | Content |
|---|---|
| `fig4_ratio_by_sex.png` | Ratio by donor sex |
| `fig4_ratio_by_ancestry.png` | Ratio by ancestry (Asian/European) |
| `fig4_ratio_by_bmi.png` | Ratio by BMI category |
| `fig4_ratio_by_age_band.png` | Ratio by age band A–F |
| `fig4_ratio_by_cmv_faceted.png` | Ratio by CMV status, faceted by dataset |
| `fig4_cmv_population_screen.csv` | Wilcoxon screen: population % change by CMV (p, p.adj) |
| `fig4_pca_cmv.png` | CD8 PCA colored by CMV status |

## Provenance

Ported from the authors' `Figures_pipeline/figures_pipeline.R`
(repo `teresho4/Pan_atlas_immune_aging_TemK_TemB_ratio`), section
"Figure 4, S4" (lines 962–1290). Ratio definition, BMI parsing
(`gsub(".0", "", BMI, fixed = TRUE)` then category cut), age-band
edges, the Young/Old CMV screen (SoundLife + KCL, <40 / >50,
`p.adj = p × 55`), and palettes follow the authors' version; the data
source is this compendium's MD5-pinned raw per-cell table via
`fig_support.R`. Donor-level covariates (age, sex, ethnicity, CMV, BMI)
are merged onto the percentage table from the same raw file (the
authors' pipeline carries them on every per-cell row upstream).

## Run

```sh
./run_fig4.sh
```

Runtime ~2–3 min (dominated by the raw-table read).
