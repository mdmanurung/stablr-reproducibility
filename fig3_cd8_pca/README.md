# Figure 3 / S3 — CD8 T-compartment composition PCA

Reproduction of **Figure 3** of Terekhova et al. (2026), *Immunity*,
DOI 10.1016/j.immuni.2026.09.012: principal-component structure of the
CD8 T-cell compartment across donors and its age gradient.

## What is reproduced

- CD8 minor-population percentages per donor → PCA (`prcomp`,
  authors' scaling).
- Scree plot (`factoextra::fviz_eig`), biplot, and PC2 contribution bars.
- PC1/PC2 colored by the age gradient (Naive-population gradient panel).
- CD8 compartment summary scatter (PC space, dataset colors).
- Tem GZMK+/GZMB+ ratio coloring per dataset.

## Outputs

| Output | Content |
|---|---|
| `fig3_pca_scree.png` | Scree / eigenvalue plot |
| `fig3_pca_biplot.png` | PC1–PC2 biplot with population loadings |
| `fig3_pca_pc2_contribution.png` | PC2 variable contributions |
| `fig3_pca_age_gradient.png` | Donors in PC space colored by age / naive gradient |
| `fig3_pca_dataset_scatter.png` | Donors in PC space colored by dataset |
| `fig3_pca_tem_ratio.png` | PC space colored by Tem GZMK+/GZMB+ ratio |

## Provenance

Ported from the authors' `Figures_pipeline/figures_pipeline.R`
(repo `teresho4/Pan_atlas_immune_aging_TemK_TemB_ratio`), section
"Figure 3, S3" (lines 843–961, plus the PCA panels at lines 880–960).
The donor × CD8-population percentage matrix is rebuilt from this
compendium's MD5-pinned raw per-cell table via `fig_support.R`
(`paper_minor_percentages`), matching the authors' percentage semantics.

## Run

```sh
./run_fig3.sh
```

Runtime ~2–3 min (dominated by the raw-table read).
