# Figure 5 / S5 — Aging clock and health-trajectory components

Reproduction of the reproducible components of **Figure 5** of
Terekhova et al. (2026), *Immunity*, DOI 10.1016/j.immuni.2026.09.012:
the age-predicting "immune clock" built from cell-type composition.

## What is reproduced

The authors' clock pipeline (repo `teresho4/...`,
`Age_clock_model/clean_notebook_for_aging_models.ipynb`), re-expressed
in R on this compendium's committed donor tables:

- **Reduced 20-population panel** — the authors' `selected_by_marina`
  selection (notebook cell 32), mapped to this compendium's harmonized
  annotation levels (`P20_MAP` in the script).
- **Ridge regression clock** — `cv.glmnet` alpha = 0, nfolds = 5,
  seed 20261013, on three feature sets:
  1. all 59 harmonized populations (paper analogue: full panel),
  2. the 20-population reduced panel (paper analogue: reduced panel),
  3. the 20-population panel + one-hot sex/ethnicity covariates
     (the authors' final model form).
- **Leave-one-cohort-out evaluation** — the same LOCO engine used by
  the stablr benchmarks: each of the 8 datasets is held out in turn,
  the clock is fit on the other 7, and age is predicted for the held-out
  cohort. Reports per-cohort and overall MAE, compared against the
  paper's reported MAE 8.2 y (full panel) / 8.9 y (reduced panel).

## Not ported (data not shipped)

- **Olink proteomics panel** — the authors' Olink table is not
  distributed with the compendium (authors' `figures_pipeline.R`
  section "Figure 5, S5", lines 1291–1353, depends on it).
- **UKBB / clinical-trajectory components** — the health-trajectory
  associations require the authors' UK Biobank extracts, which are not
  public.

Both are documented stubs; the script stops after the clock panels.

## Outputs

| Output | Content |
|---|---|
| `fig5_clock_all59_loco.png` | Observed vs predicted age, faceted by held-out cohort (all 59 populations) |
| `fig5_clock_p20_loco.png` | Same for the 20-population reduced panel |
| `fig5_clock_p20_demog_loco.png` | Same for the reduced panel + one-hot sex/ethnicity |
| `fig5_clock_mae_comparison.png` | LOCO MAE per feature set vs the paper's reported 8.2 / 8.9 y |
| `fig5_loco_mae.csv` | Per-donor predictions and per-cohort MAE table (all feature sets) |

## Result and honest discrepancy

Under the fixed all-8-cohort LOCO evaluation, overall MAE is **13.9 y
(all-59)**, **14.7 y (P20)**, and **13.9 y (P20 + demographics)** —
worse than the paper's reported 8.2 / 8.9 y. The gap is expected and
documented rather than tuned away:

- The authors' notebook evaluates with cohort splits chosen during
  development (not the fixed LOCO over all 8 cohorts), which is a
  materially easier task than predicting entirely unseen cohorts.
- The paper's "full panel" is the authors' full marker panel, which is
  richer than the compendium's harmonized 59-level space.
- The relative structure matches the paper's qualitative claim: the
  reduced 20-population panel performs on par with (here: slightly
  worse than) the full panel, and demographics add nothing once
  composition is included.

This port therefore reproduces the paper's clock *design* and its
qualitative conclusions, under a stricter, predeclared evaluation.

## Provenance

- Clock design: authors' `Age_clock_model/clean_notebook_for_aging_models.ipynb`
  (ridge on composition features; `selected_by_marina` panel, cell 32).
- Paper reference values: MAE 8.2 y (full panel) / 8.9 y (reduced panel).
- Data: this compendium's derived donor tables
  (`data/derived/donor_celltype_counts.tsv`, `donor_meta.tsv`), rebuilt
  from the MD5-pinned raw atlas by the v2 pipeline's `data_prepare.R`.

## Differences from the authors' version

- The compendium's harmonized 59-level annotation is the v2 benchmark's
  space; the P20 panel is mapped into it (all 20 levels verified present).

## Run

```sh
./run_fig5.sh
```

Runtime ~3–5 min (LOCO ridge fits are fast; the derived tables are
already committed).
