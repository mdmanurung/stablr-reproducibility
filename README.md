# stablr Aging Case Study — Reproduction Compendium

Self-contained compendium reproducing the immune-aging application case study
of the **stablr** R package (stability-based linear model selection), from raw
data download to rendered report:

<https://mdmanurung.github.io/stablr/articles/stablr-aging-case-study.html>

It bundles the full package source, all experiment/benchmark code, the raw and
derived data, committed result tables, and every figure, so the entire chain —
download, analysis, plots, report — runs on one machine with no other inputs.

## Layout

| Path | Contents |
|---|---|
| `stablr/` | Full stablr package source (v0.3.0), including `benchmarks/` (all experiment code, predeclared protocols, committed results) and `vignettes/` (9 Rmd reports + figures) |
| `data/raw/` | Raw per-cell atlas table, `all_cells_metadata.csv.gz` (13.1M cells, MD5-pinned) |
| `data/derived/` | Donor-level tables consumed by the v1 showcase |
| `figures/` | Copies of every committed figure (vignette + benchmark) for quick browsing |
| `outputs/` | Target directory for re-run results and rendered reports (empty initially) |
| `00_download_data.sh` | Fetch the raw atlas file from the source portal and verify its MD5 |
| `run_v1_showcase.sh` | v1 exploratory showcase: analysis, figures, recalibration, audits |
| `run_v2_pipeline.sh` | v2 predeclared benchmark: data prep, LOCO benchmark, null, biology, figures |
| `render_reports.sh` | Render the case-study vignette (or all vignettes) to `outputs/reports/` |
| `run_all.sh` | The whole chain in order |
| `sessionInfo.txt` | R version and package versions used for the committed results |

## Prerequisites

- **R >= 4.2** (committed results produced with R 4.4.2; see `sessionInfo.txt`)
- A C++ toolchain (the package compiles C++17 sources)
- R packages: `pkgload`, `rmarkdown`; the analysis scripts additionally use
  `glmnet`, `ranger`, `knockoff`, `data.table`, `ggplot2` (installed with the
  package or on first use)
- The stablr package itself is loaded from the bundled source via
  `pkgload::load_all()` — no separate installation is needed. To install it
  system-wide instead: `R CMD INSTALL stablr` or
  `remotes::install_github("mdmanurung/stablr")`.

## Quick path: render the report from committed results (~minutes)

The committed result tables in `stablr/benchmarks/` are exactly what the
published vignette reads. To reproduce the report without refitting anything:

```sh
./render_reports.sh                    # aging case study
RENDER_ALL_VIGNETTES=1 ./render_reports.sh   # all nine vignettes
```

Rendered HTML lands in `outputs/reports/`.

## Full path: from data download to report (many hours)

```sh
./run_all.sh                # everything: download -> v1 -> v2 -> reports
SKIP_V1=1 ./run_all.sh      # skip the exploratory v1 showcase
SKIP_V2=1 ./run_all.sh      # skip the v2 benchmark
```

Or stage by stage:

1. **Data** — `./00_download_data.sh` downloads the raw atlas file (~93 MB)
   from the Artyomov lab portal and verifies MD5
   `50651968855bf6ca3efcdf199e2a9fd5`. Skipped if the file is already present.
2. **v1 showcase** — `./run_v1_showcase.sh` (~1–2 h): the exploratory STABL
   application on 59 donor-level cell-type proportions, LOCO calibration,
   recalibration, selection-isolation, and the post-hoc audit scripts.
   Outputs to `outputs/v1/`.
3. **v2 benchmark** — `run_v2_pipeline.sh` (several hours): the predeclared,
   leakage-proof revision. `data_prepare.R` rebuilds donor tables from the raw
   file (enforcing the `Bad quality` exclusion), then the six selectors
   (M0–M5) run through the 8-fold LOCO engine in four method groups, followed
   by the age-permutation null, secondary biology analyses, and figures.
   Writes into `stablr/benchmarks/ws8_v2/results/` in place — restore the
   committed copies afterwards with
   `cd stablr && git checkout -- benchmarks/ws8_v2/results`.
4. **Report** — `./render_reports.sh` renders the vignette, which reads the
   result tables and regenerates every figure in the published article.

## What the analysis shows (headline claims)

- **v2 primary endpoint** — environment-aware STABL (M4/M5) replicates the
  source atlas' published aging associations at matched budget (K = 10) while
  simpler selectors do not; see `stablr/benchmarks/ws8_v2/claims_ledger.md`
  for the claim-to-artifact trace (V1–V15).
- **Null probe** — within-cohort age permutation saturates selection
  (50/50 features, coverage 0.942), the decisive adverse finding for
  selection-count interpretations.
- **S1 simulations** — M4/M5 select 59/59 native features under all four truth
  regimes (structural saturation); M1 stays calibrated under known truth.
- **S2 sensitivity** — M4q25 > M5 > M4 aggregation ordering (+0.0758
  [0.0245, 0.1272], p = 0.0234); orderings robust to SPAC removal.

All scripts carry the standing caveat: **application case study, not
calibration evidence** — see `stablr/benchmarks/ws8_v2/protocol.md` (the
locked scientific contract) and `stablr/RELEASE_STATUS.md`.

## Data provenance and citation

See `data/README.md`. Source: Terekhova et al. (2026), *Immunity*,
DOI 10.1016/j.immuni.2026.09.012, Artyomov lab single-cell portal, dataset
`pan_atlas_healthy_aging_all_celtypes`. Cite the source publication if you use
the data.

## Notes for pushing to a git remote

- The raw data file is ~93 MB (under GitHub's 100 MB limit, above the 50 MB
  recommendation). Prefer Git LFS for `data/raw/all_cells_metadata.csv.gz`, or
  gitignore it and let collaborators run `00_download_data.sh`.
- `stablr/` is a plain source snapshot (no `.git`); its canonical history
  lives at <https://github.com/mdmanurung/stablr>.

## License

The bundled package source is GPL (>= 2), see `stablr/LICENSE.md`. The
benchmark scripts and this compendium inherit the same license. The raw data
is public per the source publication; cite Terekhova et al. (2026).
