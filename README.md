# stablr Aging Case Study — Reproduction Compendium

[![render-quick-path](https://github.com/mdmanurung/stablr-reproducibility/actions/workflows/render-quick-path.yml/badge.svg)](https://github.com/mdmanurung/stablr-reproducibility/actions/workflows/render-quick-path.yml)

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
| `stablr/` | Git submodule pinned at `mdmanurung/stablr@6aa5a6e` (v0.3.0): full package source, `benchmarks/` (all experiment code, predeclared protocols, committed results) and `vignettes/` (9 Rmd reports + figures) |
| `results/ws8_v2/` | Canonical, committed copy of the v2 benchmark outputs (the pipeline syncs its scratch outputs here) |
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

The committed result tables in `results/ws8_v2/` (and the other committed
tables under the `stablr` submodule's `benchmarks/`) are exactly what the
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
   Writes into the `stablr` submodule's working tree
   (`stablr/benchmarks/ws8_v2/results/`) and syncs the canonical copy to
   `results/ws8_v2/` (committed). A fresh verification run (2026-10-17)
   reproduced every result file bit-for-bit; only embedded run timestamps and
   commit hashes in `run_manifest.json` / `data_manifest.csv` differ.
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
- `stablr/` is a git submodule pinned at `6aa5a6e`; its canonical history
  lives at <https://github.com/mdmanurung/stablr>. Clone with
  `git clone --recurse-submodules` (or run `git submodule update --init`
  after a plain clone). The pinned commit is exactly the package state used
  for the committed results and the bit-for-bit verification run.

## Relationship to current stablr main

**The pin is now `fdbf7f1`** (the 1.0.0 API cut, bumped 2026-10-22; the
committed results and the 2026-10-17 bit-for-bit verification evidence in
`verification/` were produced at the earlier pin `6aa5a6e`, v0.3.0). The
`6aa5a6e..fdbf7f1` span is the Phase 13 remediation (RNG hygiene,
1.0.0-contract closure, reliability diagnostics, null validation,
replication profile, gaussian validated path, lean public API), the
post-review fix batches E-1..E-4 and F-1..F-8, the Phase 14 API
streamlining (54 -> 36 exports), and documentation (the package grew from
9 to 14 vignettes).

**Sealed full v1/v2 re-run at the new pin (2026-10-22).** The full chain
was re-run from the MD5-pinned raw atlas with the environment matched to
`sessionInfo.txt` (R 4.4.2, ranger 0.18.0, glmnet 5.0, knockoff 0.3.6).
Outcome: **claims-level reproducibility PASS, bit-level FAIL** — every
re-tested claims-ledger conclusion (V1-V9) holds with the same direction
and significance (the mixed/negative headline is unchanged in structure;
null saturation 59/59 in 50/50 permutations is identical), while 81 of
the 108 regenerated v2 artifacts and all regenerated v1 tables differ
byte-level, with the delta scoped to the RNG-stream-sensitive M4/M5 STABL
paths and the exploratory v1 layer (M0-M3 are bit-exact). The committed
result tables in `results/ws8_v2/` and the submodule `benchmarks/` remain
the originally published values; the re-run record, checksum evidence,
and the full claims-impact assessment live in `verification/`
(`VERDICT.md`, `sealed_rerun_2026-10-22/`). Re-basing the committed
results (adopting the fresh numbers and updating the vignette's inlined
values) is a separate, explicitly human-gated decision. The pipeline's
final sync step has a known path bug (`run_v2_pipeline.sh` copies a
`stablr/…` relative path after `cd stablr`); the canonical copy was
historically placed manually.

Procedure for future pin moves (as executed for this bump):

1. `cd stablr && git fetch origin && git checkout <new-sha> && cd .. &&
   git add stablr && git commit` (record the bump).
2. Re-run the full chain: `./run_all.sh` (many hours; the v2 LOCO
   benchmark dominates).
3. Compare every committed result artifact against the re-run
   (checksums + claims-ledger impact) and record the evidence in
   `verification/`.
4. Refresh `sessionInfo.txt` and commit code, results, and evidence
   together.

The package's other vignettes — including the single-cell integrations
(`stablr-cebmf`, `stablr-metacells`) added after the pin — are
self-contained: their data ship as fixtures inside the package itself
(`stablr/inst/extdata/`, with generation scripts in `stablr/data-raw/`),
so they need no compendium coverage.

## License

The bundled package source is GPL (>= 2), see `stablr/LICENSE.md`. The
benchmark scripts and this compendium inherit the same license. The raw data
is public per the source publication; cite Terekhova et al. (2026).
