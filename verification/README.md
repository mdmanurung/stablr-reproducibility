# Bit-for-bit verification checksums

Evidence for the v2 reproduction verification described in the README:
MD5 checksums of every committed result artifact were taken before and
after an independent end-to-end rerun of `run_v2_pipeline.sh`, then
compared.

- `pre_run_checksums.txt` — 76 result artifacts (biology/, figures/,
  m1/, null/, run/, s1 inputs) checksummed before the rerun.
- `post_run_checksums.txt` — the same artifacts checksummed after the
  rerun.
- `checksum_diff.txt` — the complete diff: exactly one file differs,
  `m1/data_manifest.csv`, which embeds a preparation timestamp. All 75
  remaining artifacts are byte-identical across the two independent runs.

Scope note: this set covers the 76 top-level result files; per-run
`run_manifest.json` files (which embed the run timestamp and commit hash)
were verified separately and differ only in the same metadata fields.

## Submodule pin bump (2026-10-18)

The `stablr/` submodule was re-pinned from `5cf603b` (v0.3.0) to `9e5c6a8`
(v0.3.0.9000) after the documentation-audit and Phase 7 commits. The
committed verification evidence above was generated at `5cf603b`; the bump
is justified by code identity between the two revisions:

- `git diff 5cf603b..9e5c6a8 -- R/ src/` contains exactly one commit
  (`a5d1470`, the experimental sPLS-DA base learner): one new file
  (`R/learner_splsda.R`) plus additive wiring in `R/stabl_fit.R` (a new
  dispatch branch, a new defaulted formal `splsda_ncomp`, and a new switch
  case). No line executed by the v2 pipeline (lasso/rf learners,
  knockoffs, FDP+, environment stability) is modified.
- Empirical bit-identity probe (lasso + rf `stabl_fit` with
  `knockoff_equi`, environment-stability aggregation and pairwise
  similarity, fixed seeds): all outputs `identical()` between the two
  checkouts — PASS on every path.

The pin therefore tracks the current package main without invalidating the
bit-for-bit evidence; a future full rerun at the new pin would regenerate
the same checksums by the identity argument above.

## Submodule pin bump (2026-10-18, final pre-publication review)

Re-pinned `9e5c6a8` → `6aa5a6e` (the final pre-publication per-function
review commit: 1 P1 fix — rf auto-grid NULL-seed failure, 5 fail-closed
robustness fixes, 9 doc corrections, the advanced-vignette Jaccard-heatmap
axis fix, audit records; see stablr `docs/audit/2026-10-17-adversarial/
final-review.md`). The bump is justified by code identity on the pipeline
paths:

- `git diff 9e5c6a8..6aa5a6e -- R/ src/` touches 11 R files (75 insertions,
  14 deletions) and no `src/` file. Every behavioural change is fail-closed
  validation that only fires on inputs the v2 pipeline never produces
  (NULL rf grid seed — the pipeline passes explicit seeds; out-of-range
  thresholds; length-mismatched plot inputs; lone `y_valid`; mismatched
  cooperative prediction rows) or is documentation/vignette-only. The
  fixed-seed numeric paths executed by the pipeline are untouched.
- Empirical bit-identity probe (lasso + rf `stabl_fit` with
  `knockoff_equi`, `auto_rf_threshold_grid`, environment-stability
  aggregation and pairwise similarity, fixed seeds): all outputs
  `identical()` between the `9e5c6a8` and `6aa5a6e` checkouts — PASS on
  every path.
- `R CMD check` at `6aa5a6e`: 0 errors, 0 warnings, 2 benign notes; full
  test suite 0 failures.
