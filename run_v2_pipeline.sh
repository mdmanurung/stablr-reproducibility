#!/usr/bin/env bash
# Reproduce the v2 (prospective, predeclared) WS-8 benchmark end to end:
# data preparation, the leakage-proof LOCO replication benchmark for all six
# selectors (M0-M5), the age-permutation null, the real-data reference fit,
# secondary biology, and figures.
#
# Writes into the exact directories the assembly scripts consume
# (run_m012b/run_m3b/run_m4b/run_m5b — the Amendment-1 amended runs), then
# syncs a canonical copy to results/ws8_v2/ (committed at the compendium
# root; stablr is a git submodule, so its working tree is scratch).
#
# Runtime: several hours on 16 cores (the M4/M5 STABL groups dominate).
# The 50-permutation null runs as one invocation; each permutation uses its
# own seed stream (20261013 + perm_id), matching the original chunked runs
# bit-for-bit.
set -euo pipefail
cd "$(dirname "$0")"

ROOT="$(pwd)"
export WS8_DATA_DIR="$ROOT/data/raw"

cd stablr

echo "== 1/10 data_prepare.R (donor tables from the raw portal file) =="
Rscript -e 'pkgload::load_all(".", compile = FALSE, quiet = TRUE); source(commandArgs(TRUE))' benchmarks/ws8_v2/data_prepare.R

echo "== 2/10 data_checks.R (gate G1) =="
Rscript -e 'pkgload::load_all(".", compile = FALSE, quiet = TRUE); source(commandArgs(TRUE))' benchmarks/ws8_v2/data_checks.R

echo "== 3/10 folds_manifest.R (fixed outer folds) =="
Rscript -e 'pkgload::load_all(".", compile = FALSE, quiet = TRUE); source(commandArgs(TRUE))' benchmarks/ws8_v2/folds_manifest.R

echo "== 4/10 run.R: primary + secondary benchmark, all selectors =="
for GROUP in "M0,M1,M2" "M3" "M4" "M5"; do
    case "$GROUP" in
        "M0,M1,M2") RUNDIR=run_m012b ;;
        "M3")       RUNDIR=run_m3b ;;
        "M4")       RUNDIR=run_m4b ;;
        "M5")       RUNDIR=run_m5b ;;
    esac
    echo "   method group: $GROUP -> $RUNDIR"
    WS8_METHODS="$GROUP" WS8_V2_OUT="benchmarks/ws8_v2/results/$RUNDIR" \
        Rscript -e 'pkgload::load_all(".", compile = FALSE, quiet = TRUE); source(commandArgs(TRUE))' benchmarks/ws8_v2/run.R
done

echo "== 5/10 assemble_primary.R =="
Rscript -e 'pkgload::load_all(".", compile = FALSE, quiet = TRUE); source(commandArgs(TRUE))' benchmarks/ws8_v2/assemble_primary.R

echo "== 6/10 assemble_prediction.R (secondary prediction tables) =="
Rscript -e 'pkgload::load_all(".", compile = FALSE, quiet = TRUE); source(commandArgs(TRUE))' benchmarks/ws8_v2/assemble_prediction.R

echo "== 7/10 null_probe.R (age-permutation null, 50 permutations) =="
WS8_N_PERM=50 WS8_PERM_START=1 \
    Rscript -e 'pkgload::load_all(".", compile = FALSE, quiet = TRUE); source(commandArgs(TRUE))' benchmarks/ws8_v2/null_probe.R

echo "== 8/10 real_reference.R (full-atlas reference fit) =="
Rscript -e 'pkgload::load_all(".", compile = FALSE, quiet = TRUE); source(commandArgs(TRUE))' benchmarks/ws8_v2/real_reference.R

echo "== 9/10 biology.R (reference populations, GZMK/GZMB, heterogeneity) =="
Rscript -e 'pkgload::load_all(".", compile = FALSE, quiet = TRUE); source(commandArgs(TRUE))' benchmarks/ws8_v2/biology.R

echo "== 10/10 figures.R =="
Rscript -e 'pkgload::load_all(".", compile = FALSE, quiet = TRUE); source(commandArgs(TRUE))' benchmarks/ws8_v2/figures.R

echo "== sync: canonical copy at results/ws8_v2/ =="
mkdir -p "$ROOT/results/ws8_v2"
cp -r stablr/benchmarks/ws8_v2/results/. "$ROOT/results/ws8_v2/"

echo "v2 pipeline complete. Canonical outputs: results/ws8_v2/"
