#!/usr/bin/env bash
# Reproduce the v2 (prospective, predeclared) WS-8 benchmark end to end:
# data preparation, the leakage-proof LOCO replication benchmark for all six
# selectors (M0-M5), the age-permutation null, secondary biology, and figures.
#
# Follows the run order in stablr/benchmarks/ws8_v2/README.md. run.R sources
# metrics.R, loco.R and selectors.R itself, so they are not invoked separately.
#
# NOTE: writes into stablr/benchmarks/ws8_v2/results/ in place (the scripts'
# fixed contract). The committed results shipped in this compendium are the
# expected values; restore them afterwards with:
#   cd stablr && git checkout -- benchmarks/ws8_v2/results
#
# Runtime: several hours on 16 cores (four run.R method groups dominate).
set -euo pipefail
cd "$(dirname "$0")"

ROOT="$(pwd)"
export WS8_DATA_DIR="$ROOT/data/raw"

cd stablr

echo "== 1/7 data_prepare.R (donor tables from the raw portal file) =="
Rscript -e 'pkgload::load_all(".", compile = FALSE, quiet = TRUE); source(commandArgs(TRUE))' benchmarks/ws8_v2/data_prepare.R

echo "== 2/7 data_checks.R =="
Rscript -e 'pkgload::load_all(".", compile = FALSE, quiet = TRUE); source(commandArgs(TRUE))' benchmarks/ws8_v2/data_checks.R

echo "== 3/7 run.R: primary + secondary benchmark, all selectors =="
for GROUP in "M0,M1,M2" "M3" "M4" "M5"; do
    echo "   method group: $GROUP"
    WS8_METHODS="$GROUP" WS8_V2_OUT="benchmarks/ws8_v2/results/run_${GROUP//,/}" \
        Rscript -e 'pkgload::load_all(".", compile = FALSE, quiet = TRUE); source(commandArgs(TRUE))' benchmarks/ws8_v2/run.R
done

echo "== 4/7 assemble_primary.R =="
Rscript -e 'pkgload::load_all(".", compile = FALSE, quiet = TRUE); source(commandArgs(TRUE))' benchmarks/ws8_v2/assemble_primary.R

echo "== 5/7 null_probe.R (age-permutation null) =="
Rscript -e 'pkgload::load_all(".", compile = FALSE, quiet = TRUE); source(commandArgs(TRUE))' benchmarks/ws8_v2/null_probe.R

echo "== 6/7 biology.R (reference populations, GZMK/GZMB, heterogeneity) =="
Rscript -e 'pkgload::load_all(".", compile = FALSE, quiet = TRUE); source(commandArgs(TRUE))' benchmarks/ws8_v2/biology.R

echo "== 7/7 figures.R =="
Rscript -e 'pkgload::load_all(".", compile = FALSE, quiet = TRUE); source(commandArgs(TRUE))' benchmarks/ws8_v2/figures.R

echo "v2 pipeline complete. Outputs: stablr/benchmarks/ws8_v2/results/"
