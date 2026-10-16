#!/usr/bin/env bash
# Reproduce the v1 (exploratory) WS-8 aging showcase end to end.
# Runs from the package root so the scripts' relative paths resolve, loading
# stablr from the included source via pkgload. Outputs land in outputs/v1/
# (the committed copies in stablr/benchmarks/ are left untouched).
#
# Runtime: roughly 1-2 hours on 16 cores (RF fitting + LOCO loops).
set -euo pipefail
cd "$(dirname "$0")"

ROOT="$(pwd)"
export WS8_DATA_DIR="$ROOT/data/derived"
export WS8_OUT_DIR="$ROOT/outputs/v1"
export WS8_OBJECTS="$ROOT/outputs/v1/ws8_objects.rds"
mkdir -p "$WS8_OUT_DIR"

cd stablr
Rscript -e 'pkgload::load_all(".", compile = FALSE, quiet = TRUE); source(commandArgs(TRUE))' benchmarks/ws8_aging_showcase.R
Rscript -e 'pkgload::load_all(".", compile = FALSE, quiet = TRUE); source(commandArgs(TRUE))' benchmarks/ws8_figures.R
Rscript -e 'pkgload::load_all(".", compile = FALSE, quiet = TRUE); source(commandArgs(TRUE))' benchmarks/ws8_recalibration.R
Rscript -e 'pkgload::load_all(".", compile = FALSE, quiet = TRUE); source(commandArgs(TRUE))' benchmarks/ws8_loco_selection_isolated.R
Rscript -e 'pkgload::load_all(".", compile = FALSE, quiet = TRUE); source(commandArgs(TRUE))' benchmarks/ws8_review_fdp.R
Rscript -e 'pkgload::load_all(".", compile = FALSE, quiet = TRUE); source(commandArgs(TRUE))' benchmarks/ws8_review_paired_mae.R
Rscript -e 'pkgload::load_all(".", compile = FALSE, quiet = TRUE); source(commandArgs(TRUE))' benchmarks/ws8_review_sensitivity_nospac.R
Rscript -e 'pkgload::load_all(".", compile = FALSE, quiet = TRUE); source(commandArgs(TRUE))' benchmarks/ws8_rule_comparison.R

echo "v1 showcase complete. Outputs: outputs/v1/"
