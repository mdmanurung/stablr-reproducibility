#!/usr/bin/env bash
# Render the aging case-study vignette (and optionally all vignettes) to
# outputs/reports/. The vignettes read the committed result tables, so this
# is the fast path: minutes, no refitting.
set -euo pipefail
cd "$(dirname "$0")"

mkdir -p outputs/reports

cd stablr
if [ "${RENDER_ALL_VIGNETTES:-0}" = "1" ]; then
    for RMD in vignettes/*.Rmd; do
        echo "Rendering $RMD ..."
        Rscript -e 'pkgload::load_all(".", compile = FALSE, quiet = TRUE); rmarkdown::render(commandArgs(TRUE), output_dir = file.path("..", "outputs", "reports"))' "$RMD"
    done
else
    echo "Rendering the aging case study ..."
    Rscript -e 'pkgload::load_all(".", compile = FALSE, quiet = TRUE); rmarkdown::render(commandArgs(TRUE), output_dir = file.path("..", "outputs", "reports"))' vignettes/stablr-aging-case-study.Rmd
fi

echo "Reports written to outputs/reports/"
