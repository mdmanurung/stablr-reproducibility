#!/usr/bin/env bash
# Reproduce paper Figure 2 (Terekhova et al. 2026) end to end:
# ensures the raw atlas is present (MD5-checked), then runs the ported
# figure script. Outputs land in fig*/outputs/ (committed).
set -euo pipefail
cd "$(dirname "$0")"
./00_download_data.sh
Rscript fig2_age_associations/fig2_age_associations.R
