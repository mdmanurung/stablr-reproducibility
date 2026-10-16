#!/usr/bin/env bash
# Full reproduction chain: data download -> v1 showcase -> v2 pipeline ->
# report rendering. Each stage is also runnable on its own.
#
#   ./run_all.sh              # everything (many hours)
#   SKIP_V1=1 ./run_all.sh    # skip the v1 exploratory showcase
#   SKIP_V2=1 ./run_all.sh    # skip the v2 benchmark
set -euo pipefail
cd "$(dirname "$0")"

./00_download_data.sh

if [ "${SKIP_V1:-0}" != "1" ]; then ./run_v1_showcase.sh; fi
if [ "${SKIP_V2:-0}" != "1" ]; then ./run_v2_pipeline.sh; fi

./render_reports.sh
echo "Reproduction chain complete."
