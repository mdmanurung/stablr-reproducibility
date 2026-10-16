#!/usr/bin/env bash
# Download the raw per-cell annotation table for the Terekhova et al. (2026)
# healthy-aging PBMC atlas and verify its checksum.
#
# Source: Artyomov lab single-cell portal, dataset pan_atlas_healthy_aging_all_celtypes
# (Terekhova et al., Immunity 2026, DOI 10.1016/j.immuni.2026.09.012).
# The file shipped in data/raw/ already passes this check; re-download only if
# you want to fetch it fresh from the portal.
set -euo pipefail
cd "$(dirname "$0")"

URL="https://artyomovlab.wustl.edu/scn-m/var/datasets/curated/maksim/GZMB_GZMK/pan_atlas_healthy_aging_all_celtypes/files/all_cells_metadata.csv.gz"
EXPECTED_MD5="50651968855bf6ca3efcdf199e2a9fd5"
OUT="data/raw/all_cells_metadata.csv.gz"

mkdir -p data/raw

if [ -f "$OUT" ]; then
    echo "File already present: $OUT"
else
    echo "Downloading $(basename "$URL") (~93 MB)..."
    curl -fL --retry 3 --max-time 3600 -o "$OUT" "$URL"
fi

ACTUAL_MD5=$(md5sum "$OUT" | awk '{print $1}')
if [ "$ACTUAL_MD5" = "$EXPECTED_MD5" ]; then
    echo "MD5 OK: $ACTUAL_MD5"
else
    echo "MD5 MISMATCH: expected $EXPECTED_MD5, got $ACTUAL_MD5" >&2
    exit 1
fi
