# Data

## Raw

`raw/all_cells_metadata.csv.gz` — per-cell annotation table (13,137,775 cells,
one row per cell) from the public healthy-aging PBMC atlas of Terekhova et al.
(2026), *Immunity*, DOI 10.1016/j.immuni.2026.09.012. Distributed by the Artyomov
lab single-cell portal, dataset `pan_atlas_healthy_aging_all_celtypes`.
Downloaded 2026-10-13.

- MD5: `50651968855bf6ca3efcdf199e2a9fd5`
- Download URL (live as of 2026-10):
  `https://artyomovlab.wustl.edu/scn-m/var/datasets/curated/maksim/GZMB_GZMK/pan_atlas_healthy_aging_all_celtypes/files/all_cells_metadata.csv.gz`
- Re-download and verify: `./00_download_data.sh` from the compendium root.
- Columns include donor id, age, sex, CMV status, ancestry, cohort, and
  cell-type annotation. 2,609 donors across 8 cohorts (OneK1K 981, AIDA 625,
  CIMA 421, KCL 188, ABF300 166, SoundLife 96, UCSF 99, SPAC 33).

## Derived

`derived/donor_celltype_counts.tsv` and `derived/donor_meta.tsv` — donor-level
counts and metadata consumed by the v1 showcase script
(`stablr/benchmarks/ws8_aging_showcase.R`). Derived from the raw file before the
v1 run; the v1 derivation script was not archived (documented limitation — the
v2 pipeline rebuilds donor tables from the raw file instead, see
`stablr/benchmarks/ws8_v2/data_prepare.R`).

- MD5 `donor_celltype_counts.tsv`: `bc547f072feeec6b6f49a8bdd833107e`
- MD5 `donor_meta.tsv`: `241fbbcab1239098fb558169a927026c`

## Citation

If you use these data, cite the source publication:

> Terekhova et al. (2026). *Immunity*. DOI 10.1016/j.immuni.2026.09.012

## Note on git

The raw file is ~93 MB — under GitHub's 100 MB hard limit but above the
recommended 50 MB. If you push this compendium to GitHub, either use Git LFS
for `data/raw/all_cells_metadata.csv.gz` or keep the file out of git
(`.gitignore` it) and let collaborators fetch it with `00_download_data.sh`.
