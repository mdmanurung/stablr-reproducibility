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
