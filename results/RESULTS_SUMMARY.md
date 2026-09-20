# Corrected pipeline — independent rebuild with verified gyrA/qnrB fixes

Full runnable code: `corrected_pipeline.py`. Raw inputs: the 10 files in `data/raw/`.
Complete console output: `pipeline_output.txt`. Nothing below is hand-typed —
every number is copied directly from that output.

## Headline comparison vs. this repo's reported corrected numbers

| | This rebuild | Repo (table3_4_full_bootstrap_CIs.txt) |
|---|---|---|
| UK AUROC | 0.991 | 0.992 |
| Norway AUROC | 0.968 | 0.968 |
| South Asia AUROC | 0.969 | 0.970 |
| Norway calib. intercept | -0.507 | -0.488 to -0.494 |
| South Asia calib. intercept | +0.918 | +0.927 to +0.973 |

Agreement is within 0.001–0.003 AUROC and ~0.02–0.06 intercept everywhere,
despite this rebuild using 347 aac6_Ib_cr genomes vs. the ~489 reported
elsewhere (see `independent_verification.md`) — the discrepancy does not
appear to materially change any headline conclusion.

## Mechanism composition — exact match, including the new Norway finding
UK 98.5% / 0.1% / 1.4%, Norway 88.1% / 1.3% / **10.6%**, South Asia 94.3% / 4.4% / 1.3%
(parC-mutation / no-parC-acquired-gene / neither) — reproduces
`section3_4_mechanism_corrected.txt` exactly, including the 33-genome
Norway "neither" group independently.

## One new discrepancy found in this rebuild
This run found **0/300** extreme recalibration reps (|intercept|>10) at
South Asia N=200, vs. the 17/300 (up to -38.76) reported in
`section3_5_recalibration_corrected.txt`. Plausibly downstream of the
aac6_Ib_cr difference (a smaller feature-positive count changes which
genomes are separable in a 38-genome held-out fold) — flagged here rather
than silently reproducing a smoother result than what was originally found.
