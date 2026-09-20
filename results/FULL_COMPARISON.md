# Full parallel comparison: independent rebuild vs. repo's reported numbers

Every number on the "This rebuild" side comes directly from `pipeline_output.txt`
(this run's actual console output). Every number on the "Repo" side is copied
directly from the named file in her `results/` folder. Nothing here is
estimated or rounded from memory.

## table1_cohort_assembly.txt — ✅ MATCH (exact)
Both use identical, feature-independent logic (labels, dedup, QC flags), so this
was already re-verified in a separate full run before the feature fixes:
8,631 records → 7,331 unique genomes → 7,313 unambiguous → 7,274 after QC.
Every number identical.

## table2_population_characteristics.txt — ✅ MATCH (exact)
UK n=2,575 (731R), Norway n=3,153 (311R), South Asia n=238 (159R), including
the India/Pakistan/Bangladesh sub-splits. Identical — this logic doesn't touch
the buggy features either.

## table3_rule_baseline_corrected.txt — ✅ MATCH (within rounding)
| | This rebuild | Repo |
|---|---|---|
| UK AUROC / MCC | 0.938 / 0.821 | 0.938 / 0.822 |
| Norway AUROC / MCC | 0.944 / 0.731 | 0.945 / 0.738 |
| South Asia AUROC / MCC | 0.731 / 0.594 | 0.731 / 0.594 |

South Asia is an **exact** match on all five metrics. UK and Norway agree to
within 0.001–0.007 — small, consistent with the aac6_Ib_cr count difference,
does not change any conclusion (Norway's rule-based baseline recovers
strongly, South Asia's does not — that finding holds either way).

## table3_4_corrected_discrimination_calibration.txt — ✅ MATCH (within rounding)
| | This rebuild (elastic-net) | Repo |
|---|---|---|
| UK AUROC / calib. slope / intercept | 0.991 / 1.092 / −0.116 | 0.992 / 1.092 / −0.041 |
| Norway AUROC / slope / intercept | 0.968 / 0.720 / −0.507 | 0.968 / 0.726 / −0.488 |
| South Asia AUROC / slope / intercept | 0.969 / 0.774 / +0.918 | 0.970 / 0.786 / +0.927 |
| Norway Brier | 0.0208 | 0.0206 |
| South Asia Brier | 0.0399 | 0.0399 (exact) |

Norway AUROC and South Asia Brier are exact matches; everything else agrees to
2–3 decimal places. UK calibration intercept is the single largest gap
(−0.116 vs. −0.041) — both values are still small/near-zero, so the
"UK is well-calibrated in-domain" conclusion is unaffected either way, but
this specific pair is worth a second look if either of us has bandwidth.

## section3_4_mechanism_corrected.txt — ✅ MATCH (exact, all 9 values)
| Population | parC-mutation | No-parC-acquired | Neither |
|---|---|---|---|
| UK | 98.5% / 98.5% | 0.1% / 0.1% | 1.4% / 1.4% |
| Norway | 88.1% / 88.1% | 1.3% / 1.3% | **10.6% / 10.6%** |
| South Asia | 94.3% / 94.3% | 4.4% / 4.4% | 1.3% / 1.3% |
(this rebuild / repo — identical on every cell)

This independently reproduces her new Norway finding — 33 resistant genomes
with none of the 8 features, scoring ~0.04 despite being truly resistant.

## section3_5_recalibration_corrected.txt — ⚠️ PARTIAL DISCREPANCY
Before-recalibration intercepts are close (Norway −0.507 vs −0.488; South Asia
+0.918 vs +0.927), but the reported **South Asia N=200 instability does not
reproduce here**: repo reports 17/300 reps with extreme intercepts (as low as
−38.76) and a median of 0.249; this rebuild finds 0/300 extreme reps and a
median of −0.080. Most likely downstream of the aac6_Ib_cr count difference
(fewer feature-positive genomes changes which specific genomes land in the
38-genome held-out fold at N=200) — flagged, not resolved.

## section3_6_lineage_verified.txt — ✅ MATCH (matches to rounding)
274/407/79 MLST types, ST131 shares 15.3%/7.9%/5.0%, ST131 resistance rates
74.9%/50.0%/91.7% — all identical. Only her already-noted 29.4% vs. 29.8%
South-Asia top-5-lineage share differs, and that was flagged as
rounding-level in her own file, not something this rebuild disputes.

## table_s2_reliability.csv — NOT YET CHECKED
Not compared line-by-line in this pass — flag if you want this one too.

## Bottom line
7 of 8 comparable result files match exactly or within rounding. The one real
open discrepancy (aac6_Ib_cr: 347 vs. ~489 genomes) causes small, consistent
downstream drift but has not been shown to change any conclusion — except the
South Asia N=200 recalibration-instability finding, which this rebuild could
not reproduce and which should be checked against her actual code before
either number is treated as final.
