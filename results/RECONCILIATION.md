# Reconciliation with the second re-implementation of the pipeline

This note records how each difference between the reported numbers and the second re-implementation of the pipeline was resolved. The re-implementation is
commit `cb5ad9a` (`notebooks/corrected_pipeline.py`, `results/pipeline_output.txt`, `results/FULL_COMPARISON.md`, `results/RESULTS_SUMMARY.md`); those files are
unchanged. It reused the same feature-construction rules, so its agreement on the features checks the code, not the rules themselves. The rules are
documented and tested in Methods section 2.4.1.

**Status.** Both discrepancies flagged in `FULL_COMPARISON.md` (the AAC(6')-Ib-cr count and the South Asia N = 200 recalibration instability) are caused by
the re-implementation's code, and the reported numbers are unchanged.

## What matched

Table 1 and Table 2 counts; mechanism composition in every cell, including the 33 Norwegian genomes with neither a parC mutation nor an acquired gene;
lineage statistics; and discrimination (AUROC within 0.001 and MCC within 0.007 of the reported values, model coefficients within 0.08).

## Differences and their causes

| # | Difference | Cause | Evidence | Effect on the reported numbers |
|---|---|---|---|---|
| 1 | AAC(6')-Ib-cr in 347 genomes, against 485 in the cohort (489 in the pooled export) | The re-implementation removes duplicate rows on `(Genome ID, BRC ID, Source)` before matching. In this export the Ib-cr text is spread over separate evidence rows, so the informative rows are discarded. | Check V1 below: 489 genomes match with all rows; 347 after that de-duplication (142 genomes lost; rows 31,958 to 23,929). | None. 485 in the cohort stands. |
| 2 | No extreme South Asia N = 200 estimates and post-recalibration medians near zero at every N, against 17 of 300 extreme estimates and Table 7 | The re-implementation fits both logistic models with scikit-learn's default ridge penalty (C = 1). The reported experiment (Methods 2.9) uses C = 10^6, effectively unpenalized. | Check V2 below, on identical splits: C = 1 gives medians within 0.20 of zero at every N and no extreme South Asian estimates (Norway N = 25 median -0.120; the re-implementation reports -0.128). C = 10^6 reproduces the reported pattern (Norway N = 25 median -0.596; South Asia N = 200: 12 of 300 estimates beyond +-5, minimum -53.0). | None. Table 7 stands. A post hoc sensitivity analysis is added to Methods 2.9 and Results 3.5. |
| 3 | Calibration slopes and intercepts differ slightly (UK intercept -0.116 against -0.041; Norway -0.507 against -0.488; South Asia 0.918 against 0.927) | South Asia: the ridge penalty in the re-implementation's calibration fit. UK and Norway: not the penalty. The re-implementation's AAC(6')-Ib-cr and gyrA rules also differ, and it keeps the one UK genome without an MLST assignment (2,575 genomes, as its own group), so its model and folds differ. This was not tested separately. | `results/calibration_penalty_check.txt` refits the calibration on the reported predictions: with C = 1 the South Asia intercept is 0.918 (the re-implementation reports 0.918), while the UK (-0.039) and Norway (-0.488) intercepts barely change. The UK difference lies inside the reported interval (-0.654 to 0.306). | None. |
| 4 | gyrA rule | The re-implementation matches "quinolone" as a substring of the antibiotic class; the reported rule requires the class to equal "quinolone". The 31 quinolone/triclosan POINTP records are the difference. | Of the 31 genomes carrying them, 20 are gyrA-positive through other records (19 resistant); the other 11 are all susceptible and would gain the feature under the substring rule. | None for resistant genomes. The reported rule stands (Methods 2.4.1). |
| 5 | qnrS in 175 genomes against 176 (pooled exports) | The re-implementation ignores the Product field; the reported rule applies the allele pattern to Gene, Function and Product. | One genome. Not investigated further. | None. |
| 6 | Number of valid recalibration repetitions | The re-implementation discards repetitions with one phenotype class without replacement (281 and 297 valid at Norway N = 25 and 50); the reported experiment replaces them to reach 300. | Output of the re-implementation. | None. |
| 7 | Feature counts (gyrA 2,226; parC 1,828; qnrB 29; qnrS 175) | The re-implementation counts genomes over all exported records; the paper reports counts in the 7,274-genome cohort (qnrB 28; AAC(6')-Ib-cr 485). | Pooled and cohort counts differ by the genomes removed at quality control. | None. |

## Output of the confirmation checks

Code: `notebooks/reconciliation_checks.py`. The output below was copied from the notebook.

```
V1  AAC(6')-Ib-cr genomes, all rows: 489 | after the rebuild de-duplication: 347 | genomes lost: 142 | rows 31958 -> 23929

V2  Norway, C=1e+06
  N  reps  median    q25    q75  abs_gt5  abs_gt10     min
 25   281  -0.596 -0.751 -0.395        2         0 -7.506
 50   297  -0.336 -0.651  0.758        1         0 -3.882
100   300  -0.196 -0.631  0.723        0         0 -3.042
150   300  -0.006 -0.506  0.538        0         0 -2.413
200   300   0.051 -0.395  0.445        0         0 -1.720

V2  Norway, C=1
  N  reps  median    q25   q75  abs_gt5  abs_gt10     min
 25   281  -0.120 -0.598 0.416        1         1 -18.155
 50   297   0.191 -0.532 0.800        0         0  -3.003
100   300   0.094 -0.480 0.682        0         0  -2.844
150   300   0.071 -0.447 0.536        0         0  -2.236
200   300   0.063 -0.380 0.444        0         0  -1.607

V2  South Asia, C=1e+06
  N  reps  median    q25   q75  abs_gt5  abs_gt10     min
 25   300  -0.205 -0.800 1.111       16         3  -6.621
 50   300  -0.138 -0.598 0.385       12         5  -3.051
100   300  -0.067 -0.438 0.283        1         0  -2.094
150   300   0.152 -0.345 0.576        0         0  -4.757
200   300   0.233 -0.483 0.946       12        12 -53.015

V2  South Asia, C=1
  N  reps  median    q25   q75  abs_gt5  abs_gt10     min
 25   300   0.035 -0.616 0.821        0         0 -3.330
 50   300  -0.004 -0.527 0.429        0         0 -2.730
100   300  -0.054 -0.426 0.301        0         0 -1.957
150   300   0.040 -0.380 0.445        0         0 -2.126
200   300  -0.025 -0.568 0.417        0         0 -2.591
```
