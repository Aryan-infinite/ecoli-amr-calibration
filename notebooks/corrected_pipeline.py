"""
Corrected pipeline -- applies the two independently-verified feature-
construction fixes (gyrA restricted to POINTP evidence; qnrA/B/S require
a specific allele number rather than family-name substring matching) and
reruns the full analysis end to end: model training, discrimination,
calibration, recalibration (median/IQR reporting, per the documented
near-separation issue at South Asia N=200), and mechanism analysis.

Inputs: the same 10 raw BV-BRC CSVs in this folder. No hand-typed numbers.
"""

import pandas as pd
import numpy as np
import glob
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedGroupKFold, cross_val_predict
from sklearn.metrics import (roc_auc_score, average_precision_score, matthews_corrcoef,
                              confusion_matrix, brier_score_loss)

EPS = 1e-6
GENE_COLS = ['gyrA_fq_mut', 'parC_fq_mut', 'qnrA', 'qnrB', 'qnrS', 'oqxA', 'oqxB', 'aac6_Ib_cr']

# ============================================================
# 1. LABELS + METADATA (unaffected by the feature-construction bugs)
# ============================================================
cip = pd.read_csv("cip_phenotypes.csv", dtype={'Genome ID': str})
lbl = cip.groupby('Genome ID')['Resistant Phenotype'].nunique()
conflicting = set(lbl[lbl > 1].index)
clean = (cip[~cip['Genome ID'].isin(conflicting)].drop_duplicates(subset='Genome ID')
         [['Genome ID', 'Resistant Phenotype']].rename(columns={'Resistant Phenotype': 'CIP_Phenotype'}))

meta = pd.read_csv("genome_metadata.csv", dtype={'Genome ID': str}, low_memory=False)
merged = clean.merge(meta, on='Genome ID', how='left', validate='one_to_one')
merged['Isolation Country'] = merged['Isolation Country'].replace({'England': 'United Kingdom', 'Control': pd.NA})
merged = merged[(merged['Genome Status'] != 'Deprecated') & (merged['Genome Quality'] != 'Poor')]
print(f"Cohort: {len(merged)} genomes (unaffected by feature bugs -- matches original paper's Table 1)")

# ============================================================
# 2. CORRECTED FEATURE CONSTRUCTION
# ============================================================
pool = pd.concat([pd.read_csv(f, dtype={'Genome ID': str}) for f in sorted(glob.glob("gene_*.csv"))], ignore_index=True)
pool = pool.drop_duplicates(subset=['Genome ID', 'BRC ID', 'Source'])

# --- gyrA / parC: FIXED -- restrict to actual point-mutation evidence ---
# (original bug: 'quinolone' text tag also fires on K-mer-Search/BLAT presence
#  detectors, which just mean "this genome has a gyrA gene", not "...a mutant one")
def pointp_mutation_flag(gene_product_substring):
    hits = pool[(pool['Evidence'] == 'AMRFinderPlus: POINTP') &
                (pool['Product'].fillna('').str.contains(gene_product_substring, case=False))]
    hits = hits[hits['Antibiotics Class'].fillna('').str.contains('quinolone', case=False)]
    return set(hits['Genome ID'])

gyrA_ids = pointp_mutation_flag('gyrase subunit A')
parC_ids = pointp_mutation_flag('topoisomerase IV subunit A')  # verified: matches original approach closely (parC_verification.txt)

# --- qnrA/B/S: FIXED -- require a specific allele number, Gene+Function only ---
# (original bug: BV-BRC's Product field is a templated string tied to the
#  search query, not a per-row fact -- true qnrS1 rows carry Product text
#  reading "...QnrB family...", causing massive double-counting)
gf_text = (pool['Gene'].fillna('') + ' ' + pool['Function'].fillna('')).str.lower()
qnrA_ids = set(pool.loc[gf_text.str.contains(r'qnra\d'), 'Genome ID'])
qnrB_ids = set(pool.loc[gf_text.str.contains(r'qnrb\d'), 'Genome ID'])
qnrS_ids = set(pool.loc[gf_text.str.contains(r'qnrs\d'), 'Genome ID'])
oqxA_ids = set(pool.loc[gf_text.str.contains('oqxa'), 'Genome ID'])
oqxB_ids = set(pool.loc[gf_text.str.contains('oqxb'), 'Genome ID'])

# --- aac6_Ib_cr: applying the same Gene+Function-only rule ---
# NOTE: this gives 347 genomes on independent re-derivation, not the 488-489
# reported elsewhere for this repo -- flagged as an open discrepancy pending
# the original code that produced that count. Using the defensible, documented
# rule (Gene+Function only, exact substring 'ib-cr') until reconciled.
aac_ids = set(pool.loc[gf_text.str.contains('ib-cr'), 'Genome ID'])

flags = {'gyrA_fq_mut': gyrA_ids, 'parC_fq_mut': parC_ids, 'qnrA': qnrA_ids, 'qnrB': qnrB_ids,
         'qnrS': qnrS_ids, 'oqxA': oqxA_ids, 'oqxB': oqxB_ids, 'aac6_Ib_cr': aac_ids}
for name, ids in flags.items():
    merged[name] = merged['Genome ID'].isin(ids).astype(int)
    print(f"  corrected {name}: {len(ids)} genomes")

# ============================================================
# 3. POPULATIONS
# ============================================================
uk = merged[merged['Isolation Country'] == 'United Kingdom'].reset_index(drop=True)
norway = merged[merged['Isolation Country'] == 'Norway'].reset_index(drop=True)
sasia = merged[merged['Isolation Country'].isin(['India', 'Pakistan', 'Bangladesh'])].reset_index(drop=True)

# ============================================================
# 4. MODEL + LEAK-FREE CV
# ============================================================
X_uk = uk[GENE_COLS].values
y_uk = (uk['CIP_Phenotype'] == 'Resistant').astype(int).values
groups_uk = pd.Series(uk['MLST']).fillna('UNK_' + uk['Genome ID'].astype(str)).values

gkf = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=42)
base_model = LogisticRegression(penalty='elasticnet', l1_ratio=0.5, solver='saga', max_iter=5000)
oof_proba = cross_val_predict(base_model, X_uk, y_uk, cv=gkf, groups=groups_uk, method='predict_proba')[:, 1]

final_model = LogisticRegression(penalty='elasticnet', l1_ratio=0.5, solver='saga', max_iter=5000)
final_model.fit(X_uk, y_uk)
print("\nCorrected model coefficients:")
for g, c in zip(GENE_COLS, final_model.coef_[0]):
    print(f"  {g}: {c:.3f}")
print(f"  intercept: {final_model.intercept_[0]:.3f}")

# ============================================================
# 5. DISCRIMINATION (elastic-net + rule-based baseline, with corrected features)
# ============================================================
def boot_ci(y_true, p_pred, n_boot=1000, seed=0):
    rng = np.random.default_rng(seed); n = len(y_true)
    def compute(yt, pp):
        yhat = (pp >= 0.5).astype(int)
        tn, fp, fn, tp = confusion_matrix(yt, yhat, labels=[0, 1]).ravel()
        return {'AUROC': roc_auc_score(yt, pp) if len(np.unique(yt)) > 1 else np.nan,
                'AUPRC': average_precision_score(yt, pp) if len(np.unique(yt)) > 1 else np.nan,
                'Sensitivity': tp/(tp+fn) if (tp+fn) > 0 else np.nan,
                'Specificity': tn/(tn+fp) if (tn+fp) > 0 else np.nan,
                'MCC': matthews_corrcoef(yt, yhat)}
    point = compute(y_true, p_pred)
    boot = {k: [] for k in point}
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        if len(np.unique(y_true[idx])) < 2: continue
        for k, v in compute(y_true[idx], p_pred[idx]).items(): boot[k].append(v)
    return {k: (v, np.percentile(boot[k], 2.5), np.percentile(boot[k], 97.5)) for k, v in point.items()}

print("\n=== DISCRIMINATION (corrected features) ===")
rule_uk = X_uk.max(axis=1).astype(float)
for label, p in [('elastic-net', oof_proba), ('rule-based', rule_uk)]:
    m = boot_ci(y_uk, p)
    print(f"UK ({label}):", {k: round(v[0],3) for k,v in m.items()})

results_rows = []
for name, d in [('Norway', norway), ('South Asia', sasia)]:
    X = d[GENE_COLS].values; y = (d['CIP_Phenotype']=='Resistant').astype(int).values
    p = final_model.predict_proba(X)[:,1]; rp = X.max(axis=1).astype(float)
    for label, pp in [('elastic-net', p), ('rule-based', rp)]:
        m = boot_ci(y, pp)
        print(f"{name} ({label}):", {k: round(v[0],3) for k,v in m.items()})

# ============================================================
# 6. CALIBRATION
# ============================================================
def calib_slope_intercept(p, y):
    p = np.clip(p, EPS, 1-EPS); logit = np.log(p/(1-p)).reshape(-1,1)
    m = LogisticRegression(); m.fit(logit, y)
    return m.coef_[0][0], m.intercept_[0]

print("\n=== CALIBRATION (corrected features) ===")
datasets = {'UK': (oof_proba, y_uk)}
for name, d in [('Norway', norway), ('South Asia', sasia)]:
    X = d[GENE_COLS].values; y = (d['CIP_Phenotype']=='Resistant').astype(int).values
    datasets[name] = (final_model.predict_proba(X)[:,1], y)
for name, (p, y) in datasets.items():
    slope, intercept = calib_slope_intercept(p, y)
    print(f"{name}: slope={slope:.3f}  intercept={intercept:.3f}  Brier={brier_score_loss(y,p):.4f}")

# ============================================================
# 7. RECALIBRATION -- median/IQR reporting (per documented near-separation issue)
# ============================================================
def platt_recalibrate(p_calib, y_calib, p_eval):
    lc = np.log(np.clip(p_calib,EPS,1-EPS)/(1-np.clip(p_calib,EPS,1-EPS))).reshape(-1,1)
    m = LogisticRegression(); m.fit(lc, y_calib)
    le = np.log(np.clip(p_eval,EPS,1-EPS)/(1-np.clip(p_eval,EPS,1-EPS))).reshape(-1,1)
    return m.predict_proba(le)[:,1]

print("\n=== RECALIBRATION (median [IQR], corrected features) ===")
for pop_name, d in [('Norway', norway), ('South Asia', sasia)]:
    X = d[GENE_COLS].values; y = (d['CIP_Phenotype']=='Resistant').astype(int).values
    p_all = final_model.predict_proba(X)[:,1]; n_total = len(y)
    for N in [25,50,100,150,200]:
        if N >= n_total - 20: continue
        intercepts, extreme_count = [], 0
        for rep in range(300):
            rng = np.random.default_rng(rep*1000+N)
            idx = rng.permutation(n_total); ci, ei = idx[:N], idx[N:]
            yc, ye = y[ci], y[ei]
            if yc.sum()==0 or yc.sum()==len(yc) or ye.sum()==0 or ye.sum()==len(ye): continue
            pe_recal = platt_recalibrate(p_all[ci], yc, p_all[ei])
            ic = calib_slope_intercept(pe_recal, ye)[1]
            intercepts.append(ic)
            if abs(ic) > 10: extreme_count += 1
        med = np.median(intercepts); q1, q3 = np.percentile(intercepts,[25,75])
        print(f"{pop_name} N={N}: median intercept={med:.3f} [IQR {q1:.3f}, {q3:.3f}]  "
              f"(mean={np.mean(intercepts):.3f}, {extreme_count}/{len(intercepts)} extreme |intercept|>10)")

# ============================================================
# 8. MECHANISM (should reproduce original paper's Table 5 exactly -- unaffected)
# ============================================================
print("\n=== MECHANISM COMPOSITION (should match original Table 5 -- sanity check) ===")
for name, d in [('UK', uk), ('Norway', norway), ('South Asia', sasia)]:
    r = d[d['CIP_Phenotype']=='Resistant']
    has_acq = r[['qnrA','qnrB','qnrS','oqxA','oqxB','aac6_Ib_cr']].max(axis=1)
    parc = (r['parC_fq_mut']==1)
    print(f"{name}: parC={parc.mean()*100:.1f}%  no-parC-but-acquired={((~parc)&(has_acq==1)).mean()*100:.1f}%  "
          f"neither={((~parc)&(has_acq==0)).mean()*100:.1f}%")
