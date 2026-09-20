"""
Confirmation checks behind results/RECONCILIATION.md.

Run inside the analysis notebook session, after the model variables exist:
base (repository folder), qc_cohort, clf_full, norway, sasia, feature_cols.
"""
import glob, numpy as np, pandas as pd
from sklearn.linear_model import LogisticRegression as LR

# ---- V1. AAC(6')-Ib-cr: count with all rows versus after the re-implementation's de-duplication ----
pool = pd.concat([pd.read_csv(f, dtype={'Genome ID': str}) for f in sorted(glob.glob(f'{base}/data/raw/gene_*.csv'))], ignore_index=True)
def ibcr_genomes(df):
    t = (df['Gene'].fillna('') + ' ' + df['Function'].fillna('')).str.lower()
    return set(df.loc[t.str.contains('ib-cr'), 'Genome ID'])
dd = pool.drop_duplicates(subset=['Genome ID', 'BRC ID', 'Source'])
a, b = ibcr_genomes(pool), ibcr_genomes(dd)
print("V1  AAC(6')-Ib-cr genomes, all rows:", len(a), '| after the rebuild de-duplication:', len(b), '| genomes lost:', len(a - b), '| rows', len(pool), '->', len(dd))

# ---- V2. Recalibration: identical splits, unpenalized (C=1e6, as in Section 2.9) versus ridge (C=1, scikit-learn default) ----
EPS = 1e-6
def lg(p): p = np.clip(p, EPS, 1 - EPS); return np.log(p / (1 - p)).reshape(-1, 1)
def recal_experiment(pop, C, Ns=(25, 50, 100, 150, 200), reps=300):
    y = pop['y'].values; p = clf_full.predict_proba(pop[feature_cols].values)[:, 1]; n = len(y); rows = []
    for N in Ns:
        vals = []
        for rep in range(reps):
            rng = np.random.default_rng(rep * 1000 + N)
            idx = rng.permutation(n); ci, ei = idx[:N], idx[N:]
            yc, ye = y[ci], y[ei]
            if yc.sum() in (0, len(yc)) or ye.sum() in (0, len(ye)): continue
            m = LR(C=C, max_iter=1000).fit(lg(p[ci]), yc)
            pe = m.predict_proba(lg(p[ei]))[:, 1]
            vals.append(LR(C=C, max_iter=1000).fit(lg(pe), ye).intercept_[0])
        v = np.array(vals)
        rows.append(dict(N=N, reps=len(v), median=np.median(v), q25=np.percentile(v, 25), q75=np.percentile(v, 75),
                         abs_gt5=int((abs(v) > 5).sum()), abs_gt10=int((abs(v) > 10).sum()), min=v.min()))
    return pd.DataFrame(rows).round(3)
for name, pop in [('Norway', norway), ('South Asia', sasia)]:
    for C in (1e6, 1.0):
        print(f'\nV2  {name}, C={C:g}'); print(recal_experiment(pop, C).to_string(index=False))
