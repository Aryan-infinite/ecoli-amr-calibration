"""
Claims audit behind statements in the Methods and Results (blocks H1 to H4, H5a to H5c).

Run inside the analysis notebook session after the setup cell, so that these exist:
base (repository folder), qc_cohort, uk, norway, sasia, uk_model, feature_cols.
"""
import numpy as np, pandas as pd
from sklearn.model_selection import GroupKFold, StratifiedGroupKFold

pops = {'UK': uk, 'Norway': norway, 'South Asia': sasia}

# ---- H1. BioProject: missingness, top-accession share, overlap between populations; genome-ID overlap ----
print('=== H1: BioProject and identifier overlap ===')
bp = 'BioProject Accession'; acc = {}
for name, df in pops.items():
    s = df[bp]; vc = s.value_counts()
    print(f'{name:11s} missing {100 * s.isna().mean():5.1f}% | non-missing {int(s.notna().sum()):5d} | top accession {vc.index[0] if len(vc) else "none"} = {100 * vc.iloc[0] / max(1, int(s.notna().sum())):.1f}% of non-missing')
    acc[name] = set(s.dropna())
for a, b in [('UK', 'Norway'), ('UK', 'South Asia'), ('Norway', 'South Asia')]:
    print(f'shared BioProject accessions {a} / {b}: {len(acc[a] & acc[b])} | shared Genome IDs: {len(set(pops[a]["Genome ID"]) & set(pops[b]["Genome ID"]))}')

# ---- H2. Fold class balance when grouping by MLST ----
print('\n=== H2: % resistant per fold, UK grouped by MLST ===')
X_uk = uk_model[feature_cols].values; y_uk = uk_model['y'].values; groups_uk = uk_model['MLST'].values
for label, cv in [('GroupKFold', GroupKFold(n_splits=5)), ('StratifiedGroupKFold', StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=42))]:
    print(f'{label:22s}', [float(round(100 * y_uk[te].mean(), 1)) for tr, te in cv.split(X_uk, y_uk, groups_uk)])

# ---- H3. The qnrA example in Methods 2.4 (Gene = pipB2_2 with a QnrA1 identity in Function) ----
print('\n=== H3: qnrA rows whose Gene field contains pipB2 ===')
qa = pd.read_csv(f'{base}/data/raw/gene_qnrA.csv', dtype={'Genome ID': str})
qa['Gene'] = qa['Gene'].fillna(''); qa['Function'] = qa['Function'].fillna('')
x = qa[qa['Gene'].str.contains('pipB2', case=False)]
print(len(x), 'rows'); print(x.groupby(['Gene', 'Function']).size().head(8).to_string())

# ---- H4. gyrA carriage: initial (class field only) versus final (POINTP and class exactly "quinolone"), and the cohort-wide parC and gyrA figures ----
print('\n=== H4: gyrA and parC carriage ===')
g = pd.read_csv(f'{base}/data/raw/gene_gyrA.csv', dtype={'Genome ID': str})
cls = g['Antibiotics Class'].fillna('')
sets = {'initial': set(g.loc[cls.str.contains('quinolone', case=False), 'Genome ID']),
        'final': set(g.loc[(g['Evidence'] == 'AMRFinderPlus: POINTP') & (cls == 'quinolone'), 'Genome ID'])}
for label, ids_ in sets.items():
    for scope, df in [('cohort', qc_cohort), ('UK', uk), ('Norway', norway)]:
        f = df['Genome ID'].isin(ids_).astype(int)
        print(f'gyrA {label:8s} {scope:7s} susceptible {100 * f[df["y"] == 0].mean():6.2f}% | resistant {100 * f[df["y"] == 1].mean():6.2f}%')
for feat in ['gyrA_fq_mut', 'parC_fq_mut']:
    print(f'{feat} (saved feature) cohort: susceptible {100 * qc_cohort.loc[qc_cohort.y == 0, feat].mean():.2f}% | resistant {100 * qc_cohort.loc[qc_cohort.y == 1, feat].mean():.2f}%')

def mech(r):
    if r['parC_fq_mut'] == 1: return 'parC_mutation'
    if any(r[g] == 1 for g in ['qnrA', 'qnrB', 'qnrS', 'oqxA', 'oqxB', 'aac6_Ib_cr']): return 'acquired_no_parC'
    return 'neither'

# ---- H5a and H5b. South Asia: BioProject and collection year by phenotype ----
d = sasia.copy(); d['BioProject'] = d['BioProject Accession'].fillna('missing')
t = d.groupby('BioProject').agg(n=('y', 'size'), resistant=('y', 'sum'), first_year=('Collection Year', 'min'), last_year=('Collection Year', 'max')).sort_values('n', ascending=False)
print('=== H5a: South Asia by BioProject ===\n', t.head(10).to_string())
print('\n=== H5b: South Asia, collection year by phenotype ===\n', pd.crosstab(d['Collection Year'], d['y']).rename(columns={0: 'susceptible', 1: 'resistant'}).to_string())

# ---- H5c. Norway: mechanism category of the resistant genomes by collection year ----
r = norway[norway['y'] == 1].copy(); r['mechanism'] = r.apply(mech, axis=1)
print('\n=== H5c: Norway resistant genomes, collection year by mechanism ===\n', pd.crosstab(r['Collection Year'], r['mechanism']).to_string())
