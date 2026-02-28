#!/usr/bin/env python3
"""
Comprehensive statistical analysis for indel spectrum paper.
"""

import pandas as pd
import numpy as np
from scipy import stats
from scipy.spatial.distance import pdist, squareform, cosine
from sklearn.decomposition import PCA, NMF
from sklearn.preprocessing import StandardScaler
import warnings
warnings.filterwarnings('ignore')

RESULTS_DIR = "/home/yanlin/comp/cursor/analysis/results"

SUPERPOPS = ['AFR', 'AMR', 'EAS', 'EUR', 'SAS']
SUPERPOP_NAMES = {
    'AFR': 'African', 'AMR': 'American', 'EAS': 'East Asian',
    'EUR': 'European', 'SAS': 'South Asian',
}

print("=" * 70)
print("STATISTICAL ANALYSIS FOR GERMLINE INDEL SPECTRUM PAPER")
print("=" * 70)

# Load data
df = pd.read_parquet(f"{RESULTS_DIR}/classified_indels_enhanced.parquet")
spectrum_raw = pd.read_csv(f"{RESULTS_DIR}/perpop_spectrum_chr1_2_raw.csv", index_col=0)
spectrum_norm = pd.read_csv(f"{RESULTS_DIR}/perpop_spectrum_chr1_2_norm.csv", index_col=0)

# =========================================================================
# 1. Chi-square test: Indel spectrum independence across populations
# =========================================================================
print("\n" + "=" * 50)
print("1. CHI-SQUARE TEST OF INDEPENDENCE")
print("=" * 50)

chi2, pval, dof, expected = stats.chi2_contingency(spectrum_raw.values)
print(f"Chi-square statistic: {chi2:.2f}")
print(f"Degrees of freedom: {dof}")
print(f"P-value: {pval:.2e}")
print(f"Interpretation: {'SIGNIFICANT' if pval < 0.05 else 'Not significant'}")

# Effect size (Cramer's V)
n = spectrum_raw.values.sum()
k = min(spectrum_raw.shape) - 1
cramers_v = np.sqrt(chi2 / (n * k))
print(f"Cramer's V (effect size): {cramers_v:.6f}")

# =========================================================================
# 2. Pairwise population comparisons
# =========================================================================
print("\n" + "=" * 50)
print("2. PAIRWISE POPULATION COMPARISONS")
print("=" * 50)

# Cosine similarity
print("\nCosine similarity between population spectra:")
for i in range(5):
    for j in range(i+1, 5):
        sim = 1 - cosine(spectrum_norm.iloc[i], spectrum_norm.iloc[j])
        print(f"  {SUPERPOPS[i]} vs {SUPERPOPS[j]}: {sim:.8f}")

# Jensen-Shannon divergence
from scipy.spatial.distance import jensenshannon
print("\nJensen-Shannon divergence:")
for i in range(5):
    for j in range(i+1, 5):
        jsd = jensenshannon(spectrum_norm.iloc[i], spectrum_norm.iloc[j])
        print(f"  {SUPERPOPS[i]} vs {SUPERPOPS[j]}: {jsd:.8f}")

# =========================================================================
# 3. DEL/INS ratio analysis
# =========================================================================
print("\n" + "=" * 50)
print("3. DEL/INS RATIO ANALYSIS")
print("=" * 50)

del_counts = spectrum_raw.filter(like='DEL_').sum(axis=1)
ins_counts = spectrum_raw.filter(like='INS_').sum(axis=1)
ratios = del_counts / ins_counts

print("\nDEL/INS ratios per super-population:")
for spop in SUPERPOPS:
    print(f"  {spop}: {ratios[spop]:.4f} (DEL={del_counts[spop]:,}, INS={ins_counts[spop]:,})")

# Binomial test for each population
print("\nBinomial test (H0: DEL/INS = 1):")
for spop in SUPERPOPS:
    n_total = del_counts[spop] + ins_counts[spop]
    p_val = stats.binomtest(del_counts[spop], n_total, 0.5).pvalue
    print(f"  {spop}: p = {p_val:.2e}")

# =========================================================================
# 4. Channel-specific population differentiation
# =========================================================================
print("\n" + "=" * 50)
print("4. CHANNEL-SPECIFIC DIFFERENTIATION")
print("=" * 50)

# For each channel, test if proportions differ across populations using chi-square
channel_tests = []
for ch in spectrum_raw.columns:
    observed = spectrum_raw[ch].values
    expected_ch = observed.sum() / 5
    if expected_ch > 5:
        chi2_ch = np.sum((observed - expected_ch)**2 / expected_ch)
        p_ch = stats.chi2.sf(chi2_ch, df=4)
        channel_tests.append({
            'channel': ch,
            'chi2': chi2_ch,
            'p_value': p_ch,
            'max_prop': spectrum_norm[ch].max(),
            'min_prop': spectrum_norm[ch].min(),
            'range': spectrum_norm[ch].max() - spectrum_norm[ch].min(),
        })

channel_test_df = pd.DataFrame(channel_tests).sort_values('chi2', ascending=False)

# Bonferroni correction
n_tests = len(channel_test_df)
channel_test_df['p_adjusted'] = channel_test_df['p_value'] * n_tests
channel_test_df['significant'] = channel_test_df['p_adjusted'] < 0.05

print(f"\nSignificant channels (Bonferroni-corrected α=0.05): "
      f"{channel_test_df['significant'].sum()} / {n_tests}")
print("\nTop 15 differentiated channels:")
print(channel_test_df.head(15)[['channel', 'chi2', 'p_adjusted', 'range', 'significant']].to_string())

channel_test_df.to_csv(f"{RESULTS_DIR}/channel_differentiation_tests.csv", index=False)

# =========================================================================
# 5. Global indel statistics
# =========================================================================
print("\n" + "=" * 50)
print("5. GLOBAL INDEL STATISTICS")
print("=" * 50)

total = len(df)
n_del = (df['type'] == 'DEL').sum()
n_ins = (df['type'] == 'INS').sum()

print(f"Total indels genome-wide: {total:,}")
print(f"  Deletions: {n_del:,} ({100*n_del/total:.1f}%)")
print(f"  Insertions: {n_ins:,} ({100*n_ins/total:.1f}%)")
print(f"  DEL/INS ratio: {n_del/n_ins:.4f}")

print(f"\nSize statistics:")
print(f"  Mean size: {df['size'].mean():.2f} bp")
print(f"  Median size: {df['size'].median():.0f} bp")
print(f"  Max size: {df['size'].max()} bp")

print(f"\n1bp indels: {(df['size']==1).sum():,} ({100*(df['size']==1).sum()/total:.1f}%)")
print(f"  In homopolymer (hp≥2): {(df.loc[df['size']==1, 'homopolymer_len']>=2).sum():,}")
print(f"  Long homopolymer (hp≥6): {(df.loc[df['size']==1, 'homopolymer_len']>=6).sum():,}")

print(f"\nAllele frequency distribution:")
for af_bin in ['singleton-like', 'very_rare', 'rare', 'low_freq', 'common', 'major']:
    n = (df['af_bin'] == af_bin).sum()
    print(f"  {af_bin}: {n:,} ({100*n/total:.1f}%)")

# =========================================================================
# 6. Population-specific rare variant enrichment
# =========================================================================
print("\n" + "=" * 50)
print("6. POPULATION-SPECIFIC PATTERNS")
print("=" * 50)

# Key finding: AFR-specific enrichment in DEL_1bp_GC
# EAS-specific enrichment in INS_1bp_AT
key_channels = ['INS_1bp_AT', 'DEL_1bp_AT', 'DEL_1bp_GC', 'INS_1bp_GC',
                'DEL_4bp_AT-rich', 'DEL_3bp_mixed']

print("Key channel proportions across populations:")
for ch in key_channels:
    if ch in spectrum_norm.columns:
        vals = spectrum_norm[ch]
        max_pop = vals.idxmax()
        min_pop = vals.idxmin()
        print(f"\n  {ch}:")
        for sp in SUPERPOPS:
            indicator = " ← MAX" if sp == max_pop else (" ← MIN" if sp == min_pop else "")
            print(f"    {sp}: {vals[sp]:.5f}{indicator}")
        fold = vals[max_pop] / vals[min_pop]
        print(f"    Max/Min fold change: {fold:.3f}")

# =========================================================================
# 7. NMF signatures (detailed)
# =========================================================================
print("\n" + "=" * 50)
print("7. NMF SIGNATURE ANALYSIS")
print("=" * 50)

for k in [2, 3, 4]:
    model = NMF(n_components=k, init='random', random_state=42, max_iter=1000)
    W = model.fit_transform(spectrum_raw.values)
    H = model.components_
    err = model.reconstruction_err_
    
    W_norm = W / W.sum(axis=1, keepdims=True)
    H_norm = H / H.sum(axis=1, keepdims=True)
    
    print(f"\nk={k}: Reconstruction error = {err:.2f}")
    print("  Signature weights per population:")
    for i, spop in enumerate(SUPERPOPS):
        weights = " | ".join(f"S{j+1}={W_norm[i,j]:.3f}" for j in range(k))
        print(f"    {spop}: {weights}")
    
    print("  Top 3 channels per signature:")
    for j in range(k):
        top_idx = np.argsort(H_norm[j])[-3:][::-1]
        channels_str = ", ".join(f"{spectrum_raw.columns[idx]}({H_norm[j,idx]:.3f})"
                                 for idx in top_idx)
        print(f"    Sig{j+1}: {channels_str}")

# Save summary statistics
summary = {
    'total_indels': total,
    'n_deletions': n_del,
    'n_insertions': n_ins,
    'del_ins_ratio': n_del / n_ins,
    'mean_size': df['size'].mean(),
    'median_size': df['size'].median(),
    'pct_1bp': 100 * (df['size'] == 1).sum() / total,
    'pct_in_homopolymer': 100 * (df.loc[df['size']==1, 'homopolymer_len'] >= 2).sum() / (df['size']==1).sum(),
    'chi2_spectrum': chi2,
    'chi2_pvalue': pval,
    'cramers_v': cramers_v,
}
pd.Series(summary).to_csv(f"{RESULTS_DIR}/summary_statistics.csv")

print("\n" + "=" * 70)
print("ANALYSIS COMPLETE")
print("=" * 70)
