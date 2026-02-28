#!/usr/bin/env python3
"""
Generate publication-quality figures for the indel spectrum analysis paper.
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch
import seaborn as sns
from scipy import stats
from scipy.cluster.hierarchy import linkage, dendrogram
from sklearn.decomposition import NMF
import warnings
warnings.filterwarnings('ignore')

plt.rcParams.update({
    'font.size': 10,
    'font.family': 'sans-serif',
    'axes.linewidth': 0.8,
    'xtick.major.width': 0.8,
    'ytick.major.width': 0.8,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.1,
})

RESULTS_DIR = "/home/yanlin/comp/cursor/analysis/results"
FIG_DIR = "/home/yanlin/comp/cursor/analysis/figures"

SUPERPOP_COLORS = {
    'AFR': '#E41A1C',
    'AMR': '#FF7F00',
    'EAS': '#4DAF4A',
    'EUR': '#377EB8',
    'SAS': '#984EA3',
}
SUPERPOP_NAMES = {
    'AFR': 'African',
    'AMR': 'American',
    'EAS': 'East Asian',
    'EUR': 'European',
    'SAS': 'South Asian',
}

# =========================================================================
# Load data
# =========================================================================
print("Loading data...")
df = pd.read_parquet(f"{RESULTS_DIR}/classified_indels_enhanced.parquet")
spectrum_raw = pd.read_csv(f"{RESULTS_DIR}/perpop_spectrum_chr1_2_raw.csv", index_col=0)
spectrum_norm = pd.read_csv(f"{RESULTS_DIR}/perpop_spectrum_chr1_2_norm.csv", index_col=0)

SIZE_BIN_ORDER = ["1bp", "2bp", "3bp", "4bp", "5bp", "6-10bp", "11-20bp", "21-50bp", ">50bp"]
SUPERPOPS = ['AFR', 'AMR', 'EAS', 'EUR', 'SAS']

# =========================================================================
# FIGURE 1: Global indel landscape overview
# =========================================================================
print("Generating Figure 1: Global indel landscape...")
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# 1a: Indel type and size distribution
ax = axes[0, 0]
del_sizes = df[df['type'] == 'DEL']['size_bin'].value_counts().reindex(SIZE_BIN_ORDER).fillna(0)
ins_sizes = df[df['type'] == 'INS']['size_bin'].value_counts().reindex(SIZE_BIN_ORDER).fillna(0)
x = np.arange(len(SIZE_BIN_ORDER))
w = 0.35
ax.bar(x - w/2, del_sizes.values / 1e6, w, color='#2166AC', label='Deletions', alpha=0.85)
ax.bar(x + w/2, ins_sizes.values / 1e6, w, color='#B2182B', label='Insertions', alpha=0.85)
ax.set_xticks(x)
ax.set_xticklabels(SIZE_BIN_ORDER, rotation=45, ha='right')
ax.set_ylabel('Count (millions)')
ax.set_xlabel('Indel size')
ax.set_title('a) Genome-wide indel size distribution', fontweight='bold', loc='left')
ax.legend(frameon=False)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# 1b: Allele frequency spectrum
ax = axes[0, 1]
af_bins = ['singleton-like', 'very_rare', 'rare', 'low_freq', 'common', 'major']
af_labels = ['<0.1%', '0.1-0.5%', '0.5-1%', '1-5%', '5-50%', '>50%']
del_af = df[df['type'] == 'DEL']['af_bin'].value_counts().reindex(af_bins).fillna(0)
ins_af = df[df['type'] == 'INS']['af_bin'].value_counts().reindex(af_bins).fillna(0)
x = np.arange(len(af_bins))
ax.bar(x - w/2, del_af.values / 1e6, w, color='#2166AC', label='Deletions', alpha=0.85)
ax.bar(x + w/2, ins_af.values / 1e6, w, color='#B2182B', label='Insertions', alpha=0.85)
ax.set_xticks(x)
ax.set_xticklabels(af_labels, rotation=45, ha='right')
ax.set_ylabel('Count (millions)')
ax.set_xlabel('Allele frequency')
ax.set_title('b) Allele frequency spectrum', fontweight='bold', loc='left')
ax.legend(frameon=False)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# 1c: DEL/INS ratio by size
ax = axes[1, 0]
ratios = del_sizes / ins_sizes
x_ratio = np.arange(len(ratios))
ax.bar(x_ratio, ratios.values, color='#4A1486', alpha=0.85)
ax.axhline(y=1, color='gray', linestyle='--', linewidth=0.8)
ax.set_xticks(x_ratio)
ax.set_xticklabels(SIZE_BIN_ORDER, rotation=45, ha='right')
ax.set_ylabel('DEL/INS ratio')
ax.set_xlabel('Indel size')
ax.set_title('c) Deletion-to-insertion ratio by size', fontweight='bold', loc='left')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# 1d: Homopolymer context for 1bp indels
ax = axes[1, 1]
mask_1bp = df['size'] == 1
hp_data = df.loc[mask_1bp, 'homopolymer_len']
hp_bins = [1, 2, 3, 4, 5, 6, 7, 8, 10, 15, 20, 30, 55]
hp_hist_del = np.histogram(df.loc[mask_1bp & (df['type']=='DEL'), 'homopolymer_len'], bins=hp_bins)[0]
hp_hist_ins = np.histogram(df.loc[mask_1bp & (df['type']=='INS'), 'homopolymer_len'], bins=hp_bins)[0]
hp_labels = ['1', '2', '3', '4', '5', '6', '7', '8-9', '10-14', '15-19', '20-29', '30+']
x = np.arange(len(hp_labels))
ax.bar(x - w/2, hp_hist_del / 1e5, w, color='#2166AC', label='Deletions', alpha=0.85)
ax.bar(x + w/2, hp_hist_ins / 1e5, w, color='#B2182B', label='Insertions', alpha=0.85)
ax.set_xticks(x)
ax.set_xticklabels(hp_labels, rotation=45, ha='right')
ax.set_ylabel('Count (×10⁵)')
ax.set_xlabel('Homopolymer run length')
ax.set_title('d) 1bp indels by homopolymer context', fontweight='bold', loc='left')
ax.legend(frameon=False)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
fig.savefig(f"{FIG_DIR}/fig1_global_landscape.pdf")
fig.savefig(f"{FIG_DIR}/fig1_global_landscape.png")
plt.close()
print("  Figure 1 saved.")

# =========================================================================
# FIGURE 2: Per-population indel spectrum comparison
# =========================================================================
print("Generating Figure 2: Population spectrum comparison...")
fig, axes = plt.subplots(2, 2, figsize=(14, 11))

# 2a: Heatmap of normalized spectrum
ax = axes[0, 0]
# Select most variable channels
top_channels = spectrum_norm.var().nlargest(20).index
heatmap_data = spectrum_norm[top_channels].T
sns.heatmap(heatmap_data, cmap='YlOrRd', ax=ax, cbar_kws={'label': 'Proportion'},
            xticklabels=[SUPERPOP_NAMES[sp] for sp in SUPERPOPS],
            linewidths=0.5)
ax.set_title('a) Indel spectrum across super-populations', fontweight='bold', loc='left')
ax.set_ylabel('Indel channel')

# 2b: Population-specific DEL/INS ratio
ax = axes[0, 1]
del_counts = spectrum_raw.filter(like='DEL_').sum(axis=1)
ins_counts = spectrum_raw.filter(like='INS_').sum(axis=1)
ratios_pop = del_counts / ins_counts
bars = ax.bar(range(5), ratios_pop.values, color=[SUPERPOP_COLORS[sp] for sp in SUPERPOPS], alpha=0.85)
ax.set_xticks(range(5))
ax.set_xticklabels([SUPERPOP_NAMES[sp] for sp in SUPERPOPS], rotation=45, ha='right')
ax.set_ylabel('DEL/INS ratio')
ax.set_title('b) DEL/INS ratio by super-population', fontweight='bold', loc='left')
ax.axhline(y=1, color='gray', linestyle='--', linewidth=0.8)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# 2c: Population-specific SFS
ax = axes[1, 0]
sfs_df = pd.read_csv(f"{RESULTS_DIR}/perpop_sfs.csv", index_col=0)
sfs_norm = sfs_df.div(sfs_df.sum(axis=0), axis=1)
sfs_order = ['rare', 'low_freq', 'common', 'high_freq']
sfs_labels = ['Rare\n(<1%)', 'Low freq\n(1-5%)', 'Common\n(5-50%)', 'High freq\n(>50%)']
x = np.arange(len(sfs_order))
width = 0.15
for i, spop in enumerate(SUPERPOPS):
    vals = sfs_norm.loc[sfs_order, spop].values
    ax.bar(x + i*width - 2*width, vals, width, color=SUPERPOP_COLORS[spop],
           label=SUPERPOP_NAMES[spop], alpha=0.85)
ax.set_xticks(x)
ax.set_xticklabels(sfs_labels)
ax.set_ylabel('Proportion of indels')
ax.set_title('c) Site frequency spectrum by population', fontweight='bold', loc='left')
ax.legend(frameon=False, fontsize=8, loc='upper left')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# 2d: Hierarchical clustering of populations based on spectrum
ax = axes[1, 1]
# Use all channels for clustering
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, dendrogram
dist_matrix = pdist(spectrum_norm.values, metric='cosine')
Z = linkage(dist_matrix, method='average')
dn = dendrogram(Z, labels=[SUPERPOP_NAMES[sp] for sp in SUPERPOPS], ax=ax,
                leaf_rotation=45, leaf_font_size=10)
ax.set_title('d) Population clustering by indel spectrum', fontweight='bold', loc='left')
ax.set_ylabel('Cosine distance')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
fig.savefig(f"{FIG_DIR}/fig2_population_comparison.pdf")
fig.savefig(f"{FIG_DIR}/fig2_population_comparison.png")
plt.close()
print("  Figure 2 saved.")

# =========================================================================
# FIGURE 3: NMF signature decomposition
# =========================================================================
print("Generating Figure 3: NMF decomposition...")

# Build spectrum matrix for NMF (populations x channels)
# Use enhanced channels for richer decomposition
enhanced_spectrum = {}
for spop in SUPERPOPS:
    enhanced_spectrum[spop] = {}

# For NMF, we need count data per population
# Using chr1 data from the per-pop spectrum
# Build spectrum from raw counts
nmf_matrix = spectrum_raw.copy()
nmf_matrix = nmf_matrix.loc[:, nmf_matrix.sum() > 0]

# Try multiple ranks (max k = min(n_samples, n_features) - 1)
max_k = min(nmf_matrix.shape[0], nmf_matrix.shape[1])
reconstruction_errors = []
k_range = range(2, max_k)
for k in k_range:
    model = NMF(n_components=k, init='random', random_state=42, max_iter=500)
    W = model.fit_transform(nmf_matrix.values)
    reconstruction_errors.append(model.reconstruction_err_)

fig, axes = plt.subplots(2, 2, figsize=(14, 11))

# 3a: Reconstruction error vs k
ax = axes[0, 0]
ax.plot(list(k_range), reconstruction_errors, 'o-', color='#333333', markersize=6)
ax.set_xlabel('Number of signatures (k)')
ax.set_ylabel('Reconstruction error')
ax.set_title('a) NMF model selection', fontweight='bold', loc='left')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Use k=3 (good balance)
k_best = 3
model = NMF(n_components=k_best, init='random', random_state=42, max_iter=1000)
W = model.fit_transform(nmf_matrix.values)
H = model.components_

# Normalize signatures
H_norm = H / H.sum(axis=1, keepdims=True)
W_norm = W / W.sum(axis=1, keepdims=True)

# 3b: Signature profiles
ax = axes[0, 1]
channels = nmf_matrix.columns
n_ch = len(channels)
x = np.arange(n_ch)
sig_colors = ['#E41A1C', '#377EB8', '#4DAF4A']
for i in range(k_best):
    ax.plot(x, H_norm[i], color=sig_colors[i], alpha=0.8, linewidth=1.5,
            label=f'Signature {i+1}')
ax.set_xlabel('Indel channel')
ax.set_ylabel('Proportion')
ax.set_title('b) Extracted germline indel signatures', fontweight='bold', loc='left')
ax.legend(frameon=False, fontsize=8)
ax.set_xticks(x[::5])
ax.set_xticklabels([channels[i][:15] for i in range(0, n_ch, 5)], rotation=90, fontsize=6)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# 3c: Signature contributions per population (stacked bar)
ax = axes[1, 0]
bottom = np.zeros(5)
for i in range(k_best):
    ax.bar(range(5), W_norm[:, i], bottom=bottom, color=sig_colors[i],
           label=f'Signature {i+1}', alpha=0.85)
    bottom += W_norm[:, i]
ax.set_xticks(range(5))
ax.set_xticklabels([SUPERPOP_NAMES[sp] for sp in SUPERPOPS], rotation=45, ha='right')
ax.set_ylabel('Signature contribution')
ax.set_title('c) Population-specific signature weights', fontweight='bold', loc='left')
ax.legend(frameon=False, fontsize=8)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# 3d: Top channels per signature
ax = axes[1, 1]
sig_data = []
for i in range(k_best):
    top_idx = np.argsort(H_norm[i])[-5:][::-1]
    for idx in top_idx:
        sig_data.append({
            'Signature': f'Sig {i+1}',
            'Channel': channels[idx],
            'Weight': H_norm[i, idx]
        })
sig_df = pd.DataFrame(sig_data)

# Bar plot of top channels
for i, sig in enumerate([f'Sig {j+1}' for j in range(k_best)]):
    subset = sig_df[sig_df['Signature'] == sig]
    y_pos = np.arange(5) + i * 6
    ax.barh(y_pos, subset['Weight'], color=sig_colors[i], alpha=0.85)
    for j, (_, row) in enumerate(subset.iterrows()):
        ax.text(row['Weight'] + 0.005, y_pos[j], row['Channel'], fontsize=7, va='center')
    ax.text(-0.02, y_pos[2], sig, fontsize=9, fontweight='bold', va='center', ha='right',
            color=sig_colors[i])

ax.set_xlabel('Weight')
ax.set_title('d) Top channels per signature', fontweight='bold', loc='left')
ax.set_yticks([])
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)

plt.tight_layout()
fig.savefig(f"{FIG_DIR}/fig3_nmf_signatures.pdf")
fig.savefig(f"{FIG_DIR}/fig3_nmf_signatures.png")
plt.close()
print("  Figure 3 saved.")

# =========================================================================
# FIGURE 4: Chromosome-level analysis 
# =========================================================================
print("Generating Figure 4: Chromosomal patterns...")
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# 4a: Indel density per chromosome
ax = axes[0, 0]
chrom_counts = df.groupby('chrom').size().reset_index(name='count')
chrom_counts['chrom_num'] = chrom_counts['chrom'].str.replace('chr', '').astype(int)
chrom_counts = chrom_counts.sort_values('chrom_num')

# Get chromosome sizes (approximate)
chrom_sizes = {
    1: 248, 2: 242, 3: 198, 4: 190, 5: 181, 6: 170, 7: 159, 8: 145,
    9: 138, 10: 133, 11: 135, 12: 133, 13: 114, 14: 107, 15: 101,
    16: 90, 17: 83, 18: 80, 19: 58, 20: 64, 21: 46, 22: 50
}
chrom_counts['size_mb'] = chrom_counts['chrom_num'].map(chrom_sizes)
chrom_counts['density'] = chrom_counts['count'] / chrom_counts['size_mb']

ax.bar(range(22), chrom_counts['density'], color='#2166AC', alpha=0.85)
ax.set_xticks(range(22))
ax.set_xticklabels([str(i) for i in range(1, 23)], fontsize=8)
ax.set_xlabel('Chromosome')
ax.set_ylabel('Indels per Mb')
ax.set_title('a) Indel density per chromosome', fontweight='bold', loc='left')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# 4b: DEL/INS ratio per chromosome
ax = axes[0, 1]
chrom_type = df.groupby(['chrom', 'type']).size().unstack(fill_value=0)
chrom_type['chrom_num'] = chrom_type.index.str.replace('chr', '').astype(int)
chrom_type = chrom_type.sort_values('chrom_num')
chrom_ratio = chrom_type['DEL'] / chrom_type['INS']

ax.bar(range(22), chrom_ratio.values, color='#4A1486', alpha=0.85)
ax.axhline(y=chrom_ratio.mean(), color='red', linestyle='--', linewidth=0.8, label='Mean')
ax.set_xticks(range(22))
ax.set_xticklabels([str(i) for i in range(1, 23)], fontsize=8)
ax.set_xlabel('Chromosome')
ax.set_ylabel('DEL/INS ratio')
ax.set_title('b) DEL/INS ratio per chromosome', fontweight='bold', loc='left')
ax.legend(frameon=False, fontsize=8)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# 4c: Allele frequency spectrum by indel size
ax = axes[1, 0]
size_af = df.groupby(['size_bin', 'af_bin']).size().unstack(fill_value=0)
size_af = size_af.reindex(SIZE_BIN_ORDER)
af_order = ['singleton-like', 'very_rare', 'rare', 'low_freq', 'common', 'major']
size_af = size_af.reindex(columns=af_order)
size_af_norm = size_af.div(size_af.sum(axis=1), axis=0)

af_colors = ['#FDE725', '#7AD151', '#22A884', '#2A788E', '#414487', '#440154']
bottom = np.zeros(len(SIZE_BIN_ORDER))
x = np.arange(len(SIZE_BIN_ORDER))
for i, af_bin in enumerate(af_order):
    vals = size_af_norm[af_bin].values
    ax.bar(x, vals, bottom=bottom, color=af_colors[i], label=af_bin, alpha=0.85)
    bottom += vals
ax.set_xticks(x)
ax.set_xticklabels(SIZE_BIN_ORDER, rotation=45, ha='right')
ax.set_ylabel('Proportion')
ax.set_xlabel('Indel size')
ax.set_title('c) Allele frequency spectrum by indel size', fontweight='bold', loc='left')
ax.legend(frameon=False, fontsize=7, loc='center left', bbox_to_anchor=(0.0, 0.5),
          title='AF category')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# 4d: Base composition by size
ax = axes[1, 1]
size_base = df.groupby(['size_bin', 'base_comp']).size().unstack(fill_value=0)
size_base = size_base.reindex(SIZE_BIN_ORDER)
base_comps = ['AT', 'GC', 'AT-rich', 'mixed', 'GC-rich']
size_base = size_base.reindex(columns=[c for c in base_comps if c in size_base.columns])
size_base_norm = size_base.div(size_base.sum(axis=1), axis=0)

base_colors = {'AT': '#E41A1C', 'GC': '#377EB8', 'AT-rich': '#FF7F00',
               'mixed': '#999999', 'GC-rich': '#4DAF4A'}
bottom = np.zeros(len(SIZE_BIN_ORDER))
x = np.arange(len(SIZE_BIN_ORDER))
for col in size_base_norm.columns:
    vals = size_base_norm[col].values
    ax.bar(x, vals, bottom=bottom, color=base_colors.get(col, '#333333'),
           label=col, alpha=0.85)
    bottom += vals
ax.set_xticks(x)
ax.set_xticklabels(SIZE_BIN_ORDER, rotation=45, ha='right')
ax.set_ylabel('Proportion')
ax.set_xlabel('Indel size')
ax.set_title('d) Base composition by indel size', fontweight='bold', loc='left')
ax.legend(frameon=False, fontsize=7)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
fig.savefig(f"{FIG_DIR}/fig4_chromosomal_patterns.pdf")
fig.savefig(f"{FIG_DIR}/fig4_chromosomal_patterns.png")
plt.close()
print("  Figure 4 saved.")

# =========================================================================
# FIGURE 5: Population differentiation detail
# =========================================================================
print("Generating Figure 5: Population differentiation...")
fig, axes = plt.subplots(2, 2, figsize=(14, 11))

# 5a: Radar plot of spectrum proportions for key channels
ax = axes[0, 0]
key_channels = ['DEL_1bp_AT', 'DEL_1bp_GC', 'INS_1bp_AT', 'INS_1bp_GC',
                'DEL_2bp_AT-rich', 'INS_2bp_AT-rich', 'DEL_4bp_AT-rich',
                'DEL_6-10bp_AT-rich', 'DEL_11-20bp_mixed']
key_data = spectrum_norm[key_channels]
# Grouped bar chart instead of radar
x = np.arange(len(key_channels))
width = 0.15
for i, spop in enumerate(SUPERPOPS):
    ax.bar(x + i*width - 2*width, key_data.loc[spop].values, width,
           color=SUPERPOP_COLORS[spop], label=SUPERPOP_NAMES[spop], alpha=0.85)
ax.set_xticks(x)
ax.set_xticklabels(key_channels, rotation=60, ha='right', fontsize=7)
ax.set_ylabel('Proportion')
ax.set_title('a) Key channel proportions by population', fontweight='bold', loc='left')
ax.legend(frameon=False, fontsize=7)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# 5b: Population-private indel size distribution
ax = axes[0, 1]
# Compute from spectrum: private variants are those where only one pop has AC > 0
# Use proportional differences as proxy
diff_from_mean = spectrum_norm.subtract(spectrum_norm.mean(axis=0), axis=1)
top_enriched = diff_from_mean.abs().max(axis=0).nlargest(10).index

for i, spop in enumerate(SUPERPOPS):
    vals = diff_from_mean.loc[spop, top_enriched].values * 100
    ax.barh(np.arange(len(top_enriched)) + i*0.15 - 0.3, vals, 0.15,
            color=SUPERPOP_COLORS[spop], label=SUPERPOP_NAMES[spop], alpha=0.85)
ax.set_yticks(range(len(top_enriched)))
ax.set_yticklabels(top_enriched, fontsize=7)
ax.axvline(x=0, color='gray', linewidth=0.8)
ax.set_xlabel('Deviation from mean (%)')
ax.set_title('b) Population deviation from global mean', fontweight='bold', loc='left')
ax.legend(frameon=False, fontsize=7)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# 5c: Cosine similarity heatmap between populations
ax = axes[1, 0]
from scipy.spatial.distance import cosine
cos_sim = np.zeros((5, 5))
for i in range(5):
    for j in range(5):
        cos_sim[i, j] = 1 - cosine(spectrum_norm.iloc[i], spectrum_norm.iloc[j])
cos_df = pd.DataFrame(cos_sim, index=[SUPERPOP_NAMES[sp] for sp in SUPERPOPS],
                       columns=[SUPERPOP_NAMES[sp] for sp in SUPERPOPS])
sns.heatmap(cos_df, annot=True, fmt='.6f', cmap='RdYlGn', ax=ax,
            vmin=cos_sim.min() - 0.0001, vmax=1.0, linewidths=0.5)
ax.set_title('c) Cosine similarity of indel spectra', fontweight='bold', loc='left')

# 5d: Size-dependent DEL/INS ratio across populations
ax = axes[1, 1]
for spop in SUPERPOPS:
    del_by_size = {}
    ins_by_size = {}
    for ch in spectrum_raw.columns:
        parts = ch.split('_')
        itype = parts[0]
        sbin = parts[1]
        if itype == 'DEL':
            del_by_size[sbin] = del_by_size.get(sbin, 0) + spectrum_raw.loc[spop, ch]
        else:
            ins_by_size[sbin] = ins_by_size.get(sbin, 0) + spectrum_raw.loc[spop, ch]
    
    ratios = []
    for sb in SIZE_BIN_ORDER:
        d = del_by_size.get(sb, 0)
        i = ins_by_size.get(sb, 0)
        ratios.append(d / i if i > 0 else np.nan)
    
    ax.plot(range(len(SIZE_BIN_ORDER)), ratios, 'o-', color=SUPERPOP_COLORS[spop],
            label=SUPERPOP_NAMES[spop], markersize=5, alpha=0.85)

ax.axhline(y=1, color='gray', linestyle='--', linewidth=0.8)
ax.set_xticks(range(len(SIZE_BIN_ORDER)))
ax.set_xticklabels(SIZE_BIN_ORDER, rotation=45, ha='right')
ax.set_ylabel('DEL/INS ratio')
ax.set_xlabel('Indel size')
ax.set_title('d) Size-dependent DEL/INS ratio by population', fontweight='bold', loc='left')
ax.legend(frameon=False, fontsize=7)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
fig.savefig(f"{FIG_DIR}/fig5_population_differentiation.pdf")
fig.savefig(f"{FIG_DIR}/fig5_population_differentiation.png")
plt.close()
print("  Figure 5 saved.")

print("\nAll figures generated successfully!")
