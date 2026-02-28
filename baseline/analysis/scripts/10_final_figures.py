#!/usr/bin/env python3
"""
Generate final publication-quality figures with all available data.
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
from scipy import stats
from scipy.spatial.distance import pdist, squareform, cosine, jensenshannon
from scipy.cluster.hierarchy import linkage, dendrogram
from sklearn.decomposition import NMF
import warnings
warnings.filterwarnings('ignore')

plt.rcParams.update({
    'font.size': 10,
    'font.family': 'DejaVu Sans',
    'axes.linewidth': 0.8,
    'xtick.major.width': 0.8,
    'ytick.major.width': 0.8,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.15,
})

RESULTS_DIR = "/home/yanlin/comp/cursor/analysis/results"
FIG_DIR = "/home/yanlin/comp/cursor/analysis/figures"

SPOP_C = {'AFR': '#D62728', 'AMR': '#FF7F0E', 'EAS': '#2CA02C', 'EUR': '#1F77B4', 'SAS': '#9467BD'}
SPOP_N = {'AFR': 'African', 'AMR': 'American', 'EAS': 'East Asian', 'EUR': 'European', 'SAS': 'South Asian'}
SUPERPOPS = ['AFR', 'AMR', 'EAS', 'EUR', 'SAS']
SIZE_BIN_ORDER = ["1bp", "2bp", "3bp", "4bp", "5bp", "6-10bp", "11-20bp", "21-50bp", ">50bp"]

# Load data
print("Loading data...")
df = pd.read_parquet(f"{RESULTS_DIR}/classified_indels_enhanced.parquet")

# Use full spectrum if available, else chr1
try:
    spec_raw = pd.read_csv(f"{RESULTS_DIR}/perpop_spectrum_full_raw.csv", index_col=0)
    spec_norm = pd.read_csv(f"{RESULTS_DIR}/perpop_spectrum_full_norm.csv", index_col=0)
    sfs = pd.read_csv(f"{RESULTS_DIR}/perpop_sfs_full.csv", index_col=0)
    print("Using full genome per-pop data")
except:
    spec_raw = pd.read_csv(f"{RESULTS_DIR}/perpop_spectrum_chr1_2_raw.csv", index_col=0)
    spec_norm = pd.read_csv(f"{RESULTS_DIR}/perpop_spectrum_chr1_2_norm.csv", index_col=0)
    sfs = pd.read_csv(f"{RESULTS_DIR}/perpop_sfs.csv", index_col=0)
    print("Using chr1 per-pop data")

# =========================================================================
# FIGURE 1: Global indel landscape (4 panels)
# =========================================================================
print("Figure 1: Global landscape...")
fig = plt.figure(figsize=(13, 10))
gs = gridspec.GridSpec(2, 2, hspace=0.35, wspace=0.3)

# 1a: Size distribution
ax = fig.add_subplot(gs[0, 0])
del_sizes = df[df['type']=='DEL']['size_bin'].value_counts().reindex(SIZE_BIN_ORDER).fillna(0)
ins_sizes = df[df['type']=='INS']['size_bin'].value_counts().reindex(SIZE_BIN_ORDER).fillna(0)
x = np.arange(len(SIZE_BIN_ORDER))
w = 0.35
bars1 = ax.bar(x - w/2, del_sizes.values/1e6, w, color='#2166AC', label='Deletions', edgecolor='white', linewidth=0.5)
bars2 = ax.bar(x + w/2, ins_sizes.values/1e6, w, color='#B2182B', label='Insertions', edgecolor='white', linewidth=0.5)
ax.set_xticks(x)
ax.set_xticklabels(SIZE_BIN_ORDER, rotation=45, ha='right', fontsize=8)
ax.set_ylabel('Count (millions)', fontsize=9)
ax.set_title('A', fontweight='bold', fontsize=12, loc='left')
ax.legend(frameon=False, fontsize=8)
ax.spines[['top','right']].set_visible(False)

# 1b: Allele frequency spectrum
ax = fig.add_subplot(gs[0, 1])
af_bins = ['singleton-like', 'very_rare', 'rare', 'low_freq', 'common', 'major']
af_labels = ['<0.1%', '0.1-0.5%', '0.5-1%', '1-5%', '5-50%', '>50%']
for i, (abin, alab) in enumerate(zip(af_bins, af_labels)):
    n_del = ((df['type']=='DEL') & (df['af_bin']==abin)).sum()
    n_ins = ((df['type']=='INS') & (df['af_bin']==abin)).sum()
    ax.bar(i - w/2, n_del/1e6, w, color='#2166AC', edgecolor='white', linewidth=0.5)
    ax.bar(i + w/2, n_ins/1e6, w, color='#B2182B', edgecolor='white', linewidth=0.5)
ax.set_xticks(range(len(af_labels)))
ax.set_xticklabels(af_labels, rotation=45, ha='right', fontsize=8)
ax.set_ylabel('Count (millions)', fontsize=9)
ax.set_title('B', fontweight='bold', fontsize=12, loc='left')
ax.spines[['top','right']].set_visible(False)

# 1c: DEL/INS ratio by size
ax = fig.add_subplot(gs[1, 0])
ratios = del_sizes / ins_sizes
colors = ['#4A1486' if r > 1 else '#1B7837' for r in ratios]
ax.bar(x, ratios.values, color=colors, edgecolor='white', linewidth=0.5, alpha=0.9)
ax.axhline(y=1, color='#666666', linestyle='--', linewidth=0.8, zorder=0)
ax.set_xticks(x)
ax.set_xticklabels(SIZE_BIN_ORDER, rotation=45, ha='right', fontsize=8)
ax.set_ylabel('DEL/INS ratio', fontsize=9)
ax.set_title('C', fontweight='bold', fontsize=12, loc='left')
ax.spines[['top','right']].set_visible(False)
for i, v in enumerate(ratios):
    ax.text(i, v + 0.02, f'{v:.2f}', ha='center', fontsize=7, color='#333')

# 1d: Homopolymer context
ax = fig.add_subplot(gs[1, 1])
mask_1bp = df['size'] == 1
hp_bins_edges = [1, 2, 3, 4, 5, 6, 7, 8, 10, 15, 20, 55]
hp_labels = ['1', '2', '3', '4', '5', '6', '7', '8-9', '10-14', '15-19', '20+']
hp_del = np.histogram(df.loc[mask_1bp & (df['type']=='DEL'), 'homopolymer_len'], bins=hp_bins_edges)[0]
hp_ins = np.histogram(df.loc[mask_1bp & (df['type']=='INS'), 'homopolymer_len'], bins=hp_bins_edges)[0]
x_hp = np.arange(len(hp_labels))
ax.bar(x_hp - w/2, hp_del/1e5, w, color='#2166AC', label='Del', edgecolor='white', linewidth=0.5)
ax.bar(x_hp + w/2, hp_ins/1e5, w, color='#B2182B', label='Ins', edgecolor='white', linewidth=0.5)
ax.set_xticks(x_hp)
ax.set_xticklabels(hp_labels, fontsize=8)
ax.set_xlabel('Homopolymer run length', fontsize=9)
ax.set_ylabel('Count (×10⁵)', fontsize=9)
ax.set_title('D', fontweight='bold', fontsize=12, loc='left')
ax.legend(frameon=False, fontsize=8)
ax.spines[['top','right']].set_visible(False)

fig.savefig(f"{FIG_DIR}/fig1_global_landscape.pdf")
fig.savefig(f"{FIG_DIR}/fig1_global_landscape.png")
plt.close()
print("  Done.")

# =========================================================================
# FIGURE 2: Population comparison (4 panels)
# =========================================================================
print("Figure 2: Population comparison...")
fig = plt.figure(figsize=(14, 11))
gs = gridspec.GridSpec(2, 2, hspace=0.4, wspace=0.35)

# 2a: DEL/INS ratio by population
ax = fig.add_subplot(gs[0, 0])
del_c = spec_raw.filter(like='DEL_').sum(axis=1)
ins_c = spec_raw.filter(like='INS_').sum(axis=1)
pop_ratios = del_c / ins_c
bars = ax.bar(range(5), pop_ratios.values, color=[SPOP_C[sp] for sp in SUPERPOPS],
              edgecolor='white', linewidth=0.5, alpha=0.9)
ax.axhline(y=1, color='#666', linestyle='--', linewidth=0.8, zorder=0)
ax.set_xticks(range(5))
ax.set_xticklabels([SPOP_N[sp] for sp in SUPERPOPS], rotation=30, ha='right', fontsize=9)
ax.set_ylabel('DEL/INS ratio', fontsize=9)
ax.set_title('A', fontweight='bold', fontsize=12, loc='left')
ax.spines[['top','right']].set_visible(False)
for i, v in enumerate(pop_ratios):
    ax.text(i, v + 0.002, f'{v:.3f}', ha='center', fontsize=8, fontweight='bold')

# 2b: Key channel comparison
ax = fig.add_subplot(gs[0, 1])
key_ch = ['DEL_1bp_GC', 'INS_1bp_AT', 'DEL_3bp_mixed', 'DEL_4bp_AT-rich', 'INS_2bp_AT-rich']
key_ch = [c for c in key_ch if c in spec_norm.columns]
ch_labels = [c.replace('_', '\n', 1).replace('_', ' ') for c in key_ch]
x_kc = np.arange(len(key_ch))
width = 0.15
for i, spop in enumerate(SUPERPOPS):
    vals = spec_norm.loc[spop, key_ch].values * 100
    ax.bar(x_kc + i*width - 2*width, vals, width, color=SPOP_C[spop],
           label=SPOP_N[spop], edgecolor='white', linewidth=0.3, alpha=0.9)
ax.set_xticks(x_kc)
ax.set_xticklabels(ch_labels, fontsize=7, ha='center')
ax.set_ylabel('Proportion (%)', fontsize=9)
ax.set_title('B', fontweight='bold', fontsize=12, loc='left')
ax.legend(frameon=False, fontsize=7, ncol=3, loc='upper right')
ax.spines[['top','right']].set_visible(False)

# 2c: SFS by population
ax = fig.add_subplot(gs[1, 0])
sfs_norm = sfs.div(sfs.sum(axis=0), axis=1)
sfs_order = ['rare', 'low_freq', 'common', 'high_freq']
sfs_labels = ['Rare\n(<1%)', 'Low freq\n(1-5%)', 'Common\n(5-50%)', 'High freq\n(>50%)']
x_sfs = np.arange(len(sfs_order))
width = 0.14
for i, spop in enumerate(SUPERPOPS):
    if spop in sfs_norm.columns:
        vals = sfs_norm.loc[sfs_order, spop].values
        ax.bar(x_sfs + i*width - 2*width, vals, width, color=SPOP_C[spop],
               label=SPOP_N[spop], edgecolor='white', linewidth=0.3, alpha=0.9)
ax.set_xticks(x_sfs)
ax.set_xticklabels(sfs_labels, fontsize=8)
ax.set_ylabel('Proportion of segregating indels', fontsize=9)
ax.set_title('C', fontweight='bold', fontsize=12, loc='left')
ax.legend(frameon=False, fontsize=7)
ax.spines[['top','right']].set_visible(False)

# 2d: JSD heatmap
ax = fig.add_subplot(gs[1, 1])
jsd_mat = np.zeros((5, 5))
for i in range(5):
    for j in range(5):
        if i != j:
            jsd_mat[i, j] = jensenshannon(spec_norm.iloc[i], spec_norm.iloc[j])
jsd_df = pd.DataFrame(jsd_mat, index=[SPOP_N[sp] for sp in SUPERPOPS],
                       columns=[SPOP_N[sp] for sp in SUPERPOPS])
mask = np.triu(np.ones_like(jsd_mat, dtype=bool), k=1)
sns.heatmap(jsd_df, annot=True, fmt='.4f', cmap='YlOrRd', ax=ax,
            mask=mask.T, linewidths=0.5, square=True,
            cbar_kws={'label': 'Jensen-Shannon\nDivergence', 'shrink': 0.8})
ax.set_title('D', fontweight='bold', fontsize=12, loc='left')

fig.savefig(f"{FIG_DIR}/fig2_population_comparison.pdf")
fig.savefig(f"{FIG_DIR}/fig2_population_comparison.png")
plt.close()
print("  Done.")

# =========================================================================
# FIGURE 3: NMF signatures (4 panels)
# =========================================================================
print("Figure 3: NMF signatures...")
fig = plt.figure(figsize=(14, 11))
gs = gridspec.GridSpec(2, 2, hspace=0.4, wspace=0.35)

nmf_matrix = spec_raw.copy()
nmf_matrix = nmf_matrix.loc[:, nmf_matrix.sum() > 0]
max_k = min(nmf_matrix.shape)

# 3a: Model selection
ax = fig.add_subplot(gs[0, 0])
errs = []
k_range = range(2, max_k)
for k in k_range:
    m = NMF(n_components=k, init='random', random_state=42, max_iter=1000)
    m.fit_transform(nmf_matrix.values)
    errs.append(m.reconstruction_err_)
ax.plot(list(k_range), errs, 'o-', color='#333', markersize=7, linewidth=1.5)
ax.set_xlabel('Number of signatures (k)', fontsize=9)
ax.set_ylabel('Reconstruction error', fontsize=9)
ax.set_title('A', fontweight='bold', fontsize=12, loc='left')
ax.spines[['top','right']].set_visible(False)

# Fit k=2
model = NMF(n_components=2, init='random', random_state=42, max_iter=1000)
W = model.fit_transform(nmf_matrix.values)
H = model.components_
H_norm = H / H.sum(axis=1, keepdims=True)
W_norm = W / W.sum(axis=1, keepdims=True)
channels = nmf_matrix.columns.tolist()

# 3b: Signature profiles
ax = fig.add_subplot(gs[0, 1])
sig_colors = ['#E41A1C', '#377EB8']
n_ch = len(channels)
x_ch = np.arange(n_ch)

# Group channels by type for coloring
del_mask = np.array([c.startswith('DEL') for c in channels])
ins_mask = ~del_mask

for i in range(2):
    ax.bar(x_ch[del_mask], H_norm[i, del_mask], alpha=0.7 if i == 0 else 0.5,
           color=sig_colors[i], label=f'Signature {i+1}', width=0.8)

ax.set_xlabel('Indel channel', fontsize=9)
ax.set_ylabel('Weight', fontsize=9)
ax.set_title('B', fontweight='bold', fontsize=12, loc='left')

# Only label every 5th channel
tick_positions = list(range(0, n_ch, 4))
ax.set_xticks(tick_positions)
ax.set_xticklabels([channels[i][:12] for i in tick_positions], rotation=90, fontsize=5)
ax.legend(frameon=False, fontsize=8)
ax.spines[['top','right']].set_visible(False)

# 3c: Population signature weights
ax = fig.add_subplot(gs[1, 0])
bottom = np.zeros(5)
for i in range(2):
    ax.bar(range(5), W_norm[:, i], bottom=bottom, color=sig_colors[i],
           label=f'Signature {i+1}', edgecolor='white', linewidth=0.5, alpha=0.9)
    # Add percentage labels
    for j in range(5):
        y_pos = bottom[j] + W_norm[j, i] / 2
        if W_norm[j, i] > 0.1:
            ax.text(j, y_pos, f'{W_norm[j, i]*100:.1f}%', ha='center', va='center',
                    fontsize=8, fontweight='bold', color='white')
    bottom += W_norm[:, i]
ax.set_xticks(range(5))
ax.set_xticklabels([SPOP_N[sp] for sp in SUPERPOPS], rotation=30, ha='right', fontsize=9)
ax.set_ylabel('Signature contribution', fontsize=9)
ax.set_title('C', fontweight='bold', fontsize=12, loc='left')
ax.legend(frameon=False, fontsize=8)
ax.spines[['top','right']].set_visible(False)

# 3d: Top channels per signature (horizontal bars)
ax = fig.add_subplot(gs[1, 1])
n_top = 7
for i in range(2):
    top_idx = np.argsort(H_norm[i])[-n_top:][::-1]
    y_pos = np.arange(n_top) + i * (n_top + 1)
    ax.barh(y_pos, H_norm[i, top_idx], color=sig_colors[i], alpha=0.85,
            edgecolor='white', linewidth=0.3)
    for j, idx in enumerate(top_idx):
        ax.text(H_norm[i, idx] + 0.003, y_pos[j], channels[idx], fontsize=7, va='center')
    ax.text(-0.015, y_pos[n_top//2], f'Sig {i+1}', fontsize=10, fontweight='bold',
            va='center', ha='right', color=sig_colors[i])

ax.set_xlabel('Weight', fontsize=9)
ax.set_title('D', fontweight='bold', fontsize=12, loc='left')
ax.set_yticks([])
ax.spines[['top','right','left']].set_visible(False)
ax.set_xlim(-0.02, H_norm.max() * 1.4)

fig.savefig(f"{FIG_DIR}/fig3_nmf_signatures.pdf")
fig.savefig(f"{FIG_DIR}/fig3_nmf_signatures.png")
plt.close()
print("  Done.")

# =========================================================================
# FIGURE 4: Chromosomal and size-dependent patterns
# =========================================================================
print("Figure 4: Chromosomal patterns...")
fig = plt.figure(figsize=(14, 10))
gs = gridspec.GridSpec(2, 2, hspace=0.4, wspace=0.3)

# 4a: Indel density per chromosome
ax = fig.add_subplot(gs[0, 0])
chrom_counts = df.groupby('chrom').size().reset_index(name='count')
chrom_counts['num'] = chrom_counts['chrom'].str.replace('chr', '').astype(int)
chrom_counts = chrom_counts.sort_values('num')
chrom_sizes_mb = {1:248,2:242,3:198,4:190,5:181,6:170,7:159,8:145,9:138,10:133,
                  11:135,12:133,13:114,14:107,15:101,16:90,17:83,18:80,19:58,20:64,21:46,22:50}
chrom_counts['density'] = chrom_counts.apply(lambda r: r['count']/chrom_sizes_mb[r['num']], axis=1)
colors_chr = ['#1F77B4' if d < chrom_counts['density'].median() else '#D62728' 
              for d in chrom_counts['density']]
ax.bar(range(22), chrom_counts['density'], color=colors_chr, edgecolor='white', linewidth=0.3, alpha=0.85)
ax.axhline(y=chrom_counts['density'].median(), color='#666', linestyle='--', linewidth=0.8)
ax.set_xticks(range(22))
ax.set_xticklabels([str(i) for i in range(1, 23)], fontsize=7)
ax.set_xlabel('Chromosome', fontsize=9)
ax.set_ylabel('Indels per Mb', fontsize=9)
ax.set_title('A', fontweight='bold', fontsize=12, loc='left')
ax.spines[['top','right']].set_visible(False)

# 4b: DEL/INS ratio per chromosome
ax = fig.add_subplot(gs[0, 1])
ct = df.groupby(['chrom','type']).size().unstack(fill_value=0)
ct['num'] = ct.index.str.replace('chr','').astype(int)
ct = ct.sort_values('num')
chr_ratio = ct['DEL'] / ct['INS']
ax.bar(range(22), chr_ratio.values, color='#4A1486', edgecolor='white', linewidth=0.3, alpha=0.85)
ax.axhline(y=chr_ratio.mean(), color='#E41A1C', linestyle='--', linewidth=0.8, label=f'Mean={chr_ratio.mean():.3f}')
ax.set_xticks(range(22))
ax.set_xticklabels([str(i) for i in range(1, 23)], fontsize=7)
ax.set_xlabel('Chromosome', fontsize=9)
ax.set_ylabel('DEL/INS ratio', fontsize=9)
ax.set_title('B', fontweight='bold', fontsize=12, loc='left')
ax.legend(frameon=False, fontsize=8)
ax.spines[['top','right']].set_visible(False)

# 4c: Size-dependent DEL/INS ratio across populations
ax = fig.add_subplot(gs[1, 0])
for spop in SUPERPOPS:
    del_by_size, ins_by_size = {}, {}
    for ch in spec_raw.columns:
        parts = ch.split('_')
        itype, sbin = parts[0], parts[1]
        if sbin not in SIZE_BIN_ORDER:
            sbin = '_'.join(parts[1:2])
        if itype == 'DEL':
            del_by_size[sbin] = del_by_size.get(sbin, 0) + spec_raw.loc[spop, ch]
        else:
            ins_by_size[sbin] = ins_by_size.get(sbin, 0) + spec_raw.loc[spop, ch]
    pop_size_ratios = [del_by_size.get(sb, 0) / ins_by_size.get(sb, 1) for sb in SIZE_BIN_ORDER]
    ax.plot(range(len(SIZE_BIN_ORDER)), pop_size_ratios, 'o-', color=SPOP_C[spop],
            label=SPOP_N[spop], markersize=5, linewidth=1.5, alpha=0.85)
ax.axhline(y=1, color='#666', linestyle='--', linewidth=0.8)
ax.set_xticks(range(len(SIZE_BIN_ORDER)))
ax.set_xticklabels(SIZE_BIN_ORDER, rotation=45, ha='right', fontsize=8)
ax.set_ylabel('DEL/INS ratio', fontsize=9)
ax.set_title('C', fontweight='bold', fontsize=12, loc='left')
ax.legend(frameon=False, fontsize=7, loc='upper left')
ax.spines[['top','right']].set_visible(False)

# 4d: Population deviation from mean spectrum
ax = fig.add_subplot(gs[1, 1])
diff = spec_norm.subtract(spec_norm.mean(axis=0), axis=1) * 100
top_var = diff.abs().max(axis=0).nlargest(8).index
for i, spop in enumerate(SUPERPOPS):
    y_pos = np.arange(len(top_var)) + i * 0.15 - 0.3
    ax.barh(y_pos, diff.loc[spop, top_var].values, 0.14, color=SPOP_C[spop],
            label=SPOP_N[spop], alpha=0.85, edgecolor='white', linewidth=0.2)
ax.axvline(x=0, color='#333', linewidth=0.8)
ax.set_yticks(range(len(top_var)))
ax.set_yticklabels(top_var, fontsize=7)
ax.set_xlabel('Deviation from global mean (%)', fontsize=9)
ax.set_title('D', fontweight='bold', fontsize=12, loc='left')
ax.legend(frameon=False, fontsize=6, loc='lower right')
ax.spines[['top','right']].set_visible(False)

fig.savefig(f"{FIG_DIR}/fig4_chromosomal_patterns.pdf")
fig.savefig(f"{FIG_DIR}/fig4_chromosomal_patterns.png")
plt.close()
print("  Done.")

# =========================================================================
# FIGURE 5: Detailed population differentiation
# =========================================================================
print("Figure 5: Population differentiation detail...")
fig = plt.figure(figsize=(14, 10))
gs = gridspec.GridSpec(2, 2, hspace=0.4, wspace=0.35)

# 5a: Heatmap of normalized spectrum (top variable channels)
ax = fig.add_subplot(gs[0, 0])
top20 = spec_norm.var().nlargest(20).index
hm_data = spec_norm[top20].T
hm_data.columns = [SPOP_N[sp] for sp in SUPERPOPS]
sns.heatmap(hm_data * 100, cmap='YlOrRd', ax=ax, cbar_kws={'label': 'Proportion (%)', 'shrink': 0.8},
            linewidths=0.5, fmt='.2f', annot=True, annot_kws={'fontsize': 6})
ax.set_title('A', fontweight='bold', fontsize=12, loc='left')
ax.set_yticklabels(ax.get_yticklabels(), fontsize=6)

# 5b: Cosine similarity heatmap
ax = fig.add_subplot(gs[0, 1])
cos_sim = np.zeros((5, 5))
for i in range(5):
    for j in range(5):
        cos_sim[i, j] = 1 - cosine(spec_norm.iloc[i], spec_norm.iloc[j])
cos_df = pd.DataFrame(cos_sim, index=[SPOP_N[sp] for sp in SUPERPOPS],
                       columns=[SPOP_N[sp] for sp in SUPERPOPS])
sns.heatmap(cos_df, annot=True, fmt='.6f', cmap='RdYlGn', ax=ax,
            vmin=cos_sim.min() - 0.0005, vmax=1.0, linewidths=0.5, square=True,
            cbar_kws={'label': 'Cosine similarity', 'shrink': 0.8})
ax.set_title('B', fontweight='bold', fontsize=12, loc='left')

# 5c: Hierarchical clustering
ax = fig.add_subplot(gs[1, 0])
dist = pdist(spec_norm.values, metric='cosine')
Z = linkage(dist, method='average')
dn = dendrogram(Z, labels=[SPOP_N[sp] for sp in SUPERPOPS], ax=ax,
                leaf_rotation=30, leaf_font_size=9,
                above_threshold_color='#333')
ax.set_ylabel('Cosine distance', fontsize=9)
ax.set_title('C', fontweight='bold', fontsize=12, loc='left')
ax.spines[['top','right']].set_visible(False)

# 5d: DEL/INS ratio by population and type-size combination
ax = fig.add_subplot(gs[1, 1])
type_size_combos = ['1bp', '2bp', '3bp', '4-5bp', '6-10bp', '11-50bp']
for spop in SUPERPOPS:
    combo_ratios = []
    for combo in type_size_combos:
        if combo == '4-5bp':
            d = sum(spec_raw.loc[spop, c] for c in spec_raw.columns if c.startswith('DEL_4bp') or c.startswith('DEL_5bp'))
            i = sum(spec_raw.loc[spop, c] for c in spec_raw.columns if c.startswith('INS_4bp') or c.startswith('INS_5bp'))
        elif combo == '11-50bp':
            d = sum(spec_raw.loc[spop, c] for c in spec_raw.columns if c.startswith('DEL_11-20bp') or c.startswith('DEL_21-50bp'))
            i = sum(spec_raw.loc[spop, c] for c in spec_raw.columns if c.startswith('INS_11-20bp') or c.startswith('INS_21-50bp'))
        else:
            d = sum(spec_raw.loc[spop, c] for c in spec_raw.columns if c.startswith(f'DEL_{combo}'))
            i = sum(spec_raw.loc[spop, c] for c in spec_raw.columns if c.startswith(f'INS_{combo}'))
        combo_ratios.append(d / i if i > 0 else np.nan)
    ax.plot(range(len(type_size_combos)), combo_ratios, 'o-', color=SPOP_C[spop],
            label=SPOP_N[spop], markersize=5, linewidth=1.5, alpha=0.85)

ax.axhline(y=1, color='#666', linestyle='--', linewidth=0.8)
ax.set_xticks(range(len(type_size_combos)))
ax.set_xticklabels(type_size_combos, fontsize=8)
ax.set_xlabel('Indel size class', fontsize=9)
ax.set_ylabel('DEL/INS ratio', fontsize=9)
ax.set_title('D', fontweight='bold', fontsize=12, loc='left')
ax.legend(frameon=False, fontsize=7, loc='upper left')
ax.spines[['top','right']].set_visible(False)

fig.savefig(f"{FIG_DIR}/fig5_population_differentiation.pdf")
fig.savefig(f"{FIG_DIR}/fig5_population_differentiation.png")
plt.close()
print("  Done.")

print("\nAll final figures generated!")
