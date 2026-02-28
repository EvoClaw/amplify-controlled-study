#!/usr/bin/env python3
"""
Multi-variant-type population genomics analysis.
Compares population structure, differentiation, and SFS across SNVs, INDELs, SVs.
"""
import os, sys, itertools, json, warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.spatial import procrustes
from scipy.spatial.distance import pdist, squareform
from pathlib import Path

warnings.filterwarnings('ignore')
np.random.seed(42)

BASE = Path("/home/yanlin/comp/amplify")
RESULTS = BASE / "results" / "full"
FIGURES = BASE / "figures"
FIGURES.mkdir(parents=True, exist_ok=True)

SUPERPOP_MAP = {
    'YRI':'AFR','LWK':'AFR','GWD':'AFR','MSL':'AFR','ESN':'AFR','ACB':'AFR','ASW':'AFR',
    'CEU':'EUR','TSI':'EUR','FIN':'EUR','GBR':'EUR','IBS':'EUR',
    'CHB':'EAS','JPT':'EAS','CHS':'EAS','CDX':'EAS','KHV':'EAS',
    'GIH':'SAS','PJL':'SAS','BEB':'SAS','STU':'SAS','ITU':'SAS',
    'MXL':'AMR','PUR':'AMR','CLM':'AMR','PEL':'AMR'
}
SUPERPOP_COLORS = {'AFR':'#E41A1C','EUR':'#377EB8','EAS':'#4DAF4A','SAS':'#984EA3','AMR':'#FF7F00'}
SUPERPOP_ORDER = ['AFR','EUR','EAS','SAS','AMR']
VTYPE_COLORS = {'snv':'#377EB8','indel':'#E41A1C','sv':'#4DAF4A'}
VTYPE_LABELS = {'snv':'SNV','indel':'INDEL','sv':'SV'}

plt.rcParams.update({
    'font.size': 10, 'axes.titlesize': 12, 'axes.labelsize': 11,
    'xtick.labelsize': 9, 'ytick.labelsize': 9, 'legend.fontsize': 9,
    'figure.dpi': 300, 'savefig.dpi': 300, 'savefig.bbox': 'tight',
    'font.family': 'sans-serif',
})

def load_sample_info():
    df = pd.read_csv(BASE / "1000GP" / "samples.info", sep='\t', header=None,
                     names=['sample','pop','sex'])
    df['superpop'] = df['pop'].map(SUPERPOP_MAP)
    return df

def load_pca(vtype):
    if vtype == 'sv':
        f = RESULTS / vtype / "all_autosomes_pca.eigenvec"
        fv = RESULTS / vtype / "all_autosomes_pca.eigenval"
    else:
        f = RESULTS / vtype / "all_autosomes_pruned_pca.eigenvec"
        fv = RESULTS / vtype / "all_autosomes_pruned_pca.eigenval"
    if not f.exists():
        return None
    ev = pd.read_csv(f, sep='\t')
    eigenval = pd.read_csv(fv, header=None)[0].values
    return {'eigenvec': ev, 'eigenval': eigenval}

def compute_fst(sample_info):
    """Compute Hudson FST for all population pairs x 3 variant types."""
    pops = sorted(sample_info['pop'].unique())
    pairs = list(itertools.combinations(pops, 2))
    fst_results = {}

    for vtype in ['snv','indel','sv']:
        print(f"  FST {vtype}...", flush=True)
        pop_freqs = {}
        pop_ns = {}
        for pop in pops:
            ff = RESULTS / vtype / f"freq_{pop}.afreq"
            if not ff.exists():
                continue
            df = pd.read_csv(ff, sep='\t')
            pop_freqs[pop] = df['ALT_FREQS'].values
            pop_ns[pop] = df['OBS_CT'].values / 2
        
        fst_vals = {}
        for p1, p2 in pairs:
            if p1 not in pop_freqs or p2 not in pop_freqs:
                continue
            f1, f2 = pop_freqs[p1], pop_freqs[p2]
            n1, n2 = pop_ns[p1], pop_ns[p2]
            num = (f1 - f2)**2 - f1*(1-f1)/(n1-1) - f2*(1-f2)/(n2-1)
            denom = f1*(1-f2) + f2*(1-f1)
            mask = denom > 0
            fst = np.sum(num[mask]) / np.sum(denom[mask])
            fst_vals[(p1, p2)] = max(0, fst)
        
        fst_results[vtype] = fst_vals
        print(f"    {len(fst_vals)} pairs, mean FST={np.mean(list(fst_vals.values())):.4f}", flush=True)
    
    return fst_results

def compute_sfs(sample_info):
    """Compute folded SFS per variant type per super-population."""
    sfs_results = {}
    for vtype in ['snv','indel','sv']:
        sfs_results[vtype] = {}
        for spop in SUPERPOP_ORDER:
            ff = RESULTS / vtype / f"sfs_{spop}.afreq"
            if not ff.exists():
                continue
            df = pd.read_csv(ff, sep='\t')
            afs = df['ALT_FREQS'].values
            maf = np.minimum(afs, 1 - afs)
            bins = np.linspace(0, 0.5, 21)
            hist, _ = np.histogram(maf, bins=bins)
            sfs_results[vtype][spop] = {
                'hist': hist, 'bins': bins, 'maf': maf,
                'n_variants': len(maf), 'mean_maf': float(np.mean(maf)),
                'median_maf': float(np.median(maf)),
                'prop_rare': float(np.mean(maf < 0.05)),
                'skewness': float(stats.skew(maf)),
            }
    return sfs_results

def procrustes_analysis(pca_data, sample_info):
    """Procrustes analysis comparing PCA spaces across variant types."""
    types = [t for t in ['snv','indel','sv'] if pca_data.get(t)]
    if len(types) < 2:
        return {}
    
    sample_sets = []
    for vtype in types:
        ev = pca_data[vtype]['eigenvec']
        sample_sets.append(set(ev['#IID'].astype(str)))
    common = set.intersection(*sample_sets)
    print(f"  Common samples: {len(common)}", flush=True)
    
    proc_results = {}
    for t1, t2 in itertools.combinations(types, 2):
        ev1 = pca_data[t1]['eigenvec'].copy()
        ev2 = pca_data[t2]['eigenvec'].copy()
        ev1['#IID'] = ev1['#IID'].astype(str)
        ev2['#IID'] = ev2['#IID'].astype(str)
        ev1 = ev1[ev1['#IID'].isin(common)].sort_values('#IID')
        ev2 = ev2[ev2['#IID'].isin(common)].sort_values('#IID')
        
        pc_cols = [f'PC{i}' for i in range(1,11)]
        X = ev1[pc_cols].values
        Y = ev2[pc_cols].values
        X = (X - X.mean(0)) / (X.std(0) + 1e-10)
        Y = (Y - Y.mean(0)) / (Y.std(0) + 1e-10)
        
        _, _, disparity = procrustes(X, Y)
        correlation = 1 - disparity
        
        n_perm = 1000
        perm_disps = np.zeros(n_perm)
        for i in range(n_perm):
            Yp = Y[np.random.permutation(len(Y))]
            _, _, d = procrustes(X, Yp)
            perm_disps[i] = d
        p_value = float(np.mean(perm_disps <= disparity))
        
        d1 = pdist(X)
        d2 = pdist(Y)
        mantel_r, mantel_p = stats.pearsonr(d1, d2)
        
        proc_results[f"{t1}_vs_{t2}"] = {
            'disparity': float(disparity), 'correlation': float(correlation),
            'p_value': p_value, 'mantel_r': float(mantel_r),
            'mantel_p': float(mantel_p), 'n_samples': len(common),
        }
        print(f"    {t1} vs {t2}: Procrustes corr={correlation:.4f} (p={p_value}), "
              f"Mantel r={mantel_r:.4f}", flush=True)
    
    return proc_results

def compute_nj_tree_distances(fst_results, sample_info):
    """Compute Robinson-Foulds-like distances between population trees from different variant types."""
    from scipy.cluster.hierarchy import linkage, to_tree
    
    types = [t for t in ['snv','indel','sv'] if t in fst_results and len(fst_results[t]) > 0]
    pops = sorted(sample_info['pop'].unique())
    pop_to_spop = dict(zip(sample_info['pop'], sample_info['superpop']))
    
    tree_results = {}
    linkage_results = {}
    
    for vtype in types:
        mat = pd.DataFrame(0.0, index=pops, columns=pops)
        for (p1, p2), val in fst_results[vtype].items():
            mat.loc[p1, p2] = val
            mat.loc[p2, p1] = val
        
        condensed = squareform(mat.values)
        condensed = np.maximum(condensed, 0)
        Z = linkage(condensed, method='average')
        linkage_results[vtype] = {'linkage': Z, 'labels': pops, 'fst_matrix': mat}
    
    for t1, t2 in itertools.combinations(types, 2):
        if t1 in linkage_results and t2 in linkage_results:
            mat1 = linkage_results[t1]['fst_matrix']
            mat2 = linkage_results[t2]['fst_matrix']
            d1 = squareform(mat1.values)
            d2 = squareform(mat2.values)
            r, p = stats.pearsonr(d1, d2)
            tree_results[f"{t1}_vs_{t2}"] = {'matrix_corr': float(r), 'matrix_p': float(p)}
    
    return linkage_results, tree_results

# ===== PLOTTING =====

def plot_fig1_pca(pca_data, sample_info):
    """Figure 1: PCA comparison across variant types."""
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    
    variant_counts = {'snv': '2.29M', 'indel': '1.19M', 'sv': '18.5K'}
    
    for i, vtype in enumerate(['snv','indel','sv']):
        if not pca_data.get(vtype):
            continue
        ev = pca_data[vtype]['eigenvec'].copy()
        eigenval = pca_data[vtype]['eigenval']
        ev = ev.merge(sample_info[['sample','superpop']], left_on='#IID', right_on='sample', how='left')
        var_exp = eigenval / eigenval.sum() * 100
        
        for spop in SUPERPOP_ORDER:
            mask = ev['superpop'] == spop
            axes[i].scatter(ev.loc[mask, 'PC1'], ev.loc[mask, 'PC2'],
                          c=SUPERPOP_COLORS[spop], label=spop, s=2, alpha=0.4, rasterized=True)
        axes[i].set_xlabel(f"PC1 ({var_exp[0]:.1f}%)")
        axes[i].set_ylabel(f"PC2 ({var_exp[1]:.1f}%)")
        axes[i].set_title(f"{VTYPE_LABELS[vtype]} (n={variant_counts[vtype]})")
        axes[i].legend(markerscale=4, framealpha=0.9, loc='best')
    
    plt.tight_layout()
    plt.savefig(FIGURES / "fig1_pca_comparison.pdf", bbox_inches='tight')
    plt.savefig(FIGURES / "fig1_pca_comparison.png", bbox_inches='tight')
    plt.close()
    print("  Saved fig1_pca_comparison", flush=True)

def plot_fig2_procrustes_overlay(pca_data, sample_info, proc_results):
    """Figure 2: Procrustes-transformed PCA overlays."""
    types = [t for t in ['snv','indel','sv'] if pca_data.get(t)]
    pairs = list(itertools.combinations(types, 2))
    
    fig, axes = plt.subplots(1, len(pairs), figsize=(5*len(pairs), 5))
    if len(pairs) == 1:
        axes = [axes]
    
    sample_sets = [set(pca_data[t]['eigenvec']['#IID'].astype(str)) for t in types]
    common = set.intersection(*sample_sets)
    
    for idx, (t1, t2) in enumerate(pairs):
        ev1 = pca_data[t1]['eigenvec'].copy()
        ev2 = pca_data[t2]['eigenvec'].copy()
        ev1['#IID'] = ev1['#IID'].astype(str)
        ev2['#IID'] = ev2['#IID'].astype(str)
        ev1 = ev1[ev1['#IID'].isin(common)].sort_values('#IID')
        ev2 = ev2[ev2['#IID'].isin(common)].sort_values('#IID')
        
        X = ev1[['PC1','PC2']].values
        Y = ev2[['PC1','PC2']].values
        X = (X - X.mean(0)) / X.std(0)
        Y = (Y - Y.mean(0)) / Y.std(0)
        mtx1, mtx2, _ = procrustes(X, Y)
        
        samples_info_sorted = ev1.merge(sample_info[['sample','superpop']], left_on='#IID', right_on='sample')
        
        for spop in SUPERPOP_ORDER:
            mask = (samples_info_sorted['superpop'] == spop).values
            axes[idx].scatter(mtx1[mask, 0], mtx1[mask, 1], c=VTYPE_COLORS[t1],
                            marker='o', s=4, alpha=0.3, rasterized=True)
            axes[idx].scatter(mtx2[mask, 0], mtx2[mask, 1], c=VTYPE_COLORS[t2],
                            marker='x', s=4, alpha=0.3, rasterized=True)
        
        key = f"{t1}_vs_{t2}"
        corr = proc_results.get(key, {}).get('correlation', 0)
        axes[idx].set_title(f"{VTYPE_LABELS[t1]} vs {VTYPE_LABELS[t2]}\n(Procrustes corr = {corr:.4f})")
        axes[idx].legend([VTYPE_LABELS[t1], VTYPE_LABELS[t2]], loc='upper right')
        axes[idx].set_xlabel("Procrustes dim 1")
        axes[idx].set_ylabel("Procrustes dim 2")
    
    plt.tight_layout()
    plt.savefig(FIGURES / "fig2_procrustes_overlay.pdf", bbox_inches='tight')
    plt.savefig(FIGURES / "fig2_procrustes_overlay.png", bbox_inches='tight')
    plt.close()
    print("  Saved fig2_procrustes_overlay", flush=True)

def plot_fig3_fst_heatmaps(fst_results, sample_info):
    """Figure 3: FST heatmaps for 3 variant types."""
    pops = sorted(sample_info['pop'].unique())
    pop_to_spop = dict(zip(sample_info['pop'], sample_info['superpop']))
    pops_sorted = sorted(pops, key=lambda p: (SUPERPOP_ORDER.index(pop_to_spop.get(p,'AMR')), p))
    
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))
    
    for i, vtype in enumerate(['snv','indel','sv']):
        if vtype not in fst_results:
            continue
        mat = pd.DataFrame(0.0, index=pops_sorted, columns=pops_sorted)
        for (p1, p2), val in fst_results[vtype].items():
            if p1 in pops_sorted and p2 in pops_sorted:
                mat.loc[p1, p2] = val
                mat.loc[p2, p1] = val
        
        mask_upper = np.triu(np.ones_like(mat, dtype=bool), k=0)
        sns.heatmap(mat, mask=mask_upper, ax=axes[i], cmap='YlOrRd', vmin=0, vmax=0.25,
                   xticklabels=True, yticklabels=True, cbar_kws={'shrink': 0.8, 'label': 'FST'},
                   square=True)
        axes[i].set_title(f"FST — {VTYPE_LABELS[vtype]}")
        axes[i].tick_params(labelsize=6, rotation=90)
        
        spop_borders = []
        prev_spop = None
        for j, p in enumerate(pops_sorted):
            sp = pop_to_spop.get(p)
            if sp != prev_spop and prev_spop is not None:
                spop_borders.append(j)
            prev_spop = sp
        for b in spop_borders:
            axes[i].axhline(y=b, color='black', linewidth=0.5)
            axes[i].axvline(x=b, color='black', linewidth=0.5)
    
    plt.tight_layout()
    plt.savefig(FIGURES / "fig3_fst_heatmaps.pdf", bbox_inches='tight')
    plt.savefig(FIGURES / "fig3_fst_heatmaps.png", bbox_inches='tight')
    plt.close()
    print("  Saved fig3_fst_heatmaps", flush=True)

def plot_fig4_fst_crosstype(fst_results):
    """Figure 4: Cross-type FST scatter plots."""
    types = [t for t in ['snv','indel','sv'] if t in fst_results and len(fst_results[t]) > 0]
    pairs = list(itertools.combinations(types, 2))
    
    fig, axes = plt.subplots(1, len(pairs), figsize=(5*len(pairs), 5))
    if len(pairs) == 1:
        axes = [axes]
    
    fst_corr_results = {}
    for i, (t1, t2) in enumerate(pairs):
        common_keys = set(fst_results[t1].keys()) & set(fst_results[t2].keys())
        if len(common_keys) < 3:
            continue
        x = np.array([fst_results[t1][k] for k in common_keys])
        y = np.array([fst_results[t2][k] for k in common_keys])
        r, p_r = stats.pearsonr(x, y)
        rho, p_rho = stats.spearmanr(x, y)
        
        axes[i].scatter(x, y, s=8, alpha=0.5, c='steelblue', edgecolors='none', rasterized=True)
        max_val = max(x.max(), y.max()) * 1.05
        axes[i].plot([0, max_val], [0, max_val], 'k--', alpha=0.3, linewidth=0.8)
        
        slope, intercept = np.polyfit(x, y, 1)
        x_fit = np.linspace(0, max_val, 100)
        axes[i].plot(x_fit, slope*x_fit + intercept, 'r-', alpha=0.5, linewidth=1)
        
        axes[i].set_xlabel(f"FST ({VTYPE_LABELS[t1]})")
        axes[i].set_ylabel(f"FST ({VTYPE_LABELS[t2]})")
        axes[i].set_title(f"r = {r:.3f}, ρ = {rho:.3f}")
        axes[i].text(0.05, 0.95, f"n = {len(common_keys)} pairs\nPearson p = {p_r:.2e}\nSpearman p = {p_rho:.2e}",
                    transform=axes[i].transAxes, va='top', fontsize=8,
                    bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        axes[i].set_xlim(0, max_val)
        axes[i].set_ylim(0, max_val)
        
        fst_corr_results[f"{t1}_vs_{t2}"] = {
            'pearson_r': float(r), 'pearson_p': float(p_r),
            'spearman_rho': float(rho), 'spearman_p': float(p_rho),
            'n_pairs': len(common_keys), 'slope': float(slope),
        }
    
    plt.tight_layout()
    plt.savefig(FIGURES / "fig4_fst_crosstype.pdf", bbox_inches='tight')
    plt.savefig(FIGURES / "fig4_fst_crosstype.png", bbox_inches='tight')
    plt.close()
    print("  Saved fig4_fst_crosstype", flush=True)
    return fst_corr_results

def plot_fig5_sfs(sfs_results):
    """Figure 5: SFS comparison across variant types per super-population."""
    fig, axes = plt.subplots(1, 5, figsize=(20, 4))
    
    for j, spop in enumerate(SUPERPOP_ORDER):
        for vtype in ['snv','indel','sv']:
            if spop in sfs_results.get(vtype, {}):
                data = sfs_results[vtype][spop]
                bins = data['bins']
                hist = data['hist'] / data['hist'].sum()
                centers = (bins[:-1] + bins[1:]) / 2
                axes[j].plot(centers, hist, '-o', markersize=3, color=VTYPE_COLORS[vtype],
                           label=VTYPE_LABELS[vtype], alpha=0.8, linewidth=1.5)
        axes[j].set_xlabel("Minor Allele Frequency")
        if j == 0:
            axes[j].set_ylabel("Proportion of variants")
        axes[j].set_title(spop, fontweight='bold')
        axes[j].legend(fontsize=8)
        axes[j].set_yscale('log')
        axes[j].set_ylim(1e-4, 1)
    
    plt.tight_layout()
    plt.savefig(FIGURES / "fig5_sfs_comparison.pdf", bbox_inches='tight')
    plt.savefig(FIGURES / "fig5_sfs_comparison.png", bbox_inches='tight')
    plt.close()
    print("  Saved fig5_sfs_comparison", flush=True)

def plot_fig6_rare_variant_heatmap(sfs_results):
    """Figure 6: Rare variant proportion heatmap."""
    data = []
    for vtype in ['snv','indel','sv']:
        for spop in SUPERPOP_ORDER:
            if spop in sfs_results.get(vtype, {}):
                data.append({
                    'Variant Type': VTYPE_LABELS[vtype],
                    'Super-population': spop,
                    'Rare Proportion': sfs_results[vtype][spop]['prop_rare'],
                    'Mean MAF': sfs_results[vtype][spop]['mean_maf'],
                })
    
    df = pd.DataFrame(data)
    
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    
    rare_pivot = df.pivot(index='Variant Type', columns='Super-population', values='Rare Proportion')
    rare_pivot = rare_pivot[SUPERPOP_ORDER]
    sns.heatmap(rare_pivot, ax=axes[0], cmap='YlOrRd', annot=True, fmt='.3f',
               cbar_kws={'label': 'Proportion MAF < 0.05'})
    axes[0].set_title("Rare Variant Proportion (MAF < 0.05)")
    
    maf_pivot = df.pivot(index='Variant Type', columns='Super-population', values='Mean MAF')
    maf_pivot = maf_pivot[SUPERPOP_ORDER]
    sns.heatmap(maf_pivot, ax=axes[1], cmap='YlGnBu', annot=True, fmt='.4f',
               cbar_kws={'label': 'Mean MAF'})
    axes[1].set_title("Mean Minor Allele Frequency")
    
    plt.tight_layout()
    plt.savefig(FIGURES / "fig6_sfs_heatmap.pdf", bbox_inches='tight')
    plt.savefig(FIGURES / "fig6_sfs_heatmap.png", bbox_inches='tight')
    plt.close()
    print("  Saved fig6_sfs_heatmap", flush=True)

def plot_fig7_population_trees(linkage_results, sample_info):
    """Figure 7: Population trees (dendrograms) from FST distances."""
    from scipy.cluster.hierarchy import dendrogram
    
    pop_to_spop = dict(zip(sample_info['pop'], sample_info['superpop']))
    types = [t for t in ['snv','indel','sv'] if t in linkage_results]
    
    fig, axes = plt.subplots(1, len(types), figsize=(6*len(types), 6))
    if len(types) == 1:
        axes = [axes]
    
    for i, vtype in enumerate(types):
        Z = linkage_results[vtype]['linkage']
        labels = linkage_results[vtype]['labels']
        
        leaf_colors = {}
        for j, lab in enumerate(labels):
            sp = pop_to_spop.get(lab, 'AMR')
            leaf_colors[lab] = SUPERPOP_COLORS[sp]
        
        with plt.rc_context({'lines.linewidth': 1.5}):
            dend = dendrogram(Z, labels=labels, ax=axes[i], leaf_rotation=90,
                            leaf_font_size=7, color_threshold=0)
        
        xlbls = axes[i].get_xmajorticklabels()
        for lbl in xlbls:
            sp = pop_to_spop.get(lbl.get_text(), 'AMR')
            lbl.set_color(SUPERPOP_COLORS[sp])
        
        axes[i].set_title(f"{VTYPE_LABELS[vtype]}")
        axes[i].set_ylabel("FST distance")
    
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=SUPERPOP_COLORS[sp], label=sp) for sp in SUPERPOP_ORDER]
    axes[-1].legend(handles=legend_elements, loc='upper right', fontsize=8)
    
    plt.tight_layout()
    plt.savefig(FIGURES / "fig7_population_trees.pdf", bbox_inches='tight')
    plt.savefig(FIGURES / "fig7_population_trees.png", bbox_inches='tight')
    plt.close()
    print("  Saved fig7_population_trees", flush=True)

def plot_fig8_concordance_summary(proc_results, fst_corr_results):
    """Figure 8: Summary bar chart of concordance metrics."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    
    if proc_results:
        labels = [k.replace('_', ' ').replace('vs', '↔').upper() for k in proc_results.keys()]
        raw_labels = list(proc_results.keys())
        proc_vals = [proc_results[k]['correlation'] for k in raw_labels]
        mantel_vals = [proc_results[k]['mantel_r'] for k in raw_labels]
        
        x = np.arange(len(labels))
        w = 0.35
        axes[0].bar(x - w/2, proc_vals, w, label='Procrustes correlation', color='steelblue')
        axes[0].bar(x + w/2, mantel_vals, w, label='Mantel r', color='coral')
        axes[0].set_xticks(x)
        axes[0].set_xticklabels(labels, fontsize=9)
        axes[0].set_ylabel("Correlation")
        axes[0].set_title("PCA Space Concordance")
        axes[0].axhline(0.90, color='red', linestyle='--', alpha=0.5, linewidth=0.8)
        axes[0].set_ylim(0, 1.05)
        axes[0].legend(fontsize=8)
        axes[0].text(0.5, 0.88, 'threshold', color='red', alpha=0.5, fontsize=8,
                    transform=axes[0].get_xaxis_transform(), ha='center')
    
    if fst_corr_results:
        labels = [k.replace('_', ' ').replace('vs', '↔').upper() for k in fst_corr_results.keys()]
        raw_labels = list(fst_corr_results.keys())
        pearson_vals = [fst_corr_results[k]['pearson_r'] for k in raw_labels]
        spearman_vals = [fst_corr_results[k]['spearman_rho'] for k in raw_labels]
        
        x = np.arange(len(labels))
        w = 0.35
        axes[1].bar(x - w/2, pearson_vals, w, label='Pearson r', color='steelblue')
        axes[1].bar(x + w/2, spearman_vals, w, label='Spearman ρ', color='coral')
        axes[1].set_xticks(x)
        axes[1].set_xticklabels(labels, fontsize=9)
        axes[1].set_ylabel("Correlation")
        axes[1].set_title("FST Concordance")
        axes[1].axhline(0.85, color='red', linestyle='--', alpha=0.5, linewidth=0.8)
        axes[1].set_ylim(0, 1.05)
        axes[1].legend(fontsize=8)
        axes[1].text(0.5, 0.83, 'threshold', color='red', alpha=0.5, fontsize=8,
                    transform=axes[1].get_xaxis_transform(), ha='center')
    
    plt.tight_layout()
    plt.savefig(FIGURES / "fig8_concordance_summary.pdf", bbox_inches='tight')
    plt.savefig(FIGURES / "fig8_concordance_summary.png", bbox_inches='tight')
    plt.close()
    print("  Saved fig8_concordance_summary", flush=True)

def compute_sfs_statistics(sfs_results):
    """Statistical tests for SFS shape differences."""
    ks_results = []
    for spop in SUPERPOP_ORDER:
        for t1, t2 in itertools.combinations(['snv','indel','sv'], 2):
            if spop not in sfs_results.get(t1, {}) or spop not in sfs_results.get(t2, {}):
                continue
            maf1 = sfs_results[t1][spop]['maf']
            maf2 = sfs_results[t2][spop]['maf']
            ks_stat, ks_p = stats.ks_2samp(maf1, maf2)
            ks_results.append({
                'super_pop': spop,
                'type1': VTYPE_LABELS[t1], 'type2': VTYPE_LABELS[t2],
                'ks_statistic': float(ks_stat), 'ks_p_value': float(ks_p),
                'n1': len(maf1), 'n2': len(maf2),
            })
    
    df = pd.DataFrame(ks_results)
    if len(df) > 0:
        from statsmodels.stats.multitest import multipletests
        _, fdr_p, _, _ = multipletests(df['ks_p_value'], method='fdr_bh')
        df['fdr_p_value'] = fdr_p
    
    return df

def plot_fig9_sfs_ks_heatmap(ks_df):
    """Figure 9: KS test heatmap for SFS shape differences."""
    if len(ks_df) == 0:
        return
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 4))
    
    ks_pivot = ks_df.pivot_table(index=['type1','type2'], columns='super_pop', values='ks_statistic')
    ks_pivot.index = [f"{a} vs {b}" for a, b in ks_pivot.index]
    ks_pivot = ks_pivot[SUPERPOP_ORDER]
    sns.heatmap(ks_pivot, ax=axes[0], cmap='YlOrRd', annot=True, fmt='.3f',
               cbar_kws={'label': 'KS statistic'})
    axes[0].set_title("KS Statistic (SFS Shape Difference)")
    
    fdr_pivot = ks_df.pivot_table(index=['type1','type2'], columns='super_pop', values='fdr_p_value')
    fdr_pivot.index = [f"{a} vs {b}" for a, b in fdr_pivot.index]
    fdr_pivot = fdr_pivot[SUPERPOP_ORDER]
    log_p = -np.log10(fdr_pivot + 1e-300)
    sns.heatmap(log_p, ax=axes[1], cmap='YlOrRd', annot=True, fmt='.1f',
               cbar_kws={'label': '-log10(FDR-adjusted p)'})
    axes[1].set_title("Statistical Significance (-log10 FDR p)")
    
    plt.tight_layout()
    plt.savefig(FIGURES / "fig9_sfs_ks_tests.pdf", bbox_inches='tight')
    plt.savefig(FIGURES / "fig9_sfs_ks_tests.png", bbox_inches='tight')
    plt.close()
    print("  Saved fig9_sfs_ks_tests", flush=True)

# ===== MAIN =====
def main():
    print("=" * 60, flush=True)
    print("Multi-Variant-Type Population Genomics Analysis", flush=True)
    print("=" * 60, flush=True)
    
    sample_info = load_sample_info()
    print(f"\nSamples: {len(sample_info)}", flush=True)
    print(f"Populations: {sample_info['pop'].nunique()}", flush=True)
    print(f"Super-populations: {dict(sample_info['superpop'].value_counts())}", flush=True)
    
    print("\n=== Loading PCA results ===", flush=True)
    pca_data = {}
    for vtype in ['snv','indel','sv']:
        pca_data[vtype] = load_pca(vtype)
        if pca_data[vtype]:
            ev = pca_data[vtype]['eigenval']
            print(f"  {vtype}: {len(pca_data[vtype]['eigenvec'])} samples, "
                  f"top eigenval={ev[0]:.2f}, var_explained_PC1={ev[0]/ev.sum()*100:.1f}%", flush=True)
    
    print("\n=== Computing FST ===", flush=True)
    fst_results = compute_fst(sample_info)
    
    print("\n=== Computing SFS ===", flush=True)
    sfs_results = compute_sfs(sample_info)
    
    print("\n=== Procrustes Analysis ===", flush=True)
    proc_results = procrustes_analysis(pca_data, sample_info)
    
    print("\n=== Population Trees ===", flush=True)
    linkage_results, tree_dist_results = compute_nj_tree_distances(fst_results, sample_info)
    for k, v in tree_dist_results.items():
        print(f"  {k}: matrix_corr={v['matrix_corr']:.4f}", flush=True)
    
    print("\n=== SFS Statistical Tests ===", flush=True)
    ks_df = compute_sfs_statistics(sfs_results)
    if len(ks_df) > 0:
        print(ks_df[['super_pop','type1','type2','ks_statistic','fdr_p_value']].to_string(index=False), flush=True)
    
    print("\n=== Generating Figures ===", flush=True)
    plot_fig1_pca(pca_data, sample_info)
    plot_fig2_procrustes_overlay(pca_data, sample_info, proc_results)
    plot_fig3_fst_heatmaps(fst_results, sample_info)
    fst_corr_results = plot_fig4_fst_crosstype(fst_results)
    plot_fig5_sfs(sfs_results)
    plot_fig6_rare_variant_heatmap(sfs_results)
    plot_fig7_population_trees(linkage_results, sample_info)
    plot_fig8_concordance_summary(proc_results, fst_corr_results)
    plot_fig9_sfs_ks_heatmap(ks_df)
    
    # Compile results summary
    summary = {
        'dataset': {
            'samples': len(sample_info),
            'populations': int(sample_info['pop'].nunique()),
            'super_populations': 5,
        },
        'variant_counts': {
            'snv_after_ld_pruning': '2,291,746',
            'indel_after_ld_pruning': '1,194,489',
            'sv_total': '18,519',
        },
        'procrustes': proc_results,
        'fst_correlations': fst_corr_results,
        'tree_distances': tree_dist_results,
        'sfs_statistics': {},
    }
    for vtype in sfs_results:
        for spop in sfs_results[vtype]:
            key = f"{vtype}_{spop}"
            summary['sfs_statistics'][key] = {
                k: v for k, v in sfs_results[vtype][spop].items()
                if k not in ('hist', 'bins', 'maf')
            }
    
    with open(RESULTS / "analysis_summary.json", 'w') as f:
        json.dump(summary, f, indent=2, default=str)
    
    if len(ks_df) > 0:
        ks_df.to_csv(RESULTS / "sfs_ks_tests.csv", index=False)
    
    fst_df_rows = []
    for vtype, pairs in fst_results.items():
        for (p1, p2), val in pairs.items():
            fst_df_rows.append({'variant_type': vtype, 'pop1': p1, 'pop2': p2, 'fst': val})
    pd.DataFrame(fst_df_rows).to_csv(RESULTS / "fst_all_pairs.csv", index=False)
    
    print("\n" + "=" * 60, flush=True)
    print("ANALYSIS COMPLETE", flush=True)
    print(f"Results: {RESULTS}", flush=True)
    print(f"Figures: {FIGURES}", flush=True)
    print("=" * 60, flush=True)
    
    return summary

if __name__ == '__main__':
    main()
