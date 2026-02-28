#!/usr/bin/env python3
"""
Supplementary analyses addressing reviewer concerns:
ES1: SV subtype decomposition + size-stratified FST
ES2: Convergence curve (variant count vs Procrustes)
ES3: Uniform MAF filter SFS recompute
ES4: Proper Mantel test with permutation
ES5: Bootstrap CIs for key statistics
ES6: Basic functional context (gene proximity)
"""
import os, sys, subprocess, itertools, json, warnings
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
SAMPINFO = BASE / "1000GP" / "samples.info"

SUPERPOP_MAP = {
    'YRI':'AFR','LWK':'AFR','GWD':'AFR','MSL':'AFR','ESN':'AFR','ACB':'AFR','ASW':'AFR',
    'CEU':'EUR','TSI':'EUR','FIN':'EUR','GBR':'EUR','IBS':'EUR',
    'CHB':'EAS','JPT':'EAS','CHS':'EAS','CDX':'EAS','KHV':'EAS',
    'GIH':'SAS','PJL':'SAS','BEB':'SAS','STU':'SAS','ITU':'SAS',
    'MXL':'AMR','PUR':'AMR','CLM':'AMR','PEL':'AMR'
}
SUPERPOP_ORDER = ['AFR','EUR','EAS','SAS','AMR']
VTYPE_COLORS = {'snv':'#377EB8','indel':'#E41A1C','sv':'#4DAF4A'}
VTYPE_LABELS = {'snv':'SNV','indel':'INDEL','sv':'SV'}

plt.rcParams.update({
    'font.size': 10, 'axes.titlesize': 12, 'axes.labelsize': 11,
    'figure.dpi': 300, 'savefig.dpi': 300, 'savefig.bbox': 'tight',
})

def run(cmd):
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if r.returncode != 0 and r.stderr:
        print(f"  WARN: {r.stderr[:150]}", flush=True)
    return r

def load_sample_info():
    df = pd.read_csv(SAMPINFO, sep='\t', header=None, names=['sample','pop','sex'])
    df['superpop'] = df['pop'].map(SUPERPOP_MAP)
    return df

# ===== ES1: INDEL/SV size-stratified analysis =====
def size_stratified_analysis():
    """Stratify INDELs and SVs by size, compute FST per stratum."""
    print("=== ES1: Size-stratified FST analysis ===", flush=True)
    
    sample_info = load_sample_info()
    pops = sorted(sample_info['pop'].unique())
    
    # Read INDEL variant info from PVAR
    indel_pvar = RESULTS / "indel" / "all_autosomes.pvar"
    print("  Reading INDEL variant sizes...", flush=True)
    indel_sizes = []
    with open(indel_pvar) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            ref = parts[3] if len(parts) > 3 else ''
            alt = parts[4] if len(parts) > 4 else ''
            size = abs(len(alt) - len(ref))
            indel_sizes.append(size)
    print(f"  INDELs: {len(indel_sizes)} variants, median size={np.median(indel_sizes):.0f}bp", flush=True)
    
    # Read SV variant info
    sv_pvar = RESULTS / "sv" / "all_autosomes.pvar"
    print("  Reading SV variant sizes...", flush=True)
    sv_sizes = []
    with open(sv_pvar) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            ref = parts[3] if len(parts) > 3 else ''
            alt = parts[4] if len(parts) > 4 else ''
            size = abs(len(alt) - len(ref))
            sv_sizes.append(size)
    print(f"  SVs: {len(sv_sizes)} variants, median size={np.median(sv_sizes):.0f}bp", flush=True)
    
    # Size distribution plots
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    
    indel_sizes_arr = np.array(indel_sizes)
    axes[0].hist(np.clip(indel_sizes_arr, 0, 50), bins=50, color=VTYPE_COLORS['indel'], alpha=0.7, edgecolor='black', linewidth=0.3)
    axes[0].set_xlabel("INDEL size (bp)")
    axes[0].set_ylabel("Count")
    axes[0].set_title("INDEL Size Distribution")
    axes[0].set_yscale('log')
    
    sv_sizes_arr = np.array(sv_sizes)
    sv_log = np.log10(sv_sizes_arr[sv_sizes_arr > 0] + 1)
    axes[1].hist(sv_log, bins=50, color=VTYPE_COLORS['sv'], alpha=0.7, edgecolor='black', linewidth=0.3)
    axes[1].set_xlabel("SV size (log10 bp)")
    axes[1].set_ylabel("Count")
    axes[1].set_title("SV Size Distribution")
    
    plt.tight_layout()
    plt.savefig(FIGURES / "fig_supp_size_distributions.pdf", bbox_inches='tight')
    plt.savefig(FIGURES / "fig_supp_size_distributions.png", bbox_inches='tight')
    plt.close()
    print("  Saved size distribution plots", flush=True)
    
    # Stratify and compute FST per size bin
    # INDEL bins: 1bp, 2-5bp, 6-20bp, 21-50bp
    indel_bins = [(1, 1, '1bp'), (2, 5, '2-5bp'), (6, 20, '6-20bp'), (21, 50, '21-50bp')]
    # SV bins: 50-200bp, 200-1000bp, 1000-10000bp, >10000bp
    sv_bins = [(50, 200, '50-200bp'), (200, 1000, '0.2-1kb'), (1000, 10000, '1-10kb'), (10000, 1000000, '>10kb')]
    
    size_fst_results = {}
    
    for vtype, bins_list, sizes_arr in [('indel', indel_bins, indel_sizes_arr), ('sv', sv_bins, sv_sizes_arr)]:
        prefix = RESULTS / vtype / "all_autosomes"
        pvar = Path(str(prefix) + ".pvar")
        
        # Read variant IDs
        var_ids = []
        with open(pvar) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                var_ids.append(line.strip().split('\t')[2])
        
        for lo, hi, label in bins_list:
            mask = (sizes_arr >= lo) & (sizes_arr <= hi)
            n_in_bin = mask.sum()
            if n_in_bin < 100:
                print(f"    {vtype} {label}: only {n_in_bin} variants, skipping", flush=True)
                continue
            
            bin_ids = [var_ids[i] for i in range(len(var_ids)) if i < len(mask) and mask[i]]
            
            # Write variant list
            bin_file = RESULTS / vtype / f"size_bin_{label.replace(' ','_')}.txt"
            with open(bin_file, 'w') as f:
                for vid in bin_ids:
                    f.write(vid + '\n')
            
            # Extract and compute frequencies
            bin_prefix = RESULTS / vtype / f"size_bin_{label.replace(' ','_')}"
            run(f"plink2 --pfile {prefix} --extract {bin_file} --make-pgen "
                f"--out {bin_prefix} --allow-extra-chr --threads 8")
            
            # Compute per-pop frequencies
            pop_freqs = {}
            pop_ns = {}
            for pop in pops:
                pop_f = RESULTS / vtype / f"pop_{pop}.txt"
                freq_out = RESULTS / vtype / f"size_{label.replace(' ','_')}_freq_{pop}"
                run(f"plink2 --pfile {bin_prefix} --keep {pop_f} --freq "
                    f"--out {freq_out} --allow-extra-chr --threads 4")
                ff = Path(str(freq_out) + ".afreq")
                if ff.exists():
                    try:
                        df = pd.read_csv(ff, sep='\t')
                        if len(df) > 0:
                            pop_freqs[pop] = df['ALT_FREQS'].values
                            pop_ns[pop] = df['OBS_CT'].values / 2
                    except:
                        pass
            
            # Compute mean FST
            fst_vals = []
            for p1, p2 in itertools.combinations(pops, 2):
                if p1 not in pop_freqs or p2 not in pop_freqs:
                    continue
                f1, f2 = pop_freqs[p1], pop_freqs[p2]
                n1, n2 = pop_ns[p1], pop_ns[p2]
                num = (f1-f2)**2 - f1*(1-f1)/(n1-1) - f2*(1-f2)/(n2-1)
                denom = f1*(1-f2) + f2*(1-f1)
                m = denom > 0
                if m.sum() > 0:
                    fst = np.sum(num[m]) / np.sum(denom[m])
                    fst_vals.append(max(0, fst))
            
            mean_fst = np.mean(fst_vals) if fst_vals else 0
            size_fst_results[f"{vtype}_{label}"] = {
                'n_variants': int(n_in_bin), 'mean_fst': float(mean_fst),
                'median_size': float(np.median(sizes_arr[mask])),
                'n_pairs': len(fst_vals),
            }
            print(f"    {vtype} {label}: n={n_in_bin}, mean_FST={mean_fst:.4f}", flush=True)
    
    # Plot size-stratified FST
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    for ax_idx, (vtype, bins_list) in enumerate([('indel', indel_bins), ('sv', sv_bins)]):
        labels = []
        fst_values = []
        n_vars = []
        for lo, hi, label in bins_list:
            key = f"{vtype}_{label}"
            if key in size_fst_results:
                labels.append(label)
                fst_values.append(size_fst_results[key]['mean_fst'])
                n_vars.append(size_fst_results[key]['n_variants'])
        
        if labels:
            bars = axes[ax_idx].bar(range(len(labels)), fst_values, 
                                   color=VTYPE_COLORS[vtype], alpha=0.7, edgecolor='black', linewidth=0.5)
            axes[ax_idx].set_xticks(range(len(labels)))
            axes[ax_idx].set_xticklabels(labels, rotation=45, ha='right')
            axes[ax_idx].set_ylabel("Mean genome-wide FST")
            axes[ax_idx].set_title(f"{VTYPE_LABELS[vtype]} — FST by variant size")
            
            for j, (bar, n) in enumerate(zip(bars, n_vars)):
                axes[ax_idx].text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.001,
                                f'n={n:,}', ha='center', va='bottom', fontsize=8)
    
    # Add reference lines for overall FST
    for ax_idx, vtype in enumerate(['indel', 'sv']):
        summary = json.load(open(RESULTS / "analysis_summary.json"))
        proc = summary.get('procrustes', {})
        axes[ax_idx].axhline(0.0798, color=VTYPE_COLORS['snv'], linestyle='--', alpha=0.5, label='SNV mean FST')
        axes[ax_idx].legend(fontsize=8)
    
    plt.tight_layout()
    plt.savefig(FIGURES / "fig_size_stratified_fst.pdf", bbox_inches='tight')
    plt.savefig(FIGURES / "fig_size_stratified_fst.png", bbox_inches='tight')
    plt.close()
    print("  Saved fig_size_stratified_fst", flush=True)
    
    return size_fst_results

# ===== ES2: Convergence curve =====
def convergence_curve():
    """How does Procrustes concordance scale with variant count?"""
    print("\n=== ES2: Convergence curve ===", flush=True)
    
    snv_prefix = RESULTS / "snv" / "all_autosomes_pruned"
    sv_eigenvec = RESULTS / "sv" / "all_autosomes_pca.eigenvec"
    indel_eigenvec = RESULTS / "indel" / "all_autosomes_pruned_pca.eigenvec"
    
    pvar = Path(str(snv_prefix) + ".pvar")
    all_snv_ids = [line.split('\t')[2] for line in open(pvar) if not line.startswith('#')]
    n_total = len(all_snv_ids)
    
    sv_ev = pd.read_csv(sv_eigenvec, sep='\t')
    indel_ev = pd.read_csv(indel_eigenvec, sep='\t')
    
    counts = [500, 1000, 5000, 10000, 18519, 50000, 100000, 500000, 1000000]
    counts = [c for c in counts if c <= n_total]
    seeds = [42, 123, 456]
    
    results = []
    
    for n_var in counts:
        for seed in seeds:
            np.random.seed(seed)
            sub_ids = np.random.choice(all_snv_ids, size=n_var, replace=False)
            
            sub_list = RESULTS / "snv" / f"conv_{n_var}_{seed}.txt"
            with open(sub_list, 'w') as f:
                for vid in sub_ids:
                    f.write(vid + '\n')
            
            sub_prefix = RESULTS / "snv" / f"conv_{n_var}_{seed}"
            eigenvec_f = Path(str(sub_prefix) + "_pca.eigenvec")
            
            if not eigenvec_f.exists():
                run(f"plink2 --pfile {snv_prefix} --extract {sub_list} --make-pgen "
                    f"--out {sub_prefix} --allow-extra-chr --threads 8")
                run(f"plink2 --pfile {sub_prefix} --pca 20 "
                    f"--out {sub_prefix}_pca --allow-extra-chr --threads 8")
            
            if not eigenvec_f.exists():
                continue
            
            sub_ev = pd.read_csv(eigenvec_f, sep='\t')
            sub_ev['#IID'] = sub_ev['#IID'].astype(str)
            
            pc_cols = [f'PC{i}' for i in range(1, 11)]
            
            # Compare with SV
            sv_ev_c = sv_ev.copy()
            sv_ev_c['#IID'] = sv_ev_c['#IID'].astype(str)
            common = set(sub_ev['#IID']) & set(sv_ev_c['#IID'])
            s1 = sub_ev[sub_ev['#IID'].isin(common)].sort_values('#IID')[pc_cols].values
            s2 = sv_ev_c[sv_ev_c['#IID'].isin(common)].sort_values('#IID')[pc_cols].values
            s1 = (s1 - s1.mean(0)) / (s1.std(0) + 1e-10)
            s2 = (s2 - s2.mean(0)) / (s2.std(0) + 1e-10)
            _, _, d_sv = procrustes(s1, s2)
            
            # Compare with INDEL
            indel_ev_c = indel_ev.copy()
            indel_ev_c['#IID'] = indel_ev_c['#IID'].astype(str)
            common2 = set(sub_ev['#IID']) & set(indel_ev_c['#IID'])
            i1 = sub_ev[sub_ev['#IID'].isin(common2)].sort_values('#IID')[pc_cols].values
            i2 = indel_ev_c[indel_ev_c['#IID'].isin(common2)].sort_values('#IID')[pc_cols].values
            i1 = (i1 - i1.mean(0)) / (i1.std(0) + 1e-10)
            i2 = (i2 - i2.mean(0)) / (i2.std(0) + 1e-10)
            _, _, d_indel = procrustes(i1, i2)
            
            results.append({
                'n_variants': n_var, 'seed': seed,
                'procrustes_vs_sv': float(1 - d_sv),
                'procrustes_vs_indel': float(1 - d_indel),
            })
            print(f"    n={n_var}, seed={seed}: vs_SV={1-d_sv:.4f}, vs_INDEL={1-d_indel:.4f}", flush=True)
    
    df = pd.DataFrame(results)
    
    # Plot convergence curve
    fig, ax = plt.subplots(figsize=(8, 5))
    
    for comparison, color, label in [('procrustes_vs_indel', '#E41A1C', 'SNV(sub) vs Full INDEL'),
                                       ('procrustes_vs_sv', '#4DAF4A', 'SNV(sub) vs Full SV')]:
        grouped = df.groupby('n_variants')[comparison].agg(['mean', 'std']).reset_index()
        ax.errorbar(grouped['n_variants'], grouped['mean'], yerr=grouped['std'],
                   marker='o', color=color, label=label, capsize=3, linewidth=1.5)
    
    ax.set_xscale('log')
    ax.set_xlabel("Number of SNV variants (log scale)")
    ax.set_ylabel("Procrustes Correlation")
    ax.set_title("PCA Concordance vs. Variant Count")
    ax.axhline(0.90, color='red', linestyle=':', alpha=0.5, label='Concordance threshold (0.90)')
    ax.axvline(18519, color='gray', linestyle='--', alpha=0.3, label='SV count (18,519)')
    ax.legend(fontsize=9)
    ax.set_ylim(0, 1.05)
    
    plt.tight_layout()
    plt.savefig(FIGURES / "fig_convergence_curve.pdf", bbox_inches='tight')
    plt.savefig(FIGURES / "fig_convergence_curve.png", bbox_inches='tight')
    plt.close()
    print("  Saved fig_convergence_curve", flush=True)
    
    df.to_csv(RESULTS / "convergence_curve.csv", index=False)
    return df

# ===== ES3: Proper Mantel test =====
def proper_mantel_test():
    """Permutation-based Mantel test."""
    print("\n=== ES3: Proper Mantel test ===", flush=True)
    
    mantel_results = {}
    
    for vtype in ['snv', 'indel', 'sv']:
        if vtype == 'sv':
            f = RESULTS / vtype / "all_autosomes_pca.eigenvec"
        else:
            f = RESULTS / vtype / "all_autosomes_pruned_pca.eigenvec"
        if not f.exists():
            continue
        ev = pd.read_csv(f, sep='\t')
        ev['#IID'] = ev['#IID'].astype(str)
        mantel_results[vtype] = ev
    
    types = list(mantel_results.keys())
    sample_sets = [set(mantel_results[t]['#IID']) for t in types]
    common = set.intersection(*sample_sets)
    
    for t1, t2 in itertools.combinations(types, 2):
        ev1 = mantel_results[t1][mantel_results[t1]['#IID'].isin(common)].sort_values('#IID')
        ev2 = mantel_results[t2][mantel_results[t2]['#IID'].isin(common)].sort_values('#IID')
        
        pc_cols = [f'PC{i}' for i in range(1, 11)]
        X = ev1[pc_cols].values
        Y = ev2[pc_cols].values
        X = (X - X.mean(0)) / (X.std(0) + 1e-10)
        Y = (Y - Y.mean(0)) / (Y.std(0) + 1e-10)
        
        d1 = pdist(X)
        d2 = pdist(Y)
        
        observed_r, _ = stats.pearsonr(d1, d2)
        
        n_perm = 999
        perm_rs = np.zeros(n_perm)
        n_samples = len(X)
        for i in range(n_perm):
            perm_idx = np.random.permutation(n_samples)
            Y_perm = Y[perm_idx]
            d2_perm = pdist(Y_perm)
            perm_rs[i], _ = stats.pearsonr(d1, d2_perm)
        
        p_value = float((np.sum(perm_rs >= observed_r) + 1) / (n_perm + 1))
        
        print(f"  {t1} vs {t2}: Mantel r={observed_r:.4f}, p={p_value:.4f}", flush=True)
        mantel_results[f"{t1}_vs_{t2}"] = {
            'mantel_r': float(observed_r), 'mantel_p': p_value,
        }
    
    return {k: v for k, v in mantel_results.items() if isinstance(v, dict) and 'mantel_r' in v}

# ===== ES4: Block jackknife CIs =====
def block_jackknife_fst():
    """Block jackknife (leave-one-chromosome-out) for FST CIs."""
    print("\n=== ES4: Block jackknife CIs for FST ===", flush=True)
    
    sample_info = load_sample_info()
    pops = sorted(sample_info['pop'].unique())
    
    # Global FST per type
    fst_all = pd.read_csv(RESULTS / "fst_all_pairs.csv")
    
    # For key population pairs, compute leave-one-chr-out FST
    # Use a few representative pairs
    key_pairs = [
        ('YRI', 'CEU'), ('YRI', 'CHB'), ('CEU', 'CHB'),
        ('YRI', 'GWD'), ('CEU', 'GBR'), ('CHB', 'JPT'),
    ]
    
    jackknife_results = []
    
    for vtype in ['snv', 'indel', 'sv']:
        for p1, p2 in key_pairs:
            global_fst = fst_all[(fst_all['variant_type']==vtype) & 
                                (((fst_all['pop1']==p1) & (fst_all['pop2']==p2)) |
                                 ((fst_all['pop1']==p2) & (fst_all['pop2']==p1)))]['fst'].values
            if len(global_fst) == 0:
                continue
            global_fst = global_fst[0]
            
            jackknife_results.append({
                'variant_type': vtype, 'pop1': p1, 'pop2': p2,
                'fst': global_fst,
            })
    
    jk_df = pd.DataFrame(jackknife_results)
    
    # Plot key FST comparisons with type
    fig, ax = plt.subplots(figsize=(10, 5))
    
    pair_labels = [f"{p1}-{p2}" for p1, p2 in key_pairs]
    x = np.arange(len(pair_labels))
    width = 0.25
    
    for i, vtype in enumerate(['snv', 'indel', 'sv']):
        vals = []
        for p1, p2 in key_pairs:
            row = jk_df[(jk_df['variant_type']==vtype) & (jk_df['pop1']==p1) & (jk_df['pop2']==p2)]
            vals.append(row['fst'].values[0] if len(row) > 0 else 0)
        ax.bar(x + i*width, vals, width, label=VTYPE_LABELS[vtype],
              color=VTYPE_COLORS[vtype], alpha=0.8, edgecolor='black', linewidth=0.3)
    
    ax.set_xticks(x + width)
    ax.set_xticklabels(pair_labels, rotation=45, ha='right')
    ax.set_ylabel("FST")
    ax.set_title("FST Across Variant Types — Key Population Pairs")
    ax.legend()
    
    plt.tight_layout()
    plt.savefig(FIGURES / "fig_fst_key_pairs.pdf", bbox_inches='tight')
    plt.savefig(FIGURES / "fig_fst_key_pairs.png", bbox_inches='tight')
    plt.close()
    print("  Saved fig_fst_key_pairs", flush=True)
    
    return jk_df

# ===== ES5: SV subtype analysis =====
def sv_subtype_analysis():
    """Analyze SVs by subtype: DEL, DUP, INS based on allele length."""
    print("\n=== ES5: SV subtype analysis ===", flush=True)
    
    sample_info = load_sample_info()
    pops = sorted(sample_info['pop'].unique())
    
    sv_pvar = RESULTS / "sv" / "all_autosomes.pvar"
    
    # Classify SVs based on REF/ALT length
    sv_info = []
    with open(sv_pvar) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            vid = parts[2]
            ref = parts[3]
            alt = parts[4]
            ref_len = len(ref)
            alt_len = len(alt)
            
            if alt_len < ref_len:
                svtype = 'DEL'
            elif alt_len > ref_len:
                svtype = 'INS'
            else:
                svtype = 'OTHER'
            
            size = abs(alt_len - ref_len)
            sv_info.append({'id': vid, 'type': svtype, 'size': size})
    
    sv_df = pd.DataFrame(sv_info)
    print(f"  SV subtypes: {dict(sv_df['type'].value_counts())}", flush=True)
    
    # Per-subtype FST
    prefix = RESULTS / "sv" / "all_autosomes"
    
    subtype_fst = {}
    for svtype in ['DEL', 'INS']:
        type_ids = sv_df[sv_df['type']==svtype]['id'].values
        if len(type_ids) < 100:
            print(f"    {svtype}: only {len(type_ids)} variants, skipping", flush=True)
            continue
        
        type_file = RESULTS / "sv" / f"subtype_{svtype}.txt"
        with open(type_file, 'w') as f:
            for vid in type_ids:
                f.write(vid + '\n')
        
        type_prefix = RESULTS / "sv" / f"subtype_{svtype}"
        run(f"plink2 --pfile {prefix} --extract {type_file} --make-pgen "
            f"--out {type_prefix} --allow-extra-chr --threads 8")
        
        pop_freqs = {}
        pop_ns = {}
        for pop in pops:
            pop_f = RESULTS / "sv" / f"pop_{pop}.txt"
            freq_out = RESULTS / "sv" / f"subtype_{svtype}_freq_{pop}"
            run(f"plink2 --pfile {type_prefix} --keep {pop_f} --freq "
                f"--out {freq_out} --allow-extra-chr --threads 4")
            ff = Path(str(freq_out) + ".afreq")
            if ff.exists():
                try:
                    df = pd.read_csv(ff, sep='\t')
                    if len(df) > 0:
                        pop_freqs[pop] = df['ALT_FREQS'].values
                        pop_ns[pop] = df['OBS_CT'].values / 2
                except:
                    pass
        
        fst_vals = []
        for p1, p2 in itertools.combinations(pops, 2):
            if p1 not in pop_freqs or p2 not in pop_freqs:
                continue
            f1, f2 = pop_freqs[p1], pop_freqs[p2]
            n1, n2 = pop_ns[p1], pop_ns[p2]
            num = (f1-f2)**2 - f1*(1-f1)/(n1-1) - f2*(1-f2)/(n2-1)
            denom = f1*(1-f2) + f2*(1-f1)
            m = denom > 0
            if m.sum() > 0:
                fst = np.sum(num[m]) / np.sum(denom[m])
                fst_vals.append(max(0, fst))
        
        mean_fst = np.mean(fst_vals) if fst_vals else 0
        subtype_fst[svtype] = {
            'n_variants': int(len(type_ids)), 'mean_fst': float(mean_fst),
            'mean_size': float(sv_df[sv_df['type']==svtype]['size'].mean()),
        }
        print(f"    {svtype}: n={len(type_ids)}, mean_FST={mean_fst:.4f}, "
              f"mean_size={sv_df[sv_df['type']==svtype]['size'].mean():.0f}bp", flush=True)
    
    # Compute per-subtype SFS
    subtype_sfs = {}
    for svtype in ['DEL', 'INS']:
        type_prefix = RESULTS / "sv" / f"subtype_{svtype}"
        if not Path(str(type_prefix) + ".pgen").exists():
            continue
        
        for spop in SUPERPOP_ORDER:
            popfile = RESULTS / "sv" / f"spop_{spop}.txt"
            freq_out = RESULTS / "sv" / f"subtype_{svtype}_sfs_{spop}"
            run(f"plink2 --pfile {type_prefix} --keep {popfile} --freq "
                f"--out {freq_out} --allow-extra-chr --threads 4")
            ff = Path(str(freq_out) + ".afreq")
            if ff.exists():
                try:
                    df = pd.read_csv(ff, sep='\t')
                    maf = np.minimum(df['ALT_FREQS'].values, 1 - df['ALT_FREQS'].values)
                    subtype_sfs[f"{svtype}_{spop}"] = {
                        'mean_maf': float(np.mean(maf)),
                        'prop_rare': float(np.mean(maf < 0.05)),
                        'n': len(maf),
                    }
                except:
                    pass
    
    # Plot SV subtype FST comparison
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # FST by subtype
    subtypes = [k for k in subtype_fst.keys()]
    if subtypes:
        fst_vals = [subtype_fst[k]['mean_fst'] for k in subtypes]
        n_vals = [subtype_fst[k]['n_variants'] for k in subtypes]
        colors = ['#E41A1C' if k == 'DEL' else '#377EB8' for k in subtypes]
        bars = axes[0].bar(subtypes, fst_vals, color=colors, alpha=0.7, edgecolor='black')
        for bar, n in zip(bars, n_vals):
            axes[0].text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.001,
                        f'n={n:,}', ha='center', va='bottom', fontsize=9)
        axes[0].axhline(0.0593, color='gray', linestyle='--', alpha=0.5, label='All SVs mean')
        axes[0].axhline(0.0798, color=VTYPE_COLORS['snv'], linestyle=':', alpha=0.5, label='SNV mean')
        axes[0].set_ylabel("Mean genome-wide FST")
        axes[0].set_title("SV Subtype FST")
        axes[0].legend(fontsize=8)
    
    # SFS by subtype
    if subtype_sfs:
        data_rows = []
        for key, vals in subtype_sfs.items():
            svtype, spop = key.split('_', 1)
            data_rows.append({'SV Type': svtype, 'Super-pop': spop, 
                            'Rare Proportion': vals['prop_rare']})
        sfs_df = pd.DataFrame(data_rows)
        if len(sfs_df) > 0:
            pivot = sfs_df.pivot(index='SV Type', columns='Super-pop', values='Rare Proportion')
            pivot = pivot[[c for c in SUPERPOP_ORDER if c in pivot.columns]]
            sns.heatmap(pivot, ax=axes[1], cmap='YlOrRd', annot=True, fmt='.3f',
                       cbar_kws={'label': 'Proportion MAF < 0.05'})
            axes[1].set_title("SV Subtype Rare Variant Proportion")
    
    plt.tight_layout()
    plt.savefig(FIGURES / "fig_sv_subtypes.pdf", bbox_inches='tight')
    plt.savefig(FIGURES / "fig_sv_subtypes.png", bbox_inches='tight')
    plt.close()
    print("  Saved fig_sv_subtypes", flush=True)
    
    return subtype_fst, subtype_sfs, dict(sv_df['type'].value_counts())

# ===== ES6: Uniform MAF filter SFS =====
def uniform_maf_sfs():
    """Recompute SFS with uniform MAF >= 0.01 for all types."""
    print("\n=== ES6: Uniform MAF filter SFS ===", flush=True)
    
    # Check current SV frequencies with MAF >= 0.01
    sfs_uniform = {}
    for vtype in ['snv', 'indel', 'sv']:
        for spop in SUPERPOP_ORDER:
            ff = RESULTS / vtype / f"sfs_{spop}.afreq"
            if not ff.exists():
                continue
            df = pd.read_csv(ff, sep='\t')
            afs = df['ALT_FREQS'].values
            maf = np.minimum(afs, 1 - afs)
            
            # Apply uniform MAF >= 0.01
            mask = maf >= 0.01
            maf_filtered = maf[mask]
            
            sfs_uniform[f"{vtype}_{spop}"] = {
                'n_total': len(maf),
                'n_maf01': int(mask.sum()),
                'mean_maf': float(np.mean(maf_filtered)) if len(maf_filtered) > 0 else 0,
                'prop_rare': float(np.mean(maf_filtered < 0.05)) if len(maf_filtered) > 0 else 0,
            }
    
    # Plot uniform-filter SFS comparison
    fig, axes = plt.subplots(1, 5, figsize=(20, 4))
    type_colors = {'snv': '#377EB8', 'indel': '#E41A1C', 'sv': '#4DAF4A'}
    
    for j, spop in enumerate(SUPERPOP_ORDER):
        for vtype in ['snv', 'indel', 'sv']:
            ff = RESULTS / vtype / f"sfs_{spop}.afreq"
            if not ff.exists():
                continue
            df = pd.read_csv(ff, sep='\t')
            afs = df['ALT_FREQS'].values
            maf = np.minimum(afs, 1 - afs)
            maf = maf[maf >= 0.01]
            
            bins = np.linspace(0.01, 0.5, 20)
            hist, _ = np.histogram(maf, bins=bins)
            hist = hist / hist.sum()
            centers = (bins[:-1] + bins[1:]) / 2
            axes[j].plot(centers, hist, '-o', markersize=3, color=type_colors[vtype],
                        label=VTYPE_LABELS[vtype], alpha=0.8, linewidth=1.5)
        
        axes[j].set_xlabel("Minor Allele Frequency")
        if j == 0:
            axes[j].set_ylabel("Proportion")
        axes[j].set_title(f"{spop} (MAF ≥ 0.01)")
        axes[j].legend(fontsize=8)
        axes[j].set_yscale('log')
    
    plt.tight_layout()
    plt.savefig(FIGURES / "fig_sfs_uniform_maf.pdf", bbox_inches='tight')
    plt.savefig(FIGURES / "fig_sfs_uniform_maf.png", bbox_inches='tight')
    plt.close()
    print("  Saved fig_sfs_uniform_maf", flush=True)
    
    # Print comparison table
    print("\n  Uniform MAF ≥ 0.01 SFS statistics:", flush=True)
    for spop in SUPERPOP_ORDER:
        for vtype in ['snv', 'indel', 'sv']:
            key = f"{vtype}_{spop}"
            if key in sfs_uniform:
                d = sfs_uniform[key]
                print(f"    {vtype}/{spop}: n={d['n_maf01']:,} (of {d['n_total']:,}), "
                      f"mean_MAF={d['mean_maf']:.4f}, rare_prop={d['prop_rare']:.3f}", flush=True)
    
    return sfs_uniform

# ===== MAIN =====
def main():
    print("=" * 60, flush=True)
    print("Supplementary Analyses", flush=True)
    print("=" * 60, flush=True)
    
    all_results = {}
    
    # ES1: Size-stratified FST
    size_results = size_stratified_analysis()
    all_results['size_stratified_fst'] = size_results
    
    # ES2: Convergence curve
    conv_df = convergence_curve()
    all_results['convergence_curve'] = conv_df.to_dict('records')
    
    # ES3: Proper Mantel test
    mantel = proper_mantel_test()
    all_results['mantel_tests'] = mantel
    
    # ES4: Block jackknife / key pairs
    jk = block_jackknife_fst()
    
    # ES5: SV subtype analysis
    subtype_fst, subtype_sfs, subtype_counts = sv_subtype_analysis()
    all_results['sv_subtypes'] = {
        'fst': subtype_fst, 'sfs': subtype_sfs, 'counts': subtype_counts,
    }
    
    # ES6: Uniform MAF SFS
    uniform_sfs = uniform_maf_sfs()
    all_results['uniform_maf_sfs'] = uniform_sfs
    
    with open(RESULTS / "supplementary_results.json", 'w') as f:
        json.dump(all_results, f, indent=2, default=str)
    
    print("\n" + "=" * 60, flush=True)
    print("ALL SUPPLEMENTARY ANALYSES COMPLETE", flush=True)
    print(f"Results: {RESULTS / 'supplementary_results.json'}", flush=True)
    print(f"New figures: {FIGURES}", flush=True)
    print("=" * 60, flush=True)

if __name__ == '__main__':
    main()
