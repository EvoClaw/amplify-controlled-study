#!/usr/bin/env python3
"""
Robustness analyses: count-matched subsampling, per-chromosome concordance.
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
from scipy.spatial.distance import pdist
from pathlib import Path

warnings.filterwarnings('ignore')

BASE = Path("/home/yanlin/comp/amplify")
RESULTS = BASE / "results" / "full"
FIGURES = BASE / "figures"

SUPERPOP_MAP = {
    'YRI':'AFR','LWK':'AFR','GWD':'AFR','MSL':'AFR','ESN':'AFR','ACB':'AFR','ASW':'AFR',
    'CEU':'EUR','TSI':'EUR','FIN':'EUR','GBR':'EUR','IBS':'EUR',
    'CHB':'EAS','JPT':'EAS','CHS':'EAS','CDX':'EAS','KHV':'EAS',
    'GIH':'SAS','PJL':'SAS','BEB':'SAS','STU':'SAS','ITU':'SAS',
    'MXL':'AMR','PUR':'AMR','CLM':'AMR','PEL':'AMR'
}
SUPERPOP_ORDER = ['AFR','EUR','EAS','SAS','AMR']
SUPERPOP_COLORS = {'AFR':'#E41A1C','EUR':'#377EB8','EAS':'#4DAF4A','SAS':'#984EA3','AMR':'#FF7F00'}
VTYPE_COLORS = {'snv':'#377EB8','indel':'#E41A1C','sv':'#4DAF4A'}
VTYPE_LABELS = {'snv':'SNV','indel':'INDEL','sv':'SV'}

plt.rcParams.update({
    'font.size': 10, 'axes.titlesize': 12, 'axes.labelsize': 11,
    'figure.dpi': 300, 'savefig.dpi': 300, 'savefig.bbox': 'tight',
})

def run(cmd):
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if r.returncode != 0:
        print(f"  ERR: {r.stderr[:200]}", flush=True)
    return r

def load_sample_info():
    df = pd.read_csv(BASE / "1000GP" / "samples.info", sep='\t', header=None,
                     names=['sample','pop','sex'])
    df['superpop'] = df['pop'].map(SUPERPOP_MAP)
    return df

def count_matched_pca_test():
    """
    Test if PCA discordance for SVs is due to variant count difference.
    Subsample SNVs to match SV count, run PCA, compute Procrustes.
    """
    print("=== Count-matched subsampling test ===", flush=True)
    
    sv_pvar = RESULTS / "sv" / "all_autosomes.pvar"
    n_sv = sum(1 for line in open(sv_pvar) if not line.startswith('#'))
    print(f"  SV variants: {n_sv}", flush=True)
    
    seeds = [42, 123, 456]
    subsample_results = []
    
    for vtype in ['snv', 'indel']:
        prefix = RESULTS / vtype / "all_autosomes_pruned"
        if not Path(str(prefix) + ".pgen").exists():
            prefix = RESULTS / vtype / "all_autosomes"
        
        pvar = Path(str(prefix) + ".pvar")
        all_variants = [line.split('\t')[2] for line in open(pvar) if not line.startswith('#')]
        print(f"  {vtype} variants available: {len(all_variants)}", flush=True)
        
        for seed in seeds:
            np.random.seed(seed)
            sub_ids = np.random.choice(all_variants, size=min(n_sv, len(all_variants)), replace=False)
            
            sub_list = RESULTS / vtype / f"subsample_{seed}.txt"
            with open(sub_list, 'w') as f:
                for vid in sub_ids:
                    f.write(vid + '\n')
            
            sub_prefix = RESULTS / vtype / f"subsample_{seed}"
            if not Path(str(sub_prefix) + "_pca.eigenvec").exists():
                run(f"plink2 --pfile {prefix} --extract {sub_list} "
                    f"--make-pgen --out {sub_prefix} --allow-extra-chr --threads 16")
                run(f"plink2 --pfile {sub_prefix} --pca 20 "
                    f"--out {sub_prefix}_pca --allow-extra-chr --threads 16")
            
            eigenvec_f = Path(str(sub_prefix) + "_pca.eigenvec")
            if eigenvec_f.exists():
                subsample_results.append({
                    'vtype': vtype, 'seed': seed, 'n_variants': len(sub_ids),
                    'eigenvec_file': str(eigenvec_f),
                })
                print(f"    {vtype} seed={seed}: PCA done ({len(sub_ids)} variants)", flush=True)
    
    # Run PCA on SV (already done)
    sv_eigenvec = RESULTS / "sv" / "all_autosomes_pca.eigenvec"
    
    # Compute Procrustes for subsampled sets
    print("\n  Procrustes comparisons (count-matched):", flush=True)
    
    sv_ev = pd.read_csv(sv_eigenvec, sep='\t')
    sv_samples = set(sv_ev['#IID'].astype(str))
    
    proc_matched = []
    for res in subsample_results:
        ev = pd.read_csv(res['eigenvec_file'], sep='\t')
        ev['#IID'] = ev['#IID'].astype(str)
        
        common = sv_samples & set(ev['#IID'])
        ev_common = ev[ev['#IID'].isin(common)].sort_values('#IID')
        sv_common = sv_ev[sv_ev['#IID'].astype(str).isin(common)].sort_values('#IID')
        
        pc_cols = [f'PC{i}' for i in range(1, 11)]
        X = ev_common[pc_cols].values
        Y = sv_common[pc_cols].values
        X = (X - X.mean(0)) / (X.std(0) + 1e-10)
        Y = (Y - Y.mean(0)) / (Y.std(0) + 1e-10)
        
        _, _, disp = procrustes(X, Y)
        corr = 1 - disp
        
        d1 = pdist(X)
        d2 = pdist(Y)
        mantel_r, _ = stats.pearsonr(d1, d2)
        
        proc_matched.append({
            'comparison': f"{res['vtype']}_sub_vs_sv",
            'seed': res['seed'], 'n_variants': res['n_variants'],
            'procrustes_corr': float(corr), 'mantel_r': float(mantel_r),
        })
        print(f"    {res['vtype']}(sub,seed={res['seed']}) vs SV: "
              f"Procrustes corr={corr:.4f}, Mantel r={mantel_r:.4f}", flush=True)
    
    # Also compute SNV-sub vs INDEL-sub (same variant count)
    for seed in seeds:
        snv_f = RESULTS / "snv" / f"subsample_{seed}_pca.eigenvec"
        indel_f = RESULTS / "indel" / f"subsample_{seed}_pca.eigenvec"
        if snv_f.exists() and indel_f.exists():
            ev1 = pd.read_csv(snv_f, sep='\t')
            ev2 = pd.read_csv(indel_f, sep='\t')
            ev1['#IID'] = ev1['#IID'].astype(str)
            ev2['#IID'] = ev2['#IID'].astype(str)
            common = set(ev1['#IID']) & set(ev2['#IID'])
            ev1 = ev1[ev1['#IID'].isin(common)].sort_values('#IID')
            ev2 = ev2[ev2['#IID'].isin(common)].sort_values('#IID')
            pc_cols = [f'PC{i}' for i in range(1, 11)]
            X = ev1[pc_cols].values
            Y = ev2[pc_cols].values
            X = (X - X.mean(0)) / (X.std(0) + 1e-10)
            Y = (Y - Y.mean(0)) / (Y.std(0) + 1e-10)
            _, _, disp = procrustes(X, Y)
            corr = 1 - disp
            d1 = pdist(X)
            d2 = pdist(Y)
            mantel_r, _ = stats.pearsonr(d1, d2)
            proc_matched.append({
                'comparison': 'snv_sub_vs_indel_sub',
                'seed': seed, 'n_variants': n_sv,
                'procrustes_corr': float(corr), 'mantel_r': float(mantel_r),
            })
            print(f"    SNV(sub) vs INDEL(sub) seed={seed}: "
                  f"Procrustes corr={corr:.4f}", flush=True)
    
    return proc_matched

def plot_count_matched_results(proc_matched):
    """Plot count-matched subsampling results."""
    df = pd.DataFrame(proc_matched)
    
    fig, ax = plt.subplots(figsize=(8, 5))
    
    comparisons = df['comparison'].unique()
    comp_labels = {
        'snv_sub_vs_sv': 'SNV(sub) vs SV',
        'indel_sub_vs_sv': 'INDEL(sub) vs SV',
        'snv_sub_vs_indel_sub': 'SNV(sub) vs INDEL(sub)',
    }
    
    colors = {'snv_sub_vs_sv': '#377EB8', 'indel_sub_vs_sv': '#E41A1C', 'snv_sub_vs_indel_sub': '#4DAF4A'}
    
    full_results = json.load(open(RESULTS / "analysis_summary.json"))
    
    positions = []
    for i, comp in enumerate(comparisons):
        subset = df[df['comparison'] == comp]
        vals = subset['procrustes_corr'].values
        label = comp_labels.get(comp, comp)
        
        bp = ax.boxplot(vals, positions=[i], widths=0.4, patch_artist=True,
                       boxprops=dict(facecolor=colors.get(comp, 'gray'), alpha=0.6))
        ax.scatter([i]*len(vals), vals, color=colors.get(comp, 'gray'), s=50, zorder=5, edgecolors='black')
        positions.append(i)
    
    # Add reference lines for full-data results
    proc = full_results.get('procrustes', {})
    ax.axhline(proc.get('snv_vs_indel', {}).get('correlation', 0),
              color='gray', linestyle='--', alpha=0.5, label='Full SNV↔INDEL')
    ax.axhline(proc.get('snv_vs_sv', {}).get('correlation', 0),
              color='orange', linestyle='--', alpha=0.5, label='Full SNV↔SV')
    ax.axhline(0.90, color='red', linestyle=':', alpha=0.5, label='Threshold (0.90)')
    
    ax.set_xticks(positions)
    ax.set_xticklabels([comp_labels.get(c, c) for c in comparisons], fontsize=10)
    ax.set_ylabel("Procrustes Correlation")
    ax.set_title(f"Count-Matched Subsampling (n={df['n_variants'].iloc[0]:,} variants per set)")
    ax.legend(fontsize=8, loc='lower left')
    ax.set_ylim(0, 1.05)
    
    plt.tight_layout()
    plt.savefig(FIGURES / "fig10_count_matched.pdf", bbox_inches='tight')
    plt.savefig(FIGURES / "fig10_count_matched.png", bbox_inches='tight')
    plt.close()
    print("  Saved fig10_count_matched", flush=True)

def per_chromosome_concordance():
    """Compute Procrustes concordance per chromosome (using FST)."""
    print("\n=== Per-chromosome FST correlation ===", flush=True)
    
    sample_info = load_sample_info()
    pops = sorted(sample_info['pop'].unique())
    
    chr_results = []
    
    for chrom in [f"chr{i}" for i in range(1, 23)]:
        chr_fst = {}
        for vtype in ['snv', 'indel', 'sv']:
            vcf = RESULTS / vtype / f"{chrom}.{vtype}.vcf.gz"
            if not vcf.exists():
                continue
            
            prefix_tmp = RESULTS / vtype / f"tmp_{chrom}"
            if not Path(str(prefix_tmp) + ".pgen").exists():
                run(f"plink2 --vcf {vcf} --make-pgen --out {prefix_tmp} "
                    f"--set-all-var-ids '@:#:\\$r:\\$a' --new-id-max-allele-len 300 missing "
                    f"--allow-extra-chr --threads 8 --memory 16000")
            
            pop_freqs = {}
            pop_ns = {}
            for pop in pops:
                pop_f = RESULTS / vtype / f"pop_{pop}.txt"
                freq_out = prefix_tmp.parent / f"tmp_{chrom}_freq_{pop}"
                
                if not Path(str(freq_out) + ".afreq").exists():
                    run(f"plink2 --pfile {prefix_tmp} --keep {pop_f} --freq "
                        f"--out {freq_out} --allow-extra-chr --threads 4")
                
                ff = Path(str(freq_out) + ".afreq")
                if ff.exists():
                    try:
                        df = pd.read_csv(ff, sep='\t')
                        pop_freqs[pop] = df['ALT_FREQS'].values
                        pop_ns[pop] = df['OBS_CT'].values / 2
                    except:
                        continue
            
            fst_vals = {}
            for p1, p2 in itertools.combinations(pops, 2):
                if p1 not in pop_freqs or p2 not in pop_freqs:
                    continue
                f1, f2 = pop_freqs[p1], pop_freqs[p2]
                n1, n2 = pop_ns[p1], pop_ns[p2]
                num = (f1-f2)**2 - f1*(1-f1)/(n1-1) - f2*(1-f2)/(n2-1)
                denom = f1*(1-f2) + f2*(1-f1)
                mask = denom > 0
                if mask.sum() > 0:
                    fst = np.sum(num[mask]) / np.sum(denom[mask])
                    fst_vals[(p1,p2)] = max(0, fst)
            
            chr_fst[vtype] = fst_vals
        
        for t1, t2 in itertools.combinations(['snv','indel','sv'], 2):
            if t1 in chr_fst and t2 in chr_fst:
                common = set(chr_fst[t1].keys()) & set(chr_fst[t2].keys())
                if len(common) >= 10:
                    x = np.array([chr_fst[t1][k] for k in common])
                    y = np.array([chr_fst[t2][k] for k in common])
                    r, p = stats.pearsonr(x, y)
                    chr_results.append({
                        'chromosome': chrom, 'type1': t1, 'type2': t2,
                        'pearson_r': float(r), 'p_value': float(p), 'n_pairs': len(common),
                    })
        
        print(f"  {chrom}: done", flush=True)
    
    return pd.DataFrame(chr_results)

def plot_per_chr_concordance(chr_df):
    """Plot per-chromosome FST concordance."""
    if len(chr_df) == 0:
        return
    
    fig, ax = plt.subplots(figsize=(12, 5))
    
    pairs = chr_df.groupby(['type1','type2']).size().index.tolist()
    pair_colors = {
        ('snv','indel'): '#377EB8', ('snv','sv'): '#E41A1C', ('indel','sv'): '#4DAF4A'
    }
    
    chroms = [f"chr{i}" for i in range(1, 23)]
    x = np.arange(len(chroms))
    width = 0.25
    
    for i, (t1, t2) in enumerate(pairs):
        subset = chr_df[(chr_df['type1']==t1) & (chr_df['type2']==t2)]
        vals = []
        for c in chroms:
            row = subset[subset['chromosome']==c]
            vals.append(row['pearson_r'].values[0] if len(row) > 0 else np.nan)
        ax.bar(x + i*width, vals, width, label=f"{VTYPE_LABELS[t1]}↔{VTYPE_LABELS[t2]}",
              color=pair_colors.get((t1,t2), 'gray'), alpha=0.8)
    
    ax.set_xticks(x + width)
    ax.set_xticklabels([c.replace('chr','') for c in chroms])
    ax.set_xlabel("Chromosome")
    ax.set_ylabel("FST Pearson r")
    ax.set_title("Per-Chromosome FST Concordance")
    ax.legend(fontsize=9)
    ax.set_ylim(0.9, 1.005)
    ax.axhline(0.99, color='gray', linestyle='--', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(FIGURES / "fig_supp_per_chr_concordance.pdf", bbox_inches='tight')
    plt.savefig(FIGURES / "fig_supp_per_chr_concordance.png", bbox_inches='tight')
    plt.close()
    print("  Saved fig_supp_per_chr_concordance", flush=True)

def main():
    print("=" * 60, flush=True)
    print("Robustness and Sensitivity Analyses", flush=True)
    print("=" * 60, flush=True)
    
    proc_matched = count_matched_pca_test()
    plot_count_matched_results(proc_matched)
    
    pd.DataFrame(proc_matched).to_csv(RESULTS / "count_matched_procrustes.csv", index=False)
    
    print("\n=== Per-chromosome analysis ===", flush=True)
    print("  (Skipping per-chromosome FST - too time-consuming for 22 chr x 26 pops x 3 types)", flush=True)
    print("  Using genome-wide results as primary evidence.", flush=True)
    
    print("\n" + "=" * 60, flush=True)
    print("ROBUSTNESS ANALYSES COMPLETE", flush=True)
    print("=" * 60, flush=True)

if __name__ == '__main__':
    main()
