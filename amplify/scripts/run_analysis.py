#!/usr/bin/env python3
"""
Multi-variant-type population genomics analysis pipeline.
Compares population structure, differentiation, and SFS across SNVs, INDELs, and SVs
in the 1000 Genomes Project 2022 high-coverage panel.
"""
import os, sys, subprocess, warnings, itertools
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.spatial import procrustes
from collections import defaultdict, Counter
from pathlib import Path

warnings.filterwarnings('ignore')
np.random.seed(42)

BASE = Path("/home/yanlin/comp/amplify")
DATA = BASE / "1000GP" / "20220422_3202_phased_SNV_INDEL_SV"
RESULTS = BASE / "results" / "full"
FIGURES = BASE / "figures"
RESULTS.mkdir(parents=True, exist_ok=True)
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

CHROMS = [f"chr{i}" for i in range(1, 23)]

def load_sample_info():
    df = pd.read_csv(BASE / "1000GP" / "samples.info", sep='\t', header=None,
                     names=['sample','pop','sex'])
    df['superpop'] = df['pop'].map(SUPERPOP_MAP)
    return df

def run(cmd, **kwargs):
    print(f"  CMD: {cmd[:120]}...")
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True, **kwargs)
    if r.returncode != 0:
        print(f"  STDERR: {r.stderr[:300]}")
    return r

# ===== STEP 1: Extract variant types from VCFs =====
def extract_variants():
    """Extract SNVs, INDELs, SVs from each chromosome VCF."""
    for vtype in ['snv', 'indel', 'sv']:
        (RESULTS / vtype).mkdir(exist_ok=True)

    counts = {'snv': 0, 'indel': 0, 'sv': 0}

    for chrom in CHROMS:
        vcf = DATA / f"1kGP_high_coverage_Illumina.{chrom}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
        if not vcf.exists():
            print(f"  Skipping {chrom} - VCF not found")
            continue

        snv_out = RESULTS / "snv" / f"{chrom}.snv.vcf.gz"
        indel_out = RESULTS / "indel" / f"{chrom}.indel.vcf.gz"
        sv_out = RESULTS / "sv" / f"{chrom}.sv.vcf.gz"

        if not snv_out.exists():
            run(f"bcftools view -v snps -m2 -M2 --min-af 0.01:minor -i 'F_MISSING<0.05' "
                f"{vcf} -Oz -o {snv_out} --threads 8")
            run(f"bcftools index {snv_out} --threads 4")

        if not indel_out.exists():
            run(f"bcftools view -v indels -m2 -M2 --min-af 0.01:minor -i 'F_MISSING<0.05' "
                f"{vcf} -Oz -o {indel_out} --threads 8")
            run(f"bcftools index {indel_out} --threads 4")

        if not sv_out.exists():
            run(f"bcftools view -m2 -M2 --min-af 0.005:minor "
                f"-i 'F_MISSING<0.10 && (STRLEN(REF)>=50 || STRLEN(ALT)>=50)' "
                f"{vcf} -Oz -o {sv_out} --threads 8")
            run(f"bcftools index {sv_out} --threads 4")

        for vtype, path in [('snv', snv_out), ('indel', indel_out), ('sv', sv_out)]:
            if path.exists():
                r = run(f"bcftools view -H {path} | wc -l")
                n = int(r.stdout.strip()) if r.stdout.strip() else 0
                counts[vtype] += n

        print(f"  {chrom} done")

    print(f"\nTotal variant counts: {counts}")
    return counts

# ===== STEP 2: Merge and prepare for PCA =====
def merge_for_pca():
    """Merge per-chrom VCFs and convert to PLINK format for PCA."""
    for vtype in ['snv', 'indel', 'sv']:
        merged = RESULTS / vtype / "all_autosomes.vcf.gz"
        if merged.exists():
            continue

        vcf_list = RESULTS / vtype / "vcf_list.txt"
        files = sorted((RESULTS / vtype).glob("chr*.*.vcf.gz"))
        files = [f for f in files if 'all_autosomes' not in f.name and f.stat().st_size > 100]
        with open(vcf_list, 'w') as f:
            for p in files:
                f.write(str(p) + '\n')

        if len(files) == 0:
            print(f"  No files for {vtype}")
            continue

        run(f"bcftools concat -f {vcf_list} -Oz -o {merged} --threads 8")
        run(f"bcftools index {merged} --threads 4")

    for vtype in ['snv', 'indel', 'sv']:
        merged = RESULTS / vtype / "all_autosomes.vcf.gz"
        plink_prefix = RESULTS / vtype / "all_autosomes"
        if not merged.exists():
            continue
        if (RESULTS / vtype / "all_autosomes.pgen").exists():
            continue

        extra = ""
        if vtype == 'snv':
            extra = "--indep-pairwise 50 5 0.2"
        elif vtype == 'indel':
            extra = "--indep-pairwise 50 5 0.2"

        run(f"plink2 --vcf {merged} --make-pgen --out {plink_prefix} "
            f"--set-all-var-ids '@:#:\\$r:\\$a' --new-id-max-allele-len 20 "
            f"--allow-extra-chr --threads 16 --memory 64000")

        if extra:
            run(f"plink2 --pfile {plink_prefix} {extra} "
                f"--out {plink_prefix}_ldprune --allow-extra-chr --threads 16")
            prune_file = RESULTS / vtype / "all_autosomes_ldprune.prune.in"
            if prune_file.exists():
                run(f"plink2 --pfile {plink_prefix} --extract {prune_file} "
                    f"--make-pgen --out {plink_prefix}_pruned "
                    f"--allow-extra-chr --threads 16")

# ===== STEP 3: Run PCA =====
def run_pca():
    """Run PCA per variant type using plink2."""
    results = {}
    for vtype in ['snv', 'indel', 'sv']:
        if vtype in ('snv', 'indel'):
            prefix = RESULTS / vtype / "all_autosomes_pruned"
        else:
            prefix = RESULTS / vtype / "all_autosomes"

        pgen = Path(str(prefix) + ".pgen")
        if not pgen.exists():
            prefix = RESULTS / vtype / "all_autosomes"
            pgen = Path(str(prefix) + ".pgen")

        if not pgen.exists():
            print(f"  No pgen for {vtype}")
            continue

        eigenvec_f = Path(str(prefix) + "_pca.eigenvec")
        eigenval_f = Path(str(prefix) + "_pca.eigenval")
        if not eigenvec_f.exists():
            run(f"plink2 --pfile {prefix} --pca 20 "
                f"--out {prefix}_pca --allow-extra-chr --threads 16")

        if eigenvec_f.exists():
            eigenvec = pd.read_csv(eigenvec_f, sep='\t')
            eigenval = pd.read_csv(eigenval_f, header=None)[0].values
            results[vtype] = {'eigenvec': eigenvec, 'eigenval': eigenval}
            nvar_r = run(f"plink2 --pfile {prefix} --write-snplist --out /tmp/count_{vtype} --allow-extra-chr")
            snplist = Path(f"/tmp/count_{vtype}.snplist")
            nvar = sum(1 for _ in open(snplist)) if snplist.exists() else 0
            results[vtype]['nvar'] = nvar
            print(f"  PCA {vtype}: {nvar} variants, {len(eigenvec)} samples, top eigenval={eigenval[0]:.2f}")

    return results

# ===== STEP 4: Compute FST =====
def compute_fst(sample_info):
    """Compute pairwise FST for all population pairs, per variant type."""
    pop_file = RESULTS / "pop_assignments.txt"
    sample_info[['sample','pop']].to_csv(pop_file, sep='\t', index=False, header=False)

    pops = sorted(sample_info['pop'].unique())
    pairs = list(itertools.combinations(pops, 2))
    fst_results = {}

    for vtype in ['snv', 'indel', 'sv']:
        prefix = RESULTS / vtype / "all_autosomes"
        pgen = Path(str(prefix) + ".pgen")
        if not pgen.exists():
            continue

        psam = pd.read_csv(str(prefix) + ".psam", sep='\t')
        sample_col = psam.columns[0]  # #FID or #IID
        psam_samples = psam[sample_col].astype(str) if '#FID' in psam.columns else psam.iloc[:,0].astype(str)

        fst_vals = {}
        for pop1, pop2 in pairs:
            s1 = sample_info[sample_info['pop']==pop1]['sample'].values
            s2 = sample_info[sample_info['pop']==pop2]['sample'].values
            if len(s1) < 5 or len(s2) < 5:
                continue

            pop1_f = RESULTS / vtype / f"pop_{pop1}.txt"
            pop2_f = RESULTS / vtype / f"pop_{pop2}.txt"

            with open(pop1_f, 'w') as f:
                for s in s1:
                    f.write(f"{s}\t{s}\n")
            with open(pop2_f, 'w') as f:
                for s in s2:
                    f.write(f"{s}\t{s}\n")

            fst_out = RESULTS / vtype / f"fst_{pop1}_{pop2}"
            r = run(f"plink2 --pfile {prefix} --fst CATPHENO "
                    f"--pheno-name POP --allow-extra-chr --threads 8 "
                    f"--keep {pop1_f} --keep {pop2_f} "
                    f"--out {fst_out} 2>/dev/null || true")

        fst_log_pattern = RESULTS / vtype / "fst_*_*.log"
        import glob
        for logf in glob.glob(str(fst_log_pattern)):
            pname = Path(logf).stem.replace("fst_","")
            with open(logf) as f:
                for line in f:
                    if "Fst estimate:" in line:
                        val = float(line.split(":")[-1].strip())
                        fst_vals[pname] = val

        fst_results[vtype] = fst_vals
        print(f"  FST {vtype}: {len(fst_vals)} population pairs computed")

    return fst_results

def compute_fst_python(sample_info):
    """Compute genome-wide weighted FST using allele frequencies from VCFs."""
    import allel

    pops = sorted(sample_info['pop'].unique())
    pop_to_idx = {}

    fst_results = {}

    for vtype in ['snv', 'indel', 'sv']:
        merged = RESULTS / vtype / "all_autosomes.vcf.gz"
        if not merged.exists():
            continue

        # Read genotypes in chunks for memory efficiency
        print(f"\n  Computing FST for {vtype}...")

        psam = None
        prefix = RESULTS / vtype / "all_autosomes"
        psam_f = Path(str(prefix) + ".psam")
        if psam_f.exists():
            psam = pd.read_csv(psam_f, sep='\t')
            samples_in_plink = list(psam.iloc[:, 0].astype(str))
        else:
            r = run(f"bcftools query -l {merged}")
            samples_in_plink = r.stdout.strip().split('\n')

        sample_to_idx = {s: i for i, s in enumerate(samples_in_plink)}

        pop_indices = {}
        for pop in pops:
            pop_samples = sample_info[sample_info['pop']==pop]['sample'].values
            idx = [sample_to_idx[s] for s in pop_samples if s in sample_to_idx]
            if len(idx) >= 5:
                pop_indices[pop] = idx

        usable_pops = sorted(pop_indices.keys())
        pairs = list(itertools.combinations(usable_pops, 2))

        # Use allele frequencies from plink
        freq_f = RESULTS / vtype / "all_autosomes.afreq"
        if not freq_f.exists():
            run(f"plink2 --pfile {prefix} --freq --out {prefix} --allow-extra-chr --threads 16")

        # Compute per-population frequencies
        pop_freq_files = {}
        for pop in usable_pops:
            pop_f = RESULTS / vtype / f"pop_{pop}.txt"
            with open(pop_f, 'w') as f:
                for s in sample_info[sample_info['pop']==pop]['sample'].values:
                    f.write(f"{s}\t{s}\n")

            pfreq = RESULTS / vtype / f"freq_{pop}"
            if not Path(str(pfreq) + ".afreq").exists():
                run(f"plink2 --pfile {prefix} --keep {pop_f} --freq "
                    f"--out {pfreq} --allow-extra-chr --threads 8")
            pop_freq_files[pop] = str(pfreq) + ".afreq"

        # Read frequencies and compute Hudson FST
        pop_freqs = {}
        pop_ns = {}
        for pop in usable_pops:
            ff = pop_freq_files[pop]
            if not os.path.exists(ff):
                continue
            df = pd.read_csv(ff, sep='\t')
            pop_freqs[pop] = df['ALT_FREQS'].values
            pop_ns[pop] = df['OBS_CT'].values / 2  # diploid
            print(f"    {pop}: n={int(pop_ns[pop][0])}, {len(pop_freqs[pop])} variants")

        fst_vals = {}
        for pop1, pop2 in pairs:
            if pop1 not in pop_freqs or pop2 not in pop_freqs:
                continue
            p1 = pop_freqs[pop1]
            p2 = pop_freqs[pop2]
            n1 = pop_ns[pop1]
            n2 = pop_ns[pop2]

            # Hudson FST estimator
            num = (p1 - p2)**2 - p1*(1-p1)/(n1-1) - p2*(1-p2)/(n2-1)
            denom = p1*(1-p2) + p2*(1-p1)
            mask = denom > 0
            fst = np.sum(num[mask]) / np.sum(denom[mask])
            fst_vals[f"{pop1}_{pop2}"] = max(0, fst)

        fst_results[vtype] = fst_vals
        print(f"  FST {vtype}: {len(fst_vals)} pairs computed")

    return fst_results

# ===== STEP 5: Compute SFS =====
def compute_sfs(sample_info):
    """Compute folded SFS per variant type per super-population."""
    sfs_results = {}

    for vtype in ['snv', 'indel', 'sv']:
        prefix = RESULTS / vtype / "all_autosomes"
        sfs_results[vtype] = {}

        for spop in SUPERPOP_ORDER:
            spop_samples = sample_info[sample_info['superpop']==spop]['sample'].values
            pop_f = RESULTS / vtype / f"spop_{spop}.txt"
            with open(pop_f, 'w') as f:
                for s in spop_samples:
                    f.write(f"{s}\t{s}\n")

            freq_out = RESULTS / vtype / f"sfs_{spop}"
            freq_file = Path(str(freq_out) + ".afreq")
            if not freq_file.exists():
                run(f"plink2 --pfile {prefix} --keep {pop_f} --freq "
                    f"--out {freq_out} --allow-extra-chr --threads 8")

            if freq_file.exists():
                df = pd.read_csv(freq_file, sep='\t')
                afs = df['ALT_FREQS'].values
                maf = np.minimum(afs, 1 - afs)
                bins = np.linspace(0, 0.5, 21)
                hist, _ = np.histogram(maf, bins=bins)
                sfs_results[vtype][spop] = {
                    'hist': hist, 'bins': bins,
                    'n_variants': len(maf),
                    'mean_maf': np.mean(maf),
                    'median_maf': np.median(maf),
                    'prop_rare': np.mean(maf < 0.05)
                }
                print(f"  SFS {vtype}/{spop}: n={len(maf)}, mean_MAF={np.mean(maf):.4f}, rare_prop={np.mean(maf<0.05):.3f}")

    return sfs_results

# ===== STEP 6: Procrustes Analysis =====
def procrustes_analysis(pca_results, sample_info):
    """Compare PCA spaces across variant types using Procrustes analysis."""
    types = [t for t in ['snv', 'indel', 'sv'] if t in pca_results]
    if len(types) < 2:
        return {}

    # Align samples across types
    sample_sets = []
    for vtype in types:
        ev = pca_results[vtype]['eigenvec']
        col = '#IID' if '#IID' in ev.columns else ev.columns[0]
        sample_sets.append(set(ev[col].astype(str)))

    common = set.intersection(*sample_sets)
    print(f"  Common samples across types: {len(common)}")

    proc_results = {}
    for t1, t2 in itertools.combinations(types, 2):
        ev1 = pca_results[t1]['eigenvec'].copy()
        ev2 = pca_results[t2]['eigenvec'].copy()
        col1 = '#IID' if '#IID' in ev1.columns else ev1.columns[0]
        col2 = '#IID' if '#IID' in ev2.columns else ev2.columns[0]
        ev1[col1] = ev1[col1].astype(str)
        ev2[col2] = ev2[col2].astype(str)
        ev1 = ev1[ev1[col1].isin(common)].sort_values(col1)
        ev2 = ev2[ev2[col2].isin(common)].sort_values(col2)

        pc_cols1 = [c for c in ev1.columns if c.startswith('PC')]
        pc_cols2 = [c for c in ev2.columns if c.startswith('PC')]
        X = ev1[pc_cols1[:10]].values
        Y = ev2[pc_cols2[:10]].values

        # Standardize
        X = (X - X.mean(0)) / X.std(0)
        Y = (Y - Y.mean(0)) / Y.std(0)

        _, _, disparity = procrustes(X, Y)
        correlation = 1 - disparity

        # Permutation test
        n_perm = 1000
        perm_disps = []
        for _ in range(n_perm):
            Y_perm = Y[np.random.permutation(len(Y))]
            _, _, d = procrustes(X, Y_perm)
            perm_disps.append(d)
        p_value = np.mean(np.array(perm_disps) <= disparity)

        # Mantel test on distance matrices
        from scipy.spatial.distance import pdist, squareform
        d1 = pdist(X)
        d2 = pdist(Y)
        mantel_r, mantel_p = stats.pearsonr(d1, d2)

        proc_results[f"{t1}_vs_{t2}"] = {
            'disparity': disparity,
            'correlation': correlation,
            'p_value': p_value,
            'mantel_r': mantel_r,
            'mantel_p': mantel_p,
            'n_samples': len(common)
        }
        print(f"  Procrustes {t1} vs {t2}: corr={correlation:.4f} (p={p_value:.4f}), "
              f"Mantel r={mantel_r:.4f} (p={mantel_p:.2e})")

    return proc_results

# ===== PLOTTING =====
def plot_pca_comparison(pca_results, sample_info):
    """Figure 1: PCA comparison across variant types."""
    types = [t for t in ['snv', 'indel', 'sv'] if t in pca_results]
    fig, axes = plt.subplots(1, len(types), figsize=(5*len(types), 5))
    if len(types) == 1:
        axes = [axes]

    for i, vtype in enumerate(types):
        ev = pca_results[vtype]['eigenvec'].copy()
        eigenval = pca_results[vtype]['eigenval']
        col = '#IID' if '#IID' in ev.columns else ev.columns[0]
        ev[col] = ev[col].astype(str)
        ev = ev.merge(sample_info[['sample','superpop']], left_on=col, right_on='sample', how='left')

        var_exp = eigenval / eigenval.sum() * 100

        pc_cols = [c for c in ev.columns if c.startswith('PC')]
        for spop in SUPERPOP_ORDER:
            mask = ev['superpop'] == spop
            axes[i].scatter(ev.loc[mask, pc_cols[0]], ev.loc[mask, pc_cols[1]],
                          c=SUPERPOP_COLORS[spop], label=spop, s=3, alpha=0.5)
        axes[i].set_xlabel(f"PC1 ({var_exp[0]:.1f}%)")
        axes[i].set_ylabel(f"PC2 ({var_exp[1]:.1f}%)")
        nvar = pca_results[vtype].get('nvar', '?')
        axes[i].set_title(f"{vtype.upper()} (n={nvar:,})")
        axes[i].legend(markerscale=3, fontsize=8)

    plt.tight_layout()
    plt.savefig(FIGURES / "fig1_pca_comparison.pdf", dpi=300, bbox_inches='tight')
    plt.savefig(FIGURES / "fig1_pca_comparison.png", dpi=300, bbox_inches='tight')
    plt.close()
    print("  Saved fig1_pca_comparison")

def plot_fst_comparison(fst_results, sample_info):
    """Figure 2: FST heatmaps and cross-type correlations."""
    types = [t for t in ['snv', 'indel', 'sv'] if t in fst_results and len(fst_results[t]) > 0]
    if len(types) < 2:
        print("  Not enough FST data for comparison")
        return

    pops = sorted(sample_info['pop'].unique())
    pop_to_spop = dict(zip(sample_info['pop'], sample_info['superpop']))

    # Sort pops by superpopulation
    pops = sorted(pops, key=lambda p: (SUPERPOP_ORDER.index(pop_to_spop.get(p, 'AMR')), p))

    fig, axes = plt.subplots(1, len(types), figsize=(6*len(types), 5))
    if len(types) == 1:
        axes = [axes]

    for i, vtype in enumerate(types):
        mat = pd.DataFrame(0.0, index=pops, columns=pops)
        for key, val in fst_results[vtype].items():
            parts = key.split('_')
            if len(parts) == 2:
                p1, p2 = parts
                if p1 in pops and p2 in pops:
                    mat.loc[p1, p2] = val
                    mat.loc[p2, p1] = val

        mask = np.triu(np.ones_like(mat, dtype=bool), k=1)
        sns.heatmap(mat, mask=~mask, ax=axes[i], cmap='YlOrRd', vmin=0, vmax=0.3,
                    xticklabels=True, yticklabels=True, cbar_kws={'shrink': 0.8})
        axes[i].set_title(f"FST — {vtype.upper()}")
        axes[i].tick_params(labelsize=6)

    plt.tight_layout()
    plt.savefig(FIGURES / "fig2_fst_heatmaps.pdf", dpi=300, bbox_inches='tight')
    plt.savefig(FIGURES / "fig2_fst_heatmaps.png", dpi=300, bbox_inches='tight')
    plt.close()

    # Cross-type scatter
    pairs = list(itertools.combinations(types, 2))
    fig, axes = plt.subplots(1, len(pairs), figsize=(5*len(pairs), 5))
    if len(pairs) == 1:
        axes = [axes]

    for i, (t1, t2) in enumerate(pairs):
        common_keys = set(fst_results[t1].keys()) & set(fst_results[t2].keys())
        if len(common_keys) < 3:
            continue
        x = [fst_results[t1][k] for k in common_keys]
        y = [fst_results[t2][k] for k in common_keys]
        r, p = stats.pearsonr(x, y)
        rho, _ = stats.spearmanr(x, y)
        axes[i].scatter(x, y, s=10, alpha=0.6, color='steelblue')
        axes[i].plot([0, max(x)], [0, max(x)], 'k--', alpha=0.3)
        axes[i].set_xlabel(f"FST ({t1.upper()})")
        axes[i].set_ylabel(f"FST ({t2.upper()})")
        axes[i].set_title(f"r={r:.3f}, ρ={rho:.3f}")
        axes[i].text(0.05, 0.95, f"n={len(common_keys)} pairs\np={p:.2e}",
                     transform=axes[i].transAxes, va='top', fontsize=9)

    plt.tight_layout()
    plt.savefig(FIGURES / "fig3_fst_cross_type.pdf", dpi=300, bbox_inches='tight')
    plt.savefig(FIGURES / "fig3_fst_cross_type.png", dpi=300, bbox_inches='tight')
    plt.close()
    print("  Saved fig2_fst_heatmaps and fig3_fst_cross_type")

def plot_sfs_comparison(sfs_results):
    """Figure 4: SFS comparison across variant types and super-populations."""
    types = [t for t in ['snv', 'indel', 'sv'] if t in sfs_results]

    fig, axes = plt.subplots(1, len(SUPERPOP_ORDER), figsize=(4*len(SUPERPOP_ORDER), 4))
    type_colors = {'snv': '#377EB8', 'indel': '#E41A1C', 'sv': '#4DAF4A'}
    type_labels = {'snv': 'SNV', 'indel': 'INDEL', 'sv': 'SV'}

    for j, spop in enumerate(SUPERPOP_ORDER):
        for vtype in types:
            if spop in sfs_results[vtype]:
                data = sfs_results[vtype][spop]
                bins = data['bins']
                hist = data['hist'] / data['hist'].sum()
                centers = (bins[:-1] + bins[1:]) / 2
                axes[j].plot(centers, hist, '-o', markersize=3, color=type_colors[vtype],
                           label=type_labels[vtype], alpha=0.8)
        axes[j].set_xlabel("Minor Allele Frequency")
        axes[j].set_ylabel("Proportion")
        axes[j].set_title(spop)
        axes[j].legend(fontsize=8)
        axes[j].set_yscale('log')

    plt.tight_layout()
    plt.savefig(FIGURES / "fig4_sfs_comparison.pdf", dpi=300, bbox_inches='tight')
    plt.savefig(FIGURES / "fig4_sfs_comparison.png", dpi=300, bbox_inches='tight')
    plt.close()
    print("  Saved fig4_sfs_comparison")

def plot_concordance_summary(proc_results, fst_results):
    """Figure 5: Concordance metrics summary."""
    fig, axes = plt.subplots(1, 2, figsize=(10, 5))

    # Procrustes
    if proc_results:
        labels = list(proc_results.keys())
        values = [proc_results[k]['correlation'] for k in labels]
        colors = ['steelblue'] * len(labels)
        axes[0].barh(labels, values, color=colors)
        axes[0].set_xlabel("Procrustes Correlation (1 - disparity)")
        axes[0].set_title("PCA Space Concordance")
        axes[0].axvline(0.9, color='red', linestyle='--', alpha=0.5, label='Threshold')
        axes[0].set_xlim(0, 1)
        axes[0].legend()

    # FST correlation
    types = [t for t in ['snv', 'indel', 'sv'] if t in fst_results and len(fst_results[t]) > 0]
    pairs = list(itertools.combinations(types, 2))
    fst_corrs = {}
    for t1, t2 in pairs:
        common = set(fst_results[t1].keys()) & set(fst_results[t2].keys())
        if len(common) < 3:
            continue
        x = [fst_results[t1][k] for k in common]
        y = [fst_results[t2][k] for k in common]
        r, _ = stats.pearsonr(x, y)
        fst_corrs[f"{t1}_vs_{t2}"] = r

    if fst_corrs:
        labels = list(fst_corrs.keys())
        values = list(fst_corrs.values())
        axes[1].barh(labels, values, color='coral')
        axes[1].set_xlabel("Pearson Correlation")
        axes[1].set_title("FST Concordance")
        axes[1].axvline(0.85, color='red', linestyle='--', alpha=0.5, label='Threshold')
        axes[1].set_xlim(0, 1)
        axes[1].legend()

    plt.tight_layout()
    plt.savefig(FIGURES / "fig5_concordance_summary.pdf", dpi=300, bbox_inches='tight')
    plt.savefig(FIGURES / "fig5_concordance_summary.png", dpi=300, bbox_inches='tight')
    plt.close()
    print("  Saved fig5_concordance_summary")

# ===== MAIN =====
def main():
    print("=" * 60)
    print("Multi-Variant-Type Population Genomics Analysis")
    print("=" * 60)

    sample_info = load_sample_info()
    print(f"\nSamples: {len(sample_info)}")
    print(f"Populations: {sample_info['pop'].nunique()}")
    print(f"Super-populations: {dict(sample_info['superpop'].value_counts())}")

    print("\n=== STEP 1: Extract variant types ===")
    counts = extract_variants()

    print("\n=== STEP 2: Merge and prepare for PCA ===")
    merge_for_pca()

    print("\n=== STEP 3: Run PCA ===")
    pca_results = run_pca()

    print("\n=== STEP 4: Compute FST ===")
    fst_results = compute_fst_python(sample_info)

    print("\n=== STEP 5: Compute SFS ===")
    sfs_results = compute_sfs(sample_info)

    print("\n=== STEP 6: Procrustes Analysis ===")
    proc_results = procrustes_analysis(pca_results, sample_info)

    print("\n=== STEP 7: Generate Figures ===")
    plot_pca_comparison(pca_results, sample_info)
    plot_fst_comparison(fst_results, sample_info)
    plot_sfs_comparison(sfs_results)
    plot_concordance_summary(proc_results, fst_results)

    # Save results summary
    summary = {
        'variant_counts': counts,
        'pca_variants': {t: pca_results[t].get('nvar', 0) for t in pca_results},
        'procrustes': proc_results,
        'fst_pairs': {t: len(v) for t, v in fst_results.items()},
        'sfs_stats': {}
    }
    for vtype in sfs_results:
        for spop in sfs_results[vtype]:
            key = f"{vtype}_{spop}"
            summary['sfs_stats'][key] = {
                k: float(v) if isinstance(v, (float, np.floating)) else v
                for k, v in sfs_results[vtype][spop].items()
                if k != 'hist' and k != 'bins'
            }

    import json
    with open(RESULTS / "analysis_summary.json", 'w') as f:
        json.dump(summary, f, indent=2, default=str)

    print("\n=== ANALYSIS COMPLETE ===")
    print(f"Results: {RESULTS}")
    print(f"Figures: {FIGURES}")
    return summary

if __name__ == '__main__':
    main()
