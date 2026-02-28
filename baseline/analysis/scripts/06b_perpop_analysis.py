#!/usr/bin/env python3
"""
Per-superpopulation indel spectrum analysis using available chromosome data.
Computes spectrum, SFS, and population-specific patterns.
"""

import pandas as pd
import numpy as np
import gzip
import os
from collections import defaultdict

DATA_DIR = "/home/yanlin/comp/cursor/analysis/data"
RESULTS_DIR = "/home/yanlin/comp/cursor/analysis/results"
PERPOP_DIR = f"{DATA_DIR}/indels_perpop"

SUPERPOPS = ['AFR', 'AMR', 'EAS', 'EUR', 'SAS']
COMPLETE_CHROMS = ['chr1']  # fully extracted for all pops

def classify_indel(ref, alt):
    ref_len, alt_len = len(ref), len(alt)
    if alt_len > ref_len:
        return "INS", alt_len - ref_len, alt[ref_len:]
    else:
        return "DEL", ref_len - alt_len, ref[alt_len:]

def get_size_bin(size):
    if size == 1: return "1bp"
    elif size == 2: return "2bp"
    elif size == 3: return "3bp"
    elif size == 4: return "4bp"
    elif size <= 5: return "5bp"
    elif size <= 10: return "6-10bp"
    elif size <= 20: return "11-20bp"
    elif size <= 50: return "21-50bp"
    else: return ">50bp"

def get_base_comp(seq):
    seq = seq.upper()
    if len(seq) == 1:
        return 'AT' if seq in ('A', 'T') else 'GC'
    at_frac = sum(1 for b in seq if b in 'AT') / len(seq) if seq else 0
    if at_frac > 0.7: return 'AT-rich'
    elif at_frac < 0.3: return 'GC-rich'
    return 'mixed'

# Load per-population data for complete chromosomes
print("Loading per-population indel data...")

# Per-variant data: for each variant, store per-pop AC and AN
variant_data = {}  # key: (chrom, pos, ref, alt) -> {spop: (ac, an)}

for chrom in COMPLETE_CHROMS:
    for spop in SUPERPOPS:
        filepath = f"{PERPOP_DIR}/{chrom}_{spop}.tsv.gz"
        with gzip.open(filepath, 'rt') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 6:
                    continue
                chr_, pos, ref_a, alt_a, ac, an = parts[:6]
                ac, an = int(ac), int(an)
                key = (chr_, int(pos), ref_a, alt_a)
                if key not in variant_data:
                    variant_data[key] = {}
                variant_data[key][spop] = (ac, an)
    print(f"  Loaded {chrom}")

print(f"Total variants with per-pop data: {len(variant_data):,}")

# Build per-population spectrum matrix
pop_spectrum = defaultdict(lambda: defaultdict(int))
pop_sfs = defaultdict(lambda: defaultdict(int))  # (spop, af_bin) -> count
pop_del_ins_ratio = defaultdict(lambda: {'DEL': 0, 'INS': 0})
pop_size_dist = defaultdict(lambda: defaultdict(int))  # (spop, size) -> count

# Compute per-population private variant counts
pop_private = defaultdict(int)  # variants only segregating in one pop
pop_shared = defaultdict(lambda: defaultdict(int))  # shared between pops

for key, pops in variant_data.items():
    chrom, pos, ref_a, alt_a = key
    indel_type, size, seq = classify_indel(ref_a, alt_a)
    size_bin = get_size_bin(size)
    base_comp = get_base_comp(seq)
    channel = f"{indel_type}_{size_bin}_{base_comp}"
    
    segregating_pops = []
    for spop in SUPERPOPS:
        if spop in pops:
            ac, an = pops[spop]
            if ac > 0 and an > 0:
                af = ac / an
                pop_spectrum[spop][channel] += 1
                pop_del_ins_ratio[spop][indel_type] += 1
                pop_size_dist[spop][size] += 1
                
                # SFS bins (within-population AF)
                if af < 0.01:
                    pop_sfs[spop]['rare'] += 1
                elif af < 0.05:
                    pop_sfs[spop]['low_freq'] += 1
                elif af < 0.50:
                    pop_sfs[spop]['common'] += 1
                else:
                    pop_sfs[spop]['high_freq'] += 1
                
                segregating_pops.append(spop)
    
    # Track population-private variants
    if len(segregating_pops) == 1:
        pop_private[segregating_pops[0]] += 1

print("\n=== SEGREGATING INDEL COUNTS (chr1+chr2) ===")
for spop in SUPERPOPS:
    total = sum(pop_spectrum[spop].values())
    print(f"  {spop}: {total:,}")

print("\n=== DEL/INS RATIO ===")
for spop in SUPERPOPS:
    d = pop_del_ins_ratio[spop]['DEL']
    i = pop_del_ins_ratio[spop]['INS']
    print(f"  {spop}: DEL={d:,} INS={i:,} ratio={d/i:.3f}")

print("\n=== POPULATION-PRIVATE VARIANTS ===")
for spop in SUPERPOPS:
    total = sum(pop_spectrum[spop].values())
    priv = pop_private[spop]
    print(f"  {spop}: {priv:,} ({100*priv/total:.1f}%)")

print("\n=== SFS DISTRIBUTION ===")
sfs_df = pd.DataFrame({spop: pop_sfs[spop] for spop in SUPERPOPS})
sfs_norm = sfs_df.div(sfs_df.sum(axis=0), axis=1)
print(sfs_norm.round(4))

# Build spectrum matrix
all_channels = sorted(set(ch for sp in pop_spectrum.values() for ch in sp.keys()))
spectrum_df = pd.DataFrame(0, index=SUPERPOPS, columns=all_channels)
for spop in SUPERPOPS:
    for ch, count in pop_spectrum[spop].items():
        spectrum_df.loc[spop, ch] = count

spectrum_norm = spectrum_df.div(spectrum_df.sum(axis=1), axis=0)

# Compute population differentiation (pairwise channel proportion differences)
print("\n=== TOP DIFFERENTIATED CHANNELS ===")
max_diff = {}
for ch in all_channels:
    props = spectrum_norm[ch]
    max_diff[ch] = props.max() - props.min()

top_diff = sorted(max_diff.items(), key=lambda x: -x[1])[:15]
for ch, diff in top_diff:
    props = spectrum_norm[ch]
    detail = " | ".join(f"{sp}={props[sp]:.4f}" for sp in SUPERPOPS)
    print(f"  {ch}: max_diff={diff:.4f} ({detail})")

# Size distribution per population
size_df = pd.DataFrame(0, index=SUPERPOPS, columns=range(1, 52))
for spop in SUPERPOPS:
    for size, count in pop_size_dist[spop].items():
        if size <= 51:
            size_df.loc[spop, size] = count
        else:
            size_df.loc[spop, 51] = size_df.loc[spop, 51] + count

# Normalize
size_norm = size_df.div(size_df.sum(axis=1), axis=0)

# Save results
spectrum_df.to_csv(f"{RESULTS_DIR}/perpop_spectrum_chr1_2_raw.csv")
spectrum_norm.to_csv(f"{RESULTS_DIR}/perpop_spectrum_chr1_2_norm.csv")
size_norm.to_csv(f"{RESULTS_DIR}/perpop_size_distribution_norm.csv")
sfs_df.to_csv(f"{RESULTS_DIR}/perpop_sfs.csv")

# Compute FST-like measure for each channel
# Using Weir-Cockerham-like approach: variance of proportions / mean proportion
print("\n=== CHANNEL DIFFERENTIATION (FST-like) ===")
channel_fst = {}
for ch in all_channels:
    props = spectrum_norm[ch].values
    p_mean = props.mean()
    if p_mean > 0:
        channel_fst[ch] = np.var(props) / (p_mean * (1 - p_mean)) if p_mean < 1 else 0
    else:
        channel_fst[ch] = 0

top_fst = sorted(channel_fst.items(), key=lambda x: -x[1])[:15]
for ch, fst in top_fst:
    print(f"  {ch}: FST-like = {fst:.6f}")

pd.DataFrame(list(channel_fst.items()), columns=['channel', 'fst_like']).to_csv(
    f"{RESULTS_DIR}/channel_fst.csv", index=False)

print("\nDone! Results saved.")
