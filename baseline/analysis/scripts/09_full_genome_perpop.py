#!/usr/bin/env python3
"""
Full genome per-population analysis using all available extracted data.
Updates spectrum, SFS, and differentiation statistics.
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

# Find all complete chromosome sets (all 5 superpops available)
print("Scanning available per-population files...")
complete_chroms = []
for chrom_num in range(1, 23):
    chrom = f"chr{chrom_num}"
    all_present = True
    for spop in SUPERPOPS:
        f = f"{PERPOP_DIR}/{chrom}_{spop}.tsv.gz"
        if not os.path.exists(f) or os.path.getsize(f) == 0:
            all_present = False
            break
        # Check if file is complete (not truncated)
        try:
            with gzip.open(f, 'rt') as fh:
                last_line = None
                for last_line in fh:
                    pass
                if last_line is None:
                    all_present = False
                    break
        except:
            all_present = False
            break
    
    if all_present:
        complete_chroms.append(chrom)

print(f"Complete chromosomes: {len(complete_chroms)} - {', '.join(complete_chroms)}")

if len(complete_chroms) == 0:
    print("No complete chromosome data available. Using chr1 only.")
    complete_chroms = ['chr1']

# Process all complete chromosomes
pop_spectrum = defaultdict(lambda: defaultdict(int))
pop_sfs = defaultdict(lambda: defaultdict(int))
pop_del_ins = defaultdict(lambda: {'DEL': 0, 'INS': 0})
pop_totals = defaultdict(int)
pop_private = defaultdict(int)

for chrom in complete_chroms:
    print(f"  Processing {chrom}...")
    
    chrom_data = {}  # (pos, ref, alt) -> {spop: (ac, an)}
    
    for spop in SUPERPOPS:
        filepath = f"{PERPOP_DIR}/{chrom}_{spop}.tsv.gz"
        with gzip.open(filepath, 'rt') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 6:
                    continue
                chr_, pos, ref_a, alt_a, ac, an = parts[:6]
                ac, an = int(ac), int(an)
                key = (pos, ref_a, alt_a)
                if key not in chrom_data:
                    chrom_data[key] = {}
                chrom_data[key][spop] = (ac, an)
    
    for key, pops in chrom_data.items():
        pos, ref_a, alt_a = key
        indel_type, size, seq = classify_indel(ref_a, alt_a)
        size_bin = get_size_bin(size)
        base_comp = get_base_comp(seq)
        channel = f"{indel_type}_{size_bin}_{base_comp}"
        
        segregating = []
        for spop in SUPERPOPS:
            if spop in pops:
                ac, an = pops[spop]
                if ac > 0 and an > 0:
                    af = ac / an
                    pop_spectrum[spop][channel] += 1
                    pop_del_ins[spop][indel_type] += 1
                    pop_totals[spop] += 1
                    
                    if af < 0.01:
                        pop_sfs[spop]['rare'] += 1
                    elif af < 0.05:
                        pop_sfs[spop]['low_freq'] += 1
                    elif af < 0.50:
                        pop_sfs[spop]['common'] += 1
                    else:
                        pop_sfs[spop]['high_freq'] += 1
                    
                    segregating.append(spop)
        
        if len(segregating) == 1:
            pop_private[segregating[0]] += 1

print(f"\n=== FULL GENOME PER-POP ANALYSIS ({len(complete_chroms)} chromosomes) ===")

print("\nSegregating indel counts:")
for spop in SUPERPOPS:
    print(f"  {spop}: {pop_totals[spop]:,}")

print("\nDEL/INS ratios:")
for spop in SUPERPOPS:
    d, i = pop_del_ins[spop]['DEL'], pop_del_ins[spop]['INS']
    print(f"  {spop}: {d/i:.4f} (DEL={d:,}, INS={i:,})")

print("\nPopulation-private variants:")
for spop in SUPERPOPS:
    priv = pop_private[spop]
    total = pop_totals[spop]
    print(f"  {spop}: {priv:,} ({100*priv/total:.1f}%)")

# Build and save updated spectrum
all_channels = sorted(set(ch for sp in pop_spectrum.values() for ch in sp.keys()))
spectrum_df = pd.DataFrame(0, index=SUPERPOPS, columns=all_channels)
for spop in SUPERPOPS:
    for ch, count in pop_spectrum[spop].items():
        spectrum_df.loc[spop, ch] = count

spectrum_norm = spectrum_df.div(spectrum_df.sum(axis=1), axis=0)

spectrum_df.to_csv(f"{RESULTS_DIR}/perpop_spectrum_full_raw.csv")
spectrum_norm.to_csv(f"{RESULTS_DIR}/perpop_spectrum_full_norm.csv")

sfs_df = pd.DataFrame({spop: pop_sfs[spop] for spop in SUPERPOPS})
sfs_df.to_csv(f"{RESULTS_DIR}/perpop_sfs_full.csv")

print(f"\nSaved updated results for {len(complete_chroms)} chromosomes.")
print("\nDone!")
