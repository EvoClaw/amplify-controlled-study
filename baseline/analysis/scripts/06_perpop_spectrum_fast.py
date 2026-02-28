#!/usr/bin/env python3
"""
Fast per-population indel spectrum computation.
Uses bcftools + per-superpopulation sample subsetting.
Processes chromosomes sequentially, merging per-pop results.
"""

import subprocess
import pandas as pd
import numpy as np
import gzip
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from collections import defaultdict

DATA_DIR = "/home/yanlin/comp/cursor/analysis/data"
VCF_DIR = "/home/yanlin/comp/cursor/1000GP/20220422_3202_phased_SNV_INDEL_SV"
RESULTS_DIR = "/home/yanlin/comp/cursor/analysis/results"

SUPERPOPS = ['AFR', 'AMR', 'EAS', 'EUR', 'SAS']
CHROMS = [f"chr{i}" for i in range(1, 23)]

def classify_indel(ref, alt):
    """Classify an indel."""
    ref_len = len(ref)
    alt_len = len(alt)
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
    else: return 'mixed'

def get_af_bin(af):
    if af < 0.001: return 'singleton-like'
    elif af < 0.005: return 'very_rare'
    elif af < 0.01: return 'rare'
    elif af < 0.05: return 'low_freq'
    elif af < 0.5: return 'common'
    else: return 'major'

def process_perpop_file(filepath, spop):
    """Process a per-superpopulation indel file and return spectrum counts."""
    channel_counts = defaultdict(int)
    af_spectrum = defaultdict(lambda: defaultdict(int))
    total = 0
    
    with gzip.open(filepath, 'rt') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 6:
                continue
            chrom, pos, ref_allele, alt_allele, ac, an = parts[:6]
            ac = int(ac)
            an = int(an)
            if an == 0 or ac == 0:
                continue
            
            af = ac / an
            indel_type, size, seq = classify_indel(ref_allele, alt_allele)
            size_bin = get_size_bin(size)
            base_comp = get_base_comp(seq)
            af_bin = get_af_bin(af)
            
            channel = f"{indel_type}_{size_bin}_{base_comp}"
            channel_counts[channel] += 1
            af_spectrum[channel][af_bin] += 1
            total += 1
    
    return spop, dict(channel_counts), dict(af_spectrum), total

print("Processing per-superpopulation indel spectra...")
print("Waiting for per-pop extraction files...")

# Check which files are available
available = {}
for chrom in CHROMS:
    for spop in SUPERPOPS:
        f = f"{DATA_DIR}/indels_perpop/{chrom}_{spop}.tsv.gz"
        if os.path.exists(f) and os.path.getsize(f) > 0:
            available.setdefault(spop, []).append(f)

for spop in SUPERPOPS:
    files = available.get(spop, [])
    print(f"  {spop}: {len(files)}/22 chromosomes available")

# Process whatever is available
pop_spectra = defaultdict(lambda: defaultdict(int))
pop_af_spectra = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
pop_totals = defaultdict(int)

for spop in SUPERPOPS:
    files = available.get(spop, [])
    for f in files:
        with gzip.open(f, 'rt') as fh:
            for line in fh:
                parts = line.strip().split('\t')
                if len(parts) < 6:
                    continue
                chrom, pos, ref_allele, alt_allele, ac, an = parts[:6]
                ac = int(ac)
                an = int(an)
                if an == 0 or ac == 0:
                    continue
                
                af = ac / an
                indel_type, size, seq = classify_indel(ref_allele, alt_allele)
                size_bin = get_size_bin(size)
                base_comp = get_base_comp(seq)
                af_bin = get_af_bin(af)
                
                channel = f"{indel_type}_{size_bin}_{base_comp}"
                pop_spectra[spop][channel] += 1
                pop_af_spectra[spop][channel][af_bin] += 1
                pop_totals[spop] += 1
    
    print(f"  {spop}: {pop_totals[spop]:,} segregating indels from {len(files)} chroms")

# Build spectrum matrix
all_channels = sorted(set(ch for sp in pop_spectra.values() for ch in sp.keys()))
spectrum_matrix = pd.DataFrame(0, index=SUPERPOPS, columns=all_channels)
for spop in SUPERPOPS:
    for ch, count in pop_spectra[spop].items():
        spectrum_matrix.loc[spop, ch] = count

# Normalize to proportions
spectrum_norm = spectrum_matrix.div(spectrum_matrix.sum(axis=1), axis=0)

print("\n=== RAW SPECTRUM MATRIX (top channels) ===")
top_channels = spectrum_matrix.sum().nlargest(15).index
print(spectrum_matrix[top_channels])

print("\n=== NORMALIZED SPECTRUM (proportions, top channels) ===")
print(spectrum_norm[top_channels].round(4))

# Save
spectrum_matrix.to_csv(f"{RESULTS_DIR}/perpop_spectrum_raw.csv")
spectrum_norm.to_csv(f"{RESULTS_DIR}/perpop_spectrum_norm.csv")

# Build AF spectrum per population per channel
af_bins = ['singleton-like', 'very_rare', 'rare', 'low_freq', 'common', 'major']
af_data = []
for spop in SUPERPOPS:
    for channel in all_channels:
        for af_bin in af_bins:
            count = pop_af_spectra[spop].get(channel, {}).get(af_bin, 0)
            af_data.append({
                'superpop': spop,
                'channel': channel,
                'af_bin': af_bin,
                'count': count
            })

af_df = pd.DataFrame(af_data)
af_df.to_csv(f"{RESULTS_DIR}/perpop_af_spectrum.csv", index=False)

print(f"\nSaved per-population spectra to {RESULTS_DIR}/")
print(f"Channels: {len(all_channels)}")
