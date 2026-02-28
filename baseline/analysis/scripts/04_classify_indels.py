#!/usr/bin/env python3
"""
Classify all indels genome-wide by type, size, and basic properties.
Produce a unified classification for indel spectrum analysis.
"""

import pandas as pd
import numpy as np
import gzip
import os
from collections import Counter

DATA_DIR = "/home/yanlin/comp/cursor/analysis/data"
INDEL_DIR = f"{DATA_DIR}/indels"
OUT_DIR = "/home/yanlin/comp/cursor/analysis/results"
os.makedirs(OUT_DIR, exist_ok=True)

def classify_indel(ref, alt):
    """Classify an indel by type and size."""
    ref_len = len(ref)
    alt_len = len(alt)
    
    if alt_len > ref_len:
        indel_type = "INS"
        size = alt_len - ref_len
        seq = alt[ref_len:]  # inserted sequence
    else:
        indel_type = "DEL"
        size = ref_len - alt_len
        seq = ref[alt_len:]  # deleted sequence
    
    return indel_type, size, seq

def get_size_bin(size):
    """Bin indel sizes."""
    if size == 1:
        return "1bp"
    elif size == 2:
        return "2bp"
    elif size == 3:
        return "3bp"
    elif size == 4:
        return "4bp"
    elif size <= 5:
        return "5bp"
    elif size <= 10:
        return "6-10bp"
    elif size <= 20:
        return "11-20bp"
    elif size <= 50:
        return "21-50bp"
    else:
        return ">50bp"

def get_repeat_unit(seq):
    """Check if sequence is a repeat of a shorter unit."""
    seq = seq.upper()
    for unit_len in range(1, len(seq)//2 + 1):
        unit = seq[:unit_len]
        if unit * (len(seq) // unit_len) == seq[:unit_len * (len(seq) // unit_len)] and len(seq) % unit_len == 0:
            return unit_len, unit
    return len(seq), seq

def get_base_composition(seq):
    """Get dominant base for 1bp indels, AT/GC content for longer."""
    seq = seq.upper()
    if len(seq) == 1:
        if seq in ('A', 'T'):
            return 'AT'
        else:
            return 'GC'
    at_frac = sum(1 for b in seq if b in 'AT') / len(seq) if seq else 0
    if at_frac > 0.7:
        return 'AT-rich'
    elif at_frac < 0.3:
        return 'GC-rich'
    else:
        return 'mixed'

def get_af_bin(af):
    """Bin allele frequency."""
    if af < 0.001:
        return 'singleton-like'
    elif af < 0.005:
        return 'very_rare'
    elif af < 0.01:
        return 'rare'
    elif af < 0.05:
        return 'low_freq'
    elif af < 0.5:
        return 'common'
    else:
        return 'major'

SIZE_BIN_ORDER = ["1bp", "2bp", "3bp", "4bp", "5bp", "6-10bp", "11-20bp", "21-50bp", ">50bp"]

print("Loading and classifying indels from all chromosomes...")

all_records = []
chrom_stats = {}

for chrom_num in list(range(1, 23)):
    chrom = f"chr{chrom_num}"
    filepath = f"{INDEL_DIR}/{chrom}_indels.tsv.gz"
    
    chrom_records = []
    with gzip.open(filepath, 'rt') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 7:
                continue
            chr_, pos, ref, alt, ac, an, af = parts
            pos = int(pos)
            ac = int(ac)
            an = int(an)
            af = float(af)
            
            indel_type, size, seq = classify_indel(ref, alt)
            size_bin = get_size_bin(size)
            base_comp = get_base_composition(seq)
            af_bin = get_af_bin(af)
            
            chrom_records.append({
                'chrom': chr_,
                'pos': pos,
                'ref': ref,
                'alt': alt,
                'ac': ac,
                'an': an,
                'af': af,
                'type': indel_type,
                'size': size,
                'size_bin': size_bin,
                'base_comp': base_comp,
                'af_bin': af_bin,
                'seq': seq[:20],  # truncate long seqs
            })
    
    n = len(chrom_records)
    chrom_stats[chrom] = n
    all_records.extend(chrom_records)
    print(f"  {chrom}: {n:,} indels")

print(f"\nTotal: {len(all_records):,} indels")

df = pd.DataFrame(all_records)

# Create combined classification channel
df['channel'] = df['type'] + '_' + df['size_bin'] + '_' + df['base_comp']

# Summary statistics
print("\n=== INDEL TYPE DISTRIBUTION ===")
print(df['type'].value_counts())

print("\n=== SIZE DISTRIBUTION ===")
type_size = pd.crosstab(df['type'], df['size_bin'])
type_size = type_size.reindex(columns=[s for s in SIZE_BIN_ORDER if s in type_size.columns])
print(type_size)

print("\n=== ALLELE FREQUENCY DISTRIBUTION ===")
print(df['af_bin'].value_counts())

print("\n=== CHANNEL COUNTS (top 30) ===")
print(df['channel'].value_counts().head(30))

# Save classified indels
outpath = f"{OUT_DIR}/classified_indels.parquet"
df.to_parquet(outpath, index=False)
print(f"\nSaved classified indels to {outpath}")

# Also save summary tables
type_size.to_csv(f"{OUT_DIR}/indel_type_size_crosstab.csv")

# Save channel definitions
channels = sorted(df['channel'].unique())
pd.DataFrame({'channel': channels}).to_csv(f"{OUT_DIR}/indel_channels.csv", index=False)
print(f"Defined {len(channels)} indel channels")
