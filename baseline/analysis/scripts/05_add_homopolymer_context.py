#!/usr/bin/env python3
"""
Add homopolymer context to 1bp indels using the reference genome.
Optimized: loads full chromosome sequence and does vectorized lookups.
"""

import pandas as pd
import numpy as np
import pysam
import os
from collections import Counter

REF_GENOME = "/home/yanlin/public/referenceGenomes/hg38.fa.gz"
RESULTS_DIR = "/home/yanlin/comp/cursor/analysis/results"

print("Loading classified indels...")
df = pd.read_parquet(f"{RESULTS_DIR}/classified_indels.parquet")
print(f"Total indels: {len(df):,}")

ref = pysam.FastaFile(REF_GENOME)

def get_homopolymer_lengths_batch(chrom, positions, bases, ref_seq):
    """Get homopolymer lengths for a batch of positions efficiently."""
    results = np.ones(len(positions), dtype=np.int32)
    seq_len = len(ref_seq)
    
    for i in range(len(positions)):
        pos = positions[i] - 1  # 0-based
        base = bases[i].upper()
        count = 1
        
        # Look upstream
        j = pos - 1
        while j >= 0 and ref_seq[j] == base:
            count += 1
            j -= 1
        
        # Look downstream  
        j = pos + 1
        while j < seq_len and ref_seq[j] == base:
            count += 1
            j += 1
        
        results[i] = count
    
    return results

mask_1bp = df['size'] == 1
df_1bp_idx = df.index[mask_1bp]
print(f"1bp indels to process: {len(df_1bp_idx):,}")

hp_lengths = np.ones(len(df), dtype=np.int32)

for chrom_name in sorted(df['chrom'].unique()):
    chrom_mask = (df['chrom'] == chrom_name) & mask_1bp
    if chrom_mask.sum() == 0:
        continue
    
    print(f"  Loading {chrom_name}...")
    chrom_seq = ref.fetch(chrom_name).upper()
    
    subset = df[chrom_mask]
    positions = subset['pos'].values
    bases = [s[0] if s else 'N' for s in subset['seq'].values]
    
    print(f"    Processing {len(positions):,} 1bp indels...")
    hp_lens = get_homopolymer_lengths_batch(chrom_name, positions, bases, chrom_seq)
    hp_lengths[chrom_mask.values] = hp_lens
    
    del chrom_seq

df['homopolymer_len'] = hp_lengths

# Only meaningful for 1bp indels
df.loc[~mask_1bp, 'homopolymer_len'] = 0

# Create enhanced channel classification
def make_enhanced_channel(row_type, row_size, row_seq, row_hp_len, row_base_comp, row_size_bin):
    if row_size == 1:
        base = row_seq[0].upper() if row_seq else 'N'
        base_class = 'T' if base in ('A', 'T') else 'C'
        hp = min(row_hp_len, 6)
        hp_str = '6+' if row_hp_len >= 6 else str(row_hp_len)
        return f"{row_type}_1bp_{base_class}_hp{hp_str}"
    elif row_size <= 5:
        return f"{row_type}_{row_size_bin}_{row_base_comp}"
    else:
        return f"{row_type}_{row_size_bin}_{row_base_comp}"

print("\nCreating enhanced channels...")
df['enhanced_channel'] = [
    make_enhanced_channel(t, s, sq, hp, bc, sb)
    for t, s, sq, hp, bc, sb in zip(
        df['type'], df['size'], df['seq'], 
        df['homopolymer_len'], df['base_comp'], df['size_bin']
    )
]

print("\n=== ENHANCED CHANNEL DISTRIBUTION (1bp) ===")
mask = df['size'] == 1
print(df.loc[mask, 'enhanced_channel'].value_counts().head(20))

print(f"\n=== HOMOPOLYMER LENGTH STATS (1bp indels) ===")
hp_1bp = df.loc[mask, 'homopolymer_len']
print(f"  Mean: {hp_1bp.mean():.2f}")
print(f"  Median: {hp_1bp.median():.0f}")
print(f"  Max: {hp_1bp.max()}")
print(f"  Distribution:")
for length in range(1, 8):
    label = f'{length}' if length < 6 else ('6+' if length == 6 else '')
    if length <= 6:
        if length < 6:
            n = (hp_1bp == length).sum()
        else:
            n = (hp_1bp >= 6).sum()
        print(f"    hp={label}: {n:,} ({100*n/len(hp_1bp):.1f}%)")

# Total enhanced channels
all_channels = sorted(df['enhanced_channel'].unique())
print(f"\nTotal enhanced channels: {len(all_channels)}")

# Save
df.to_parquet(f"{RESULTS_DIR}/classified_indels_enhanced.parquet", index=False)
print(f"Saved to {RESULTS_DIR}/classified_indels_enhanced.parquet")

pd.DataFrame({'channel': all_channels}).to_csv(
    f"{RESULTS_DIR}/enhanced_channels.csv", index=False)

ref.close()
