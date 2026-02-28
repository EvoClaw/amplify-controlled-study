#!/usr/bin/env python3
"""Set up population metadata and create sample lists per population and super-population."""

import pandas as pd
import os

DATA_DIR = "/home/yanlin/comp/cursor/1000GP"
OUT_DIR = "/home/yanlin/comp/cursor/analysis/data"
os.makedirs(OUT_DIR, exist_ok=True)

POP_TO_SUPERPOP = {
    'GBR': 'EUR', 'FIN': 'EUR', 'IBS': 'EUR', 'TSI': 'EUR', 'CEU': 'EUR',
    'CHB': 'EAS', 'JPT': 'EAS', 'CHS': 'EAS', 'CDX': 'EAS', 'KHV': 'EAS', 'CHD': 'EAS',
    'YRI': 'AFR', 'LWK': 'AFR', 'GWD': 'AFR', 'MSL': 'AFR', 'ESN': 'AFR', 'ACB': 'AFR', 'ASW': 'AFR',
    'MXL': 'AMR', 'PUR': 'AMR', 'CLM': 'AMR', 'PEL': 'AMR',
    'GIH': 'SAS', 'PJL': 'SAS', 'BEB': 'SAS', 'STU': 'SAS', 'ITU': 'SAS',
}

samples = pd.read_csv(f"{DATA_DIR}/samples.info", sep='\t', header=None,
                       names=['sample', 'pop', 'sex'])
samples['superpop'] = samples['pop'].map(POP_TO_SUPERPOP)

# Drop samples without superpop mapping (if any)
unmapped = samples[samples['superpop'].isna()]
if len(unmapped) > 0:
    print(f"WARNING: {len(unmapped)} samples without superpop mapping:")
    print(unmapped['pop'].value_counts())
    samples = samples.dropna(subset=['superpop'])

samples.to_csv(f"{OUT_DIR}/samples_metadata.tsv", sep='\t', index=False)

# Create sample lists per super-population
for spop in sorted(samples['superpop'].unique()):
    subset = samples[samples['superpop'] == spop]
    subset['sample'].to_csv(f"{OUT_DIR}/samples_{spop}.txt", index=False, header=False)
    print(f"{spop}: {len(subset)} samples")

# Create sample lists per population
for pop in sorted(samples['pop'].unique()):
    subset = samples[samples['pop'] == pop]
    subset['sample'].to_csv(f"{OUT_DIR}/samples_{pop}.txt", index=False, header=False)

# Create a cluster file for plink2 (sample -> superpop)
clusters = samples[['sample', 'sample', 'superpop']].copy()
clusters.columns = ['FID', 'IID', 'CLUSTER']
clusters.to_csv(f"{OUT_DIR}/superpop_clusters.txt", sep='\t', index=False, header=False)

# Population-level cluster file
pop_clusters = samples[['sample', 'sample', 'pop']].copy()
pop_clusters.columns = ['FID', 'IID', 'CLUSTER']
pop_clusters.to_csv(f"{OUT_DIR}/pop_clusters.txt", sep='\t', index=False, header=False)

print(f"\nTotal samples: {len(samples)}")
print(f"\nSuper-population counts:")
print(samples['superpop'].value_counts().sort_index())
print(f"\nPopulation counts:")
print(samples['pop'].value_counts().sort_index())

# Save population -> superpop mapping
pd.DataFrame(list(POP_TO_SUPERPOP.items()), columns=['pop', 'superpop']).to_csv(
    f"{OUT_DIR}/pop_superpop_map.tsv", sep='\t', index=False)
