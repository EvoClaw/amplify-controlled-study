# Literature Review — Deep Reading Table

## Core Dataset

| Paper | Access | Problem | Method | Datasets | Key Results | Limitations |
|-------|--------|---------|--------|----------|-------------|-------------|
| Byrska-Bishop 2022 | full_text | Expand 1kGP to high-coverage | 30x Illumina WGS + phasing | 3,202 samples, 26 pops | 73.5M variants (64M SNV, 9.5M INDEL, 102K SV); improved sensitivity vs Phase 3 | SV calling from short reads limited; no systematic population analysis of SVs |
| Sirén et al. 2025 | abstract | Long-read SV characterization | ONT + linear & graph pangenome | 1,019 of 1kGP samples | 167K SVs; population stratification; L1/SVA mediated transductions | Only ~1/3 of samples; long-read only; doesn't compare with short-read SVs |

## Multi-Variant-Type Population Studies

| Paper | Access | Problem | Method | Datasets | Key Results | Limitations |
|-------|--------|---------|--------|----------|-------------|-------------|
| Hansson et al. 2024 | full_text | SV/INDEL/SNP differentiation in salmon | FST analysis per variant type | Atlantic salmon populations | SVs show different differentiation patterns than SNPs; some SVs under divergent selection | Non-human; limited populations; no systematic SFS comparison |
| Galtier et al. 2024 | full_text | SV complexity in pangenomes | Comparative pangenomics | Multiple species | SVs show unexpected complexity; fitness effects vary | Not human-focused; pangenome-based |

## Selection and Frequency Analysis

| Paper | Access | Problem | Method | Datasets | Key Results | Limitations |
|-------|--------|---------|--------|----------|-------------|-------------|
| Johri et al. 2024 | full_text | Separate selection from demography | SFRatios method | Drosophila | SFS ratios avoid demographic confounds | Not applied to variant-type comparison in humans |
| Stolyarova et al. 2025 | abstract | Deleterious variant distribution | Mutation-drift-selection models | gnomAD | Highly pathogenic variants similar across ancestries | Focus on LoF only, not multi-type |

## Mutation Spectrum

| Paper | Access | Problem | Method | Datasets | Key Results | Limitations |
|-------|--------|---------|--------|----------|-------------|-------------|
| Harris 2015 | full_text | Population-specific mutation rate | 3-mer spectrum analysis | 1kGP Phase 3 | TCC→TTC enriched in Europeans ~50% | Phase 3 low-coverage; SNVs only |
| Mathieson 2017 | full_text | Rare variant spectrum variation | SFS in trinucleotide context | 1kGP Phase 3 | CCG→CTG enriched in South Asians | Low-coverage data; no INDEL/SV comparison |

## Key Observations

1. **No paper compares population structure (PCA, FST, clustering) across variant types (SNV, INDEL, SV) in the same human dataset**
2. The Atlantic salmon paper (Hansson 2024) is the closest precedent but in a non-human species
3. The 2022 high-coverage panel provides the first opportunity to do this systematically in humans with short-read SVs
4. Mutation spectrum variation has been studied for SNVs but not extended to INDELs or SVs
5. Selection pressures differ by variant type (SVs have larger fitness effects) but comparative analysis is lacking
