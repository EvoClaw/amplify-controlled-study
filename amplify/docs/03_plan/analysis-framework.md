# Analysis Framework — Type D (Revised after Multi-Agent Deliberation)

## Main Story Line
**"Different variant types are independent witnesses to population history. Their concordance reveals shared demography; their discordance reveals type-specific evolutionary forces."**

## Execution Pipeline

### Phase A: Data Preparation & Quality Control (Pilot: chr22, then full genome)

**A1. Variant extraction and type-specific filtering**
- Extract SNVs, INDELs, SVs separately from VCF files using bcftools
- Remove related individuals (keep unrelated set using PED)
- Type-specific quality filters:
  - SNVs: FILTER=PASS, biallelic, MAF≥0.01, missingness<5%
  - INDELs: FILTER=PASS, biallelic, MAF≥0.01, missingness<5%
  - SVs: FILTER=PASS, biallelic, MAF≥0.005, missingness<10%
- SV encoding: standard diploid genotype (0/1/2 allele count)
- Handle multiallelic SVs: decompose or filter to biallelic
- Report post-filter counts per type, per chromosome

**A2. LD pruning (type-specific parameters)**
- SNVs: --indep-pairwise 50 5 0.2
- INDELs: --indep-pairwise 50 5 0.2
- SVs: skip LD pruning (already sparse) or use relaxed threshold (r²<0.5, window 500kb)
- Document and justify all parameters

**A3. Count-matched subsampling**
- Random subsample of SNVs to match INDEL count (N~500K after LD pruning)
- Random subsample of SNVs to match SV count (N~50-80K after filtering)
- Seeds: 42, 123, 456 (3 replicates for robustness)
- Run analyses on BOTH full and subsampled sets

**A4. Quality checks**
- Verify no sequencing batch effects (PCA colored by sequencing batch if info available)
- Check SV genotyping quality distribution
- Confirm population assignments match metadata
- **Table 1: Dataset summary (variants per type, post-filter, per-population)**

### Phase B: Core Population Structure Analysis

**B1. PCA per variant type**
- scikit-allel: Patterson PCA on genotype arrays
- Extract PC1-PC20 per variant type
- Run on: (a) full variant sets, (b) count-matched subsampled sets
- Color by super-population (AFR, AMR, EAS, EUR, SAS) and by population
- **Figure 1: Side-by-side PCA plots (3 panels) + Procrustes overlay**

**B2. Procrustes analysis for PCA concordance**
- scipy.spatial.procrustes on PC1-PC10 spaces
- Compute Procrustes dissimilarity (d²) for all 3 pairwise comparisons
- Permutation test (1000 permutations) for significance
- Bootstrap CIs (1000 replicates) for Procrustes dissimilarity
- **Effect size threshold: Procrustes correlation <0.90 = meaningful discordance**

**B3. ADMIXTURE analysis per variant type**
- plink2 → PLINK BED format per variant type
- Run ADMIXTURE for K=2 to K=8 (or use sNMF/LEA in R)
- Compare ancestry proportions across variant types
- Quantify admixture concordance per population
- **Figure 2: ADMIXTURE bar plots (K=5) side-by-side for 3 variant types**

### Phase C: Differentiation Analysis

**C1. Pairwise FST**
- Weir-Cockerham FST for all 325 population pairs × 3 variant types
- Report mean genome-wide FST per pair per type
- FDR correction (Benjamini-Hochberg) for all pairwise comparisons
- Bootstrap CIs for FST estimates (1000 replicates, block jackknife)
- **Figure 3: FST heatmaps (3 panels) + scatter plots of cross-type FST**
- **Table 2: Concordance metrics (Procrustes d², FST Pearson/Spearman r, RF distance)**

**C2. FST cross-type correlation**
- Pearson and Spearman correlation of FST values across types
- Permutation-based p-values (10,000 permutations)
- Identify population pairs with maximal discordance (|FST_SNV - FST_SV| > threshold)
- **Effect size threshold: FST correlation <0.85 = meaningful discordance**

**C3. Population trees**
- NJ trees from FST distance matrices per variant type
- Robinson-Foulds distance for topology comparison
- Bootstrap support (1000 replicates)
- **Figure 4: Side-by-side population trees (3 panels)**

### Phase D: Allele Frequency Spectrum Comparison

**D1. Folded SFS computation**
- Compute folded SFS per variant type per super-population (5 groups)
- Normalize by total variant count for cross-type comparison
- Compute summary statistics: proportion rare (MAF<0.05), mean, variance, skewness
- **Figure 5: SFS overlay plots per super-population (5 panels)**

**D2. SFS shape comparison**
- KS tests for SFS shape differences (variant type pairs × super-populations)
- FDR correction across all 30 KS tests (3 pairs × 5 super-populations × 2 sets)
- Moment-based comparison (mean MAF, variance, skewness per type)
- **Figure 6: Rare variant proportion heatmap (super-population × variant type)**
- **Table 3: SFS shape statistics per super-population × variant type**

### Phase E: Variant Property Analysis

**E1. Per-variant FST distributions**
- Compute per-variant FST for all variants per type
- Compare FST distributions across types (violin plots)
- KS test for distribution differences
- **Figure 7: FST distribution violin plots per variant type**

**E2. Variant size effects (INDELs and SVs)**
- Partition INDELs by size: 1-5bp, 5-20bp, 20-50bp
- Partition SVs by size: <1kb, 1-10kb, >10kb
- Partition SVs by type: DEL, DUP, INS, INV
- Compare FST distributions and PCA concordance per partition
- Regression: FST ~ log(variant_size) + variant_type + population_pair
- **Figure 8: Size-dependent FST patterns**

**E3. Per-chromosome concordance**
- Compute Procrustes dissimilarity and FST correlation per chromosome (chr1-22)
- Identify chromosomes with maximum cross-type discordance
- Exclude chrX from main analysis (analyze separately)
- FDR correction across 22 chromosomes
- **Table 4: Per-chromosome concordance metrics**

### Phase F: Mechanism Investigation (conditional on discordance)

**F1. Functional enrichment of high-FST variants**
- For top 1% FST variants per type: annotate functional categories
- Categories: coding, intronic, intergenic, regulatory (if annotation available)
- Compare enrichment patterns across variant types
- Fisher's exact test with FDR correction

**F2. Selection signal overlap**
- Identify FST outlier loci (top 1%) per variant type
- Window-based approach (100kb windows) for SV/INDEL overlap with SNV outliers
- Compute Jaccard index of selected regions across types
- **Figure 9: Selection signal overlap (Venn or upset plot)**

### Phase G: Robustness & Sensitivity

**G1. Filter sensitivity**
- Repeat core analyses with strict (QUAL>30) and relaxed (QUAL>10) SV filters
- Compare results across filter thresholds

**G2. Count matching sensitivity**
- Compare full vs. count-matched results (3 random subsamples)
- Report variance across subsamples

**G3. Population size sensitivity**
- Flag analyses where population n<80
- Bootstrap CIs to assess stability for small populations

**G4. Reproducibility protocol**
- All random seeds documented: 42, 123, 456
- Software versions logged
- Complete pipeline as Python scripts

## Statistical Framework

**Multiple comparisons:**
- FDR (Benjamini-Hochberg, α=0.05) for all hypothesis tests
- Permutation-based p-values for Procrustes, Mantel, and correlation tests
- Block jackknife for FST CIs (blocks = chromosomes)

**Effect size thresholds:**
- Procrustes correlation <0.90: meaningful PCA discordance
- FST cross-type r <0.85: meaningful FST discordance
- Robinson-Foulds distance >20% of maximum: meaningful tree discordance

**Power considerations:**
- SVs (~50-80K post-filter): sufficient for population-level PCA (variants >> samples)
- Small populations (n~60-80): report with wider CIs, flag as lower power
- SV-specific FST: higher variance expected; address via bootstrap and jackknife

## Planned Outputs
- **Figures:** 9 main + supplementary
- **Tables:** 4 main + supplementary (full FST matrices, SFS statistics)
- **Content points:** 9 (PCA concordance, admixture, FST, trees, SFS, variant properties, size effects, chromosomes, robustness)
