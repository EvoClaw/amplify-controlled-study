# Methods Evaluation: Population Genomics Analysis Framework

**Evaluator Role:** Methods Expert in Population Genomics  
**Date:** 2026-02-28  
**Target Venue:** PLOS Genetics  
**Research Type:** Type D (Discovery)

---

## EXECUTIVE SUMMARY

**VERDICT: CONDITIONAL PASS**

The proposed framework is methodologically sound and uses appropriate tools, but requires **critical modifications** to handle variant-type-specific challenges before execution. Key gaps: SV encoding strategy, LD pruning for sparse variants, count disparity handling, and filtering threshold justification.

---

## 1. METHOD APPROPRIATENESS

### ✅ STRENGTHS

**Data Preparation Pipeline:**
- ✅ **Type-specific VCF filtering**: Correct approach. SNVs, INDELs, and SVs have different quality metrics and require different thresholds (e.g., SV QUAL scores vs SNV GQ).
- ✅ **LD pruning per type**: Essential for fair comparison. LD structure differs dramatically between SNVs (high LD blocks) and SVs (sparse, often independent).
- ✅ **Relatedness removal using PED**: Appropriate for population-level analysis. Prevents family structure from confounding population signals.

**Population Structure Analysis:**
- ✅ **PCA using scikit-allel**: Appropriate tool. `GenotypeArray` → PCA is standard and handles missing data well.
- ✅ **Per-variant-type PCA**: Correct design choice. Enables direct comparison of structure signals.

**Differentiation Analysis:**
- ✅ **Weir-Cockerham FST**: Gold standard estimator. Handles unequal sample sizes and accounts for sampling variance.
- ✅ **Per-population-pair FST**: Standard approach for multi-population studies.

**Phylogenetic Analysis:**
- ✅ **NJ trees from FST distance matrices**: Appropriate for visualization. FST is a valid distance metric for tree construction.

**Frequency Spectrum Analysis:**
- ✅ **Folded SFS**: Correct choice for discovery analysis. Avoids phasing errors and is robust to ancestral state uncertainty.

**Concordance Metrics:**
- ✅ **Procrustes analysis**: Standard for comparing PCA spaces. Tests whether variant types capture the same population structure.
- ✅ **Mantel test**: Appropriate for distance matrix correlation. Accounts for non-independence of pairwise distances.
- ✅ **Robinson-Foulds distance**: Standard metric for tree topology comparison.

**Statistical Testing:**
- ✅ **KS tests for SFS comparison**: Appropriate for shape comparison. Non-parametric, robust to outliers.

### ⚠️ CONCERNS

**1. SV Encoding Strategy (CRITICAL)**
- **Issue**: Binary (0/1) vs. allele count (0/1/2) encoding fundamentally changes PCA interpretation.
- **Binary encoding**: Treats SVs as presence/absence markers. Appropriate for deletions/duplications where copy number is ambiguous.
- **Allele count encoding**: Assumes diploid genotyping (0/1/2). Appropriate if SV calls are phased and copy-number resolved.
- **Recommendation**: 
  - **For deletions**: Use binary encoding (0=absent, 1=present) if copy number is uncertain.
  - **For duplications/CNVs**: Use allele count if copy-number resolved; otherwise binary.
  - **For inversions/translocations**: Binary is appropriate (presence/absence).
  - **Document the encoding choice** and justify based on SV call set characteristics.
  - **Sensitivity analysis**: Run PCA with both encodings and report if results differ.

**2. LD Pruning for SVs (CRITICAL)**
- **Issue**: Standard LD pruning (e.g., `--indep-pairwise 50 5 0.2` in PLINK) assumes dense markers. SVs are sparse (~102K vs 64M SNVs).
- **Problem**: With sparse variants, LD pruning may remove too many variants or fail to detect LD properly.
- **Recommendation**:
  - **For SVs**: Use a relaxed LD threshold (r² < 0.8 or 0.9) or skip LD pruning entirely if variants are sufficiently sparse.
  - **For SNVs**: Standard LD pruning (r² < 0.2) is appropriate.
  - **For INDELs**: Intermediate approach (r² < 0.5) may be appropriate.
  - **Justify the choice** based on variant density and LD decay curves.

**3. Multiallelic SV Handling (CRITICAL)**
- **Issue**: Some SVs may be multiallelic (e.g., multiple deletion alleles at the same locus).
- **Current framework**: Unclear how multiallelic variants are handled.
- **Recommendation**:
  - **Option A**: Split multiallelic variants into biallelic (bcftools `--split-multi`), then treat each as independent.
  - **Option B**: Use multiallelic-aware PCA (scikit-allel supports this via `GenotypeArray`).
  - **Option C**: Filter to biallelic only (may lose information).
  - **Document the approach** and report how many multiallelic variants were handled.

**4. Filtering Threshold Justification (MODERATE)**
- **Issue**: Framework states "type-specific filtering" but doesn't specify thresholds.
- **Recommendation**:
  - **SNVs**: Standard thresholds (QUAL > 20, GQ > 20, DP > 10, MAF > 0.01).
  - **INDELs**: Similar to SNVs but may need stricter QUAL (INDELs are harder to call).
  - **SVs**: SV-specific metrics (e.g., SVLEN, SVTYPE, SUPPORT, PE, SR). Use thresholds from Byrska-Bishop 2022 paper or SV call set documentation.
  - **Document all thresholds** and justify based on call set quality metrics.

---

## 2. STATE-OF-THE-ART PRACTICES

### ✅ CURRENT BEST PRACTICES USED

1. **Weir-Cockerham FST**: Current gold standard (preferred over Nei's FST or Hudson's FST for sample-size correction).
2. **scikit-allel**: Modern, well-maintained Python library. Preferred over older tools (e.g., EIGENSOFT) for integration.
3. **Procrustes analysis**: Standard for comparing ordination spaces (used in microbiome, ecology, genomics).
4. **Folded SFS**: Appropriate for discovery analysis without ancestral state.

### ⚠️ MISSING STATE-OF-THE-ART ELEMENTS

**1. Admixture Analysis (MODERATE GAP)**
- **Missing**: ADMIXTURE or STRUCTURE analysis per variant type.
- **Why important**: PCA captures continuous structure; admixture quantifies discrete ancestry components. Both are needed for comprehensive structure comparison.
- **Recommendation**: Add ADMIXTURE analysis (K=2 to K=10) per variant type. Compare ancestry proportions across types.

**2. LD Decay Analysis (MINOR GAP)**
- **Missing**: Quantification of LD decay per variant type.
- **Why important**: LD structure differences could explain PCA/FST differences. SVs may have different LD patterns than SNVs.
- **Recommendation**: Calculate r² vs. distance for each variant type. Report LD decay curves.

**3. Allele Frequency Spectrum Ratios (MINOR GAP)**
- **Missing**: SFRatios method (Johri et al. 2024) to separate selection from demography.
- **Why important**: SFS shape differences could reflect selection vs. demography. SFRatios can help distinguish.
- **Recommendation**: Consider adding SFRatios analysis if SFS discordances are found.

**4. Permutation Testing for Concordance (MODERATE GAP)**
- **Missing**: Null distribution for Procrustes/Mantel statistics.
- **Why important**: Need to know if observed concordance is significantly higher than random.
- **Recommendation**: 
  - **Procrustes**: Permute sample labels, recompute PCA, calculate Procrustes correlation. Repeat 1000x to get null distribution.
  - **Mantel test**: Built-in permutation in scipy.stats.mantel (use this instead of correlation alone).

**5. Sex Chromosome Handling (MINOR GAP)**
- **Issue**: chrX included but framework doesn't specify how it's handled.
- **Recommendation**: 
  - **Option A**: Analyze chrX separately (different ploidy, different effective population size).
  - **Option B**: Exclude chrX from main analysis, analyze separately.
  - **Document the choice**.

---

## 3. METHODOLOGICAL GAPS

### CRITICAL GAPS

**1. Count Disparity Handling (102K SVs vs 64M SNVs)**
- **Problem**: PCA is sensitive to the number of variants. More variants = more variance explained, potentially biasing comparison.
- **Current approach**: Unclear how this is handled.
- **Recommendations**:
  - **Option A**: Subsample SNVs to match SV count (random sampling or LD-pruned subset). Run PCA on matched sets.
  - **Option B**: Report variance explained per variant type. If SV PCA explains similar variance with fewer variants, that's biologically interesting.
  - **Option C**: Use variance-standardized PCA (normalize by number of variants).
  - **Recommendation**: Use **Option A + Option B**. Subsample SNVs to 102K (after LD pruning) for fair comparison, but also report full SNV PCA for context.

**2. Missing Data Handling**
- **Issue**: SVs may have more missing data than SNVs (harder to call).
- **Current approach**: scikit-allel handles missing data, but strategy should be documented.
- **Recommendation**:
  - **Report missingness per variant type** (mean, distribution).
  - **Filter variants with >10% missing** (or justify threshold).
  - **Document imputation strategy** (if any). For PCA, mean imputation is often acceptable.

**3. Sample Size Per Population**
- **Issue**: 26 populations with 3,202 samples = ~123 samples/population on average, but distribution is uneven.
- **Problem**: Small populations may have unreliable FST estimates.
- **Recommendation**:
  - **Report sample sizes per population**.
  - **Filter population pairs** with <10 samples each (or justify minimum).
  - **Use Weir-Cockerham FST** (already planned) — it accounts for unequal sample sizes.

**4. Multiple Testing Correction**
- **Issue**: Multiple comparisons across variant types, population pairs, chromosomes.
- **Current approach**: Not specified.
- **Recommendation**:
  - **FST comparisons**: FDR correction across all population pairs × variant types.
  - **SFS KS tests**: FDR correction across populations × variant type pairs.
  - **Procrustes/Mantel**: Single test per comparison (less concern, but still report p-values).

### MODERATE GAPS

**5. Chromosome-Level Analysis**
- **Issue**: Framework doesn't specify if analysis is genome-wide or per-chromosome.
- **Recommendation**:
  - **Primary analysis**: Genome-wide (all autosomes).
  - **Sensitivity analysis**: Per-chromosome PCA/FST to check for chromosome-specific effects.
  - **chrX**: Analyze separately (different ploidy).

**6. Batch Effect Control**
- **Issue**: 1kGP samples were sequenced in batches. Batch effects could confound population structure.
- **Recommendation**:
  - **Check for batch effects**: PCA colored by sequencing batch. If batch explains variance, include as covariate.
  - **Document sequencing batches** from metadata.

**7. Relatedness Threshold**
- **Issue**: PED file used to remove related individuals, but threshold not specified.
- **Recommendation**:
  - **Standard approach**: Remove one individual from each pair with PI_HAT > 0.2 (second-degree relatives or closer).
  - **Report**: How many individuals removed, which relationships.

---

## 4. ROBUSTNESS & CROSS-VALIDATION

### ✅ STRENGTHS

1. **Multiple concordance metrics**: Procrustes (PCA), Mantel (FST), Robinson-Foulds (trees), KS (SFS). Good triangulation.
2. **Per-variant-type analysis**: Enables direct comparison.
3. **Standard statistical tests**: KS tests, Mantel tests are well-established.

### ⚠️ MISSING ROBUSTNESS CHECKS

**1. Bootstrap/Jackknife Resampling (CRITICAL)**
- **Missing**: Confidence intervals for FST estimates, PCA loadings, tree topologies.
- **Recommendation**:
  - **FST**: Bootstrap samples (resample with replacement) 1000x, recompute FST, report 95% CI.
  - **PCA**: Jackknife resampling (leave-one-population-out) to check stability.
  - **Trees**: Bootstrap support values (use PHYLIP or RAxML bootstrap).

**2. Sensitivity Analysis (MODERATE)**
- **Missing**: How sensitive are results to filtering thresholds, LD pruning parameters, encoding choices?
- **Recommendation**:
  - **Filtering**: Run with relaxed/strict thresholds, report if conclusions change.
  - **LD pruning**: Compare with/without LD pruning for SVs.
  - **SV encoding**: Compare binary vs. allele count (if applicable).

**3. Replication Across Chromosomes (MINOR)**
- **Missing**: Are concordance patterns consistent across chromosomes?
- **Recommendation**: Report per-chromosome concordance metrics. If some chromosomes show discordance, investigate (e.g., selection, structural variation).

**4. Alternative Distance Metrics (MINOR)**
- **Missing**: FST is one distance metric. Others exist (Nei's D, Reynolds' distance).
- **Recommendation**: Consider reporting correlation between FST and alternative metrics per variant type. If they disagree, investigate.

---

## 5. VERDICT & RECOMMENDATIONS

### VERDICT: **CONDITIONAL PASS**

The framework is methodologically sound and uses appropriate tools, but requires **critical modifications** before execution:

### REQUIRED MODIFICATIONS (Must Address)

1. **SV Encoding Strategy**: Document and justify binary vs. allele count encoding. Run sensitivity analysis.
2. **LD Pruning for SVs**: Specify relaxed thresholds or skip LD pruning. Justify choice.
3. **Count Disparity**: Subsample SNVs to match SV count for fair PCA comparison. Report both matched and full analyses.
4. **Multiallelic Handling**: Document how multiallelic SVs are handled (split vs. filter vs. multiallelic-aware).
5. **Filtering Thresholds**: Specify all thresholds per variant type. Justify based on call set quality.
6. **Permutation Testing**: Add null distributions for Procrustes/Mantel statistics.
7. **Bootstrap/Jackknife**: Add confidence intervals for FST, PCA stability checks.

### RECOMMENDED ADDITIONS (Should Address)

8. **Admixture Analysis**: Add ADMIXTURE per variant type for comprehensive structure comparison.
9. **LD Decay Analysis**: Quantify LD decay per variant type to explain potential differences.
10. **Missing Data Reporting**: Report missingness per variant type, filter high-missingness variants.
11. **Batch Effect Check**: Verify sequencing batch doesn't confound population structure.
12. **Sensitivity Analysis**: Test robustness to filtering thresholds and parameter choices.

### OPTIONAL ENHANCEMENTS (Nice to Have)

13. **SFRatios Analysis**: If SFS discordances found, use SFRatios to separate selection from demography.
14. **Per-Chromosome Analysis**: Check if concordance patterns are consistent across chromosomes.
15. **Alternative Distance Metrics**: Compare FST with Nei's D, Reynolds' distance.

---

## IMPLEMENTATION PRIORITY

**Phase 1 (Before Execution):**
- Address items 1-7 (Required Modifications)

**Phase 2 (During Execution):**
- Address items 8-12 (Recommended Additions)

**Phase 3 (If Time Permits):**
- Address items 13-15 (Optional Enhancements)

---

## FINAL ASSESSMENT

**Method Appropriateness:** ✅ PASS (with modifications)  
**State-of-the-Art:** ⚠️ CONDITIONAL (missing admixture, permutation testing)  
**Methodological Gaps:** ⚠️ CONDITIONAL (count disparity, SV encoding, LD pruning need specification)  
**Robustness:** ⚠️ CONDITIONAL (needs bootstrap/jackknife, sensitivity analysis)

**OVERALL VERDICT: CONDITIONAL PASS**

The framework is **scientifically sound** and uses **appropriate tools**, but requires **critical methodological specifications** before execution. Once modified, this framework will produce **publication-quality results** suitable for PLOS Genetics.

---

## REFERENCES FOR BEST PRACTICES

1. **Weir-Cockerham FST**: Weir & Cockerham (1984) Evolution
2. **Procrustes Analysis**: Peres-Neto & Jackson (2001) Oikos
3. **Mantel Test**: Mantel (1967) Cancer Research
4. **SFRatios**: Johri et al. (2024) Genetics
5. **SV Population Genomics**: Hansson et al. (2024) Molecular Ecology (salmon example)
6. **1kGP 2022 Panel**: Byrska-Bishop et al. (2022) Cell
7. **scikit-allel**: Miles et al. (2020) Bioinformatics

---

*Evaluation completed: 2026-02-28*
