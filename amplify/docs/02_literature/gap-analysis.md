# Gap Analysis

## Gap 1: No Systematic Multi-Variant-Type Population Genomics in Humans
- **Status:** Not solved
- **Why not solved:** Previous 1kGP studies focused on SNVs; SVs were too few/unreliable from low-coverage data; the 2022 panel is the first to provide reliable short-read SV calls alongside SNVs and INDELs in the same cohort
- **Approach:** Systematic comparison of population structure, differentiation, and selection signals across SNV/INDEL/SV using the 2022 panel
- **Confidence:** high (confirmed through extensive search — only non-human analogues exist)

## Gap 2: Variant-Type-Specific Allele Frequency Spectrum Signatures
- **Status:** Not solved for human multi-type comparison
- **Why not solved:** Most SFS analyses focus on SNVs; INDEL and SV SFS require different handling due to size/complexity
- **Approach:** Joint SFS analysis partitioned by variant type and population, comparing shapes and population-specific deviations
- **Confidence:** high (SFRatios and other methods exist but haven't been applied cross-type)

## Gap 3: Population Structure Concordance/Discordance Across Variant Types
- **Status:** Not solved
- **Why not solved:** Nobody has asked whether PCA on SVs alone gives the same population tree as PCA on SNVs
- **Approach:** Independent PCA, admixture analysis, and phylogenetic trees per variant type; quantify concordance
- **Confidence:** high (direct gap, novel question)

## Gap 4: Variant-Type-Specific Selection Signatures
- **Status:** Partially explored (individual type studies exist)
- **Why not solved:** FST-based outlier analysis has been done for SNPs and separately for SVs, but not compared within the same framework
- **Approach:** Parallel selection scans across variant types; identify loci where different types give conflicting signals
- **Confidence:** medium (some SNP-based selection scans exist, but cross-type comparison is novel)

## Gap 5: INDEL Length Spectrum and Population-Specific Patterns
- **Status:** Not solved
- **Why not solved:** INDELs are often grouped as a monolithic category; length distribution across populations is unexplored
- **Approach:** Analyze INDEL length spectrum by population; test for population-specific length preferences (possibly related to repair pathway differences)
- **Confidence:** medium (less clear whether findings will be interesting)

## Gap 6: Functional Enrichment of Population-Differentiated Variants by Type
- **Status:** Partially explored
- **Why not solved:** Functional annotation of SVs is harder than SNVs; systematic cross-type functional enrichment comparison is lacking
- **Approach:** Annotate high-FST variants by type; compare functional category enrichment
- **Confidence:** medium (depends on functional annotation quality for SVs)
