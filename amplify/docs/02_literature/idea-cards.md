# Candidate Research Ideas

═══════════════════════════════════════
Idea 1: Discordant Evolutionary Signals — How SNVs, INDELs, and SVs Paint Different Pictures of Human Population History
═══════════════════════════════════════
Core question:   Do different variant types (SNVs, INDELs, SVs) produce concordant or discordant population structure, differentiation, and selection signals in human populations, and what do the discordances reveal about type-specific evolutionary forces?
Novelty source:  Contradiction Mining (Insight 1) + Assumption Challenging (Insight 2) + Trend Extrapolation (Insight 6)
Why it matters:  If variant types tell different population stories, it challenges the universal use of SNP-derived population structure as a proxy for all genomic variation — with implications for GWAS correction, imputation, and ancestry inference.
Feasibility:     high — Data available (1kGP 2022 panel has all 3 types); standard pop-gen tools (PCA, FST, SFS, admixture); no DL needed; 128 CPUs + 1TB RAM sufficient
Risk:            Low — even if concordance is high, quantifying it is novel. If discordance exists, it's a strong finding.
Competition:     low — Only precedent is Hansson 2024 in salmon. No human study exists.
Estimated scope: journal (PLOS Genetics / Genome Research / MBE level)

═══════════════════════════════════════
Idea 2: Variant-Type-Specific Selection Landscapes Across Human Populations
═══════════════════════════════════════
Core question:   Do genomic regions under population-specific selection show different signals when assayed with SNVs vs. INDELs vs. SVs, and can cross-type discordances identify novel selection targets missed by SNP-only scans?
Novelty source:  Assumption Challenging (Insight 7) + Counterfactual Reasoning (Insight 5)
Why it matters:  Selection scans are overwhelmingly SNP-based. SVs and INDELs have larger functional effects on average, meaning SV/INDEL-based selection scans might reveal adaptive loci invisible to SNP scans.
Feasibility:     medium — Requires careful FST outlier analysis and possibly PBS (population branch statistics) for multiple variant types; SV annotation needed
Risk:            Medium — SV counts may be too low for robust outlier analysis in some chromosomal regions
Competition:     low-medium — Some SV FST work exists (Asian catalog 2024) but no systematic cross-type comparison
Estimated scope: journal (Genome Research / MBE level)

═══════════════════════════════════════
Idea 3: Updated Germline Mutation Spectrum Analysis with High-Coverage Data and Extension to INDELs
═══════════════════════════════════════
Core question:   Does the high-coverage 1kGP data confirm, refine, or revise the known population-specific mutation spectrum (TCC→TTC in Europeans, CCG→CTG in South Asians), and do INDELs show analogous population-specific length/context signatures?
Novelty source:  Trend Extrapolation (Insight 11) + Cross-Domain Transfer (Insight 3)
Why it matters:  Mutation spectrum analysis is foundational to population genetics. Updating with 30x data and extending to INDELs could reveal new biology about DNA repair pathway variation across populations.
Feasibility:     high — Trinucleotide context extraction is straightforward with bcftools/Python; INDEL length analysis is simple
Risk:            Low-medium — Risk of merely confirming known results without new findings
Competition:     medium — Harris group and others actively publish in this area
Estimated scope: conference/journal (MBE / Genome Biology level)

═══════════════════════════════════════
Idea 4: Population-Specific Genetic Load Across Variant Types
═══════════════════════════════════════
Core question:   Does the distribution of deleterious variants across human ancestry groups differ between SNVs, INDELs, and SVs, and does the Stolyarova et al. 2025 "similar distribution" finding hold for large-effect structural variants?
Novelty source:  Contradiction Mining (Insight 10)
Why it matters:  If SVs show different load patterns than SNVs, it would challenge the mutation-drift-selection model for large-effect variants and have implications for disease risk assessment across populations.
Feasibility:     medium — Requires functional annotation of SVs (harder than SNVs); CADD scores available for SNVs/INDELs but not all SVs
Risk:            Medium-high — SV functional annotation quality may limit conclusions
Competition:     medium — Stolyarova 2025 is very recent; extending to SVs is natural follow-up
Estimated scope: journal (AJHG / PNAS level if findings are strong)

═══════════════════════════════════════
Idea 5: Allele Frequency Spectrum Shapes as Population-Specific Evolutionary Fingerprints Across Variant Types
═══════════════════════════════════════
Core question:   How do the shapes of allele frequency spectra differ across variant types (SNV, INDEL, SV) and across populations, and can population × variant-type SFS interactions reveal distinct evolutionary regimes (drift-dominated vs. selection-dominated)?
Novelty source:  Counterfactual Reasoning (Insight 5) + Cross-Domain Transfer (applying SFRatios framework cross-type)
Why it matters:  SFS shape is a fundamental summary of evolutionary forces. Comparing SFS across variant types within the same populations controls for demography, isolating mutation/selection effects.
Feasibility:     high — SFS computation is standard; comparison framework exists (SFRatios, dadi, moments); data available
Risk:            Low — SFS will definitely differ; the question is whether differences are biologically interesting
Competition:     low — Nobody has done systematic cross-type SFS comparison in humans
Estimated scope: journal (Genetics / MBE level)
