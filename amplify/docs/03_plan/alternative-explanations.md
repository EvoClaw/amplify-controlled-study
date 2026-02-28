# Alternative Explanations

| Expected Finding | Alternative Explanation | How to Rule Out |
|-----------------|----------------------|-----------------|
| SVs show different PCA from SNVs | SV calling artifacts (not biology) | Sensitivity analysis with strict/relaxed filters; subsample SNVs to match SV count; check if SV PCA improves with higher-quality subset |
| INDELs show intermediate concordance | Variant count difference (statistical power) | Random subsample of SNVs to INDEL count; bootstrap CIs |
| FST discordance at specific loci | Ascertainment bias in SV calling | Check SV calling rates per population; normalize by call rate |
| SFS shape differences across types | Different mutation rates (not selection) | Compare within neutral regions vs. functional regions; use SFRatio-like framework |
| Some chromosomes show high discordance | Inversions or other structural features | Annotate high-discordance regions for known inversions/CNVs |
| Larger SVs show higher FST | Fewer large SVs → sampling noise | Report counts per size bin; bootstrap CIs; permutation tests |
| Admixed populations show more discordance | Population size effects | Control for sample size; bootstrap within admixed groups |
