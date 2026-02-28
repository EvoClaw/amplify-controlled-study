# Multi-Agent Brainstorming Results

## Deliberation Summary (1 round — convergence achieved)

### Agent Verdicts

| Idea | Visionary | Pragmatic | Scout | Consensus |
|------|-----------|-----------|-------|-----------|
| 1. Discordant Signals | **STRONG** | **STRONG** | **STRONG** | **STRONG** ✓ |
| 2. Selection Landscapes | VIABLE | VIABLE | VIABLE | VIABLE |
| 3. Mutation Spectrum | VIABLE | WEAK | WEAK | WEAK |
| 4. Genetic Load | WEAK | WEAK | WEAK | KILL |
| 5. SFS Fingerprints | VIABLE | STRONG | VIABLE | VIABLE+ |
| 6. Multi-Witness Framework (new) | **STRONG** | — | — | STRONG (Visionary proposed) |

### Key Debates and Resolutions

- **Idea 1 vs 6:** Visionary proposed Idea 6 (Multi-Witness Framework) as a refined version of Idea 1. All agents agree these should merge — the "multi-witness" framing adds conceptual depth.
- **Idea 5 as component:** Pragmatic and Visionary both recommend incorporating SFS analysis (Idea 5) as a core component of Idea 1, not as a standalone paper. This strengthens the analytical depth.
- **Ideas 3 and 4 eliminated:** Scout flagged very high competition risk for Idea 3 (Harris/Mathieson groups likely already updating) and high derivative risk for Idea 4 (Stolyarova/Coop/Przeworski). All agents agree to kill these.
- **Idea 2 as optional extension:** Selection landscapes analysis could be included as a secondary analysis within the merged framework.

### Ideas Merged
- Idea 1 + Idea 5 + Idea 6 → **Merged Idea: "Discordant Evolutionary Signals — A Multi-Witness Framework for Human Population Genomics"**

### Ideas Eliminated
- Idea 3: Too incremental; high competition risk (Harris/Mathieson groups)
- Idea 4: Derivative of Stolyarova 2025; SV annotation quality limits conclusions

---

## FINAL SELECTED DIRECTION

═══════════════════════════════════════
**Discordant Evolutionary Signals: A Multi-Witness Framework Reveals Variant-Type-Specific Evolutionary Forces in Human Populations**
═══════════════════════════════════════

**Core question:** Do SNVs, INDELs, and SVs produce concordant or discordant signals of population structure, differentiation, and evolutionary history in human populations, and what do the discordances reveal about type-specific evolutionary forces?

**Framework:** Each variant type is treated as an independent "witness" to population history. Concordance = shared demographic signal. Discordance = evidence for type-specific evolutionary forces (selection, mutation rate differences, structural constraints).

**Key analyses:**
1. Population structure (PCA, UMAP) separately for each variant type → concordance quantification
2. Pairwise FST across all population pairs, separately per variant type → FST concordance/discordance
3. Allele frequency spectrum (SFS) shapes per variant type per population → evolutionary fingerprints
4. Admixture proportions per variant type → test if admixture estimates are type-dependent
5. Where discordant: mechanistic investigation (functional enrichment, size effects, selection signals)

**Novelty:** First systematic multi-variant-type population genomics study in humans. Unique "multi-witness" conceptual framework. Challenges the universal SNP-as-proxy assumption.

**Feasibility:** High — data available (1kGP 2022), standard tools, ~2-3 weeks for core analyses.

**Competition:** Low — only non-human precedent (salmon). 12-18 month window.

**Target venues:** PLOS Genetics / Genome Research / MBE

**Expected output:** 
- Quantification of concordance/discordance across variant types
- Identification of genomic regions/categories driving discordance
- Mechanistic hypotheses for type-specific evolutionary forces
- Implications for ancestry inference and GWAS methodology
