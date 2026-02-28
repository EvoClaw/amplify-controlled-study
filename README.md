# Amplify Controlled Study

> **Controlled evaluation of AI-assisted research automation using the 1000 Genomes Project dataset**

---

## Key Findings

This experiment was run on **Claude claude-4.6-sonnet-medium** (opus 4.6). Both setups were instructed to operate as autonomously as possible — from data exploration through to a completed research paper. Human intervention was minimal throughout: the researcher occasionally issued simple continuation prompts (e.g., "continue") but did not guide the scientific direction, correct methods, or edit the writing.

The central result is straightforward: **opus 4.6 can already produce scientific research papers independently.** Both agents — with and without Amplify — completed the full research cycle and delivered manuscripts. This is a demonstration of large language model capability, not a limitation.

**Whether Amplify is used determines the quality ceiling, not whether a paper gets written.** With Amplify, the agent followed a structured research methodology — literature review, multi-agent problem validation, hypothesis-driven design, convergence analysis — and produced a more coherent, theoretically grounded population genetics paper. Without Amplify, the agent still produced a complete, technically solid paper, but defaulted to a more computational, descriptive approach with a narrower scope.

**More human involvement would improve both outputs.** A researcher actively engaged at key decision points — refining the research question, connecting disparate analyses, selecting which figures to include — would elevate the quality of either setup. This is entirely expected: the autonomous experiments here represent a lower bound on what the human–AI combination can achieve.

---

## ⚠️ Important Disclaimer

> **All scientific content in this repository — experiments, analyses, figures, and papers — was generated entirely by AI agents operating autonomously. This repository exists solely as a system performance benchmark for AI research automation tools. The papers have NOT been peer-reviewed. Do not treat any scientific claims herein as validated research. Please exercise critical judgment when reading.**

---

## What Is This?

This repository contains a head-to-head comparison of two AI research automation setups on an identical task:

- **`amplify/`** — an AI agent augmented with [Amplify](https://evoclaw.github.io/amplify/), a structured research automation framework with 24 skills, 7 phases, 4 mandatory human-checkpoint gates, and multi-agent deliberation panels.
- **`baseline/`** (originally `cursor/`) — a baseline AI agent without Amplify, given a similar starting prompt.

Both agents were given the same task: **use the publicly available 1000 Genomes Project high-coverage dataset (3,202 individuals, 26 populations) to conduct original scientific discovery — without training deep learning models — and produce a full research paper.** The computer used has a GPU available (used for non-DL computation), bcftools, PLINK2, Python, and related bioinformatics tools.

---

## 📄 Output Papers

Both agents produced complete research manuscripts. PDFs are included in this repository.

| Setup | Paper Title | PDF |
|---|---|---|
| **With Amplify** | *Size-dependent attenuation of population differentiation across the mutational spectrum in 3,202 human genomes* | [`amplify/paper/main.pdf`](amplify/paper/main.pdf) |
| **Baseline (no Amplify)** | *Population-Stratified Germline Indel Mutation Spectra Reveal Differential Mutagenic Process Contributions Across Diverse Human Populations* | [`baseline/analysis/manuscript/main.pdf`](baseline/analysis/manuscript/main.pdf) |

> ⚠️ Both PDFs carry an explicit AI-generated content disclaimer on the first page.

---

## Evaluation: Amplify vs. Baseline

### Amplify — *Cross-variant-type population structure comparison*

**Research focus:** Systematic comparison of population structure and differentiation signals across SNVs (12.7M), INDELs (3.3M), and SVs (18.5K) — a unified analysis of all three variant classes.

**What the agent did well:**
- Formulated a clear, theoretically-grounded research question: do different variant types encode the same or different population signals, and if so, why?
- Designed a convergence analysis framework to distinguish statistical power effects from genuine biological discordance — a methodologically sound and original contribution.
- Found a continuous, monotonic FST gradient with variant size (SNV > INDEL > SV), consistent with size-dependent purifying selection.
- Identified asymmetric selective constraint between insertions and deletions (insertions ~25% lower FST than deletions).
- Used Procrustes analysis to show SV-SNV concordance plateaus at r ≈ 0.60 regardless of variant count — establishing the discordance as biology, not sampling noise.
- Writing style is clearly population genetics: hypothesis-driven, with explicit reference to classical population genetics theory (Wright's FST, purifying selection, bottleneck demography).
- Ran extensive analyses and generated 16 distinct figures in total.

**Limitations observed:**
- SV-based analyses are inherently limited by the smaller SV catalog (18.5K vs. millions of SNVs/INDELs).
- The insertion-deletion asymmetry finding has a caveat (differential short-read ascertainment), which the paper acknowledges.
- **Figure integration gap:** The agent produced 16 distinct figures across the `figures/` directory, but only 4 were incorporated into the final paper (Fig. 1 PCA comparison, convergence curve, size-stratified FST, population trees). Many analytical results — FST heatmaps, SFS comparisons, cross-type FST, concordance summaries, SFS KS tests, and more — were computed and plotted but never written up. This suggests a disconnect between the execution phase and the paper-writing phase: the agent ran more analyses than it could coherently narrate.

---

### Baseline (no Amplify) — *Germline indel mutation spectrum*

**Research focus:** Population-stratified analysis of germline indel mutation spectra using NMF decomposition — focused exclusively on indels, borrowing a computational framework from cancer mutational signature analysis.

**What the agent did well:**
- Produced a thorough quantitative analysis with solid statistical testing (chi-square, JSD, Bonferroni correction).
- Classified 9M indels into a 52-channel system and demonstrated statistically significant population differences across all channels.
- Applied NMF decomposition (borrowed from COSMIC cancer signature methodology) to decompose spectra into two principal signatures.
- Correctly identified that African populations show elevated GC-deletion enrichment and East Asian populations show AT-homopolymer insertion enrichment.
- Writing is technically correct and references relevant literature.

**Where it diverges from pure population genetics:**
- The NMF signature approach introduces an algorithmic/machine learning component (from cancer genomics) — the paper has a noticeably more computational flavor.
- The connection to mutational mechanisms (oxidative damage, polymerase slippage) is speculative and more weakly supported than the population structure claims.
- The research question is narrower — restricted to indels only, missing the cross-variant comparison that would reveal size-dependent selection.
- The "discovery" (populations differ in their indel spectra) is less novel than the Amplify agent's convergence framework.
- The paper lacks a unifying biological framework — it describes patterns but the mechanistic story is fragmented.

---

### Head-to-Head Summary

| Criterion | Amplify | Baseline |
|---|---|---|
| Research scope | SNVs + INDELs + SVs (all variant classes) | INDELs only |
| Research question clarity | Strong theoretical framing | Descriptive / exploratory |
| Methodological novelty | Convergence analysis framework (original) | NMF borrowed from cancer genomics |
| Scientific cohesion | Unified narrative, clear takeaway | Multiple observations, fragmented story |
| Population genetics depth | High — classical theory well-integrated | Moderate — computation-heavy |
| Figures produced | 16 distinct figures generated | 5 figures generated |
| Figures in paper | 4 (large gap vs. produced) | 5 (fully integrated) |
| Algorithmic content | Minimal (bioinformatics pipeline) | Moderate (NMF decomposition) |
| Paper quality | Closer to a population genetics journal paper | Closer to a computational genomics report |

**Overall assessment:** The Amplify-augmented agent produced a more cohesive, hypothesis-driven population genetics paper with a novel analytical framework (convergence analysis) and a broader scope. The baseline agent produced competent but narrower work with stronger computational tool dependence and weaker biological narrative integration — reflecting the tendency of unstructured AI coding assistants to reach for algorithmic methods (NMF) rather than building domain-specific scientific frameworks.

---

## Repository Structure

```
amplify/
├── paper/          # LaTeX source and compiled PDF
├── figures/        # All publication figures (PDF/PNG)
├── scripts/        # Analysis scripts
└── results/        # (large files, not tracked in git)

baseline/           # Originally named "cursor/", renamed to remove branding
├── analysis/
│   ├── manuscript/ # LaTeX source and compiled PDF
│   ├── figures/    # All publication figures
│   └── scripts/    # Analysis scripts
└── results/        # (large files, not tracked in git)
```

Data (`1000GP/`) and large intermediate results are not tracked in this repository. The 1000 Genomes Project high-coverage dataset is publicly available at [https://www.internationalgenome.org/](https://www.internationalgenome.org/).

---

## About Amplify

[Amplify](https://evoclaw.github.io/amplify/) is an open-source, skills-based framework that teaches AI coding assistants how to conduct scientific research with rigor — from literature review and multi-agent problem validation through experiment execution and paper writing. It enforces mandatory human checkpoints, metric locking, anti-cherry-picking, and multi-agent deliberation at every stage.

- Website: [https://evoclaw.github.io/amplify/](https://evoclaw.github.io/amplify/)
- GitHub: [https://github.com/EvoClaw/amplify](https://github.com/EvoClaw/amplify)

---

*This study was designed and supervised by a human researcher. All AI agent actions were run autonomously within the structured Amplify framework (or without it, for the baseline).*
