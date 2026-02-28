# Size-dependent attenuation of population differentiation across the mutational spectrum in 3,202 human genomes

## Paper Structure

```
paper/
├── main.tex              # Main document (inputs all sections)
├── preamble.tex          # LaTeX packages and macros
├── references.bib        # Bibliography (34 references)
├── sections/
│   ├── abstract.tex      # Abstract + Author Summary
│   ├── introduction.tex  # Introduction
│   ├── results.tex       # Results (6 subsections, 2 tables)
│   ├── discussion.tex    # Discussion (5 subsections)
│   └── methods.tex       # Materials and Methods
├── supplementary/
│   └── supplementary.tex # Supporting Information descriptions
└── figures/              # -> symlink to ../figures/
```

## Figures

### Main Text
- **Fig 1**: PCA comparison + Procrustes heatmap + FST scatter (fig1_pca_comparison, fig2_procrustes_overlay, fig4_fst_crosstype)
- **Fig 2**: Convergence curve (fig_convergence_curve)
- **Fig 3**: Size-stratified FST + SV subtypes (fig_size_stratified_fst, fig_sv_subtypes)
- **Fig 4**: Population trees + FST key pairs (fig7_population_trees, fig_fst_key_pairs)

### Supplementary
- S1: SFS per super-population (fig5_sfs_comparison, fig6_sfs_heatmap)
- S2: Count-matched PCA (fig10_count_matched)
- S3: SV subtype SFS (fig_sv_subtypes)
- S4: Size distributions (fig_supp_size_distributions)
- S5: Uniform MAF SFS (fig_sfs_uniform_maf)

## Compilation

```bash
cd paper/
pdflatex main
bibtex main
pdflatex main
pdflatex main
```

## Data

1000 Genomes Project 2022 high-coverage release:
- Source: https://www.internationalgenome.org/data-portal/data-collection/30x-grch38
- 3,202 samples, 26 populations, chromosomes 1-22

## Analysis Scripts

```
scripts/
├── 01_extract_parallel.sh    # Variant extraction (SNV, INDEL, SV)
├── 04_frequencies.sh         # Per-population allele frequencies
├── 05_analysis.py            # Core analyses (PCA, FST, SFS, Procrustes)
├── 06_robustness.py          # Count-matched subsampling
└── 07_supplements.py         # Size-stratified, convergence, SV subtypes
```

## Key Results

| Metric | SNV | INDEL | SV |
|--------|-----|-------|-----|
| Variants | 12.67M | 3.34M | 18,519 |
| Mean FST | 0.080 | 0.072 | 0.059 |
| Procrustes vs SNV | — | 0.993 | 0.598 |
