# ChromaTyper

Identification of whole-gene chromatin configurations (**chromatypes**) from epigenomic features.

A chromatype is a recurrent chromatin configuration of one gene, shared by a group of
epigenomes. ChromaTyper reads the chromatin state of every 200 bp bin along a gene in
many epigenomes, and groups the epigenomes whose gene-body patterns are alike. Like a
genotype, a chromatype is discrete, so it can be counted, shared and tested for
enrichment; unlike a genotype, it varies between cell types.

The method is unsupervised: no disease, expression or cell-type labels are used to
define chromatypes.

## Method

For each gene, independently:

1. **States.** Each epigenome has one ChromHMM state per 200 bp bin along the gene
   (developed on the IHEC EpiAtlas 18-state model from six histone marks: H3K4me3,
   H3K4me1, H3K27ac, H3K36me3, H3K27me3, H3K9me3).
2. **State distance.** Two states are as far apart as their emission profiles
   (the probability of each mark in each state):
   `d(i, j) = sqrt(mean((emission_i - emission_j)^2))`.
   States with similar marks (e.g. weak transcription and quiescent) are close;
   distinct ones (active promoter vs transcribed body) are far.
3. **Sample distance.** For each pair of epigenomes, the state distances are combined
   bin by bin along the gene as a weighted root mean square:
   `D(A, B) = sqrt(sum(w * d^2) / sum(w))`, with `w` the bin's overlap with the gene in bp.
4. **Chromatypes.** Average-linkage hierarchical clustering of the epigenomes,
   cut adaptively with `dynamicTreeCut` into groups of at least `minClusterSize`
   epigenomes. Epigenomes that do not fit a group are left unassigned.

## Input

| File | Content |
|---|---|
| Emission matrix (CSV) | One row per ChromHMM state, one column per histone mark; first column the state ID. |
| Gene bins (long table, CSV or RDS) | `gene_id`, `sample_id`, `bin_id`, `state`, optional `bin_weight` (bp overlap; 1 if absent). Within a gene, `bin_id` must refer to the same position in every sample. |

## Output

- A gene x epigenome matrix of chromatype labels (`Chromotype_1`, `Chromotype_2`, ... or `Unassigned`).
- A per-gene summary (number of epigenomes, bins and chromatypes) and a list of genes that could not be processed.

## Main parameters

| Parameter | Default | Meaning |
|---|---|---|
| `deepSplit` | 0 | `dynamicTreeCut` sensitivity (0 = fewest, largest groups) |
| `minClusterSize` | 5 | Minimum epigenomes per chromatype |
| `pamRespectsDendro` | FALSE | PAM reassignment step of `dynamicTreeCut` |
| aggregation | `rms` | How bin distances are combined (`rms` or `mean`) |
| linkage | `average` | Hierarchical clustering linkage |
| min valid fraction | 0.8 | Minimum fraction of shared, valid bins for a sample pair |

Labels are per gene: `Chromotype_4` of one gene has no relation to `Chromotype_4` of another.

## Requirements

R (>= 4.2) with `data.table`, `dynamicTreeCut`, `fastcluster` (or base `hclust`) and
`parallel`.

## Usage

From the repository root:

```bash
Rscript chromatyper.R --emissions example/emissions.csv --bins example/gene_bins.csv --out results/toy
```

This writes `results/toy_chromatypes.csv` (epigenome x gene labels), `results/toy_summary.csv`
and `results/toy_failures.csv`. Optional flags: `--deepSplit`, `--minClusterSize`,
`--pamRespectsDendro`, `--aggregation`, `--minValid` and `--cores` (see Main parameters).

From R:

```r
library(data.table)
source("R/chromatyper.R")
res <- chromatyper(read_gene_bins("example/gene_bins.csv"),
                   read_emissions("example/emissions.csv"),
                   deepSplit = 0, minClusterSize = 5)
res$matrix[1:5, ]
```

`example/emissions.csv` is the emission matrix of the IHEC EpiAtlas 18-state model.
`example/gene_bins.csv` is a toy table (two genes, 30 epigenomes, three chromatin patterns
with noise), made by `example/make_example.R`.

## Citation

Ochoteco-Asensio J, Rico D. *Stolen, Not Invented? Cancer-associated chromatin states
shared with healthy cells.* Manuscript in preparation. Centro Andaluz de Biología
Molecular y Medicina Regenerativa (CABIMER), Seville, Spain.

## License

GPL-3.0 (see `LICENSE`).
