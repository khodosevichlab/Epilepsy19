# Epilepsy19
Analysis for the Epilepsy paper

To see the compiled notebooks, please visit [the website](https://khodosevichlab.github.io/Epilepsy19/).

- Pre-processing
  - [Raw data filtering](./analysis/prep_filtration.Rmd)
  - [Alignment of the datasets](./analysis/prep_alignment.Rmd)
- Analyses
  - [Overview of the datasets](./analysis/fig_overview.Rmd)
  - [Distances between cell types](./analysis/fig_type_distance.Rmd)
  - [Differential Expression and Gene Ontology analyses](./analysis/fig_go.Rmd)
  - [Cell Type prioritization summary](./analysis/fig_summary.Rmd)
  - [Comparison to Smart-Seq2 and Allen Brain Institute data](./analysis/fig_smart_seq.Rmd)
  - [Gene module analysis](./gene_modules/)

<!--
## Prepare data

Before running the analyses, create the necessary folders

```
mkdir -p output/overview output/fig_type_diff output/fig_go/gene_scattermaps output/fig_go/pathway_clustering output/total_ranking output/fig_neun output/fig_smart_seq
```
-->