# Epilepsy19
Analysis for the [Identification of epilepsy-associated neuronal subtypes and gene expression underlying epileptogenesis](https://doi.org/10.1038/s41467-020-18752-7) paper.

To see the compiled notebooks, please visit [the website](https://khodosevichlab.github.io/Epilepsy19/).

## Interactive data exploration

- [Conos alignment](http://pklab.med.harvard.edu/viktor/pagodaURL/index.html?fileURL=http://pklab.med.harvard.edu/viktor/publications/Epilepsy19/con_filt_cells.bin) (*con_filt_cells*)
- [Conos alignment with bad samples removed](http://pklab.med.harvard.edu/viktor/pagodaURL/index.html?fileURL=http://pklab.med.harvard.edu/viktor/publications/Epilepsy19/con_filt_samples.bin) (*con_filt_samples*)
- By sample
  - Control: 
[C1](http://pklab.med.harvard.edu/viktor/pagodaURL/index.html?embedding=UMAP&fileURL=http://pklab.med.harvard.edu/viktor/publications/Epilepsy19/C1.bin), [C2](http://pklab.med.harvard.edu/viktor/pagodaURL/index.html?embedding=UMAP&fileURL=http://pklab.med.harvard.edu/viktor/publications/Epilepsy19/C2.bin), [C3](http://pklab.med.harvard.edu/viktor/pagodaURL/index.html?embedding=UMAP&fileURL=http://pklab.med.harvard.edu/viktor/publications/Epilepsy19/C3.bin), [C4](http://pklab.med.harvard.edu/viktor/pagodaURL/index.html?embedding=UMAP&fileURL=http://pklab.med.harvard.edu/viktor/publications/Epilepsy19/C4.bin), [C5](http://pklab.med.harvard.edu/viktor/pagodaURL/index.html?embedding=UMAP&fileURL=http://pklab.med.harvard.edu/viktor/publications/Epilepsy19/C5.bin), [C6](http://pklab.med.harvard.edu/viktor/pagodaURL/index.html?embedding=UMAP&fileURL=http://pklab.med.harvard.edu/viktor/publications/Epilepsy19/C6.bin), [C7](http://pklab.med.harvard.edu/viktor/pagodaURL/index.html?embedding=UMAP&fileURL=http://pklab.med.harvard.edu/viktor/publications/Epilepsy19/C7.bin), [C8](http://pklab.med.harvard.edu/viktor/pagodaURL/index.html?embedding=UMAP&fileURL=http://pklab.med.harvard.edu/viktor/publications/Epilepsy19/C8.bin), [C9](http://pklab.med.harvard.edu/viktor/pagodaURL/index.html?embedding=UMAP&fileURL=http://pklab.med.harvard.edu/viktor/publications/Epilepsy19/C9.bin), [C10](http://pklab.med.harvard.edu/viktor/pagodaURL/index.html?embedding=UMAP&fileURL=http://pklab.med.harvard.edu/viktor/publications/Epilepsy19/C10.bin)
  - Epilepsy: 
[E1](http://pklab.med.harvard.edu/viktor/pagodaURL/index.html?embedding=UMAP&fileURL=http://pklab.med.harvard.edu/viktor/publications/Epilepsy19/E1.bin), [E2](http://pklab.med.harvard.edu/viktor/pagodaURL/index.html?embedding=UMAP&fileURL=http://pklab.med.harvard.edu/viktor/publications/Epilepsy19/E2.bin), [E3](http://pklab.med.harvard.edu/viktor/pagodaURL/index.html?embedding=UMAP&fileURL=http://pklab.med.harvard.edu/viktor/publications/Epilepsy19/E3.bin), [E4](http://pklab.med.harvard.edu/viktor/pagodaURL/index.html?embedding=UMAP&fileURL=http://pklab.med.harvard.edu/viktor/publications/Epilepsy19/E4.bin), [E5](http://pklab.med.harvard.edu/viktor/pagodaURL/index.html?embedding=UMAP&fileURL=http://pklab.med.harvard.edu/viktor/publications/Epilepsy19/E5.bin), [E6](http://pklab.med.harvard.edu/viktor/pagodaURL/index.html?embedding=UMAP&fileURL=http://pklab.med.harvard.edu/viktor/publications/Epilepsy19/E6.bin), [E7](http://pklab.med.harvard.edu/viktor/pagodaURL/index.html?embedding=UMAP&fileURL=http://pklab.med.harvard.edu/viktor/publications/Epilepsy19/E7.bin), [E8](http://pklab.med.harvard.edu/viktor/pagodaURL/index.html?embedding=UMAP&fileURL=http://pklab.med.harvard.edu/viktor/publications/Epilepsy19/E8.bin), [E9](http://pklab.med.harvard.edu/viktor/pagodaURL/index.html?embedding=UMAP&fileURL=http://pklab.med.harvard.edu/viktor/publications/Epilepsy19/E9.bin)

*\*See the demonstration for [Pagoda 2 web apps](https://www.youtube.com/watch?v=j6PmRtOBTRM)*

## The raw notebooks
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

## Data, used for analysis

You can find all data files with download guides here:
http://kkh.bric.ku.dk/Epilepsy19/ 

Conos objects and filtered count matrices (must be saved into the [`cache`](./cache) folder):
- [con_filt_cells.rds] alignment of all neuronal cells
- [con_filt_samples.rds]: 3 bad samples removed from the analysis
- [count_matrices.rds]: filtered count matrices

Additional data:
- [cms_raw.rds]: raw count matrices (as they were loaded in chunks [2](https://github.com/khodosevichlab/Epilepsy19/blob/611a4eb5c893d273258c56717d0682747db65d2f/analysis/prep_filtration.Rmd#L23) and [5](https://github.com/khodosevichlab/Epilepsy19/blob/611a4eb5c893d273258c56717d0682747db65d2f/analysis/prep_filtration.Rmd#L51) of [the pre-processing notebook](https://khodosevichlab.github.io/Epilepsy19/prep_filtration.html))

Relevant meta-data:
- [CellAnnotatoR-compatible list of marker genes](./metadata/neuron_markers.txt)
- [Cell type annotation](./metadata/annotation.csv)
- [Sample information](./metadata/sample_info.csv)

<!--
## Prepare data

Before running the analyses, create the necessary folders

```
mkdir -p output/overview output/fig_type_diff output/fig_go/gene_scattermaps output/fig_go/pathway_clustering output/total_ranking output/fig_neun output/fig_smart_seq
```
-->

## Citation

If you used this code for your analysis, please cite the paper:

> Pfisterer, U., Petukhov, V., Demharter, S. et al. Identification of epilepsy-associated neuronal subtypes and gene expression underlying epileptogenesis. Nat Commun 11, 5038 (2020). https://doi.org/10.1038/s41467-020-18752-7
