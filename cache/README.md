# Epilepsy19
Analysis for the Epilepsy paper

To reproduce the analyses you need first download the count matrices from **TODO** and pre-process them using **TODO**.
It stores the pre-processed matrices in `cache/count_matrices.rds`. Then, you need to run [alignment.Rmd](**TODO**) to get 
`cache/con_filt_cells_10x.rds` and `cache/con_filt_samples_10x.rds`.

Alternatively, you can download the ready Conos objects from the pklab server:

```
cd cache
wget http://pklab.med.harvard.edu/viktor/publications/Epilepsy19/count_matrices.rds.gz
wget http://pklab.med.harvard.edu/viktor/publications/Epilepsy19/con_filt_cells_10x.rds.gz
wget http://pklab.med.harvard.edu/viktor/publications/Epilepsy19/con_filt_samples_10x.rds.gz
gunzip ./*.gz
```

To investigate data interactively, you can access the Pagoda 2 Web App:
- [All datasets]()
- [Filtered samples]()
- Epilepsy
  - [E1]()
  - [E2]()
  - [E3]()
  - [E4]()
  - [E5]()
  - [E6]()
  - [E7]()
  - [E8]()
  - [E9]()
- Control
  - [C1]()
  - [C2]()
  - [C3]()
  - [C4]()
  - [C5]()
  - [C6]()
  - [C7]()
  - [C8]()
  - [C9]()
  - [C10]()