# rWGCNA gene module analyses

## Quick start

0. clone this repo: `git clone https://github.com/khodosevichlab/Epilepsy19.git`
1. install packages
    1. `cd Epilepsy19`
    2. start an R session
    3. install renv: `install.packages("renv")`
    4. `renv::restore(lockfile="./gene_modules/renv.lock")`
2. preproces.Rmd
3. rWGCNA
    1. clone the rWGCNA repo: `git clone --recurse-submodules https://github.com/perslab/wgcna-toolbox.git`
    2. run rWGCNA after updating paths: `Rscript ./rwgcna_main_seurat3.0.R --prefixData ep_10x_v2 --dirTmp <TMP_FILES_DIR> --prefixRun wgcna3 --pathDatExpr <DATA_DIR>/ep_10x_v2_logcounts.csv.gz   --pathMetadata <DATA_DIR>/ep_10x_v2_annotation_filter.csv   --dirProject <WGCNA_OUTPUT_DIR>  --colIdents l2  --minGeneCells 20 --minCellClusterSize 50 --featuresUse PCLoading --nFeatures 5000 --nPC 120 --corFnc cor --networkType "c('signed hybrid')" --nRepTOM 100 --hclustMethod average --minClusterSize 15L --deepSplit 2 --moduleMergeCutHeight "c(.15)" --pamStage "c(FALSE)" --kMReassign T --kMSignifFilter T --nRepJackStraw 0 --fuzzyModMembership kIM --RAMGbMax 200 &> log_ep_10x_v2_wgcna3.txt` 
4. prune overlapping gene modules
    1. clone the gene network scripts repo: `git clone https://github.com/perslab/geneNetworkScripts.git`
    2. run the module merging/pruning script after updating paths: `Rscript ./gene_module_merge2.R --path_df_NWA <DIR_PROJECT>/tables/ep_10x_v2_wgcna3_geneMod.csv.gz --colGeneWeights pkMs --colGeneNames genes --minPropIntersect 0.75 --minWeightedCor 0.75 --mergeOrPrune prune --dirOut $DIR_PROJECT --prefixOut _merged &> log_wgcna3_prune.txt`
5. filter_PMI_injury.Rmd
6. compute_module_embeddings.Rmd 
7. linear_models.Rmd
8. geneset_enrichment.Rmd 
9. GO.Rmd
10. Compute module preservation in level 4 cell subtypes
    1. `cd geneNetworkScripts`
    2. run module preservation after updating paths: `Rscript ./gene_module_preservation.R --path_dt_geneMod <PROJECT_DIR>/ep_10x_v2_wgcna3_geneMod_merged_filterPMIsurg.csv.gz --colGeneNames genes --colGeneWeights pkMs --colMod module_filter_PMI_surgery --colCellClust cell_cluster_filter_PMI_surgery --vec_pathsDatExpr "c('ep_10x_v2_l2'='<DATA_DIR>/ep_10x_v2_logcounts.csv.gz','ep_10x_v2_l4'='<DATA_DIR>/ep_10x_v2_logcounts.csv.gz')" --vec_pathsMetadata "c('ep_10x_v2_l2'='<DATA_DIR>/ep_10x_v2_cell_sample_annot_filter.csv', 'ep_10x_v2_l4'=<DATA_DIR>/ep_10x_v2_cell_sample_annot_filter.csv')" --vec_metadataIdentCols "c('ep_10x_v2_l2'='l2', 'ep_10x_v2_l4'='l4')" --list_list_vec_identLvlsMatch 'list("L2_3_Cux2"=list("ep_10x_v2_l4"=c("L2_3_Cux2_Frem3_Unc5d","L2_Cux2_Lamp5", "L3_Cux2_Prss12")), "L4_Rorb"=list("ep_10x_v2_l4"=c("L4_Rorb_Schlap1_Met","L4_Rorb_Schlap1_Mme", "L4_Rorb_Schlap1_Arhgap15")), "L5_6_Fezf2"=list("ep_10x_v2_l4"=c("L5_6_Fezf2_Lrrk1_Sema3e","L5_6_Fezf2_Lrrk1_Pcp4", "L5_6_Fezf2_Tle4_Abo","L5_6_Fezf2_Tle4_Htr2c","L5_6_Fezf2_Tle4_Scube1")),    "L5_6_Themis"=list("ep_10x_v2_l4"=c("L5_6_Themis_Sema3a","L5_6_Themis_Ntng2")), "Pvalb"=list("ep_10x_v2_l4"=c("Pvalb_Crh", "Pvalb_Lgr5", "Pvalb_Nos1", "Pvalb_Mepe")))' --minCellClusterSize 10 --minpropModGenesInTestLvl 0.5 --minPropModGenesNon0inTestLvl 0.1 --dirOut <PROJECT_DIR> --prefixOut ep_10x_v2_l2_vs_l4 --dirTmp <TMP_FILES_DIR> --networkType 'c("signed hybrid")' --corFnc cor 
--RAMGbMax 200`
11. makeplots.Rmd   
