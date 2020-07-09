#' @export
GetPagoda <- function (cm, n.cores = 30, clustering.type = "leiden", embeding.type = "tSNE", verbose=TRUE,
                       n.pcs=100, distance="cosine", trim=5, n.odgenes=1000, graph.k=30,
                       od.genes=NULL, clustering.resolution=2, build.graph=TRUE, var.scale=TRUE, min.dist=0.5, spread=1.5, ...) {
  r <- pagoda2::Pagoda2$new(cm, trim=trim, n.cores=n.cores, verbose=verbose, ...)

  if (var.scale) {
    r$adjustVariance(plot = F, do.par = F, gam.k = 10, verbose = verbose)
  }

  if (n.pcs > 0) {
    r$calculatePcaReduction(nPcs = n.pcs, n.odgenes = n.odgenes, odgenes=od.genes, maxit = 1000, verbose=verbose, var.scale=var.scale)
  }

  if (!build.graph)
    return(r)

  r$makeKnnGraph(k = graph.k, type = "PCA", center = T, distance = distance,
                 weight.type = "none", verbose = verbose)

  for (ct in clustering.type) {
    switch (ct,
            infomap = r$getKnnClusters(method = igraph::infomap.community, type = "PCA", name = "infomap"),
            multilevel = r$getKnnClusters(method = igraph::multilevel.community, type = "PCA", name = "multilevel"),
            leiden = r$getKnnClusters(method = conos::leiden.community, type = "PCA", name = "leiden", resolution=clustering.resolution),
            stop("Unknown clustering type: ", ct)
    )
  }

  for (et in embeding.type) {
    r$getEmbedding(type = "PCA", embeddingType = et, distance=distance, min_dist=min.dist, spread=spread)
  }

  return(r)
}

#' @export
GetPagodaWebApp <- function(p2, clusters, organism=NULL, additional.metadata=list(), verbose=T, go.sets=NULL, go.env=NULL, test.pathways=TRUE) {
  if (is.null(go.env)) {
    if (is.null(organism))
      stop("Either organism or go.env must be provided")

    if (verbose) cat("Generate go environment\n")
    go.env <- pagoda2::p2.generate.go(p2, organism=organism)
  }

  if (is.null(go.sets)) {
    if (verbose) cat("Generate genesets\n")

    go.sets <- ExtractGoSets(go.env, verbose=verbose)
  }

  if (verbose) cat("Generate de geneset\n")
  de.sets <- pagoda2::get.de.geneset(p2, groups = clusters, prefix = 'de_')

  if (test.pathways) {
    if (verbose) cat("Test pathway overdispersion\n")
    p2$testPathwayOverdispersion(setenv = go.env, verbose = verbose, correlation.distance.threshold = 0.8,
                                 recalculate.pca = F, min.pathway.size = 50, max.pathway.size = 1000)
  }

  if (verbose) cat("Create app\n")
  p2.web <- p2 %>%
    pagoda2::make.p2.app(
      dendrogramCellGroups = as.factor(clusters),
      geneSets = c(go.sets, de.sets),
      additionalMetadata=additional.metadata,
      show.clusters = T);

  if (verbose) cat("All done!\n")

  return(p2.web)
}

ExtractGoSets <- function(go.env) {
  names(go.env) %>% setNames(., .) %>% sccore:::plapply(function(x)
    list(properties = list(locked=T, genesetname=x, shortdescription=GO.db::GOTERM[[x]]@Term),
         genes = c(go.env[[x]])))
}

#' @export
Pagoda2FromConos <- function(con, embedding.type="tSNE", annotation=NULL, ...) {
  p2 <- con$getJointCountMatrix(raw=T) %>% Matrix::t() %>% GetPagoda(build.graph=F, ...)
  p2$graphs$conos <- con$graph

  if (!is.null(con$embedding)) {
    p2$embeddings$PCA[[embedding.type]] <- con$embedding
  }

  if (length(con$clusters) > 0) {
    p2$clusters %<>% c(con$clusters %>% setNames(paste0("conos_", names(.))) %>%
                         lapply(`[[`, "groups") %>% lapply(as.factor))
  }

  p2$clusters$dataset <- con$getDatasetPerCell()

  if (!is.null(annotation)) {
    p2$clusters$annotation <- as.factor(annotation[names(p2$clusters$dataset)])
  }

  return(p2)
}

#' @export
ConvertMetadataToPagoda2Format <- function(...) {
  metadata.list <- list(...)# lapply(list(...), as.factor)
  metadata.list <- metadata.list[!sapply(metadata.list, is.null)] %>% lapply(as.factor)
  return(mapply(pagoda2::p2.metadata.from.factor, metadata.list, names(metadata.list), SIMPLIFY=F))
}

#' @export
PagodaWebAppFromConos <- function(con, clustering.type=1, output.path=NULL,
                                  go.env=NULL, go.sets=NULL, ...) {
  p2 <- Pagoda2FromConos(con, count.matrices, umap.embedding, annotation=annotation, ...)

  metadata <- ConvertMetadataToPagoda2Format(Annotation=p2$clusters$annotation, Dataset=p2$clusters$dataset)

  if (is.null(con$clusters[[clustering.type]])) {
    stop("Conos object doesn't have clustering of type ", clustering.type)
  }

  p2 <- GetPagodaWebApp(p2, con$clusters[[clustering.type]]$groups, additional.metadata=metadata,
                        go.sets=go.sets, go.env=go.env)

  if (!is.null(output.path)) {
    p2$serializeToStaticFast(binary.filename=output.path);
  }

  return(p2)
}

#' @export
GetScrubletScores <- function(mat, python_path, min.molecules.per.gene=10) {
  tf.in <- tempfile()
  tf.out <- tempfile()
  dt <- mat[Matrix::rowSums(mat)>=min.molecules.per.gene,] %>% Matrix::t() %>% as.matrix() %>% data.table::data.table()
  data.table::fwrite(dt, file=tf.in)

  cmd <- paste0(python_path, " -c 'import sys; import pandas; import scrublet; ",
                "df = pandas.read_csv(\"", tf.in, "\"); scrub = scrublet.Scrublet(df); ",
                "doublet_scores, predicted_doublets = scrub.scrub_doublets();",
                "pandas.DataFrame(dict(score=doublet_scores, is_doublet=predicted_doublets)).to_csv(\"",tf.out,"\");'",
                sep='')
  system(cmd, intern=F)
  x <- data.table::fread(tf.out,sep=',')[,2:3] %>% as.data.frame() %>% as.list() %>%
    lapply(`names<-`, colnames(mat))

  return(x)
}