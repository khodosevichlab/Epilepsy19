#' @import dplyr
NULL

fTest <- function(de.genes, n.all.genes, gene.set, ...) { # Where k is randomly drawn sample
  # See https://www.pathwaycommons.org/guide/primers/statistics/fishers_exact_test/
  n <- length(intersect(de.genes, gene.set))

  n.de.in.gs <- length(intersect(de.genes, gene.set))
  c.table <- matrix(c(n.all.genes - length(union(de.genes, gene.set)), length(de.genes) - n.de.in.gs,
                      length(gene.set) - n.de.in.gs, n.de.in.gs), nrow = 2)
  return(fisher.test(c.table, ...))
}

fTestPerDe <- function(de.gene.list, gene.set, expressed.genes, ...) {
  names(de.gene.list) %>%
    lapply(function(n) fTest(de.gene.list[[n]], length(expressed.genes), intersect(gene.set, expressed.genes), ...)) %>%
    lapply(function(x) data.frame(stat=x$estimate, stat_min=x$conf.int[1],stat_max=x$conf.int[2])) %>%
    Reduce(rbind, .) %>% dplyr::mutate(Type=names(de.gene.list))
}

plotFTestResults <- function(test.res, neuron.type.per.type, y.max, text.angle=90, line.width=0.3, y.lab="F-test statistics") {
  ggplot2::ggplot(test.res) +
    ggplot2::geom_bar(ggplot2::aes(x=Type, y=stat, fill=neuron.type.per.type[as.character(Type)]), stat="identity") +
    ggplot2::geom_hline(ggplot2::aes(yintercept=1), color="darkred") +
    ggplot2::geom_errorbar(ggplot2::aes(x=Type, ymin=stat_min, ymax=pmin(stat_max, y.max)), size=line.width, width=.2, position=ggplot2::position_dodge(.9)) +
    ggplot2::scale_y_continuous(limits=c(0, y.max), expand=c(0, 0)) +
    ggplot2::labs(x="", y=y.lab) +
    ggplot2::theme(axis.text.x=ggplot2::element_text(angle = text.angle, hjust = 1, vjust=0.5), legend.position="none",
                   panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank())
}

distanceBetweenTerms <- function(go.df) {
  genes.per.go <- sapply(go.df$geneID, strsplit, "/") %>% setNames(go.df$Description)
  all.go.genes <- unique(unlist(genes.per.go))
  all.gos <- unique(go.df$Description)

  genes.per.go.mat <- matrix(0, length(all.go.genes), length(all.gos)) %>%
    `colnames<-`(all.gos) %>% `rownames<-`(all.go.genes)

  for (i in 1:length(genes.per.go)) {
    genes.per.go.mat[genes.per.go[[i]], go.df$Description[[i]]] <- 1
  }

  return(dist(t(genes.per.go.mat), method="binary"))
}

# Normally, clusterProfiler should cache results, so each request to OrgDb uses local copy.
# But by some reason it didn't work properly, so I cached it by hands
enrichGOOpt <- function (gene, OrgDb, goData, keyType = "ENTREZID", ont = "MF", pvalueCutoff = 0.05,
                         pAdjustMethod = "BH", universe=NULL, qvalueCutoff = 0.2, minGSSize = 10,
                         maxGSSize = 500, readable = FALSE, pool = FALSE) {
  ont %<>% toupper %>% match.arg(c("BP", "CC", "MF"))

  res <- clusterProfiler:::enricher_internal(gene, pvalueCutoff = pvalueCutoff,
                                             pAdjustMethod = pAdjustMethod, universe = universe,
                                             qvalueCutoff = qvalueCutoff, minGSSize = minGSSize,
                                             maxGSSize = maxGSSize, USER_DATA = goData)
  if (is.null(res))
    return(res)

  res@keytype <- keyType
  res@organism <- clusterProfiler:::get_organism(OrgDb)
  if (readable) {
    res <- DOSE::setReadable(res, OrgDb)
  }
  res@ontology <- ont

  return(res)
}

### Embedding

getPathwayClustGenesPerType <- function(pathway.clust.name, go) {
  go %>% dplyr::filter(GOClustName == pathway.clust.name) %$% split(geneID, as.character(Type)) %>%
    lapply(strsplit, "/") %>% lapply(unlist)
}

getGenePathwayMatrix <- function(pathway.clust.name, go, de.info=NULL) {
  genes.per.type <- getPathwayClustGenesPerType(pathway.clust.name, go)
  genes.all <- unlist(genes.per.type) %>% unique()

  gpt.mat <- matrix(0, length(genes.per.type), length(genes.all), dimnames=list(names(genes.per.type), genes.all))
  for (n in names(genes.per.type)) {
    gpt.mat[n, genes.per.type[[n]]] <- if (is.null(de.info)) 1 else de.info[[n]][genes.per.type[[n]],]$stat
  }

  return(gpt.mat)
}

embedPathwayTypesUmap <- function(pathways, gos, k=10, sep="!", return.all=F, ...) {
  gpt.per.path <- pathways %>% setNames(., .) %>% lapply(getPathwayClustGenesPerType, gos)
  # print(names(gpt.per.path) %>% lapply(function(n) gpt.per.path[[n]]) %>% sapply(length))
  # print(names(gpt.per.path) %>% lapply(function(n) paste0(n, sep, names(gpt.per.path[[n]]))) %>% sapply(length))

  genes.per.tp <- names(gpt.per.path) %>% lapply(function(n)
    setNames(gpt.per.path[[n]], paste0(n, sep, names(gpt.per.path[[n]])))) %>% unlist(recursive=F)

  tp.dists <- 1 - sapply(genes.per.tp, function(s1)
    sapply(genes.per.tp, function(s2) length(intersect(s1, s2)) / length(union(s1, s2))))

  k <- min(k, nrow(tp.dists))
  idx <- apply(tp.dists, 1, function(row) order(row)[1:k]) %>% t()
  dist <- apply(tp.dists, 1, function(row) sort(row)[1:k]) %>% t()

  emb <- uwot::umap(NULL, nn_method=list(idx=idx, dist=dist), ...) %>%
    `rownames<-`(rownames(tp.dists))
  # emb <- Rtsne::Rtsne(tp.dists, is_distance=TRUE, perplexity=30, num_threads=1)$Y %>% set_rownames(rownames(tp.dists))

  if (return.all)
    return(list(emb=emb, dists=tp.dists))

  return(emb)
}

estimateHeatmapClusters <- function(heatmap.df, cut.h, min.rows.per.clust=5, distance="minkowski",
                                    method="average", neuron.type.per.type=NULL) {
  cl.row <- heatmap.df %>% `/`(rowSums(.)) %>% dist(distance, p=1) %>% hclust(method=method)
  clusters <- cutree(cl.row, h=cut.h)

  ann.df <- tibble::as_tibble(clusters[cl.row$order], rownames="Type") %>%
    dplyr::mutate(Type=factor(Type, levels=Type)) %>%
    dplyr::mutate(value = ifelse(value %in% names(which(table(value) >= min.rows.per.clust)), value, NA))

  seq.info <- rle(ann.df$value)
  rect.df = data.frame(x1=0.5, x2=ncol(heatmap.df) + 0.5, y2=cumsum(seq.info$lengths)[!is.na(seq.info$values)]) %>%
    dplyr::mutate(y1=y2 - seq.info$lengths[!is.na(seq.info$values)], cluster=seq.info$values[!is.na(seq.info$values)])

  # if (!is.null(neuron.type.per.type)) {
  #   types.per.nt <- colnames(heatmap.df) %>% split(neuron.type.per.type[.])
  #   for (n in names(types.per.nt)) {
  #     rect.df[[paste0("Frac", n)]] <- ann.df %$% split(Type, value) %>% lapply(as.character) %>%
  #       sapply(function(gs) mean(abs(heatmap.df[gs, types.per.nt[[n]]]) > 1e-10))
  #   }
  # }

  return(list(rect=rect.df, ann=ann.df, order=cl.row$labels[cl.row$order]))
}

clusterIndividualGOs <- function(gos, cut.h) {
  gos.split <- gos %>% split(.$Type)

  clusts.per.go <- lapply(gos.split, distanceBetweenTerms) %>% lapply(function(ld)
    if (ncol(as.matrix(ld)) == 1) 1 else {hclust(ld) %>% cutree(h=cut.h)})

  clust.df <- names(gos.split) %>%
    lapply(function(n) mutate(gos.split[[n]], Cluster=clusts.per.go[[n]])) %>%
    Reduce(rbind, .) %>%
    mutate(Cluster=factor(Cluster, levels=c(0, unique(Cluster)))) %>%
    select(Type, Cluster, Description) %>%
    tidyr::spread(Type, Cluster) %>% as.data.frame() %>%
    set_rownames(.$Description) %>% .[, 2:ncol(.)]

  return(clust.df)
}

## Cell type coordination

estimatePathwayWeights <- function(genes.per.pathway) {
  gene.vecs <- genes.per.pathway %>% strsplit("/")
  p.dists <- sapply(gene.vecs, function(gv1) sapply(gene.vecs, function(gv2) jaccardDist(gv1, gv2)))
  return(1 / rowSums(p.dists))
}

jaccardDist <- function(v1, v2){
  length(intersect(v1, v2)) / length(union(v1, v2))
}

jaccardDistWeighted <- function(v1, v2){
  sum(pmin(v1, v2)) / sum(pmax(v1, v2))
}

estimateTypeSimilarityOverPathways <- function(gos.dfs) {
  go.df.weighted <- names(gos.dfs) %>% lapply(function(ng)
    names(gos.dfs[[ng]]) %>% .[sapply(gos.dfs[[ng]], nrow) > 1] %>%
      lapply(function(nt) mutate(gos.dfs[[ng]][[nt]], PathwayType=ng, Weight=estimatePathwayWeights(geneID), Type=nt)) %>%
      dplyr::bind_rows()) %>%
    dplyr::bind_rows()

  path.weights <- go.df.weighted %>% dplyr::select(Description, Type, Weight) %>%
    tidyr::spread(Description, Weight) %>% as.data.frame() %>% magrittr::set_rownames(.$Type) %>% .[, 2:ncol(.)] %>%
    as.matrix()
  path.weights[is.na(path.weights)] <- 0
  cor.mat <- path.weights %>% apply(1, function(r1) apply(., 1, function(r2) jaccardDistWeighted(r1, r2)))

  return(cor.mat)
}

getPathwayGeneExpression <- function(de, cm.collapsed, pathway.name=NULL, go=NULL, genes=NULL, neuron.type.per.type=NULL, type.order=NULL, max.size=4) {
  if (is.null(genes)) {
    if (is.null(go) || is.null(pathway.name)) stop("Either genes or both go and pathway.name must be provided")
    genes <- lapply(go, function(df) strsplit(df$geneID[df$Description == pathway.name], "/")) %>% unlist() %>% unique()
  }

  if (is.null(type.order)) {
    type.order <- names(de)
  } else {
    type.order %<>% intersect(names(de))
  }
  expr.mat <- matrix(0.0, nrow=length(genes), ncol=length(type.order), dimnames=list(genes, type.order))
  stat.mat <- expr.mat
  for (n in type.order) {
    df <- de[[n]][genes, ] %>% as_tibble(rownames="Gene")
    expr.mat[df$Gene, n] <- log10(t(cm.collapsed)[df$Gene, n]) %>% `-`(min(.)) %>% pmin(2.5) %>% `+`(0.01)
    stat.mat[df$Gene, n] <- ifelse(is.na(df$stat), 0.0, df$stat)
  }

  p.df <- as_tibble(expr.mat, rownames="Gene") %>% reshape2::melt(variable.name="Type", value.name="Expression") %>%
    mutate(Z=reshape2::melt(as_tibble(stat.mat, rownames="Gene"))$value,
           ZBin=ifelse(abs(Z) > 3, ">=3", ifelse(abs(Z) > 2, "1.5-3", "<1.5")) %>% factor(levels=c("<1.5", "1.5-3", ">=3")))

  if (length(genes) > 2) {
    gene.order <- as.dist(1 - cor(t(stat.mat))) %>% hclust() %>% .$order
    p.df %<>% mutate(Gene=factor(Gene, levels=rownames(stat.mat)[gene.order]))
  }

  gg <- ggplot(p.df, mapping=aes(x=Type, y=Gene)) +
    geom_point(aes(color=Z, size=Expression, alpha=ZBin)) +
    scale_color_distiller(palette="RdYlBu", limits=c(-7, 7), oob=scales::squish) +
    scale_size_continuous(range=c(0.1, max.size), limits=c(0, 2.52), breaks=c(0.0, 1.25, 2.5)) +
    scale_alpha_manual(values=c(0.05, 0.2, 1.0)) +
    guides(alpha=guide_legend(title="|Z|", override.aes=list(size=2))) +
    theme(panel.background=element_rect(fill="gray"), panel.grid=element_blank(),
          axis.text.x=element_text(angle=60, hjust=1), axis.title=element_blank(),
          legend.key.height=unit(0.1, "in"), legend.margin=margin())

  if (!is.null(neuron.type.per.type)) {
    gg <- gg + gNtVline(neuron.type.per.type)
  }

  return(list(expr=expr.mat, stat=stat.mat, gg=gg))
}
