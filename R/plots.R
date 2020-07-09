#' @import ggplot2

plotSimilarityStat <- function(plot.df, y.lims, fill.color, text.angle=60, size=0.3, width=0.25) {
  ggplot(dplyr::filter(plot.df, !SameCondition), aes(x=Type, y=DiffStat)) +
    geom_boxplot(fill=fill.color, outlier.alpha=0, show.legend=F) +
    geom_hline(aes(yintercept=split(DiffStat, Type) %>% sapply(median) %>% median()), color="darkred", size=1) +
    geom_hline(aes(yintercept=0), color="darkgreen", size=1) +
    ggbeeswarm::geom_quasirandom(size=size, width=width) +
    labs(x="", y="Similarity score") + scale_y_continuous(limits=y.lims, expand=c(0, 0)) +
    ggrastr::theme_pdf() +
    theme(axis.text.x=element_text(angle=text.angle, vjust=ifelse(text.angle < 89, 1, 0.5), hjust=1, size=9))
}

plotCellFracPerType <- function(cell.num.df, y.max, legend.pos=c(1, 1), text.angle=90) {
  cell.num.df$Condition %<>% as.factor()
  ggplot(cell.num.df, aes(x=Type, fill=Condition, y=Frac * 100)) +
    geom_boxplot(outlier.alpha=0.0, alpha=0.5) +
    ggbeeswarm::geom_beeswarm(dodge.width=1, size=0.1) +
    ggbeeswarm::geom_beeswarm(aes(alpha=Condition), data=filter(cell.num.df, (Sample == "C1") | (Sample == "E1")),
                  color="#d95f02", size=1, shape=17, show.legend=F, dodge.width=1) +
    scale_alpha_manual(values=c(1, 0)) +
    scale_fill_manual("", values=c("#0571b0", "#ca0020")) +
    scale_y_continuous(limits=c(0, y.max), expand=c(0, 0.1)) +
    labs(x="", y="% of nuclei") +
    theme(legend.position=c(1, 1.05), legend.justification=c(1, 1),
          legend.background=element_blank(), legend.title=element_blank(),
          axis.text.x=element_text(angle=text.angle, hjust=1, vjust=0.5, size=9),
          panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(),
          panel.border = element_blank(), axis.line = element_line())
}

### Gene plots

prepareExpressionDfs <- function(genes, samples, annotation, neuron.types, condition.per.sample) {
  genes %>% setNames(., .) %>%
    lapply(function(gn) {
      names(samples) %>% lapply(function(n)
        data.frame(Expression=samples[[n]]$counts[,gn]) %>% tibble::as_tibble(rownames="Cell") %>%
          dplyr::mutate(Type=as.character(annotation[Cell]), Sample=n)) %>%
        dplyr::bind_rows() %>%
        dplyr::mutate(NeuronType=neuron.types[Type], Condition=condition.per.sample[Sample])
    })
}

averageGeneDfBySamples <- function(gene.df, trim=0.3) {
  gene.df %>% dplyr::group_by(Sample, Type, NeuronType, Condition) %>%
    dplyr::summarise(Med=median(Expression), UQ=quantile(Expression, 0.75), LQ=quantile(Expression, 0.25),
                     Q90=quantile(Expression, 0.9), Q10=quantile(Expression, 0.1)) %>%
    dplyr::group_by(Type, NeuronType, Condition) %>%
    dplyr::summarise(LQ=mean(LQ, trim=trim), UQ=mean(UQ, trim=trim), Med=mean(Med, trim=trim),
                     Q90=mean(Q90, trim=trim), Q10=mean(Q10, trim=trim))
}

plotExpression <- function(gene.df, geom=geom_violin, text.size=8, ...) {
  elem.text <- element_text(size=text.size)
  c("Inhibitory", "Excitatory") %>% setNames(., .) %>% lapply(function(n) {
    ggplot(dplyr::filter(gene.df, NeuronType==n)) +
      geom(aes(x=Type, y=Expression, fill=Condition), ...) +
      labs(x="", y="Expression") +
      theme(axis.text.x=element_text(angle=45, hjust=1)) +
      ggrastr::theme_pdf(legend.pos=c(0, 1)) +
      theme(axis.text=elem.text, axis.title=elem.text, legend.text=elem.text, legend.title=elem.text)
  })
}

plotExpressionAveraged <- function(gene.df.sum, text.size=8, text.angle=90, text.vjust=0.5, ...) {
  elem.text <- element_text(size=text.size)
  c("Inhibitory", "Excitatory") %>% setNames(., .) %>% lapply(function(n) {
    ggplot(dplyr::filter(gene.df.sum, NeuronType==n), aes(x=Type, group=Condition, color=Condition)) +
      geom_point(aes(y=Med), position=position_dodge(width=0.5)) +
      geom_errorbar(aes(ymin=LQ, ymax=UQ), position=position_dodge(width=0.5), width=0.5) +
      labs(x="", y="Expression") +
      theme(axis.text.x=element_text(angle=text.angle, hjust=1, vjust=text.vjust)) +
      ggrastr::theme_pdf(legend.pos=c(0, 1)) +
      theme(axis.text=elem.text, axis.title=elem.text, legend.text=elem.text, legend.title=elem.text, ...)
  })
}

## GO

#' Plot GO Heatmap
#'
#' @param df columns are cell types, rows are pathways and values are plotted with colors
plotGOHeatmap <- function(df, color.per.type=NULL, row.order=NULL, col.order=F,
                          legend.position="right", legend.justification=c(1, 1), legend.title.size=12,
                          legend.key.width=unit(8, "pt"), legend.title="-log10(p-value)") {
  if (is.null(row.order)) {
    row.order <- rownames(df)[dist(df) %>% hclust() %>% .$order]
  } else if (is.logical(row.order) && row.order) {
    row.order <- rownames(df)
  }

  if (is.null(col.order)) {
    col.order <- colnames(df)[dist(df) %>% hclust() %>% .$order]
  } else if (is.logical(col.order) && col.order) {
    col.order <- colnames(df)
  }

  df %<>% as_tibble(rownames="GO") %>%
    melt(id.vars="GO", variable.name="Type", value.name="p.value")

  if (!is.logical(row.order)) {
    df %<>% mutate(GO=factor(GO, levels=row.order))
  }

  if (!is.logical(col.order)) {
    df %<>% mutate(Type=factor(Type, levels=col.order))
  }

  if (is.null(color.per.type)) {
    color.per.type <- "black"
  } else {
    color.per.type <- color.per.type[levels(df$Type)]
  }

  ggplot(df) + geom_tile(aes(x=Type, y=GO, fill=pmin(p.value, 10)), colour = "grey50") +
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, color=color.per.type),
          axis.text=element_text(size=8), axis.ticks=element_blank(), axis.title=element_blank()) +
    scale_fill_distiller(palette="RdYlBu", limits=c(0, 10)) +
    guides(fill=guide_colorbar(title=legend.title, title.position="left", title.theme=element_text(angle=90, hjust=0.5, size=legend.title.size))) +
    scale_y_discrete(position="right", expand=c(0, 0)) +
    scale_x_discrete(expand=c(0, 0)) +
    theme(legend.position=legend.position, legend.justification=legend.justification,
          legend.key.width=legend.key.width, legend.background=element_blank())
}

prepareGOHeatmapClustered <- function(gos, type.order, color.per.type="red", coloring.go.threshold=20, ...) {
  go.df <- gos %>% mutate(p.adjust=-log10(p.adjust)) %>%
    dplyr::select(Type, p.adjust, Description) %>%
    tidyr::spread(Type, p.adjust) %>% as.data.frame() %>%
    magrittr::set_rownames(.$Description) %>% .[, 2:ncol(.)] %>%
    .[, type.order[type.order %in% colnames(.)]]

  go.df[is.na(go.df)] <- 0

  plot.df <- go.df %>% .[rowSums(. > 0.1) > 1,] %>% .[colSums(.) > 0]

  clust.info <- estimateHeatmapClusters(plot.df, ...)
  ann.tile <- ggplot(clust.info$ann) +
    geom_tile(aes(x="1", y=Type, fill=as.character(value))) +
    theme(legend.position="none", axis.text=element_blank(), axis.title=element_blank(), axis.ticks=element_blank()) +
    scale_x_discrete(expand=c(0, 0))

  if (length(color.per.type) > 0) {
    color.per.type[colnames(plot.df)[colSums(plot.df > 0) < coloring.go.threshold]] <- "black"
  }

  gg.go.heatmap <- plotGOHeatmap(plot.df, legend.position=c(1.9, 1), row.order=clust.info$order, color.per.type=color.per.type) +
    gNtVline(neuron_type_per_type) +
    geom_rect(aes(xmin=x1, xmax=x2, ymin=y1 + 0.5, ymax=y2 + 0.5), clust.info$rect, fill="transparent", color="black") +
    theme(plot.margin=margin())

  clust.info$gg <- cowplot::plot_grid(
    ann.tile + theme(plot.margin=margin()),
    gg.go.heatmap,
    ncol=2, align="h", rel_widths=c(0.02, 1)
  )

  clust.info$df <- plot.df
  clust.info$color.per.type <- color.per.type

  return(clust.info)
}

plotNumberOfClustersPerHeight <- function(clusters, step=0.01) {
  h.vals <- seq(0.0, 1.0, step)
  n.clusts.per.h <- sapply(h.vals, function(h) length(unique(cutree(clusters, h=h))))
  ggplot(data.frame(Height=h.vals, NClusters=n.clusts.per.h), aes(x=Height, y=NClusters)) +
    geom_line() +
    scale_y_continuous("Num. of clusters", sec.axis = sec_axis(~ . * 100 / length(clusters$labels), name = "% of clusters"))
}

plotTypeSimilarityOverPathways <- function(cor.mat, cut.height, color.per.type=NULL,
                                            max.val=max(cor.mat[upper.tri(cor.mat)]), min.rows.per.clust=3) {
  if (min.rows.per.clust < 3)
    warning("The visualization is unstable for min.rows.per.clust < 3")

  clusts <- dist(cor.mat) %>% hclust()
  type.order <- clusts %$% labels[order]
  clust.labels <- cutree(clusts, h=cut.height)[type.order]

  clust.labels[clust.labels %in% names(which(table(clust.labels) < min.rows.per.clust))] <- max(clust.labels) + 1

  if (is.null(color.per.type)) {
    color.per.type <- "black"
  } else {
    color.per.type <- color.per.type[colnames(cor.mat)][type.order]
  }

  clust.lengths <- rle(clust.labels)$lengths
  line.df <- data.frame(x=cumsum(clust.lengths)[clust.lengths > 1] + 0.5)
  clust.lengths <- rev(rle(clust.labels)$lengths)
  line.df2 <- data.frame(x=cumsum(clust.lengths)[clust.lengths > 1] + 0.5)

  plotGOHeatmap(cor.mat, color.per.type=color.per.type, row.order=type.order, col.order=rev(type.order),
                legend.title="GO similarity") +
    scale_fill_distiller(palette="RdYlBu", limits=c(0, max.val)) +
    geom_vline(aes(xintercept=x), line.df2) +
    geom_hline(aes(yintercept=x), line.df) +
    theme(axis.text=element_text(color=color.per.type))
}

plotIndividualClustersPerGO <- function(clust.df, row.order) {
  plot.df <- clust.df %>% as_tibble(rownames="GO") %>%
    melt(id.vars="GO", variable.name="Type", value.name="Cluster") %>%
    mutate(Cluster=as.factor(Cluster), GO=factor(GO, levels=row.order)) %>%
    filter(!is.na(GO))

  max.clust <- as.matrix(clust.df) %>% .[!is.na(.)] %>% as.integer() %>% max()
  palette <- colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(max.clust) %>%
    setNames(sample(1:max.clust))

  ggplot(plot.df) + geom_tile(aes(x=Type, y=GO, fill=Cluster)) +
    scale_fill_manual(values=palette) + theme(legend.position="nothing") +
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
          axis.text=element_text(size=8), axis.ticks=element_blank(),
          axis.title=element_blank())
}

gNtVline <- function(neuron.type.per.type, color=alpha("black", 0.6), linetype=2, size=0.5) {
  geom_vline(aes(xintercept=sum(neuron.type.per.type[unique(as.character(Type))] == "Excitatory") + 0.5),
             color=color, linetype=linetype, size=size)
}
