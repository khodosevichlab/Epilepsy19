plotGGHeatmap <- function(df, legend.title="") {
  df %<>% tibble::as_tibble(rownames=".id.") %>%
    reshape2::melt(id.vars=".id.") %>%
    dplyr::mutate(.id.=factor(.id., levels=rev(rownames(df))), variable=factor(variable, levels=colnames(df)))

  ggplot(df) + geom_tile(aes(x=variable, y=.id., fill=value), colour = "grey50") +
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
          axis.text=element_text(size=8), axis.ticks=element_blank(), axis.title=element_blank()) +
    guides(fill=guide_colorbar(title=legend.title, title.position="left", title.theme=element_text(angle=90, hjust=0.5))) +
    scale_y_discrete(position="right", expand=c(0, 0)) +
    scale_x_discrete(expand=c(0, 0))
}

estimateFoldChanges <- function(cm, markers, annotation) {
  annotation %<>% .[intersect(names(.), rownames(cm))]
  cm %<>% .[names(annotation), markers]
  mean.mtx <- split(1:length(annotation), annotation) %>%
    sapply(function(ids) Matrix::colMeans(cm[ids,]))

  fold.change.mtx <- 1:ncol(mean.mtx) %>%
    sapply(function(i) log2(1e-10 + mean.mtx[,i]) - log2(1e-10 + rowMeans(mean.mtx[,-i]))) %>%
    `colnames<-`(colnames(mean.mtx))

  return(fold.change.mtx)
}

getNGenes <- function(cm) {
  cm@x <- (cm@x > 0) * 1
  return(Matrix::colSums(Matrix::t(cm)))
}