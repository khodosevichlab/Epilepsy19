#' @importFrom magrittr %<>% %$% %>%

#' @export
Read10x <- function(path, feature.name="features") {
  cm <- Matrix::readMM(file.path(path, "matrix.mtx")) %>% as("dgCMatrix")
  colnames(cm) <- data.table::fread(file.path(path, "barcodes.tsv"), header=F)[[1]]
  rownames(cm) <- data.table::fread(file.path(path, paste0(feature.name, ".tsv")), header=F)[[1]]
  return(cm)
}

#' @export
getColumnwiseCorrelations <- function(m1, m2) {
  m1 <- m1[, intersect(colnames(m1), colnames(m2))]
  m2 <- m2[, intersect(colnames(m1), colnames(m2))]

  return(lineup::corbetw2mat(m1, m2))
}

#' @export
prepareDiffStatDf <- function(pca, sample.per.cell, annotation, neuron.type.per.type, condition.per.sample) {
  mean.pc.per.samp.per.type <- split(names(sample.per.cell), sample.per.cell) %>%
    lapply(function(nsa) split(nsa, annotation[nsa]) %>% .[sapply(., length) > 0] %>%
             sapply(function(nss) colMeans(pca[nss,,drop=F])))

  res.df <- combn(names(mean.pc.per.samp.per.type), 2) %>% apply(2, function(ns)
    data.frame(S1=ns[1], S2=ns[2], value=getColumnwiseCorrelations(mean.pc.per.samp.per.type[[ns[1]]], mean.pc.per.samp.per.type[[ns[2]]]), stringsAsFactors=F) %>%
      tibble::as_tibble(rownames="Type")) %>% Reduce(rbind, .) %>%
    dplyr::mutate(NeuronType=neuron.type.per.type[Type],
                  SameCondition=(condition.per.sample[S1] == condition.per.sample[S2]),
                  Condition=ifelse(!SameCondition, "Between", condition.per.sample[S1]))

  res.df %<>% split(.$Type) %>%
    lapply(function(df) dplyr::mutate(df, DiffStat=(value - mean(value[Condition == "control"], trim=0.4)) / mad(value[Condition == "control"]))) %>%
    Reduce(rbind, .)

  clust.size.per.samp <- split(sample.per.cell, annotation[names(annotation)]) %>%
    lapply(table)

  res.df$Size1 <-  res.df %>%
    apply(1, function(row) clust.size.per.samp[[row["Type"]]][[row["S1"]]])

  res.df$Size2 <-  res.df %>%
    apply(1, function(row) clust.size.per.samp[[row["Type"]]][[row["S2"]]])

  res.df %<>% dplyr::mutate(MinSize=pmin(Size1, Size2), MaxSize=pmax(Size1, Size2), SumSize=Size1 + Size2)
  return(res.df)
}
