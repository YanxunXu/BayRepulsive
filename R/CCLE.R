#' CCLE
#'
#' This dataset was generated from the pure cell line expression data from CCLE .
#' @docType data
#' @name CCLE
#' @usage data(CCLE)
#' @format This is a data frame with three components: (a) the pure cell line expression, (b) the sample proportion, (c) the mixed data.
#' Z: the pure cell line expression of three cancer cell lines -- NCIH524_LUNG, NCIH209_LUNG, SBC5_LUNG. Each line is for one cancer cell line. We selected top 100 differentially expressed gene.
#' W: the sample proportion. Each row is the proportion of one sample. We used this sample porportion to mix 24 mixed samples.
#' DATA: the mixed data. Each line is the gene expression for one sample. DATA = WZ.
#'
#' @examples
#' # import the data
#' data(CCLE)
#' # get the mixed data
#' CCLE$DATA
#' # get the gene expression level of pure cell line
#' CCLE$Z
#' # get the proportion
#' CCLE$W
#' @keywords datasets
NULL
