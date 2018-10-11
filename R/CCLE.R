#' The CCLE Dataset
#'
#' This dataset is the pure cell line expression from \emph{Cancer Cell Line Encyclopedia} (CCLE) \insertCite{barretina2012cancer}{BayRepulsive}.
#' @docType data
#' @name CCLE
#' @usage data(CCLE)
#' @format From the CCLE dataset, we randomly chose three lung-related cell lines, NCIH524_LUNG, NCIH209_LUNG and SBC5_LUNG.
#' We then selected the top 100 differentially expressed genes by ranking the standard deviations of genes across pure samples. The expression levels of these 100 genes in the selected three cell lines compose our simluated \code{Z} matrix.
#' This is a data frame with one components: the pure cell line expression, \code{Z}.
#'
#' @examples
#' # import the data
#' data(CCLE)

#' # get the gene expression level of pure cell line
#' CCLE$Z
#' # get the name of cell lines included in the dataset
#' colnames(CCLE$Z)
#' @keywords datasets
#'
#' @references
#' \insertRef{barretina2012cancer}{BayRepulsive}
NULL
