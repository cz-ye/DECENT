#' Simulated data
#'
#' A simulated dataset for demonstration. 1,000 genes and 200 cells, with 100 cells in each group.
#'
#' @docType data
#'
#' @format An list of five objects:
#' \describe{
#'   \item{data.obs}{
#'     a numeric matrix with 1,000 genes and 200 cells: observed count matrix for endogeneous genes
#'   }
#'   \item{cell.type}{
#'     a numeric vector of length 200: cell-type information
#'   }
#'   \item{ercc.obs}{
#'     a numeric matrix with 50 spike-ins and 200 cells: observed count matrix for spike-in molecules
#'   }
#'   \item{ercc.true}{
#'     a numeric vector of length 50: true count of spike-ins added into each cell
#'   }
#'   \item{DE.gene}{
#'     a logical vector of length 1,000: whether the gene is actually DE or not.
#'   }
#' }
#' 
"sim"