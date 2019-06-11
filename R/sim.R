#' Simulated data
#'
#' A simulated dataset for demonstration. 2,944 genes and 500 cells from two groups (224 vs 276).
#'
#' @docType data
#'
#' @format A list containing five objects:
#' \describe{
#'   \item{\code{data.obs}}{
#'     a numeric matrix with 2,944 genes and 500 cells: the observed count matrix for endogeneous genes
#'   }
#'   \item{\code{cell.type}}{
#'     a numeric vector of length 500: the binary indicator of cellular groups
#'   }
#'   \item{\code{sp.obs}}{
#'     a numeric matrix of 50 spike-ins and 500 cells: the observed count matrix for spike-in molecules
#'   }
#'   \item{\code{sp.true}}{
#'     a numeric vector of length 50: the true count of spike-ins added into each cell
#'   }
#'   \item{\code{DE.gene}}{
#'     a logical vector of length 2,944: whether the gene is actually DE or not.
#'   }
#' }
#'
"sim"
