#' DECENT
#'
#' Differential Expression with Capture Efficiency adjustmeNT
#'
#' @param data.obs Observed count matrix for endogeneous genes, rows represent genes, columns represent cells.
#' @param spike Observed count matrix for spike-ins, rows represent spike-ins, columns represent cells
#' (ONLY if spikes = \code{TRUE}).
#' @param spike.conc A vector of theoretical count for each spike-in in one cell (ONLY if spikes = \code{TRUE}).
#' @param CE.range A two-element vector of the lower limit and upper limit for the estimated
#' capture efficiencies (ONLY if spikes = \code{FALSE}, default [0.02, 0.10]).
#' @param cell.type A factor or a integer/numeric vector starting from 1 providing cell-type labels.
#' @param use.spikes If \code{TRUE}, use spike-ins to estimate capture efficiencies.
#' @param normalize Normalization method, either 'ML' (maximum likelihood) or 'TMM' (Robinson et al., 2010).
#' @param GQ.approx If \code{TRUE}, use GQ approximation to speed up E-step.
#' @param parallel If \code{TRUE}, run DECENT in parallel.
#' @param n.cores Number of CPU cores to use, default is all (ONLY if \code{TRUE}).
#' @param imputed If \code{TRUE}, include imputed data matrix in the output (Not available now).
#' @param dir Output directory for the DE model estimates, no-DE model estimates as well as the LRT output.
#'
#' @return A data frame containing the result of differential expression analysis.
#'
#' @import parallel
#' @import foreach
#' @import doParallel
#'
#' @export
decent <- function (data.obs, cell.type, spikes = NULL, spike.conc = NULL, CE.range = c(0.02, 0.1),
                     use.spikes = FALSE, normalize = 'ML', GQ.approx = TRUE, maxit = 30,
                     parallel = T, n.cores = 0, imputed = F, dir = './') {

  if (length(unique(cell.type)) != 2) stop('Number of groups (cell types) should be two.')
  if (CE.range[1] < 0 | CE.range[1] > CE.range[2] | CE.range[2] > 1) stop('CE.range invalid.')

  if (parallel) {
    t.cores <- detectCores() - 1
    if (n.cores > 0 & n.cores < t.cores) {
      cl <- makeCluster(n.cores)
    } else {
      cl <- makeCluster(t.cores, outfile=paste0(dir,'/out.txt'))
      n.cores <- t.cores
    }
    registerDoParallel(cl)
    print(paste0('Running in parallel, ', n.cores, ' cores used.'))
  }

  # Fit DE model
  out.DE <- fitDE(data.obs, cell.type, spikes, spike.conc, CE.range, use.spikes, normalize,
                  GQ.approx, maxit, parallel)
  saveRDS(out.DE, paste0(dir, '/decent.DE.rds'))

  # Fit no-DE model
  out.noDE <- fitNoDE(data.obs, CE = out.DE$CE, normalize, GQ.approx, maxit, parallel)
  saveRDS(out.noDE, paste0(dir, '/decent.noDE.rds'))

  # Likelihood-ratio test
  out <- lrTest(data.obs, out = out.DE, out2 = out.noDE, cell.type, parallel)
  saveRDS(out, paste0(dir, '/decent.lrt.rds'))

  # get imputed data
  # TODO
  if (imputed) {
  }

  # Generate table
  de.table <- data.frame(gene = rownames(out.DE$est.mu),
                         logfc = out$par.DE[,3],
                         pvalue = out$pval,
                         stat = sign(out$par.DE[,3])*sqrt(out$stat))

  return(de.table)
}
