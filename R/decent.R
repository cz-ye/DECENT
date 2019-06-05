#' DECENT
#'
#' Differential Expression with Capture Efficiency adjustmeNT
#'
#' @param data.obs Observed count matrix for endogeneous genes, rows represent genes, columns represent cells.
#' @param spike Observed count matrix for spike-ins, rows represent spike-ins, columns represent cells. Only needed if spikes = \code{TRUE}).
#' @param spike.conc A vector of theoretical count for each spike-in in one cell (ONLY needed if spikes = \code{TRUE}).
#' @param CE.range A two-element vector of the lower limit and upper limit for the estimated range of
#' capture efficiencies (ONLY needed if spikes = \code{FALSE}, default [0.02, 0.10]).
#' @param tau.init initial estimates (intcp,slope) that link Beta-Binomial dispersion parameter to the mean expression.
#' @param tau.global whether to use the same tau parameters across cell. Default TRUE
#' @param tau.est Methods to estimate tau parameters. The default 'endo' corresponds to using endogeneous genes. Other options
#' are 'spikes' that corresponds to using spike-ins and 'none', which means tau.init is not further estimated.
#' @param X An R model formula with covariates of interest (cell-type) as factor on the righthand side.
#' @param W An R model formula with other covariates to adjust DE analysis on the righthand side. Default NULL
#' @param use.spikes If \code{TRUE}, use spike-ins to estimate capture efficiencies.
#' @param normalize Method for estimating size factors, either 'ML' (maximum likelihood, Ye et al., 2017) or 'TMM' (Robinson et al., 2010).
#' @param GQ.approx If \code{TRUE}, use Gaussian-Quadrature approximation to speed up E-step.
#' @param parallel If \code{TRUE}, run DECENT in parallel.
#' @param n.cores Number of CPU cores to use, default is all (ONLY if parallel=\code{TRUE}).
#' @param s.imputed If \code{TRUE}, save the single imputed data matrix under the output diretory.
#' @param E.imputed If \code{TRUE}, save the mean imputed data matrix under the output diretory.
#' @param dir Directory to save all outputs, including EM algorithm estimates of no-DE model and the LRT output.
#'
#' @return A list containing the result of differential expression analysis, with the following components: stat,pval,par.DE and par.noDE.
#'
#' @import parallel
#' @import foreach
#' @import doParallel
#'
#' @export
decent <- function (data.obs, X, W=NULL, spikes = NULL, spike.conc = NULL, CE.range = c(0.02, 0.1), tau.init = c(-5,0),
                    tau.global = T, tau.est = "endo", use.spikes = FALSE, normalize = 'ML', GQ.approx = TRUE, 
                    maxit = 30, parallel = T, n.cores = 0, s.imputed = F, E.imputed = F, dir = './') {

  X <- model.matrix(X) 
  if(!is.null(W)) 
   W <- model.matrix(W)[,-1, drop = F]

  #if (length(unique(cell.type)) != 2) stop('Number of groups (cell types) should be two.')
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
    message('Running in parallel, ', n.cores, ' cores used.')
  }

  # Fit no-DE model
  out.noDE <- fitNoDE(data.obs, spikes, spike.conc, use.spikes, CE.range, tau.init, tau.global, tau.est, 
                      normalize, GQ.approx, maxit, parallel)
  saveRDS(out.noDE, paste0(dir, '/decent.noDE.rds'))


  # Likelihood-ratio test
  out <- lrTest(data.obs, out = out.noDE, X=X, W=W, tau= cbind(out.noDE$tau0, out.noDE$tau1), parallel)
  #
  saveRDS(out, paste0(dir, '/decent.lrt.rds'))

  # get imputed data
  if (E.imputed) {
    data.simp <- getImputed(data.obs, out$par.DE, out.noDE$est.sf, out.noDE$CE, cbind(X, W), tau = cbind(out.noDE$tau0,out.noDE$tau1), parallel)
    saveRDS(data.simp, paste0(dir, '/mean.imputed.rds'))
  }
  if (s.imputed) {
    data.eimp <- getSImputed(data.obs, out$par.DE, out.noDE$est.sf, out.noDE$CE, cbind(X, W), tau = cbind(out.noDE$tau0,out.noDE$tau1), parallel)
    saveRDS(data.eimp, paste0(dir, '/single.imputed.rds'))
  }

  # Generate table
  de.table <- data.frame(gene = rownames(out.noDE$est.mu),
                         logfc = out$par.DE[,3],
                         pvalue = out$pval,
                         stat = sign(out$par.DE[,3])*sqrt(out$stat))

  return(de.table)
}
