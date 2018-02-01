#' Likelihood-ratio test
#'
#' Likelihood-ratio test for DE analysis of scRNA-seq data
#'
#' @param data.obs Observed count matrix for endogeneous genes, rows represent genes, columns represent cells
#' @param out Output of fitDE, it contains EM algorithm estimates for models with DE between cell-types.
#' @param out2 Output of fitNoDE, it contains EM algorithm estimates for models without DE between cell-types
#' @param cell.type A factor or a integer/numeric vector starting from 1 giving cell-type labels
#' @param parallel If \code{TRUE}, run in parallel
#'
#' @return A list containing statistics, p-values and parameter estimates for models with and without DE.
#'
#' @import RcppNumerical
#' @useDynLib DECENT
#' @importFrom Rcpp sourceCpp
#'
#' @export
#'
lrTest <- function(data.obs, out, out2, cell.type, parallel) {

  message('Likelihood ratio test started at ', Sys.time())

  if (class(cell.type) %in% c('factor', 'numeric')) {
    cell.type.names <- levels(cell.type)
    cell.type <- as.integer(cell.type)
  } else {
    cell.type.names <- NULL
  }

  ngene <- nrow(out2$est.mu)
  ncelltype <- length(unique(cell.type))
  ncell <- length(out2$est.sf)

  DO.par <- matrix(0,ncell,2)
  DO.par[,1] <- log(out2$CE/(1-out2$CE))

  par1 <- matrix(0, ngene, 2+ncelltype)
  par2 <- matrix(0, ngene, 1+ncelltype)
  logl1<- rep(0, ngene)
  logl2<- rep(0, ngene)

  if (parallel) {
    temp <- foreach (i = 1:ngene, .combine = 'rbind', .packages = c('DECENT')) %dopar% {
      y <- 1 : max(100, qnbinom(0.9999, mu = max(out2$est.mu[i]*out2$est.sf), size = 1/out2$est.disp[i]))

      res2 <- tryCatch(optimLRTCpp(p = c(log(out2$est.pi0[i,1]/(1-out2$est.pi0[i,1])), log(out2$est.mu[i]), -2), y = y,
                                   sf = out2$est.sf, ct = rep(1, ncell), DO_par = DO.par, z = data.obs[i,]),
                       error = function(e) {
                         warning("Numerical problem in noDE model for gene ", i);
                         NA
                       })

      res1 <- tryCatch(optimLRTCpp(p = c(res2$p[1:2], 0, res2$p[3]), y = y,
                                   sf = out2$est.sf, ct = cell.type, DO_par = DO.par, z = data.obs[i, ]),
                       error = function(e) {
                         warning("Numerical problem in DE model for gene ", i);
                         NA
                       })
      if (is.na(res1) | is.na(res2)) {
        return(rep(0, 5+2*ncelltype))
      } else {
        if (res1$status < 0) {
          warning("DE model failed to converge for gene ", i)
        }
        if (res2$status < 0) {
          warning("noDE model failed to converge for gene ", i)
        }
        return(c(res1$par, res2$par, -res1$fopt, -res2$fopt))
      }
    }
    par1 <- temp[, 1:(2+ncelltype)]
    par2 <- temp[, (3+ncelltype):(ncol(temp)-2)]
    logl1 <- temp[, ncol(temp)-1]
    logl2 <- temp[, ncol(temp)]

  } else {
    for(i in 1:ngene) {
      y <- 1 : max(100, qnbinom(0.9999, mu = max(out2$est.mu[i]*out2$est.sf), size = 1/out2$est.disp[i]))
      res2 <- tryCatch(optimLRTCpp(p = c(log(out2$est.pi0[i,1]/(1-out2$est.pi0[i,1])), log(out2$est.mu[i]), -2), y = y,
                                   sf = out2$est.sf, ct = rep(1, ncell), DO_par = DO.par, z = data.obs[i,]),
                       error = function(e) {
                         warning("numerical problem in noDE model for gene ", i);
                         NA
                       })
      res1 <- tryCatch(optimLRTCpp(p = c(res2$p[1:2], 0, res2$p[3]), y = y,
                                   sf = out2$est.sf, ct = cell.type, DO_par = DO.par, z = data.obs[i, ]),
                       error = function(e) {
                         warning("numerical problem in DE model for gene ", i);
                         NA
                       })
      if (is.na(res1) | is.na(res2)) {
      } else {
        if (res1$status < 0) {
          warning("DE model failed to converge for gene ", i)
        }
        if (res2$status < 0) {
          warning("noDE model failed to converge for gene ", i)
        }
        par1[i, ] <- res1$par
        par2[i, ] <- res2$par
        logl1[i] <- -res1$fopt
        logl2[i] <- -res2$fopt
      }
    }
  }
  message('Likelihood ratio test finished at ', Sys.time())

  output <- list()
  lrt.stat <- 2*(logl1 - logl2)
  lrt.stat <- ifelse(lrt.stat<0, 0, lrt.stat)
  pval <- exp(pchisq(lrt.stat, df = 1, lower.tail = FALSE, log = TRUE))
  names(lrt.stat) <- rownames(data.obs)
  names(pval) <- rownames(data.obs)
  rownames(par1) <- rownames(data.obs)
  rownames(par2) <- rownames(data.obs)

  output[['stat']] <- lrt.stat
  output[['pval']] <- pval
  output[['par.DE']] <- par1
  output[['par.noDE']] <- par2

  return(output)
}
