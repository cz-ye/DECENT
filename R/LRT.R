#' Likelihood-ratio test
#'
#' Likelihood-ratio test
#'
#' @param data.obs Observed count matrix for endogeneous genes, rows represent genes, columns represent cells
#' @param out Output object from EM for DE model
#' @param out2 Output object from EM for no-DE model
#' @param cell.type A factor or a integer/numeric vector starting from 1 giving cell-type labels
#' @param parallel If \code{TRUE}, run in parallel
#'
#' @return A list including statistics, p-values and parameters for likelihood-ratio test.
#'
#' @useDynLib DECENT
#' @importFrom Rcpp sourceCpp
#'
#' @export
#'
lrTest <- function(data.obs, out, out2, cell.type, parallel) {

  print(paste0('Likelihood ratio test started at ', Sys.time()))

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
      print(paste0(i, ' started'))
      y <- 1 : max(100, qnbinom(0.9999, mu = max(out2$est.mu[i]*out2$est.sf), size = 1/out2$est.disp[i]))

      res1 <- tryCatch(optim(p = c(log(out2$est.pi0[i,1]/(1-out2$est.pi0[i,1])), log(out2$est.mu[i]), 0, -2), fn = loglIBBCpp, y = y,
                       sf = out$est.sf, ct = cell.type, DO_par = DO.par, z = data.obs[i,], lower = -30),
                       error = function(e) {
                         print(paste("numerical problem in DE model for gene", i));
                         NA
                       })

      res2 <- tryCatch(optim(p = c(log(out2$est.pi0[i,1]/(1-out2$est.pi0[i,1])), log(out2$est.mu[i]), -2), fn = loglIBBCpp, y = y,
                    sf = out2$est.sf, ct = rep(1, ncell), DO_par = DO.par, z = data.obs[i, ], lower = -30),
                       error = function(e) {
                         print(paste("numerical problem in noDE model for gene", i));
                         NA
                       })
      print(paste0(i, ' ended'))
      if (is.na(res1) | is.na(res2)) {
        return(rep(0, 5+2*ncelltype))
      } else {
        return(c(res1$p, res2$p, -res1$v, -res2$v))
      }
    }
    par1 <- temp[, 1:(2+ncelltype)]
    par2 <- temp[, (2+ncelltype):(3+2*ncelltype)]
    logl1 <- temp[, 4+2*ncelltype]
    logl2 <- temp[, 5+2*ncelltype]

  } else {
    for(i in 1:ngene) {
      y <- 1 : max(100, qnbinom(0.9999, mu = max(out2$est.mu[i]*out2$est.sf), size = 1/out2$est.disp[i]))
      res1 <- tryCatch(optim(p = c(log(out2$est.pi0[i,1]/(1-out2$est.pi0[i,1])), log(out2$est.mu[i]), 0, -2), fn = loglIBBCpp, y = y,
                       sf = out$est.sf, ct = cell.type, DO_par = DO.par, z = data.obs[i,], lower = -30),
                       error = function(e) {
                         print(paste("numerical problem in DE model for gene", i));
                         NA
                       })
      res2 <- tryCatch(optim(p = c(log(out2$est.pi0[i,1]/(1-out2$est.pi0[i,1])), log(out2$est.mu[i]), -2), fn = loglIBBCpp, y = y,
                    sf = out2$est.sf, ct = rep(1, ncell), DO_par = DO.par, z = data.obs[i, ], lower = -30),
                       error = function(e) {
                         print(paste("numerical problem in noDE model for gene", i));
                         NA
                       })
      if (is.na(res1) | is.na(res2)) {
      } else {
        par1[i, ] <- res1$p
        par2[i, ] <- res2$p
        logl1[i] <- -res1$v
        logl2[i] <- -res2$v
      }
    }
  }
  print(paste0('Likelihood ratio test finished at ', Sys.time()))

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



