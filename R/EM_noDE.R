#' Fitting no-DE model with EM algorithm
#'
#' Fit the DECENT model with no-DE assumption
#'
#' @param data.obs Observed count matrix for endogeneous genes, rows represent genes, columns represent cells.
#' @param CE A  vector containing capture efficiencies for all cells.
#' @param normalize Method for estimating size factors, either 'ML' (maximum likelihood, Ye et al., 2017) or 'TMM' (Robinson et al., 2010).
#' @param GQ.approx If \code{TRUE}, use Gaussian-Quadrature approximation to speed up E-step.
#' @param maxit maximum number of iterations for EM algorithm.
#' @param parallel If \code{TRUE}, run DECENT in parallel.
#'
#' @return A list of no DE model estimates
#' @examples
#'
#' @import MASS
#' @import ZIM
#' @import statmod
#' @import edgeR
#'
#' @export
fitNoDE <- function (data.obs, CE, normalize, GQ.approx, maxit, parallel) {

  ncell <- ncol(data.obs)
  ngene <- nrow(data.obs)

  cell.type <- rep(1, ncell)
  cell.type.names <- NULL
  ncelltype <- length(unique(cell.type))


  DO.coef <- matrix(0, ncell, 2)
  DO.coef[, 1] <- log(CE/(1-CE))

  # Initialize size factor
  data.obs.adj <- data.obs %*% diag(1/CE)
  est.sf <- apply(data.obs.adj, 2, mean, trim = 0.025)
  est.sf <- est.sf/mean(est.sf)

  # Initialize other ZINB parameters
  est.mu <- matrix(0, ngene, ncelltype)
  for (K in 1:ncelltype) {
    est.mu[, K] <- 2 + apply(data.obs[, cell.type == K], 1, quantile, prob = 0.9)
  }
  est.disp  <- rbeta(ngene, 0.1 ,0.6)
  est.pi0   <- matrix(0, ngene, ncelltype)
  est.pi0[, 1]   <-  rbeta(ngene, 3, 15)
  if(ncelltype>1) {
    for (K in 2:ncelltype) {
       est.pi0[, K] <- est.pi0[, 1]
    }
  }
  est.dmu <- rgamma(ngene, 5, 5)

  # Initialize other variables
  loglik.vec <- rep(0, maxit)
  data.imp <- data.obs
  PE <- matrix(0, ngene, ncell)
  iter <- 1
  converge <- FALSE
  if (GQ.approx) gq <- gauss.quad(64, kind = 'legendre') else gq <- NULL

  message('No-DE model fitting started at ', Sys.time())
  # Begin EM algorithm
  for (iter in 1:maxit) {

    # E-step gene by gene
    if (parallel) {
      if (!GQ.approx) {
        temp <- foreach (i = 1:ngene, .combine = 'rbind', .packages = c('DECENT')) %dopar% {
          out <- EstepByGene(par = DO.coef, z = data.obs[i, ], z.ind = data.obs[i, ] > 0, sf = est.sf,
                              pi0 = est.pi0[i, cell.type], mu = est.mu[i, cell.type], disp = est.disp[i])
          return(c(ifelse(is.na(out$EYZ0E1),data.obs[i, ],out$EYZ0E1), 1 - out$PE0Z0))
        }
      } else {
        temp <- foreach (i = 1:ngene, .combine = 'rbind', .packages = c('MASS','ZIM', 'DECENT')) %dopar% {
          out <- Estep2ByGene(par = DO.coef,z = data.obs[i, ], z.ind = data.obs[i, ]>0, sf = est.sf,
                               pi0 = est.pi0[i, cell.type], mu = est.mu[i, cell.type], disp = est.disp[i], GQ.object = gq)
          return(c(ifelse(is.na(out$EYZ0E1),data.obs[i, ],out$EYZ0E1), 1 - out$PE0Z0))
        }
      }
      data.imp <- temp[, 1:ncell]
      PE <- temp[, (ncell+1):(2*ncell)]

    } else {
      if (!GQ.approx) {
        for (i in 1:ngene) {
          # use E-step with expected value evaluated using GQ integral
          out <- EstepByGene(par = DO.coef, z = data.obs[i, ], z.ind = data.obs[i, ] > 0, sf = est.sf,
                             pi0 = est.pi0[i, cell.type], mu = est.mu[i, cell.type], disp = est.disp[i])
          data.imp[i, ] <- ifelse(is.na(out$EYZ0E1),data.obs[i, ],out$EYZ0E1)
          PE[i, ]<- 1 - out$PE0Z0
        }
      } else {
        for (i in 1:ngene) {
          out <- Estep2ByGene(par = DO.coef,z = data.obs[i, ], z.ind = data.obs[i, ]>0, sf = est.sf,
                              pi0 = est.pi0[i, cell.type], mu = est.mu[i, cell.type], disp = est.disp[i], GQ.object = gq)
          data.imp[i, ] <- ifelse(is.na(out$EYZ0E1), data.obs[i, ], out$EYZ0E1)
          PE[i, ] <- 1 - out$PE0Z0
        }
      }
    }

    # M-step 1: Update SF
    data.imp = as.matrix(data.imp)
    data.imp2 = data.imp*PE

    # M-step 2: Estimate SF by maximum-likelihood
    if (normalize == 'ML') {
      for (i in 1:ncell) {
        p0 <- est.pi0[, cell.type[i]] +
          (1 - est.pi0[, cell.type[i]])*dnbinom(0, mu = est.sf[i]*est.mu[, cell.type[i]], size = 1/est.disp)
        w  <- ((p0 - est.pi0[, cell.type[i]]*(1-PE[, i]))*(1 - est.pi0[, cell.type[i]]))/p0
        est.sf[i]  <- sum(data.imp2[, i], na.rm=T)/sum(w, na.rm=T)
      }
    } else if (normalize == 'TMM') {
      tmm <- calcNormFactors(data.imp2)
      est.sf <- colSums(data.imp2)*tmm
    } else {
      stop('Normalization method should either be "ML" or "TMM"')
    }
    est.sf <- est.sf/mean(est.sf)

    # M-step 3: Update pi_0, mu and phi, gene-by-gene
    loglik <- rep(0, ngene)
    if (parallel) {
      temp <- foreach (i = 1:ngene, .combine = 'rbind', .packages = c('ZIM', 'DECENT')) %dopar% {
        if (sum(data.imp[i, ])>sum(data.obs[i, ])) {
          prop0 <- ifelse(est.pi0[i, 1] < 0.01,
                          0.025, ifelse(est.pi0[i, 1] > 0.99, 0.975, est.pi0[i,1]))
          out <- optim(par = c(log(prop0/(1-prop0)), log(mean(data.imp[i, ], na.rm=T)), rep(0, ncelltype-1), -2),
                       fn = MstepNB, y = data.imp[i, ], sf = est.sf, status = PE[i, ], ct = cell.type,lower=-30)#,
                       #gr = zinbGrad, method = 'L-BFGS-B')
          new.pi0 <- rep(1/(1 + exp(-out$p[1])), ncelltype)
          new.mu <- exp(out$p[2])
          new.disp <- exp(out$p[length(out$p)])
          if(!GQ.approx){
            new.loglik <- -loglI(p = out$p, sf = est.sf, ct = cell.type, DO.par = DO.coef, z = data.obs[i, ])
          } else {
            new.loglik <- -loglI2(p = out$p, sf = est.sf, ct = cell.type, DO.par = DO.coef, z = data.obs[i, ],
                                    GQ.object = gq)
          }
          return(c(new.pi0, new.mu, new.disp, new.loglik))
        } else {
          return(c(est.pi0[i, ], est.mu[i, ], est.disp[i], loglik[i]))
        }
      }
      est.pi0 <- as.matrix(temp[, 1:ncelltype])
      est.mu <- as.matrix(temp[, (ncelltype+1):(2*ncelltype)])
      est.disp <- temp[, 2*ncelltype+1]
      loglik <- temp[, 2*ncelltype+2]

    } else {
      for (i in 1:ngene) {
        if (sum(data.imp[i, ])>sum(data.obs[i, ])) {
          prop0 <- ifelse(est.pi0[i, 1] < 0.01,
                          0.025, ifelse(est.pi0[i, 1] > 0.99, 0.975, est.pi0[i,1]))
          out <- optim(par = c(log(prop0/(1-prop0)), log(mean(data.imp[i, ], na.rm=T)), rep(0, ncelltype-1), -2),
                       fn = MstepNB, y = data.imp[i, ], sf = est.sf, status = PE[i, ], ct = cell.type,lower=-30)#,
                       #gr = zinbGrad, method = 'L-BFGS-B')
          est.pi0[i, ] <- rep(1/(1 + exp(-out$p[1])), ncelltype)
          est.mu[i, ]  <- exp(out$p[2])
          est.disp[i] <- exp(out$p[length(out$p)])
          if (!GQ.approx) {
            loglik[i] <- -loglI(p = out$p, sf = est.sf, ct = cell.type, DO.par = DO.coef, z = data.obs[i, ])
          } else {
            loglik[i] <- -loglI2(p = out$p, sf = est.sf, ct = cell.type, DO.par = DO.coef, z = data.obs[i, ],
                                    GQ.object = gq)
          }
        }
      }
    }

    loglik.vec[iter] <- sum(loglik)
    message('EM iteration ', iter, ' finished at ', Sys.time(), '  Log-likelihood: ', loglik.vec[iter])

    if (iter > 5) {
      if ( (loglik.vec[iter] - loglik.vec[iter-1])/abs(loglik.vec[iter-1]) < 1e-03 | iter == maxit ) converge <- TRUE
    }

    if (converge) {
  # NOT CALCULATING SE
      break
    }
  } # end of EM loop

  message('No-DE model fitting finished at ', Sys.time())

  # Output
  if(ncelltype>1) {
   rownames(est.mu) <- rownames(data.obs)
   rownames(est.pi0) <- rownames(data.obs)
  }
  if(ncelltype==1) {
   names(est.mu) <- rownames(data.obs)
   names(est.pi0) <- rownames(data.obs)
  }
  names(est.disp) <- rownames(data.obs)
  names(est.sf) <- colnames(data.obs)

  if (!is.null(cell.type.names)){
    colnames(est.mu) <- cell.type.names
    colnames(est.pi0) <- cell.type.names
  }

  output <- list()
  output[['est.pi0']] <- est.pi0
  output[['est.mu']] <- est.mu
  output[['est.disp']] <- est.disp
  output[['est.sf']] <- est.sf
  #output[['var.est']] <- EM.var
  output[['CE']] <- CE
  output[['loglik']] <- loglik.vec[1:iter]
  output[['logl']] <- loglik
  output[['GQ']] <- gq

  return(output)
}
