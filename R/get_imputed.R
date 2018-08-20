#' Get imputed data
#'
#' Get mean imputed data by running E-step once without GQ approximation using the LRT estimates.
#'
#' @param data.obs Observed count matrix for endogeneous genes, rows represent genes, columns represent cells
#' @param par.est LRT estimates, using par.DE by default
#' @param sf A vector containing size factors for all cells
#' @param CE A  vector containing capture efficiencies for all cells
#' @param XW A matrix containing information for X (cell-type) and W (other covariates)
#' @param tau cell-specific estimates (intcp,slope) that link Beta-Binomial dispersion parameter to the mean expression.
#' @param parallel If \code{TRUE}, run in parallel
#'
#' @return A matrix containing the mean imputed data.
#'
#' @export
#'
getImputed <- function(data.obs, par.est, sf, CE, XW, tau, parallel = T) {
  ngene <- nrow(data.obs)
  ncell <- ncol(data.obs)

  DO.coef <- matrix(0, ncell, 2)
  DO.coef[, 1] <- log(CE/(1-CE))
  gq <- gauss.quad(16)
  if (parallel) {
    temp <- foreach (i = 1:ngene, .combine = 'rbind', .packages = c('VGAM', 'DECENT2')) %dopar% {
      est.pi0 <- exp(par.est[i,1])/(1 + exp(par.est[i,1]))
      est.mu  <- c(exp( XW %*% as.matrix(par.est[i,-c(1,ncol(par.est))]) ) )
      est.disp<- exp(par.est[i,ncol(par.est)])
      out <- Estep2ByGene(par = DO.coef, z = data.obs[i, ], sf = sf,
                         pi0 = rep(est.pi0,ncell), mu = est.mu, disp = est.disp,
                         k = tau[,2], b = tau[,1], GQ.object=gq)
      return(c(ifelse(is.na(out$EYZ0E1), data.obs[i, ]/CE, out$EYZ0E1), 1 - out$PE0Z0))
    }
    data.imp <- temp[, 1:ncell]
    PE <- temp[, (ncell+1):(2*ncell)]

  } else {
    data.imp <- data.obs
    PE <- matrix(0, ngene, ncell)
    for (i in 1:ngene) {
      est.pi0 <- exp(par.est[i,1])/(1 + exp(par.est[i,1]))
      est.mu  <- c(exp( XW %*% as.matrix(par.est[i,-c(1,ncol(par.est))]) ) )
      est.disp<- exp(par.est[i,ncol(par.est)])
      out <- Estep2ByGene(par = DO.coef, z = data.obs[i, ], sf = sf,
                         pi0 = rep(est.pi0,ncell), mu = est.mu, disp = est.disp,
                         k = tau[,2], b = tau[,1], GQ.object=gq)
      data.imp[i, ] <- ifelse(is.na(out$EYZ0E1),data.obs[i, ]/CE,out$EYZ0E1)
      PE[i, ]<- 1 - out$PE0Z0
    }
  }
  data.imp = as.matrix(data.imp)
  dimnames(data.imp) <- dimnames(data.obs)
  return(data.imp*PE)
}

#' Get imputed data
#'
#' Get single imputed data by sampling from the posterior distribution of the complete data given the observed data and LRT estimates. 
#'
#' @param data.obs Observed count matrix for endogeneous genes, rows represent genes, columns represent cells
#' @param par.est LRT/EM estimates, using par.DE by default
#' @param sf A vector containing size factors for all cells
#' @param CE A  vector containing capture efficiencies for all cells
#' @param XW A matrix containing information for X (cell-type) and W (other covariates)
#' @param tau cell-specific estimates (intcp,slope) that link Beta-Binomial dispersion parameter to the mean expression.
#' @param M number of imputed values to be generated.
#' @param parallel If \code{TRUE}, run in parallel
#'
#' @return A matrix containing the single imputed data.
#'
#' @export
#'
getSImputed <- function(data.obs, par.est, sf, CE, XW, tau, parallel = T) {
  ngene <- nrow(data.obs)
  ncell <- ncol(data.obs)

  DO.coef <- matrix(0, ncell, 2)
  DO.coef[, 1] <- log(CE/(1-CE))

  if (parallel) {
    temp <- foreach (i = 1:ngene, .combine = 'rbind', .packages = c('VGAM', 'DECENT2')) %dopar% {
      est.pi0 <- exp(par.est[i,1])/(1 + exp(par.est[i,1]))
      est.mu  <- c(exp( XW %*% as.matrix(par.est[i,-c(1,ncol(par.est))]) ) * sf)
      est.disp<- exp(par.est[i,ncol(par.est)])
      out <- SImputeByGene(par = DO.coef, z = data.obs[i, ], sf=sf,
                         pi0 = rep(est.pi0,ncell), mu = est.mu, disp = est.disp, 
                         k = tau[,2], b = tau[,1],M=1)
      return(out)
    }
    data.imp <- temp[,1:ncell]

  } else {
    data.imp <- matrix(0,nrow(data.obs),ncol(data.obs))
    for (i in 1:ngene) {
      est.pi0 <- exp(par.est[i,1])/(1 + exp(par.est[i,1]))
      est.mu  <- c(exp( XW %*% as.matrix(par.est[i,-c(1,ncol(par.est))]) ) * sf)
      est.disp<- exp(par.est[i,ncol(par.est)])
 
      data.imp[i,] <- SImputeByGene(par = DO.coef, z = data.obs[i, ], sf=sf,
                         pi0 = rep(est.pi0,ncell), mu = est.mu, disp = est.disp, 
                         k = tau[,2], b = tau[,1],M=1)
    }
  }
  dimnames(data.imp) <- dimnames(data.obs)
  
  return(data.imp)
}

