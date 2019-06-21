#' E-step without Gauss–Legendre quadrature approximation
#'
#' Perform E-step for one gene, without GQ approximation
#'
#' @param par A matrix with two columns containing coefficient for the dropout model. The second column is currently unused
#' and set to zero but may be used in the future to model droput rates dependence on gene-specific factors.
#' @param z vector of observed count
#' @param sf Size factors for different cells
#' @param pi0 Gene-specific zero-inflated parameter
#' @param mu Gene-specific mean parameter
#' @param disp Reciprocal of gene-specific size parameter for
#'
EstepByGene <- function(par, z, sf, pi0, mu, disp, k, b) {

  rho <- 1/(1+exp(-k*log((1-pi0)*sf*mu)-b))
  y <- 1:max(2, qnbinom(0.999, mu = max(mu*sf), size = 1/disp)) #changed to 350 on 01/12/16
  DO.prob <- apply(as.matrix(cbind(z, par, rho)), 1, calcDOProb, y = y)

  NB.prob <- t(apply(as.matrix(y), 1, calcNBProb, mu = mu*sf, size = 1/disp))
  NB.prob <- NB.prob %*% diag(1-pi0)

  # calc output (PZ0S0)
  PZ0E1 <- colSums(DO.prob*NB.prob,na.rm=T)
  PE0   <- pi0 + (1-pi0)*dnbinom(0, mu = mu*sf, size = 1/disp)
  PE0Z0   <- PE0/(PE0 + PZ0E1)

  # if PE0=0 -> PE0Z0=0
  PE0Z0[is.na(PE0Z0)] <- 0.9999999

  # for genes with Z > 0, probability of 'expressed' is 1
  PE0Z0[z>0] <- 0
  # evaluate EY
  EY.wt <- (DO.prob*NB.prob) %*% diag(1/(PZ0E1))
  # impute for all obs
  EYZ0E1 <- colSums(EY.wt*y)

  out <- list()
  out[['PZ0E1']] <- PZ0E1
  out[['EYZ0E1']] <- EYZ0E1
  out[['PE0Z0']] <- PE0Z0
  return(out)
}

#' Calculate Log Likelihood
#'
loglI <- function(p, z, sf, ct, DO.par) {
  pi0 <- exp(p[1])/(1 + exp(p[1]))
  mu  <- exp(p[2] + p[(ct + 1)]*(ct > 1))*sf
  size <- exp(-p[length(p)])
  rho <- 1/(1+exp(-k*log((1-pi0)*mu)-b))

  f0  <- pi0 + (1-pi0)*dnbinom(0, mu = mu, size = size)

  y <- 1:qzinb(0.999, omega = mean(pi0), lambda = max(mu), k = size)
  DO.prob <- apply(as.matrix(cbind(z,DO.par,rho)), 1, calcDOProb, y = y)
  NB.prob <- t(apply(as.matrix(y), 1, calcNBProb, mu = mu, size = size))*(1-pi0)

  # evaluate PZ (prob of observed data)
  PZ <- colSums(DO.prob*NB.prob)
  PZ <- PZ + f0*(z == 0)
  PZ <- ifelse(PZ == 0, .Machine$double.xmin, PZ)
  return(-sum(log(PZ)))
}

#' Calculate probability for dropout model
#'
#' @import VGAM
#'
calcDOProb <- function(x, y) {
  return(dbetabinom(x[1], prob = exp(x[2] + x[3]*log(y+1))/(1 + exp(x[2] + x[3]*log(y+1))), size = y, log = TRUE, rho = x[4]))
}

#' Calculate probability for negative binomial model
#'
calcZINBProb <- function(y, pi0, mu, size) {
  return(ZIM::dzinb(y,omega=pi0,lambda = mu, k = size,log=TRUE))
}
