#' E-step with Gauss–Legendre quadrature approximation
#'
#' Perform E-step for one gene, with GQ approximation
#'
#' @param par A matrix with two columns containing coefficient for the dropout model. The second column is currently unused
#' and set to zero but may be used in the future to model droput rates dependence on gene-specific factors.
#' @param z vector of observed count
#' @param z.ind A binary indicator of (z>0) in the original obs count (before imputation)
#' @param sf Size factors for different cells
#' @param pi0 gene-specific zero-inflated parameter
#' @param mu gene-specific mean parameter
#' @param disp Reciprocal of gene-specific size parameter for
#' @param GQ.object Gauss-Legendre quadrature points and weights
#'
Estep2ByGene <- function(par, z, z.ind, sf, pi0, mu,disp, GQ.object) {

  if (disp >= exp(-30)) {
    a <- 1
    b <- max(2, qzinb(0.999, omega = mean(pi0), lambda = max(mu*sf), k = 1/disp))
    new.nodes <- (b - a)/2*GQ.object$nodes + (a + b)/2

    # DO.probs is length(new.nodes) x ncell matrix
    DO.prob <- apply(as.matrix(cbind(z,par)), 1, calc2DOProb, y = new.nodes)

    NB.prob <- (b - a)/2*diag(GQ.object$weights) %*% t(apply(as.matrix(new.nodes), 1, calc2NBProb, mu = mu*sf, size = 1/disp))
    NB.prob <- NB.prob %*% diag(1-pi0)

    # calc output (PZ0S0)
    PZ0E1 <- colSums(DO.prob*NB.prob)
    PE0 <- pi0 + (1-pi0)*dnbinom(0, mu = mu*sf, size = 1/disp)
    PE0Z0 <- PE0/(PE0 + PZ0E1)

    # for genes with Z>0, probability of 'expressed' is 1
    PE0Z0[z.ind] <- 0
    # evaluate EY
    EY.wt <- (DO.prob*NB.prob) %*% diag(1/(PZ0E1))
    # impute for all obs
    EYZ0E1 <- colSums(EY.wt*new.nodes)
  } else {
    a <- 1
    b <- max(2, qzip(0.999, omega = mean(pi0), lambda = max(mu*sf)))
    new.nodes <- (b - a)/2*GQ.object$nodes + (a + b)/2

    # DO.probs is length(new.nodes) x ncell matrix
    DO.prob <- apply(as.matrix(cbind(z,par)), 1, calc2DOProb, y = new.nodes)

    NB.prob <- (b - a)/2*diag(GQ.object$weights) %*% t(apply(as.matrix(new.nodes), 1, calc2PoisProb, mu = mu*sf))
    NB.prob <- NB.prob %*% diag(1-pi0)

    # calc output (PZ0S0)
    PZ0E1 <- colSums(DO.prob*NB.prob)
    PE0   <- pi0 + (1-pi0)*dpois(0, lambda = mu*sf)
    PE0Z0   <- PE0/(PE0 + PZ0E1)

    # for genes with Z>0, probability of 'expressed' is 1
    PE0Z0[z.ind] <- 0
    # evaluate EY
    EY.wt <- (DO.prob*NB.prob) %*% diag(1/(PZ0E1))
    # impute for all obs
    EYZ0E1 <- colSums(EY.wt*new.nodes)
  }

  out <- list()
  out[['PZ0E1']] <- PZ0E1
  out[['EYZ0E1']] <- EYZ0E1
  out[['PE0Z0']] <- PE0Z0
  return(out)
}

#' Calculate Log Likelihood
#'
#'
loglI2 <- function(p, z, sf, ct, DO.par, GQ.object) {
  pi0 <- exp(p[1])/(1 + exp(p[1]))
  mu  <- exp(p[2]+p[(ct + 1)]*(ct > 1))*sf
  size <- exp(-p[length(p)])

  # DO.probs is length(new.nodes) x ncell matrix
  if (size <= exp(30)) {
    f0  <- pi0 + (1-pi0)*dnbinom(0, mu = mu, size = size)
    a <- 1
    b <- max(2, qzinb(0.999, lambda = max(mu), k = size, omega = pi0))
    new.nodes <- (b-a)/2*GQ.object$nodes + (a+b)/2

    DO.prob <- apply(as.matrix(cbind(z,DO.par)), 1, calc2DOProb, y = new.nodes)
    NB.prob <- (b - a)/2*diag(GQ.object$weights) %*% t(apply(as.matrix(new.nodes), 1, calc2NBProb, mu = mu, size = size)) * (1-pi0)
    # evaluate PZ (prob of observed data)
    PZ <- colSums(DO.prob*NB.prob)
    PZ <- PZ + f0*(z == 0)
    PZ <- ifelse(PZ == 0, .Machine$double.xmin, PZ)
  }
  if (size > exp(30)) {
    f0  <- pi0+ (1-pi0)*dpois(0, lambda = mu)
    a <- 1
    b <- max(2, qzip(0.999, lambda = max(mu), omega = pi0))
    new.nodes <- (b - a)/2*GQ.object$nodes + (a + b)/2

    DO.prob <- apply(as.matrix(cbind(z, DO.par)), 1, calc2DOProb, y = new.nodes)
    NB.prob <- (b - a)/2*diag(GQ.object$weights) %*% t(apply(as.matrix(new.nodes), 1, calc2PoisProb, mu = mu))*(1-pi0)
    # evaluate PZ (prob of observed data)
    PZ <- colSums(DO.prob*NB.prob)
    PZ <- PZ + f0*(z == 0)
    PZ <- ifelse(PZ == 0, 1e-08, PZ)
  }
  return(-sum(log(PZ)))
}

#' Calculate probability for dropout model
#'
calc2DOProb <- function(x, y, rho = 0.2) {
  return(dbetabinom2(x[1], prob = exp(x[2] + x[3]*log(y+1))/(1 + exp(x[2] + x[3]*log(y+1))), size = y, rho=rho))
#  return(dbinom2(x[1], prob = exp(x[2] + x[3]*log(y+1))/(1 + exp(x[2] + x[3]*log(y+1))), size = y))
}

#' Calculate probability for negative binomial model
#'
calc2NBProb <- function(y, mu, size) {
  return(dnbinom2(y, mu = mu, size = size))
}

#' dnbinom function for non-integer input
#'
dnbinom2 <- function(x, mu, size) {
  C <- ifelse(x == 0, 0, -lbeta(x, size) - log(x))
  logl <- C + x*log(mu/(mu+size)) + dnbinom(0, mu = mu, size = size, log = T)
  return(exp(logl))
}

#' dbinom function for non-integer input
#'
dbinom2 <- function(x, prob, size) {
  corr.size <- ifelse(size < x, x, size)
  logl <- lgamma(corr.size + 1) - lgamma(x + 1) - lgamma(corr.size - x + 1) + x*log(prob) + (corr.size - x)*log(1-prob)
  return(exp(logl)*(x <= size))
}

#' dbetabinom function for non-integer input
#'
dbetabinom2 <- function(x,prob,size,rho) {
  a <- prob*(1-rho)/rho ; b <- (1-prob)*(1-rho)/rho
  corr.size <- ifelse(size<x,x,size)
  logl <- lgamma(size+1) - lgamma(x+1) - lgamma(corr.size - x + 1) + lbeta(x+a,corr.size-x+b)-lbeta(a,b)
  return(exp(logl)*(x <= size))
}

#' Calculate probability for Poisson model
#'
calc2PoisProb <- function(y, mu) {
  return(dpois2(y, mu = mu))
}

#' dpois function for non-integer output
#'
dpois2 <- function(x, mu) {
  logl <- -mu + x*log(mu) - lgamma(x + 1)
  return(exp(logl))
}
