#' E-step with Gauss-Legendre quadrature approximation
#'
#' Perform E-step for one gene, with GQ approximation
#'
#' @param par A matrix with two columns containing coefficient for the dropout model. The second column is currently unused
#' and set to zero but may be used in the future to model droput rates dependence on gene-specific factors.
#' @param z vector of observed count
#' @param sf Size factors for different cells
#' @param pi0 gene-specific zero-inflated parameter
#' @param mu gene-specific mean parameter
#' @param disp Reciprocal of gene-specific size parameter for
#' @param GQ.object Gauss-Legendre quadrature points and weights
#'
Estep2ByGene <- function(par, z, sf, pi0, mu, disp, k, b, GQ.object) {
    rho <- 1/(1+exp(-k*log((1-pi0)*mu)-b))
    # calc PZ0E1 and EY
    out <- dBBNB(z,pi0=pi0,mu=mu*sf,size=1/disp,CE=1/(1+exp(-par[,1])),rho=rho,GQ.object=GQ.object,EY=TRUE)
    PZ0E1 <- out$PZ
    PZ0E1 <- ifelse(PZ0E1 == 0 | is.na(PZ0E1) | is.infinite(PZ0E1), .Machine$double.xmin, PZ0E1)
    PE0   <- pi0 + (1-pi0)*dnbinom(0, mu = mu*sf, size = 1/disp)
    PE0Z0 <- PE0/(PE0 + PZ0E1)
    # for genes with Z>0, probability of 'expressed' is 1
    PE0Z0[z>0] <- 0
    # impute for all obs
    EYZ0E1 <- out$EY

  out <- list()
  out[['PZ0E1']] <- PZ0E1
  out[['EYZ0E1']] <- EYZ0E1
  out[['PE0Z0']] <- PE0Z0
  return(out)
}

#' Calculate Log Likelihood
#'
loglI2 <- function(p, z, sf, ct, DO.par, k, b, GQ.object) {
  pi0 <- exp(p[1])/(1 + exp(p[1]))
  mu  <- exp(p[2]+p[(ct + 1)]*(ct > 1))*sf
  size <- exp(-p[length(p)])
  rho <- 1/(1+exp(-k*log((1-pi0)*mu/sf)-b))

  # DO.probs is length(new.nodes) x ncell matrix
  if (size <= exp(30)) {
    f0  <- pi0 + (1-pi0)*dnbinom(0, mu = mu, size = size)
    a <- 0.5
    b <- 0.5 + max(2, qzinb(0.999, lambda = max(mu), k = size, omega = pi0))
    new.nodes <- (b-a)/2*GQ.object$nodes + (a+b)/2

    DO.prob <- apply(as.matrix(cbind(z,DO.par,rho)), 1, calc2DOProb, y = new.nodes)
    NB.prob <- (b - a)/2*diag(GQ.object$weights) %*% t(apply(as.matrix(new.nodes), 1, calc2NBProb, mu = mu, size = size)) * (1-pi0)
    # evaluate PZ (prob of observed data)
    PZ <- colSums(DO.prob*NB.prob)
    PZ <- PZ + f0*(z == 0)
    PZ <- ifelse(PZ == 0, .Machine$double.xmin, PZ)
  }
  if (size > exp(30)) {
    f0  <- pi0+ (1-pi0)*dpois(0, lambda = mu)
    a <- 0.5
    b <- 0.5 + max(2, qzip(0.999, lambda = max(mu), omega = pi0))
    new.nodes <- (b - a)/2*GQ.object$nodes + (a + b)/2

    DO.prob <- apply(as.matrix(cbind(z, DO.par, rho)), 1, calc2DOProb, y = new.nodes)
    NB.prob <- (b - a)/2*diag(GQ.object$weights) %*% t(apply(as.matrix(new.nodes), 1, calc2PoisProb, mu = mu))*(1-pi0)
    # evaluate PZ (prob of observed data)
    PZ <- colSums(DO.prob*NB.prob)
    PZ <- PZ + f0*(z == 0)
    PZ <- ifelse(PZ == 0, .Machine$double.xmin, PZ)
  }
  return(-sum(log(PZ)))
}

#' Calculate probability for dropout model
#'
calc2DOProb <- function(x, y) {
  return(dbetabinom2(x[1], prob = exp(x[2] + x[3]*log(y+1))/(1 + exp(x[2] + x[3]*log(y+1))), size = y, rho=x[4]))
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

#' Cost function for optimizing tau using endogenous genes
#'
update.rho <- function(p,x,size,CE,sf=1) {
  require(VGAM)
  rho  <- 1/(1+exp(-p))
  logl <- dbetabinom3(x,prob=CE,size=size*sf,rho=rho)
  logl[is.na(logl) | is.infinite(logl)] <- log(.Machine$double.xmin)
  -sum(logl)
}

#' Cost function for optimizing tau using endogenous genes
#'
update.rho3 <- function(p,z,size,CE) {
  require(VGAM)
  logit.rho <- p[1] + p[2]*log(size)
  rho  <- 1/(1+exp(-logit.rho))
  rho  <- ifelse(rho<1e-05 | is.na(rho),1e-05,rho)
  lik <- dbetabinom(z,prob=CE,size=size,rho=rho,log=FALSE)
  lik <- ifelse(lik==0 | is.na(lik),.Machine$double.xmin,lik)
 -sum(log(lik[lik>0]))
}

#' Cost function for optimizing tau for spike-in data
#'
cbbinom.logl <- function(p, z, y, prob, c) {
  z <- as.matrix(z)
  logl <- sapply(1:nrow(y), function(i) {
    if (1/(1+exp(-p[1]-p[2]*log(c[i]))) < 1e-4 ){
      return(dbinom(x = z[i, ], prob=prob, size=y[i, ], log = T))
    } else {
      return(VGAM::dbetabinom(x = z[i, ], prob=prob, size=y[i, ], rho=1/(1+exp(-p[1]-p[2]*log(c[i]))), log = T))
    }
  })
  logl[is.infinite(logl)] <- log(.Machine$double.xmin)
  # logl
  -sum(logl)
}


# dbetabinom function for non-integer input
dbetabinom3 <- function(x,prob,size,rho,log=TRUE) {
  a <- prob*(1-rho)/rho ; b <- (1-prob)*(1-rho)/rho
  corr.size <- ifelse(size<x,x,size)
  logl <- lgamma(size+1) - lgamma(x+1) - lgamma(corr.size - x + 1) + lbeta(x+a,corr.size-x+b)-lbeta(a,b)
  if(!log)
    logl <- exp(logl)
  return(logl)
}
