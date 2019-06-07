#' M-step
#'
#' Calculate negative complete data log-likelihood of DECENT model
#' @param p vector of the following parameters: log odds-ratio of pi0, log mean of the first cell type,
#' log fold-change for the mean parameter and log of dispersion (1/size) parameter.
#' @param y vector containing the expected number of molecules as output from the E-step.
#' @param sf vector of size factor estimates
#' @param status the estimated probability of non-zero pre-dropout count, output of the E-step.
#' @param ct A factor or a integer/numeric vector starting from 1 giving cell-type labels
#'
#' @return negative complete data log-likelihood evaluated at the parameter values.
#'
MstepNB <- function(p, y, vy, sf, status, ct) {
  #n.ct=max(ct)
  pi0 <- (1 + exp(-p[1]))^-1
  mu  <- exp(p[2] + p[(ct + 1)] * (ct > 1)) * sf
  size <- exp(-p[length(p)])
  logl <- (1 - status) * dzinb(0, lambda = mu, k = size, omega = pi0, log = TRUE) +
    status * (log(1-pi0) + dnbinom(0, mu = mu, size = size, log = TRUE)) +
    status * (-Elbeta(y, size,vy=vy)-log(y)) + y * status * log(mu/(mu+size))
  logl[which(is.na(logl) | !is.finite(logl))] <- -1e+20
  return(-sum(logl))
}


Elbeta <- function(y,size,vy) {
  lbeta(y,size) + (trigamma(y)+trigamma(size)-trigamma(y+size))*vy/2
}

#' Gradient of complete data log-likelihood function for DECENT model
#'
#' Not used currently
#'
#' @param p A vector of parameter
#' p[1] = logistic regression intercept for pi0
#' p[2] = log mu for the first cell type
#' p[3] = log fold-change (FC) for mean parameter
#' p[4] = log dispersion (1/size) parameter
#'
zinbGrad <- function(p,y,status,sf,ct) {

  nct <- max(ct)
  dmu <- exp(c(0, p[(nct + 1):(2*nct -1)]))
  mu <- exp(p[nct])*dmu
  pi0 <- rep(1/(1 + exp(-p[1])), nct)
  theta <- exp(-p[(2*nct)])

  grad.pi0 <- (1 - status)*(1 - (theta/(sf*mu[ct] + theta)) ^ theta)/
    (pi0[ct] + (1 - pi0[ct])*(theta/(sf*mu[ct] + theta)) ^ theta) +
    status*(-1/(1 - pi0[ct]))
  grad.a0  <- grad.pi0*pi0[ct]*(1 - pi0[ct])

  grad.lmu1 <- mu[1]*( (1 - status)*(-dmu[ct]*sf*(1 - pi0[ct])*(theta/(sf*mu[ct] + theta)) ^ (theta + 1)/
                                       (pi0[ct] + (1-pi0[ct])*(theta/(sf*mu[ct] + theta)) ^ theta)) +
                         status*(y/mu[1] - (y + theta)*dmu[ct]*sf/(theta + sf*mu[ct])) )

  grad.ldmu <- dmu[ct]*(ct > 1)*
               ( (1 - status)*sf*mu[1]*(-(1 - pi0[ct])*(theta/(sf*mu[ct] + theta)) ^ (theta + 1)/
                                        (pi0[ct] + (1 - pi0[ct])*(theta/(sf*mu[ct] + theta)) ^ theta)) +
                  status*(y/dmu[ct] - (y + theta)*sf*mu[1]/(theta + sf*mu[ct])) )

  grad.ltheta <- -theta*( (1 - status)*(log(theta) + 1 - log(sf*mu[ct]+theta) - theta/(sf*mu[ct]+theta))/
                           (1 + pi0[ct]/(1 - pi0[ct])*((sf*mu[ct] + theta)/theta) ^ theta) +
                           status*(digamma(y + theta) - digamma(theta) + log(theta) + 1 -
                                   log(sf*mu[ct] + theta) - (y + theta)/(sf*mu[ct] + theta)) )

  grad.vec <- -c(sum(grad.a0), sum(grad.lmu1), sum(grad.ldmu[ct == 2]), sum(grad.ltheta))

  return(grad.vec)
}
