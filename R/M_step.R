#' M-step
#'
#' Perform the last part of M-step gene by gene, update pi_0, mu and phi
#'
#'
MstepNB <- function(p, y, sf, status, ct) {
  #n.ct=max(ct)
  pi0 <- exp(p[1]) / (1 + exp(p[1]))
  mu  <- exp(p[2] + p[(ct + 1)] * (ct > 1)) * sf
  size <- exp(-p[length(p)])
  logl <- (1 - status) * dzinb(0, lambda = mu, k = size, omega = pi0, log = TRUE) +
    status * (log(1-pi0) + dnbinom(0, mu = mu, size = size, log = TRUE)) +
    status * (-lbeta(y, size) - log(y)) + y * status * log(mu/(mu+size))
  logl <- ifelse(is.na(logl) | !is.finite(logl), -1e+20, logl)
  return(-sum(logl))
}

#' Gradient of log-likelihood function for ZINB model
#'
#' Not used currently
#'
#' @param p A vector of parameter
#' p[1] = logistic regression intercept for pi0
#' p[2] = log mu1
#' p[3] = log DE par for mu
#' p[4] = log disp parameter
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
