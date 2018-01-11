#' Get imputed data
#'
#' Get imputed data by running E-step once without GQ approximation using the LRT estimates.
#'
#' @param data.obs Observed count matrix for endogeneous genes, rows represent genes, columns represent cells
#' @param p LRT estimates, using par.DE by default
#' @param sf A vector containing size factors for all cells
#' @param CE A  vector containing capture efficiencies for all cells
#' @param cell.type A factor or a integer/numeric vector starting from 1 giving cell-type labels
#' @param parallel If \code{TRUE}, run in parallel
#'
#' @return A matrix containing the imputed data assuming DE.
#'
#' @export
#'
getImputed <- function(data.obs, p, sf, CE, cell.type, parallel = T) {
  cell.type <- as.numeric(cell.type)
  ngene <- nrow(data.obs)
  ncell <- ncol(data.obs)

  DO.coef <- matrix(0, ncell, 2)
  DO.coef[, 1] <- log(CE/(1-CE))

  est.pi0 <- exp(p[,1])/(1+exp(p[,1]))
  est.pi0 <- cbind(est.pi0, est.pi0)
  est.mu <- cbind(exp(p[,2]), exp(p[,2]+p[,3]))
  est.disp <- exp(p[,4])
  if (parallel) {
    temp <- foreach (i = 1:ngene, .combine = 'rbind', .packages = c('VGAM', 'DECENT')) %dopar% {
      out <- EstepByGene(par = DO.coef, z = data.obs[i, ], z.ind = data.obs[i, ] > 0, sf = sf,
                         pi0 = est.pi0[i, cell.type], mu = est.mu[i, cell.type], disp = est.disp[i])
      return(c(ifelse(is.na(out$EYZ0E1), data.obs[i, ], out$EYZ0E1), 1 - out$PE0Z0))
    }
    data.imp <- temp[, 1:ncell]
    PE <- temp[, (ncell+1):(2*ncell)]

  } else {
    data.imp <- data.obs
    PE <- matrix(0, ngene, ncell)
    for (i in 1:ngene) {
      out <- EstepByGene(par = DO.coef, z = data.obs[i, ], z.ind = data.obs[i, ] > 0, sf = sf,
                         pi0 = est.pi0[i, cell.type], mu = est.mu[i, cell.type], disp = est.disp[i])
      data.imp[i, ] <- ifelse(is.na(out$EYZ0E1),data.obs[i, ],out$EYZ0E1)
      PE[i, ]<- 1 - out$PE0Z0
    }
  }
  data.imp = as.matrix(data.imp)
  dimnames(data.imp) <- dimnames(data.obs)
  return(data.imp*PE)
}
