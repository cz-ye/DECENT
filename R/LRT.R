#' Likelihood-ratio test
#'
#' Likelihood-ratio test for DE analysis of scRNA-seq data
#'
#' @param data.obs Observed count matrix for endogeneous genes, rows represent genes, columns represent cells
#' @param out Output of fitNoDE, it contains EM algorithm estimates for models without DE between cell-types
#' @param tau cell-specific estimates (intcp,slope) that link Beta-Binomial dispersion parameter to the mean expression.
#' @param X Data matrix containing dummy variables for covariates of interest (cell-type) 
#' @param W Data matrix containing other covariates to adjust DE analysis. Default NULL
#' @param parallel If \code{TRUE}, run in parallel
#'
#' @return A list containing statistics, p-values and parameter estimates for models without DE.
#'
#' @export
#'
lrTest <- function(data.obs, out, X, W=NULL, tau, parallel) {

  message('Likelihood ratio test started at ', Sys.time())

  ngene <- nrow(out$est.mu)
  ncelltype <- ncol(X)
  ncell <- length(out$est.sf)
  nW    <- 0
  if(!is.null(W)) 
    nW <- ncol(W)
  DO.par <- matrix(0,ncell,2)
  DO.par[,1] <- log(out$CE/(1-out$CE))

  XW <- cbind(X,W)
  XW.H0 <- cbind(X[,1],W)
  
  tau0 <- tau[,1]; tau1 <- tau[,2]

  par1 <- matrix(0, ngene, 1+ncelltype+nW+1)
  par2 <- matrix(0, ngene, 1+1+nW+1)
  logl1<- rep(0, ngene)
  logl2<- rep(0, ngene)
  gq <- gauss.quad(16,kind='legendre')
  
  if (parallel) {
    temp <- foreach (i = 1:ngene, .combine = 'rbind', .packages = c('DECENT')) %dopar% {
      out$est.pi0[i,1] <- ifelse(is.na(out$est.pi0[i,1]) | is.infinite(out$est.pi0[i,1]),0.1,out$est.pi0[i,1])
      out$est.mu[i,1] <- ifelse(is.na(out$est.mu[i,1]) | is.infinite(out$est.mu[i,1]), mean(data.obs[i,]/out$CE,na.rm=T), out$est.mu[i,1])
      p.init <- c(log(out$est.pi0[i,1]/(1-out$est.pi0[i,1])), log(out$est.mu[i,1]), rep(0,nW+1))
      res2 <- tryCatch(optim(p = p.init, fn = loglI.GQ,
                             rho = (1+exp(-tau0-tau1*log(out$est.sf*out$est.mu[i,1]*(1-out$est.pi0[i,1]))))^-1,
                             sf = out$est.sf, XW = XW.H0, DO.par = DO.par, z = data.obs[i,], GQ.object=gq, control=list(maxit=300)),
                       error = function(e) {
                         warning("Numerical problem in noDE model for gene ", i);
                         NA
                       })
      if(!is.na(res2)) {
        res2$par[res2$par < -100] <- -100
        p2.init <- c(res2$p[1:2], rep(0,ncelltype+nW-1),res2$p[ncol(par2)])
      } else {
        p2.init <- c(p.init[1:2], rep(0,ncelltype+nW))
      }
      res1 <- tryCatch(optim(p = p2.init, fn = loglI.GQ, sf = out$est.sf, XW = XW,
                             rho = (1+exp(-tau0-tau1*log(out$est.sf*out$est.mu[i,1]*(1-out$est.pi0[i,1]))))^-1,
                             DO.par = DO.par, z = data.obs[i, ], GQ.object=gq, control=list(maxit=300)),
                       error = function(e) {
                         warning("Numerical problem in DE model for gene ", i);
                         NA
                       })
      message(i)
      if (is.na(res1) & is.na(res2)) {
        return(rep(0, ncol(par1)+ncol(par2)+2))
      } else {
        if (res1$conv > 0) {
          warning("DE model failed to converge for gene ", i)
        }
        if (res2$conv > 0) {
          warning("noDE model failed to converge for gene ", i)
        }
        return(c(res1$par, res2$par, -res1$value, -res2$value))
      }
    }
    par1 <- temp[, 1:ncol(par1)]
    par2 <- temp[, (ncol(par1)+1):(ncol(par1)+ncol(par2))]
    logl1 <- temp[, ncol(temp)-1]
    logl2 <- temp[, ncol(temp)]

  } else {
    for(i in 1:ngene) {
      out$est.pi0[i,1] <- ifelse(is.na(out$est.pi0[i,1]) | is.infinite(out$est.pi0[i,1]),0.1,out$est.pi0[i,1])
      out$est.mu[i,1] <- ifelse(is.na(out$est.mu[i,1]) | is.infinite(out$est.mu[i,1]), mean(data.obs[i,]/out$CE,na.rm=T), out$est.mu[i,1])
      p.init <- c(log(out$est.pi0[i,1]/(1-out$est.pi0[i,1])), log(out$est.mu[i,1]), rep(0,nW+1))
      res2 <- tryCatch(optim(p = p.init, fn = loglI.GQ,
                             rho = (1+exp(-tau0-tau1*log(out$est.sf*out$est.mu[i,1]*(1-out$est.pi0[i,1]))))^-1,
                             sf = out$est.sf, XW=XW.H0, DO.par = DO.par, z = data.obs[i,], GQ.object=gq, control=list(maxit=300)),
                       error = function(e) {
                         warning("numerical problem in noDE model for gene ", i);
                         NA
                       })
      if(!is.na(res2)) {
        res2$par[res2$par < -100] <- -100
        p2.init <- c(res2$p[1:2], rep(0,ncelltype+nW-1),res2$p[ncol(par2)])
      } else {
        p2.init <- c(p.init[1:2], rep(0,ncelltype+nW))
      }
      res1 <- tryCatch(optim(p = p2.init, fn = loglI.GQ, sf = out$est.sf, XW = XW,
                             rho = (1+exp(-tau0-tau1*log(out$est.sf*out$est.mu[i,1]*(1-out$est.pi0[i,1]))))^-1,
                             DO.par = DO.par, z = data.obs[i, ], GQ.object=gq, control=list(maxit=300)),
                       error = function(e) {
                         warning("numerical problem in DE model for gene ", i);
                         NA
                       })
      message(i)
      if (is.na(res1) | is.na(res2)) {
      } else {
        if (res1$conv > 0) {
          warning("DE model failed to converge for gene ", i)
        }
        if (res2$conv > 0) {
          warning("noDE model failed to converge for gene ", i)
        }
        par1[i, ] <- res1$par
        par2[i, ] <- res2$par
        logl1[i] <- -res1$value
        logl2[i] <- -res2$value
      }
    }
  }
  message('Likelihood ratio test finished at ', Sys.time())
  
  output <- list()
  lrt.stat <- 2*(logl1 - logl2)
  lrt.stat <- ifelse(lrt.stat<0, 0, lrt.stat)
  pval <- exp(pchisq(lrt.stat, df = ncol(XW)-ncol(XW.H0), lower.tail = FALSE, log = TRUE))
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


#' Calculate Log Likelihood of Incomplete Data
#'
loglI.GQ <- function(p, z, sf, XW, DO.par,rho, GQ.object) {
  pi0 <- exp(p[1])/(1 + exp(p[1]))
  mu  <- c(exp( XW %*% as.matrix(p[-c(1,length(p))]) ) * sf)
  size <- exp(-p[length(p)])
  # DO.probs is length(new.nodes) x ncell matrix
  f0  <- pi0 + (1-pi0)*dnbinom(0, mu = mu, size = size)
  # evaluate PZ (prob of observed data)
  PZ <- dBBNB(z,pi0=pi0,mu=mu,size=size,CE=1/(1+exp(-DO.par[,1])),rho=rho,GQ.object=GQ.object)
  PZ <- PZ + f0*(z == 0)
  PZ <- ifelse(PZ == 0, .Machine$double.xmin, PZ)
  return(-sum(log(PZ)))
}


#' Probablity mass function of the distribution for observed data z_{ij}
#'
dBBNB <- function(z,pi0,mu,size,CE,rho,GQ.object,EY=FALSE) {
  a <- pmax.int(z+0.5,1.5)
  b <- pmax.int(a+1,qzinb(1-0.0001,omega=pi0,lambda=mu,k=size) + 0.5)
  rho[which(rho>0.999)] <- 0.999
  pi0[which(pi0< 1e-04)]<- 1e-04
  pi0[which(pi0> 0.999)]<- 0.999
  new.nodes <- outer((b-a)/2,GQ.object$nodes,'*') + (a+b)/2
  new.weight<- outer((b-a)/2,GQ.object$weights,'*')
  out.mat <- dnbinom2(new.nodes,mu=mu,size=size)*dbetabinom2(z,prob=CE,size=new.nodes,rho=rho)*new.weight
  out <- (1-pi0)*(rowSums(out.mat) + dnbinom2(a-0.5,mu=mu,size=size)*dbetabinom2(z,prob=CE,size=a-0.5,rho=rho))
  if(EY) {
    out <- list()
    out[['PZ']] <- (1-pi0)*(rowSums(out.mat) + dnbinom2(a-0.5,mu=mu,size=size)*dbetabinom2(z,prob=CE,size=a-0.5,rho=rho))
    EY.wt <- out.mat * (1/(out$PZ))
    # impute for all obs
    EY  <- (1-pi0)*(rowSums(EY.wt*new.nodes) + (a-0.5)*dnbinom2(a-0.5,mu=mu,size=size)*dbetabinom2(z,prob=CE,size=a-0.5,rho=rho)/out$PZ)
    EY[which(out$PZ==0)] <- .Machine$double.xmin
    out[['EY']] <- EY
  }
 out
}

#' Perform single imputation for one gene using the fitted model
#'
SImputeByGene <- function(par,z, pi0,mu,sf,disp,k,b,M=1) {
  rho <- 1/(1+exp(-k*log((1-pi0)*mu/sf)-b))
  pi0[which(pi0>0.999)] <- 0.999
  rho[which(rho>0.999)] <- 0.999
  ncell <- length(z)
  z.imp <- z
  for(iter in 1:ncell) {
   y <- z[iter]:pmax(z[iter]+1, ZIM::qzinb(0.999, omega=pi0[iter],lambda = mu[iter], k = 1/disp))
   DO.prob <- calcDOProb(x=c(z[iter],par[iter,],rho[iter]), y = y)
   NB.prob <- calcZINBProb(y, pi0=pi0[iter], mu = mu[iter], size = 1/disp)
   lpost.prob <- c(scale(DO.prob+NB.prob,scale=FALSE))
   post.prob  <- exp(lpost.prob)
   post.prob  <- post.prob/sum(post.prob,na.rm=T)
   y  <- y[which(! (is.na(post.prob) | is.infinite(post.prob)))]
   post.prob <- post.prob[which(! (is.na(post.prob) | is.infinite(post.prob)))]
   if(any(post.prob>0) & length(post.prob)>0) z.imp[iter]  <- sample(y,prob=post.prob,size=M)
   
 }
 z.imp
}

#' Draw one sample from the given distribution
#'
draw.sample <- function(x) {
  x <- x/sum(x,na.rm=T)
  x <- ifelse(is.na(x) | x==0, 1e-12,x)
  sample(length(x),prob=x,size=1,replace=TRUE)-1
}
