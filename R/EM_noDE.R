#' Fitting the unrestricted model with EM algorithm
#'
#' Fit the DECENT model assuming no differentially-expressed (DE) genes
#'
#' @param data.obs Observed count matrix for endogeneous genes, rows represent genes, columns represent cells.
#' @param spike Observed count matrix for spike-ins, rows represent spike-ins, columns represent cells. Only needed if spikes = \code{TRUE}).
#' @param spike.conc A vector of theoretical count for each spike-in in one cell (Only needed if spikes = \code{TRUE}).
#' @param use.spikes If \code{TRUE}, use spike-ins to estimate capture efficiencies.
#' @param CE.range A two-element vector of the lower limit and upper limit for the estimated range of
#' capture efficiencies (ONLY needed if spikes = \code{FALSE}, default [0.02, 0.10]).
#' @param tau.init initial estimates (intcp,slope) that link Beta-Binomial dispersion parameter to the mean expression.
#' @param tau.global whether to use the same tau parameters across cell. Default TRUE
#' @param tau.est Methods to estimate tau parameters. The default 'endo' corresponds to using endogeneous genes. Other options
#' are 'none' which means tau.init is not further estimated and 'spikes' corresponds to using spike-ins. 
#' @param normalize Method for estimating size factors, either 'ML' (maximum likelihood, Ye et al., 2017) or 'TMM' (Robinson et al., 2010).
#' @param GQ.approx If \code{TRUE}, use Gaussian-Quadrature approximation to speed up E-step.
#' @param maxit maximum number of iterations for EM algorithm.
#' @param parallel If \code{TRUE}, run DECENT in parallel.
#'
#' @return A list containing estimates of DE model
#' @examples
#'
#' @import MASS
#' @import ZIM
#' @import statmod
#' @import edgeR
#'
#' @export

fitNoDE <- function(data.obs, spikes, spike.conc, use.spikes, CE.range, tau.init, tau.global, tau.est,
                    normalize, GQ.approx, maxit, parallel)
 {
  ncell <- ncol(data.obs)
  ngene <- nrow(data.obs)
  cell.type <- rep(1, ncell)
  cell.type.names <- NULL
  ncelltype <- length(unique(cell.type))
  XW <- as.matrix(model.matrix(~cell.type)[,1])

  # Get capture efficiency. Calculate with spike-ins, if available;
  # If not, randomize and sort by libsize.
  if (use.spikes) {
    capeff.spike <- apply(spikes, 2, sum)/sum(spike.conc)
    DO.coef <- matrix(0, ncell, 2)
    DO.coef[,1] <- log(capeff.spike/(1-capeff.spike))
    CE <- capeff.spike
  } else {
    obs.ls  <- log10(colSums(data.obs))
    max.ls  <- max(obs.ls)
    min.ls  <- min(obs.ls)
    # generate rand.CE within CE.range but following the dist of obs.ls closely
    ls.wt <- (obs.ls-min.ls)/(max.ls-min.ls)
    rand.CE <- (1-ls.wt)*CE.range[1] + ls.wt*CE.range[2]
    CE <- rand.CE
    DO.coef <- matrix(0, ncell, 2)
    DO.coef[, 1] <- log(CE/(1-CE))
  }

  # Initialize size factor
  data.obs.adj <- sweep(data.obs,2,CE,'/')
  est.sf <- apply(data.obs.adj, 2, mean, trim = 0.025)
  est.sf <- est.sf/mean(est.sf)

  # Initialize other ZINB parameters
  est.disp  <- rbeta(ngene, 0.1 ,0.6)
  est.pi0   <- matrix(0, ngene, ncelltype)
  # start with small pi0 (close to zero)
  est.pi0[, 1]   <-  rbeta(ngene, 1,9)
  if(ncelltype>1) {
    for (K in 2:ncelltype) {
       est.pi0[, K] <- est.pi0[, 1]
    }
  }

  est.mu <- matrix(0, ngene, ncelltype)
  # start with est.mu close to method of moments estimate
  est.mu[, 1] <- rowMeans( sweep(data.obs.adj,2, est.sf, '/' ))/(1-est.pi0[,1])

  #....new block of codes, june 7 2019 
  # use MoM to get a reasonable starting value
  dlogZINB <- function(p,y,sf) {
       pi0 <- 1/(1+exp(-p[1]))
       mu  <- exp(p[2])
       size<- exp(-p[3])
       -sum(ZIM::dzinb(y,omega=pi0,lambda=sf*mu,k=size,log=TRUE))
  }
  print('Finding starting values for EM algorithm...')  
  if (parallel) {
      temp <- foreach (i = 1:ngene, .combine = 'rbind', .packages = c('ZIM', 'DECENT')) %dopar% {
          out <- tryCatch(optim(par = c(0,log(mean(data.obs.adj[i,], na.rm=T)),0),
                                fn = dlogZINB, y = data.obs[i, ], sf = est.sf*CE,lower=-30),
                          error = function(e) {
                            list(p = c(0,log(mean(data.obs.adj[i,], na.rm=T)),0))
                          })
          new.pi0 <- rep(1/(1 + exp(-out$p[1])), ncelltype)
          new.mu <- exp(out$p[2])
          new.disp <- exp(out$p[length(out$p)])
          return(c(new.pi0, new.mu, new.disp))
      } 
      est.pi0 <- as.matrix(temp[, 1:ncelltype])
      est.mu <- as.matrix(temp[, (ncelltype+1):(2*ncelltype)])
      est.disp <- temp[, 2*ncelltype+1]

    } else {
      for (i in 1:ngene) {
          out <- tryCatch(optim(par = c(0,log(mean(data.obs.adj[i,], na.rm=T)),0),
                                fn = dlogZINB, y = data.obs[i, ], sf = est.sf*CE,lower=-30),
                          error = function(e) {
                            list(p = c(0,log(mean(data.obs.adj[i,], na.rm=T)),0))
			  })
	  est.pi0[i, ] <- rep(1/(1 + exp(-out$p[1])), ncelltype)
          est.mu[i, ]  <- exp(out$p[2])
          est.disp[i] <- exp(out$p[length(out$p)])
        }
    }
  #...end of new block

  # Initialize other variables
  loglik.vec <- rep(0, maxit)
  data.imp <- data.obs
  PE <- matrix(0, ngene, ncell)

  if (tau.global) {
    tau0 <- tau.init[1]; tau1 <- tau.init[2]
    tau.old <- c(tau0, tau1)
  } else{
    if (is.null(dim(tau.init))) {
      tau0 <- rep(tau.init[1],ncell)
      tau1 <- rep(tau.init[2],ncell)
    } else {
      if(dim(tau.init) != c(ncell, 2)) {
        stop('tau.init must be either a vector of 2 or a matrix with shape (#cells, 2)')
      }
      tau0 <- tau.init[, 1]; tau1 <- tau.init[, 2]
    }
    tau.old <- cbind(tau0,tau1)
  }
  tau.conv<- FALSE
  
  if (tau.est == 'spikes') {
    if (tau.global) {
      set.seed(1)
      y.sim <- t(sapply(1:nrow(spikes), function(i) { rpois(ncell, spike.conc[i])}))
      est.p <- optim(par = c(tau0, tau1), f = cbbinom.logl, z = spikes, y = y.sim, prob = CE, c = spike.conc)$p
      tau0 <- est.p[1]; tau1 <- est.p[2]
    } else {
      est.p <- matrix(0, ncell, 2)
      for (j in 1:ncell) {
        set.seed(1)
        y.sim <- t(sapply(1:nrow(spikes), function(i) rpois(10, spike.conc[i])))
        est.p[j, ] <- optim(par = c(tau0[j], tau1[j]), f = cbbinom.logl, z = spikes[, j], y = y.sim, 
                            prob = CE[j], c = spike.conc)$p
      }
      tau0 <- est.p[, 1]; tau1 <- est.p[, 2]
    }

  }

  iter <- 1
  converge <- FALSE
  if (GQ.approx) gq <- gauss.quad(16, kind = 'legendre') else gq <- NULL

  message('No-DE model fitting started at ', Sys.time())
  # Begin EM algorithm
  for (iter in 1:maxit) {

    # E-step gene by gene
    if (parallel) {
      if (!GQ.approx) {
        temp <- foreach (i = 1:ngene, .combine = 'rbind', .packages = c('DECENT')) %dopar% {
          out <- EstepByGene(par = DO.coef, z = data.obs[i, ], sf = est.sf,
                              pi0 = est.pi0[i, cell.type], mu = est.mu[i, cell.type], disp = est.disp[i])
          return(c(ifelse(is.na(out$EYZ0E1) | is.infinite(out$EYZ0E1),data.obs[i, ],out$EYZ0E1), 1 - out$PE0Z0))
        }
      } else {
        temp <- foreach (i = 1:ngene, .combine = 'rbind', .packages = c('MASS','ZIM', 'DECENT')) %dopar% {
          out <- Estep2ByGene(par = DO.coef,z = data.obs[i, ], sf = est.sf,
                              pi0 = est.pi0[i, cell.type], mu = est.mu[i, cell.type], disp = est.disp[i],
                              k = tau1, b = tau0, GQ.object = gq)
          return(c(ifelse(is.na(out$EYZ0E1) | is.infinite(out$EYZ0E1),data.obs[i, ],out$EYZ0E1), 1 - out$PE0Z0))
        }
      }
      data.imp <- temp[, 1:ncell]
      PE <- temp[, (ncell+1):(2*ncell)]

    } else {
      if (!GQ.approx) {
        for (i in 1:ngene) {
          # use E-step with expected value evaluated using GQ integral
          out <- EstepByGene(par = DO.coef, z = data.obs[i, ], sf = est.sf,
                             pi0 = est.pi0[i, cell.type], mu = est.mu[i, cell.type], disp = est.disp[i])
          data.imp[i, ] <- ifelse(is.na(out$EYZ0E1) | is.infinite(out$EYZ0E1),data.obs[i, ],out$EYZ0E1)
          PE[i, ]<- 1 - out$PE0Z0
        }
      } else {
        for (i in 1:ngene) {
          out <- Estep2ByGene(par = DO.coef,z = data.obs[i, ], sf = est.sf,
                              pi0 = est.pi0[i, cell.type], mu = est.mu[i, cell.type], disp = est.disp[i],
                              k = tau1, b = tau0, GQ.object = gq)
          data.imp[i, ] <- ifelse(is.na(out$EYZ0E1) | is.infinite(out$EYZ0E1), data.obs[i, ], out$EYZ0E1)
          PE[i, ] <- 1 - out$PE0Z0
        }
      }
    }

    # M-step 1: Update SF
    data.imp <- as.matrix(data.imp)
    data.imp2 <- data.imp*PE

    # M-step 2: Estimate SF by maximum-likelihood
    if (normalize == 'ML') {
      for (i in 1:ncell) {
        p0 <- est.pi0[, cell.type[i]] +
          (1 - est.pi0[, cell.type[i]])*dnbinom(0, mu = est.sf[i]*est.mu[, cell.type[i]], size = 1/est.disp)
        w  <- ((p0 - est.pi0[, cell.type[i]]*(1-PE[, i]))*(1 - est.pi0[, cell.type[i]]))/p0
        est.sf[i]  <- sum(data.imp2[, i], na.rm=T)/sum(w, na.rm=T)
      }
    } else if (normalize == 'TMM') {
      tmm <- calcNormFactors(ifelse(is.na(data.imp2) | is.infinite(data.imp2), data.obs.adj, data.imp2))
      est.sf <- colSums(ifelse(is.na(data.imp2) | is.infinite(data.imp2), data.obs.adj, data.imp2))*tmm
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
          out <- tryCatch(optim(par = c(log(prop0/(1-prop0)), log(mean(data.imp[i, ], na.rm=T)), rep(0, ncelltype-1), -2),
                                fn = MstepNB, y = data.imp[i, ], sf = est.sf, status = PE[i, ], ct = cell.type,lower=-30),
                                #gr = zinbGrad, method = 'L-BFGS-B'),
                          error = function(e) {
                            list(p = c(log(prop0/(1-prop0)), log(mean(data.imp[i, ], na.rm=T)), rep(0, ncelltype-1), -2))
                          })
          new.pi0 <- rep(1/(1 + exp(-out$p[1])), ncelltype)
          new.mu <- exp(out$p[2])
          new.disp <- exp(out$p[length(out$p)])
          
          if(!GQ.approx){
            new.loglik <- -loglI(p = out$p, sf = est.sf, ct = cell.type, DO.par = DO.coef, z = data.obs[i, ])
          } else {
             rho <- 1/(1+exp(-tau0-tau1*log(est.mu[i]*(1-est.pi0[i]))))
             rho <- ifelse(rho<1e-05,1e-05,rho)
             new.loglik <- -loglI.GQ(p=out$p, z=data.obs[i,], sf = est.sf, XW=XW, DO.par=DO.coef,rho=rho, GQ.object=gq)
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
          out <- tryCatch(optim(par = c(log(prop0/(1-prop0)), log(mean(data.imp[i, ], na.rm=T)), rep(0, ncelltype-1), -2),
                                fn = MstepNB, y = data.imp[i, ], sf = est.sf, status = PE[i, ], ct = cell.type,lower=-30),
                                #gr = zinbGrad, method = 'L-BFGS-B')
                          error = function(e) {
                            list(p = c(log(prop0/(1-prop0)), log(mean(data.imp[i, ], na.rm=T)), rep(0, ncelltype-1), -2))
                          })
          est.pi0[i, ] <- rep(1/(1 + exp(-out$p[1])), ncelltype)
          est.mu[i, ]  <- exp(out$p[2])
          est.disp[i] <- exp(out$p[length(out$p)])
          if (!GQ.approx) {
            loglik[i] <- -loglI(p = out$p, sf = est.sf, ct = cell.type, DO.par = DO.coef, z = data.obs[i, ])
          } else {
            loglik[i] <- -loglI2(p = out$p, sf = est.sf, ct = cell.type, DO.par = DO.coef, k = tau1, b = tau0,
                                 z = data.obs[i, ], GQ.object = gq)
          }
        }
      }
    }

    # update tau1 and tau0 when no-spikeins: NOTE this code is not parallelized yet and will only update (k,b) until (k,b) converges
    if(!tau.conv & tau.est=='endo') {
      if(tau.global) {
        logit.rho <- foreach (i = 1:ngene, .combine = 'rbind', .packages = c('DECENT')) %dopar% {
          out.rho <- optim(p=-2,fn=update.rho,x=data.obs[i,],sf=est.sf,size=mean(data.imp2[i,],na.rm=T),CE=CE,method='Brent',upper=10,lower=-10)
          out.rho$p
        }
        data.rho <- data.frame(y=logit.rho,x=log(rowMeans(data.imp2,na.rm=T)),w=rowMeans(data.imp2,na.rm=T))
        rho.model <- glm(y ~ x,data=data.rho)
        # weighted est is better when number of genes is small
        if(ngene<=5000) 
          rho.model <- glm(y ~ x,weights=w, data=data.rho)

        tau.old    <- c(tau0, tau1)
        tau.new    <- coef(rho.model)
        tau0 <- tau.new[1] ; tau1 <- tau.new[2]
        tau.reltol <- sum( abs(tau.old-tau.new)/abs(tau.old) )
        tau.conv   <- ifelse(tau.reltol< 1e-04, TRUE,FALSE)
      } else {
        tau.old <- cbind(tau0,tau1)
        tau.new <- foreach (i = 1:ncell, .combine = 'rbind', .packages = c('DECENT')) %dopar% {
          #print(i)
          #size.bb <- rowMeans(data.imp2,na.rm=TRUE)*est.sf[i]
          size.bb <- ifelse(!is.na(data.imp2[,i]) & is.finite(data.imp2[,i]),data.imp2[,i],data.obs.adj[,i])
          out.tau <- optim(p=tau.old[i,],fn=update.rho3,z=data.obs[,i], size=size.bb,CE=rep(CE[i],ngene))
          out.tau$p
        }
        tau0 <- tau.new[, 1]; tau1 <- tau.new[, 2]
        tau.reltol <- apply( abs(tau.new-tau.old)/abs(tau.old), 2, mean)
        tau.conv   <- ifelse(any(tau.reltol > 1e-04), FALSE,TRUE)
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
  rownames(est.mu) <- rownames(data.obs)
  rownames(est.pi0) <- rownames(data.obs)
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
  output[['CE']] <- CE
  output[['loglik']] <- loglik.vec[1:iter]
  output[['logl']] <- loglik
  output[['GQ']] <- gq
  if(tau.global) {
    output[['tau1']] <- rep(tau1, ncell)
    output[['tau0']] <- rep(tau0, ncell)
  } else {
    output[['tau1']] <- tau1
    output[['tau0']] <- tau0
  }
  return(output)
}
