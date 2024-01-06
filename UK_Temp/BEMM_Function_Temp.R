######################################################################################
# Nimble code for running MCMC on our bivariate extreme mixture model using 
# Automated Factor Slice Sampling
######################################################################################


run_MCMC_parallel <- function(seed, dat, niter, nburnin, thin){
  dir.work <- '/home/pgrad2/2448355h/My_PhD_Project/Biv_Ext_Mix_Mod'
  library(nimble)
  library(mvtnorm)
  source(file.path(dir.work, "KRSW/RevExp_U_Functions.r"))
  source(file.path(dir.work, "KRSW/CommonFunctions.r"))
  
  R_pmnorm_chol <- function(lower, upper, mean, cholesky){
    sigma <- t(cholesky) %*% cholesky
    return(pmvnorm(lower=lower, upper=as.vector(upper), mean=as.vector(mean), sigma=sigma, keepAttr = F))
  }
  
  
  R_dmvnorm_chol <- function(x, mean, cholesky, log){
    sigma <- t(cholesky) %*% cholesky
    dvect <- dmvnorm(x, mean=mean, sigma=sigma, log=log)
    if (log) {
      return(sum(dvect))
    }else{
      return(prod(dvect))
    }
    
  }
  
  
  pmnorm_chol <- nimbleRcall(function(lower = double(1), upper = double(1), 
                                      mean=double(1),cholesky=double(2)){}, 
                             Rfun = 'R_pmnorm_chol',
                             returnType = double(0))
  
  
  
  dmvnorm_chol <- nimbleRcall(function(x = double(2), mean = double(1), 
                                       cholesky=double(2),log=logical(0, default = 0)){}, 
                              Rfun = 'R_dmvnorm_chol',
                              returnType = double(0))
  
  
  
  nll.powunif.GPD.1<-function(theta,x,u,a.ind,lam.ind,sig.ind,gamma.ind, lamfix=FALSE, balthresh=FALSE, marg.scale.ind,marg.shape.ind)
  { 
    x.mat.ind <- 1
    if (dim(x)[1]==1){
      d <- dim(x)[2]
      x.mat.ind <- 0
    }else{
      d<-dim(x)[2]
    }
    
    a<-theta[a.ind]
    if(length(a)==1)
    {
      a<-rep(a,d)
    }
    
    if(lamfix){
      lam<-rep(1,d)
    }else{
      lam<-c(theta[lam.ind],1)
    }
    
    if(balthresh){
      lam<-1/(1+a)
    }
    
    sig<-theta[sig.ind]
    gamma<-theta[gamma.ind]
    
    sig<-sig[marg.scale.ind]
    gamma<-gamma[marg.shape.ind]
    
    rej<-NULL
    # upper bound when xi is greater than 0
    if(x.mat.ind){
      for(j in 1:d)
      {
        rej[j]<- gamma[j]<0 && any(x[,j]>-sig[j]/gamma[j])
      }
    }else{
      for(j in 1:d)
      {
        rej[j]<- gamma[j]<0 && any(x[j]>-sig[j]/gamma[j])
      }
    }
    
    if(any(lam<0.01)||any(a<=0.01)||any(sig<=0.001)||any(rej)){return(10e7)}
    
    nll.uc <- 0
    nll.pc <- 0
    if (!x.mat.ind){
      uc <- comp.gt(x, u)
      if (uc){
        L <- fX.powunif(x=x, a=a, lam=lam, sig=sig, gamma=gamma)
        nll.uc <- -log(L)
      }else{
        L2 <- fX.powunif.cens(x=x, u=u, lam=lam, a=a, sig=sig, gamma=gamma)
        nll.pc <- -log(L2)
      }
    }else{
      uc<-apply(x,1,comp.gt,u=u)
      
      x.uc<-x[uc,]
      x.pc<-x[!uc,]
      
      L<-apply(x.uc,1,fX.powunif,a=a,lam=lam,sig=sig,gamma=gamma)
      nll.uc<--sum(log(L))
      
      if(sum(!uc)>0)
      {
        L2<-apply(x.pc,1,fX.powunif.cens,a=a,lam=lam,u=u,sig=sig,gamma=gamma)
        nll.pc<--sum(log(L2))
      }
    }
    if (is.nan(nll.uc)|is.nan(nll.pc)|(nll.uc==-Inf)|(nll.uc==Inf)){
      return(10e7)
    }
    nll<-nll.uc+nll.pc
    
    return(nll)
  }
  
  
  
  
  
  nim_nll_powunif_GPD_SG <- nimbleRcall(function(theta=double(1), x=double(1), u=double(0), a.ind=double(1),
                                                 lam.ind=double(0), sig.ind=double(1), gamma.ind=double(1), 
                                                 lamfix=logical(0, default = 0), balthresh=logical(0, default = 0), 
                                                 marg.scale.ind=double(1), marg.shape.ind=double(1)){}, 
                                        Rfun = 'nll.powunif.GPD.1',
                                        returnType = double(0))
  
  nim_nll_powunif_GPD_MAT <- nimbleRcall(function(theta=double(1), x=double(2), u=double(0), a.ind=double(1),
                                                  lam.ind=double(0), sig.ind=double(1), gamma.ind=double(1), 
                                                  lamfix=logical(0, default = 0), balthresh=logical(0, default = 0), 
                                                  marg.scale.ind=double(1), marg.shape.ind=double(1)){}, 
                                         Rfun = 'nll.powunif.GPD.1',
                                         returnType = double(0))
  
  
  
  dbiextmix_ns <- nimbleFunction(
    run = function(x=double(2), para.mg=double(1), beta.a=double(2), beta.b=double(2),
                   X=double(2), thres=double(1), mu=double(1), 
                   cholesky=double(2), D=integer(0, default=2),
                   a.ind=double(1), lam.ind=double(0), lamfix=logical(0, default = 0), 
                   sig.ind=double(1), gamma.ind=double(1),
                   log = logical(0, default = 0)) {
      returnType(double(0))
      
      K <- length(beta.a[,1])
      cond <- (x[,1]>thres[1]) | (x[,2]>thres[2])
      n.tail <- sum(cond)
      n.bulk <- sum(!cond)
      y.tail <- matrix(c(x[cond,1] - thres[1], x[cond,2] - thres[2]), ncol=D)
      y.bulk <- x[!cond,]
      
      pi <- pmnorm_chol(lower=rep(-Inf,D), upper=thres, mean=mu, cholesky = cholesky)
      
      sig <- para.mg[1:D]
      gamma <- para.mg[(D+1):(2*D)]
      eta <- -sig/gamma
      eta[which(gamma<=0)] <- -Inf
      
      dtail <- 0
      dbulk <- 0
      
      if (n.tail>0){
        y.min <- eta
        for (i in 1:D){
          y.min[i] <- min(y.tail[,i])
        }
        if (all(y.min>eta)){
          X.tail <- X[cond,]
          llt <- 0
          alpha.inv <- rep(NA, D)
          lam <- rep(NA, D)
          for (i in 1:n.tail){
            for (j in 1:D){
              alpha.inv[j] <- exp(inprod(X.tail[i,(1+(j-1)*K):(j*K)], beta.a[,j]))
              if (j<D){
                lam[j] <- exp(inprod(X.tail[i,(1+(j-1)*K):(j*K)] ,beta.b[,j]))
              }else{
                lam[j] <- 1
              }
            }
            theta <- c(alpha.inv,lam[1:(D-1)],para.mg)
            llt.sg <- -nim_nll_powunif_GPD_SG(x=y.tail[i,], theta=theta, u=min(y.tail)-0.01, a.ind=a.ind,
                                              lam.ind=lam.ind, sig.ind=sig.ind, gamma.ind=gamma.ind, 
                                              lamfix=lamfix, balthresh=FALSE, 
                                              marg.scale.ind=1:2, marg.shape.ind=1:2)
            llt <- llt + llt.sg
          }
          if (log){
            dtail <- llt
          }else{
            dtail <- exp(llt)
          }
        }else{
          if (log) dtail <- -10e10
        }
      }
      
      if (n.bulk>0){
        dbulk <- dmvnorm_chol(y.bulk, mean=mu, cholesky = cholesky, log = log )
      }
      
      if (log) {
        totalProb <- n.tail*log(1-pi) + dtail + dbulk
      }else{
        totalProb <- (1-pi)^n.tail *dtail*dbulk
      }
      return(totalProb)
    })
  
  dbiextmix <- nimbleFunction(
    run = function(x=double(2), theta=double(1), thres=double(1), mu=double(1), 
                   cholesky=double(2), D=integer(0, default=2),
                   a.ind=double(1), lam.ind=double(0), lamfix=logical(0, default = 0), 
                   sig.ind=double(1), gamma.ind=double(1),
                   log = logical(0, default = 0)) {
      returnType(double(0))

      cond <- (x[,1]>thres[1]) | (x[,2]>thres[2])
      n.tail <- sum(cond)
      n.bulk <- sum(!cond)
      y.tail <- matrix(c(x[cond,1] - thres[1], x[cond,2] - thres[2]), ncol=D)
      y.bulk <- x[!cond,]
      
      pi <- pmnorm_chol(lower=rep(-Inf,D), upper=thres, mean=mu, cholesky = cholesky)
      
      sig <- theta[sig.ind]
      gamma <- theta[gamma.ind]
      eta <- -sig/gamma
      eta[which(gamma<=0)] <- -Inf
      
      dtail <- 0
      dbulk <- 0
      
      if (n.tail>0){
        y.min <- eta
        for (i in 1:D){
          y.min[i] <- min(y.tail[,i])
        }
        if (all(y.min>eta)){
          llt <- -nim_nll_powunif_GPD_MAT(x=y.tail, theta=theta, u=min(y.tail)-0.01, a.ind=a.ind,
                                          lam.ind=lam.ind, sig.ind=sig.ind, gamma.ind=gamma.ind, 
                                          lamfix=lamfix, balthresh=FALSE, 
                                          marg.scale.ind=1:2, marg.shape.ind=1:2)
          if (log){
            dtail <- llt
          }else{
            dtail <- exp(llt)
          }
        }else{
          if (log) dtail <- -10^10
        }
      }
      
      if (n.bulk>0){
        dbulk <- dmvnorm_chol(y.bulk, mean=mu, cholesky = cholesky, log = log )
      }
      
      if (log) {
        totalProb <- n.tail*log(1-pi) + dtail + dbulk
      }else{
        totalProb <- (1-pi)^n.tail *dtail*dbulk
      }
      
      return(totalProb)
    })
  
  
  
  
  rbiextmix_ns <- nimbleFunction(
    run = function(n=integer(0), para.mg=double(1), beta.a=double(2), beta.b=double(2),
                   X=double(2), thres=double(1), mu=double(1), 
                   cholesky=double(2), D=integer(0, default=2),
                   a.ind=double(1), lam.ind=double(0), lamfix=logical(0, default = 0), 
                   sig.ind=double(1), gamma.ind=double(1)) {
      returnType(double(2))
      
      totalProb <- matrix(1,nrow=n, ncol=D)
      return(totalProb)
    })
  
  registerDistributions(list(
    dbiextmix_ns = list(
      BUGSdist = "dbiextmix_ns(para.mg, beta.a, beta.b, thres, X, mu, cholesky, D, 
                  a.ind, lam.ind, lamfix, sig.ind, gamma.ind)",
      types = c('value = double(2)', 'para.mg = double(1)', 'beta.a = double(2)', 
                'beta.b = double(2)', 'thres = double(1)', 'X = double(2)',
                'mu = double(1)', 'cholesky = double(2)', 'D = integer(0)', 'a.ind = double(1)', 
                'lam.ind = double(0)', 'lamfix = logical(0)', 'sig.ind = double(1)',
                'gamma.ind = double(1)')
    )))
  
  rbiextmix <- nimbleFunction(
    run = function(n=integer(0), theta=double(1), thres=double(1), mu=double(1), 
                   cholesky=double(2), D=integer(0, default=2),
                   a.ind=double(1), lam.ind=double(0), lamfix=logical(0, default = 0), 
                   sig.ind=double(1), gamma.ind=double(1)) {
      returnType(double(2))
      
      totalProb <- matrix(1,nrow=n, ncol=D)
      return(totalProb)
    })
  
  registerDistributions(list(
    dbiextmix = list(
      BUGSdist = "dbiextmix(theta, thres, mu, cholesky, D, a.ind, lam.ind, lamfix, sig.ind, gamma.ind)",
      types = c('value = double(2)', 'theta = double(1)', 'thres = double(1)', 
                'mu = double(1)', 'cholesky = double(2)', 'D = integer(0)', 'a.ind = double(1)', 
                'lam.ind = double(0)', 'lamfix = logical(0)', 'sig.ind = double(1)',
                'gamma.ind = double(1)')
    )))
  
  uppertri_mult_diag <- nimbleFunction(
    run = function(mat = double(2), vec = double(1)) {
      returnType(double(2))
      p <- length(vec)
      out <- matrix(nrow = p, ncol = p, init = FALSE)
      for(i in 1:p)
        out[ , i] <- mat[ , i] * vec[i]
      return(out)
    })
  
  
  
  assign('R_pmnorm_chol', R_pmnorm_chol, envir = .GlobalEnv)
  assign('R_dmvnorm_chol', R_dmvnorm_chol, envir = .GlobalEnv)
  assign('pmnorm_chol', pmnorm_chol, envir = .GlobalEnv)
  assign('dmvnorm_chol', dmvnorm_chol, envir = .GlobalEnv)
  assign('nll.powunif.GPD.1', nll.powunif.GPD.1, envir = .GlobalEnv)
  assign('nim_nll_powunif_GPD_SG', nim_nll_powunif_GPD_SG, envir = .GlobalEnv)
  assign('nim_nll_powunif_GPD_MAT', nim_nll_powunif_GPD_MAT, envir = .GlobalEnv)
  assign('dbiextmix', dbiextmix, envir = .GlobalEnv)
  assign('dbiextmix_ns', dbiextmix_ns, envir = .GlobalEnv)
  assign('rbiextmix', rbiextmix, envir = .GlobalEnv)
  assign('rbiextmix_ns', rbiextmix_ns, envir = .GlobalEnv)
  assign('uppertri_mult_diag', uppertri_mult_diag, envir = .GlobalEnv)
  
  BivExtMixcode <- nimbleCode({
    for (i in 1:D)
      sds[i] ~ dunif(0, 50)
    Ustar[1:D,1:D] ~ dlkj_corr_cholesky(1.3, D)
    U[1:D,1:D] <- uppertri_mult_diag(Ustar[1:D, 1:D], sds[1:D])
    
    for (i in 1:D)
      mu[i] ~ dnorm(0, sd=50)
    
    for (i in 1:D)
      thres[i] ~ T(dnorm(mu.thres[i], sd=sd.thres[i]),min.thres[i],max.thres[i])
    
    if (station.ind){
      # priors for a
      for (i in 1:2)
        theta[i] ~ dunif(0,50)
      # priors for sig
      for (i in 3:5)
        theta[i] ~ dunif(0,50)
      # priors for gamma 
      for (i in 6:7)
        theta[i] ~ dunif(-1,1)
      
      # constraints for a finite marginal expectation gamma > -1/a
      constraint_data1 ~ dconstraint( theta[6] + 1/theta[1] > 0 )
      constraint_data2 ~ dconstraint( theta[7] + 1/theta[2] > 0 )
        
      y[1:N,1:D] ~ dbiextmix(theta=theta[1:7], thres=thres[1:D], mu=mu[1:D], 
                             cholesky=U[1:D,1:D], D=D,
                             a.ind=a.ind[1:D], lam.ind=lam.ind, lamfix=lamfix, 
                             sig.ind=sig.ind[1:D], gamma.ind=gamma.ind[1:D])
    }else{
      for (i in 1:2)
        para.mg[i] ~ dunif(0,50)
      for (i in 3:4)
        para.mg[i] ~ dunif(-1,1)
      for (i in 1:D){
        beta.a[1:D.pred,i] ~ dmnorm(mu_beta[1:D.pred], cov=cov_beta[1:D.pred,1:D.pred])
        beta.b[1:D.pred,i] ~ dmnorm(mu_beta[1:D.pred], cov=cov_beta[1:D.pred,1:D.pred])
      }
      y[1:N,1:D] ~ dbiextmix_ns(para.mg=para.mg[1:4], beta.a=beta.a[1:D.pred,1:D],
                                beta.b=beta.b[1:D.pred,1:D], thres=thres[1:2],
                                X=X[1:N,1:(D*D.pred)], mu=mu[1:D], 
                                cholesky=U[1:D,1:D], D=D,
                                a.ind=a.ind[1:D], lam.ind=lam.ind, lamfix=lamfix, 
                                sig.ind=sig.ind[1:D], gamma.ind=gamma.ind[1:D])
    }
    
    
  })
  set.seed(seed)
  # [1] "Ustar[1, 1]" "Ustar[2, 1]" "Ustar[1, 2]" "Ustar[2, 2]"
  # [5] "mu[1]"       "mu[2]"       "sds[1]"      "sds[2]"     
  # [9] "theta[1]"    "theta[2]"    "theta[3]"    "theta[4]"   
  # [13] "theta[5]"    "theta[6]"    "theta[7]"    "thres[1]"   
  # [17] "thres[2]" 
  Inits <- list(U=matrix( c(1,0,runif(2,0,1)), nrow=2), mu= runif(2,-5,5), sds=runif(2,0.5,5),
                theta=c(runif(5,0,5), runif(2,0,1)), thres=(runif(2,2,5)))
  BivExtMixmodel <- nimbleModel(BivExtMixcode, constants = list(N = nrow(dat), 
                                                                D = 2,
                                                                station.ind = TRUE,
                                                                D.pred = 2,
                                                                X = matrix(1,nrow=nrow(dat),ncol=2),
                                                                mu.thres = apply(dat,2,quantile,0.9),#90 quantile
                                                                sd.thres = c(50,50), 
                                                                min.thres = apply(dat,2,quantile,0.8),
                                                                max.thres = apply(dat,2,quantile,0.99),
                                                                mu_beta = rep(0,2),
                                                                cov_beta = 100*diag(2),
                                                                a.ind = 1:2,
                                                                lam.ind = 3,
                                                                sig.ind = c(4,5),
                                                                gamma.ind = c(6,7),
                                                                lamfix=TRUE), check = FALSE)
  
  
  BivExtMixmodel$setData(list(y = dat, constraint_data1=1, constraint_data2=1 ))  ## Set those values as data in the model
  BivExtMixmodel$setInits(Inits)  
  cBivExtMixmodel <- compileNimble(BivExtMixmodel, showCompilerOutput = TRUE)
  
  BivExtMixconf <- configureMCMC(BivExtMixmodel,
                                 enableWAIC = TRUE, time=TRUE)
  BivExtMixconf$removeSamplers(c('theta[1:7]', 'thres[1:2]'))
  BivExtMixconf$addSampler(target = c('theta[1:7]', 'thres[1:2]'), type = 'AF_slice')
  
  BivExtMixMCMC <- buildMCMC(BivExtMixconf)
  # BivExtMixMCMC$run(1)
  cBivExtMixMCMC <- compileNimble(BivExtMixMCMC, project = BivExtMixmodel, showCompilerOutput = TRUE)
  
  results <- runMCMC(cBivExtMixMCMC, niter = niter, nburnin=nburnin,thin=thin,
                     summary = TRUE, WAIC = TRUE,setSeed = seed)
  return(results)
}

#--------------------Define the functions globally to calculate WAIC-------------------------

# 
# R_pmnorm_chol <- function(lower, upper, mean, cholesky){
#   sigma <- t(cholesky) %*% cholesky
#   return(pmvnorm(lower=lower, upper=as.vector(upper), mean=as.vector(mean), sigma=sigma, keepAttr = F))
# }
# 
# 
# R_dmvnorm_chol <- function(x, mean, cholesky, log){
#   sigma <- t(cholesky) %*% cholesky
#   dvect <- dmvnorm(x, mean=mean, sigma=sigma, log=log)
#   if (log) {
#     return(sum(dvect))
#   }else{
#     return(prod(dvect))
#   }
#   
# }
# 
# 
# pmnorm_chol <- nimbleRcall(function(lower = double(1), upper = double(1), 
#                                     mean=double(1),cholesky=double(2)){}, 
#                            Rfun = 'R_pmnorm_chol',
#                            returnType = double(0))
# 
# 
# 
# dmvnorm_chol <- nimbleRcall(function(x = double(2), mean = double(1), 
#                                      cholesky=double(2),log=logical(0, default = 0)){}, 
#                             Rfun = 'R_dmvnorm_chol',
#                             returnType = double(0))
# 
# 
# 
# nll.powunif.GPD.1<-function(theta,x,u,a.ind,lam.ind,sig.ind,gamma.ind, lamfix=FALSE, balthresh=FALSE, marg.scale.ind,marg.shape.ind)
# { 
#   x.mat.ind <- 1
#   if (is.null(dim(x))){
#     d <- length(x)
#     x.mat.ind <- 0
#   }else{
#     d<-dim(x)[2]
#   }
#   
#   a<-theta[a.ind]
#   if(length(a)==1)
#   {
#     a<-rep(a,d)
#   }
#   
#   if(lamfix){
#     lam<-rep(1,d)
#   }else{
#     lam<-c(theta[lam.ind],1)
#   }
#   
#   if(balthresh){
#     lam<-1/(1+a)
#   }
#   
#   sig<-theta[sig.ind]
#   gamma<-theta[gamma.ind]
#   
#   sig<-sig[marg.scale.ind]
#   gamma<-gamma[marg.shape.ind]
#   
#   rej<-NULL
#   # upper bound when xi is greater than 0
#   if(x.mat.ind){
#     for(j in 1:d)
#     {
#       rej[j]<- gamma[j]<0 && any(x[,j]>-sig[j]/gamma[j])
#     }
#   }else{
#     for(j in 1:d)
#     {
#       rej[j]<- gamma[j]<0 && any(x[j]>-sig[j]/gamma[j])
#     }
#   }
#   
#   if(any(lam<0.01)||any(a<=0.01)||any(sig<=0.001)||any(rej)){return(10e7)}
#   
#   nll.uc <- 0
#   nll.pc <- 0
#   if (!x.mat.ind){
#     uc <- comp.gt(x, u)
#     if (uc){
#       L <- fX.powunif(x=x, a=a, lam=lam, sig=sig, gamma=gamma)
#       nll.uc <- -log(L)
#     }else{
#       L2 <- fX.powunif.cens(x=x, u=u, lam=lam, a=a, sig=sig, gamma=gamma)
#       nll.pc <- -log(L2)
#     }
#   }else{
#     uc<-apply(x,1,comp.gt,u=u)
#     
#     x.uc<-x[uc,]
#     x.pc<-x[!uc,]
#     
#     L<-apply(x.uc,1,fX.powunif,a=a,lam=lam,sig=sig,gamma=gamma)
#     nll.uc<--sum(log(L))
#     
#     if(sum(!uc)>0)
#     {
#       L2<-apply(x.pc,1,fX.powunif.cens,a=a,lam=lam,u=u,sig=sig,gamma=gamma)
#       nll.pc<--sum(log(L2))
#     }
#   }
#   if (is.nan(nll.uc)|is.nan(nll.pc)|(nll.uc==-Inf)|(nll.uc==Inf)){
#     return(10e7)
#   }
#   nll<-nll.uc+nll.pc
#   
#   return(nll)
# }
# 
# 
# 
# 
# 
# nim_nll_powunif_GPD_SG <- nimbleRcall(function(theta=double(1), x=double(1), u=double(0), a.ind=double(1),
#                                                lam.ind=double(0), sig.ind=double(1), gamma.ind=double(1), 
#                                                lamfix=logical(0, default = 0), balthresh=logical(0, default = 0), 
#                                                marg.scale.ind=double(1), marg.shape.ind=double(1)){}, 
#                                       Rfun = 'nll.powunif.GPD.1',
#                                       returnType = double(0))
# 
# nim_nll_powunif_GPD_MAT <- nimbleRcall(function(theta=double(1), x=double(2), u=double(0), a.ind=double(1),
#                                                 lam.ind=double(0), sig.ind=double(1), gamma.ind=double(1), 
#                                                 lamfix=logical(0, default = 0), balthresh=logical(0, default = 0), 
#                                                 marg.scale.ind=double(1), marg.shape.ind=double(1)){}, 
#                                        Rfun = 'nll.powunif.GPD.1',
#                                        returnType = double(0))
# 
# 
# 
# dbiextmix_ns <- nimbleFunction(
#   run = function(x=double(2), para.mg=double(1), beta.a=double(2), beta.b=double(2),
#                  X=double(2), thres=double(1), mu=double(1), 
#                  cholesky=double(2), D=integer(0, default=2),
#                  a.ind=double(1), lam.ind=double(0), lamfix=logical(0, default = 0), 
#                  sig.ind=double(1), gamma.ind=double(1),
#                  log = logical(0, default = 0)) {
#     returnType(double(0))
#     
#     K <- length(beta.a[,1])
#     cond <- (x[,1]>thres[1]) | (x[,2]>thres[2])
#     n.tail <- sum(cond)
#     n.bulk <- sum(!cond)
#     y.tail <- matrix(c(x[cond,1] - thres[1], x[cond,2] - thres[2]), ncol=D)
#     y.bulk <- x[!cond,]
#     
#     pi <- pmnorm_chol(lower=rep(-Inf,D), upper=thres, mean=mu, cholesky = cholesky)
#     
#     sig <- para.mg[1:D]
#     gamma <- para.mg[(D+1):(2*D)]
#     eta <- -sig/gamma
#     eta[which(gamma<=0)] <- -Inf
#     
#     dtail <- 0
#     dbulk <- 0
#     
#     if (n.tail>0){
#       y.min <- eta
#       for (i in 1:D){
#         y.min[i] <- min(y.tail[,i])
#       }
#       if (all(y.min>eta)){
#         X.tail <- X[cond,]
#         llt <- 0
#         alpha.inv <- rep(NA, D)
#         lam <- rep(NA, D)
#         for (i in 1:n.tail){
#           for (j in 1:D){
#             alpha.inv[j] <- exp(inprod(X.tail[i,(1+(j-1)*K):(j*K)], beta.a[,j]))
#             if (j<D){
#               lam[j] <- exp(inprod(X.tail[i,(1+(j-1)*K):(j*K)] ,beta.b[,j]))
#             }else{
#               lam[j] <- 1
#             }
#           }
#           theta <- c(alpha.inv,lam[1:(D-1)],para.mg)
#           llt.sg <- -nim_nll_powunif_GPD_SG(x=y.tail[i,], theta=theta, u=min(y.tail)-0.01, a.ind=a.ind,
#                                             lam.ind=lam.ind, sig.ind=sig.ind, gamma.ind=gamma.ind, 
#                                             lamfix=lamfix, balthresh=FALSE, 
#                                             marg.scale.ind=1:2, marg.shape.ind=1:2)
#           llt <- llt + llt.sg
#         }
#         if (log){
#           dtail <- llt
#         }else{
#           dtail <- exp(llt)
#         }
#       }else{
#         if (log) dtail <- -10e10
#       }
#     }
#     
#     if (n.bulk>0){
#       dbulk <- dmvnorm_chol(y.bulk, mean=mu, cholesky = cholesky, log = log )
#     }
#     
#     if (log) {
#       totalProb <- n.tail*log(1-pi) + dtail + dbulk
#     }else{
#       totalProb <- (1-pi)^n.tail *dtail*dbulk
#     }
#     return(totalProb)
#   })
# 
# dbiextmix <- nimbleFunction(
#   run = function(x=double(2), theta=double(1), thres=double(1), mu=double(1), 
#                  cholesky=double(2), D=integer(0, default=2),
#                  a.ind=double(1), lam.ind=double(0), lamfix=logical(0, default = 0), 
#                  sig.ind=double(1), gamma.ind=double(1),
#                  log = logical(0, default = 0)) {
#     returnType(double(0))
#     
#     cond <- (x[,1]>thres[1]) | (x[,2]>thres[2])
#     n.tail <- sum(cond)
#     n.bulk <- sum(!cond)
#     y.tail <- matrix(c(x[cond,1] - thres[1], x[cond,2] - thres[2]), ncol=D)
#     y.bulk <- x[!cond,]
#     
#     pi <- pmnorm_chol(lower=rep(-Inf,D), upper=thres, mean=mu, cholesky = cholesky)
#     
#     sig <- theta[sig.ind]
#     gamma <- theta[gamma.ind]
#     eta <- -sig/gamma
#     eta[which(gamma<=0)] <- -Inf
#     
#     dtail <- 0
#     dbulk <- 0
#     
#     if (n.tail>0){
#       y.min <- eta
#       for (i in 1:D){
#         y.min[i] <- min(y.tail[,i])
#       }
#       if (all(y.min>eta)){
#         llt <- -nim_nll_powunif_GPD_MAT(x=y.tail, theta=theta, u=min(y.tail)-0.01, a.ind=a.ind,
#                                         lam.ind=lam.ind, sig.ind=sig.ind, gamma.ind=gamma.ind, 
#                                         lamfix=lamfix, balthresh=FALSE, 
#                                         marg.scale.ind=1:2, marg.shape.ind=1:2)
#         if (log){
#           dtail <- llt
#         }else{
#           dtail <- exp(llt)
#         }
#       }else{
#         if (log) dtail <- -10^10
#       }
#     }
#     
#     if (n.bulk>0){
#       dbulk <- dmvnorm_chol(y.bulk, mean=mu, cholesky = cholesky, log = log )
#     }
#     
#     if (log) {
#       totalProb <- n.tail*log(1-pi) + dtail + dbulk
#     }else{
#       totalProb <- (1-pi)^n.tail *dtail*dbulk
#     }
#     
#     return(totalProb)
#   })
# 
# 
# 
# 
# rbiextmix_ns <- nimbleFunction(
#   run = function(n=integer(0), para.mg=double(1), beta.a=double(2), beta.b=double(2),
#                  X=double(2), thres=double(1), mu=double(1), 
#                  cholesky=double(2), D=integer(0, default=2),
#                  a.ind=double(1), lam.ind=double(0), lamfix=logical(0, default = 0), 
#                  sig.ind=double(1), gamma.ind=double(1)) {
#     returnType(double(2))
#     
#     totalProb <- matrix(1,nrow=n, ncol=D)
#     return(totalProb)
#   })
# 
# registerDistributions(list(
#   dbiextmix_ns = list(
#     BUGSdist = "dbiextmix_ns(para.mg, beta.a, beta.b, thres, X, mu, cholesky, D, 
#                   a.ind, lam.ind, lamfix, sig.ind, gamma.ind)",
#     types = c('value = double(2)', 'para.mg = double(1)', 'beta.a = double(2)', 
#               'beta.b = double(2)', 'thres = double(1)', 'X = double(2)',
#               'mu = double(1)', 'cholesky = double(2)', 'D = integer(0)', 'a.ind = double(1)', 
#               'lam.ind = double(0)', 'lamfix = logical(0)', 'sig.ind = double(1)',
#               'gamma.ind = double(1)')
#   )))
# 
# rbiextmix <- nimbleFunction(
#   run = function(n=integer(0), theta=double(1), thres=double(1), mu=double(1), 
#                  cholesky=double(2), D=integer(0, default=2),
#                  a.ind=double(1), lam.ind=double(0), lamfix=logical(0, default = 0), 
#                  sig.ind=double(1), gamma.ind=double(1)) {
#     returnType(double(2))
#     
#     totalProb <- matrix(1,nrow=n, ncol=D)
#     return(totalProb)
#   })
# 
# registerDistributions(list(
#   dbiextmix = list(
#     BUGSdist = "dbiextmix(theta, thres, mu, cholesky, D, a.ind, lam.ind, lamfix, sig.ind, gamma.ind)",
#     types = c('value = double(2)', 'theta = double(1)', 'thres = double(1)', 
#               'mu = double(1)', 'cholesky = double(2)', 'D = integer(0)', 'a.ind = double(1)', 
#               'lam.ind = double(0)', 'lamfix = logical(0)', 'sig.ind = double(1)',
#               'gamma.ind = double(1)')
#   )))
# 
# uppertri_mult_diag <- nimbleFunction(
#   run = function(mat = double(2), vec = double(1)) {
#     returnType(double(2))
#     p <- length(vec)
#     out <- matrix(nrow = p, ncol = p, init = FALSE)
#     for(i in 1:p)
#       out[ , i] <- mat[ , i] * vec[i]
#     return(out)
#   })
# 
# BivExtMixcode <- nimbleCode({
#   for (i in 1:D)
#     sds[i] ~ dunif(0, 100)
#   Ustar[1:D,1:D] ~ dlkj_corr_cholesky(1.3, D)
#   U[1:D,1:D] <- uppertri_mult_diag(Ustar[1:D, 1:D], sds[1:D])
#   
#   for (i in 1:D)
#     mu[i] ~ dnorm(0, sd=50)
#   
#   for (i in 1:D)
#     thres[i] ~ T(dnorm(mu.thres[i], sd=sd.thres[i]),min.thres[i],max.thres[i])
#   
#   if (station.ind){
#     # priors for sig
#     for (i in 1:5)
#       theta[i] ~ dunif(0,50)
#     # priors for gamma 
#     for (i in 6:7)
#       theta[i] ~ dunif(-1,1)
#     y[1:N,1:D] ~ dbiextmix(theta=theta[1:7], thres=thres[1:D], mu=mu[1:D], 
#                            cholesky=U[1:D,1:D], D=D,
#                            a.ind=a.ind[1:D], lam.ind=lam.ind, lamfix=lamfix, 
#                            sig.ind=sig.ind[1:D], gamma.ind=gamma.ind[1:D])
#   }else{
#     for (i in 1:2)
#       para.mg[i] ~ dunif(0,50)
#     for (i in 3:4)
#       para.mg[i] ~ dunif(-1,1)
#     for (i in 1:D){
#       beta.a[1:D.pred,i] ~ dmnorm(mu_beta[1:D.pred], cov=cov_beta[1:D.pred,1:D.pred])
#       beta.b[1:D.pred,i] ~ dmnorm(mu_beta[1:D.pred], cov=cov_beta[1:D.pred,1:D.pred])
#     }
#     y[1:N,1:D] ~ dbiextmix_ns(para.mg=para.mg[1:4], beta.a=beta.a[1:D.pred,1:D],
#                               beta.b=beta.b[1:D.pred,1:D], thres=thres[1:2],
#                               X=X[1:N,1:(D*D.pred)], mu=mu[1:D], 
#                               cholesky=U[1:D,1:D], D=D,
#                               a.ind=a.ind[1:D], lam.ind=lam.ind, lamfix=lamfix, 
#                               sig.ind=sig.ind[1:D], gamma.ind=gamma.ind[1:D])
#   }
#   
#   
# })
# 
# 
# dat <- Y
# BivExtMixmodel <- nimbleModel(BivExtMixcode, constants = list(N = nrow(dat), 
#                                                               D = 2,
#                                                               station.ind = TRUE,
#                                                               D.pred = 2,
#                                                               X = matrix(1,nrow=nrow(dat),ncol=2),
#                                                               mu.thres = apply(dat,2,quantile,0.9),#90 quantile
#                                                               sd.thres = c(10,10), 
#                                                               min.thres = apply(dat,2,quantile,0.8),
#                                                               max.thres = apply(dat,2,quantile,0.99),
#                                                               mu_beta = rep(0,2),
#                                                               cov_beta = 100*diag(2),
#                                                               a.ind = 1:2,
#                                                               lam.ind = 3,
#                                                               sig.ind = c(4,5),
#                                                               gamma.ind = c(6,7),
#                                                               lamfix=FALSE),check = FALSE)
# 
# 
# BivExtMixmodel$setData(list(y = dat))  ## Set those values as data in the model
# cBivExtMixmodel <- compileNimble(BivExtMixmodel)