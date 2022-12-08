source("KRSW/RevExp_U_Functions.r")
source("KRSW/CommonFunctions.r")
library(nimble, warn.conflicts = F)
library(mvtnorm)
library(parallel)

detectCores()
this_cluster <- makeCluster(3)

dir.in <- '/home/pgrad2/2448355h/My_PhD_Project/01_Output/Biv_Ext_Mix_Mod/biv_ext_mix_mod_simdat'

load(file.path(dir.in,'simulation_data_pos_shape_fixed.RData'))


X <- X.p09
#X <- X.p05.10t

# total 5%, not marginal 5%
u.x <- u.x.p09



# dat <- X
# seed <- 1
run_MCMC_parallel <- function(seed, dat){
  library(nimble)
  library(mvtnorm)
  source("KRSW/RevExp_U_Functions.r")
  source("KRSW/CommonFunctions.r")
  
  R_pmnorm_chol <- function(lower, upper, mean, cholesky){
    sigma <- t(cholesky) %*% cholesky
    return(pmvnorm(lower=lower, upper=as.vector(upper), mean=as.vector(mean), sigma=sigma, keepAttr = F))
  }
  assign('R_pmnorm_chol', R_pmnorm_chol, envir = .GlobalEnv)
  
  R_dmvnorm_chol <- function(x, mean, cholesky, log){
    sigma <- t(cholesky) %*% cholesky
    dvect <- dmvnorm(x, mean=mean, sigma=sigma, log=log)
    if (log) {
      return(sum(dvect))
    }else{
      return(prod(dvect))
    }
    
  }
  assign('R_dmvnorm_chol', R_dmvnorm_chol, envir = .GlobalEnv)
  
  pmnorm_chol <- nimbleRcall(function(lower = double(1), upper = double(1), 
                                      mean=double(1),cholesky=double(2)){}, 
                             Rfun = 'R_pmnorm_chol',
                             returnType = double(0))
  
  assign('pmnorm_chol', pmnorm_chol, envir = .GlobalEnv)
  
  dmvnorm_chol <- nimbleRcall(function(x = double(2), mean = double(1), 
                                       cholesky=double(2),log=logical(0, default = 0)){}, 
                              Rfun = 'R_dmvnorm_chol',
                              returnType = double(0))
  
  assign('dmvnorm_chol', dmvnorm_chol, envir = .GlobalEnv)
  
  #why switching theta and x will cause error
  assign('nll.powunif.GPD', nll.powunif.GPD, envir = .GlobalEnv)
  nim_nll_powunif_GPD <- nimbleRcall(function(theta=double(1), x=double(2), u=double(0), a.ind=double(0),
                                              lam.ind=double(0), sig.ind=double(1), gamma.ind=double(1), 
                                              lamfix=logical(0, default = 0), balthresh=logical(0, default = 0), 
                                              marg.scale.ind=double(1), marg.shape.ind=double(1)){}, 
                                     Rfun = 'nll.powunif.GPD',
                                     returnType = double(0))
  assign('nim_nll_powunif_GPD', nim_nll_powunif_GPD, envir = .GlobalEnv)
  
  dbiextmix <- nimbleFunction(
    run = function(x=double(2), theta=double(1), thres=double(1), mu=double(1), 
                   cholesky=double(2), D=integer(0, default=2),
                   a.ind=double(0), lam.ind=double(0), lamfix=logical(0, default = 0), 
                   sig.ind=double(1), gamma.ind=double(1),
                   log = logical(0, default = 0)) {
      returnType(double(0))
      
      cond <- (x[,1]>thres[1]) | (x[,2]>thres[2])
      n.tail <- sum(cond)
      n.bulk <- sum(!cond)
      y.tail <- matrix(c(x[cond,1] - thres[1], x[cond,2] - thres[2]), ncol=D)
      y.bulk <- x[!cond,]
      
      pi <- pmnorm_chol(lower=rep(0,D), upper=thres, mean=mu, cholesky = cholesky)
      
      sig <- theta[sig.ind]
      gamma <- theta[gamma.ind]
      eta <- -sig/gamma
      
      dtail <- 0
      dbulk <- 0
      
      if (n.tail>0){
        y.min <- eta
        for (i in 1:D){
          y.min[i] <- min(y.tail[,i])
        }
        if (all(y.min>eta)){
          llt <- -nim_nll_powunif_GPD(x=y.tail, theta=theta, u=min(y.tail)-0.01, a.ind=a.ind,
                                      lam.ind=lam.ind, sig.ind=sig.ind, gamma.ind=gamma.ind, 
                                      lamfix=lamfix, balthresh=FALSE, 
                                      marg.scale.ind=1:2, marg.shape.ind=1:2)
          if (log){
            dtail <- llt
          }else{
            dtail <- exp(llt)
          }
        }else{
          if (log) dtail <- -10^100
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
  assign('dbiextmix', dbiextmix, envir = .GlobalEnv)
  
  rbiextmix <- nimbleFunction(
    run = function(n=integer(0), theta=double(1), thres=double(1), mu=double(1), 
                   cholesky=double(2), D=integer(0, default=2),
                   a.ind=double(0), lam.ind=double(0), lamfix=logical(0, default = 0), 
                   sig.ind=double(1), gamma.ind=double(1),
                   log = logical(0, default = 0)) {
      returnType(double(2))
      
      totalProb <- matrix(1,nrow=1000, ncol=D)
      return(totalProb)
    })
  assign('rbiextmix', rbiextmix, envir = .GlobalEnv)

  registerDistributions(list(
    dbiextmix = list(
      BUGSdist = "dbiextmix(theta, thres, mu, cholesky, a.ind, lam.ind, lamfix, sig.ind, gamma.ind)",
      types = c('value = double(2)', 'theta = double(1)', 'thres = double(1)', 
                'mu = double(1)', 'cholesky = double(2)', 'D = integer(0)', 'a.ind = double(0)', 
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
  assign('uppertri_mult_diag', uppertri_mult_diag, envir = .GlobalEnv)
  
  BivExtMixcode <- nimbleCode({
    for (i in 1:D)
      sds[i] ~ dunif(0, 100)
    Ustar[1:D,1:D] ~ dlkj_corr_cholesky(1.3, D)
    U[1:D,1:D] <- uppertri_mult_diag(Ustar[1:D, 1:D], sds[1:D])
    # S[1,1] <- 100
    # S[1,2] <- 0
    # S[2,1] <- 0
    # S[2,2] <- 100
    # U[1:D,1:D] ~ dinvwish(S=S[1:2,1:2], df=3)
    
    for (i in 1:3)
      theta[i] ~ dunif(0,50)
    for (i in 4:5)
      theta[i] ~ dunif(0,1)
    
    for (i in 1:2)
      mu[i] ~ T(dnorm(0, sd=100),0, 4.6)
    
    for (i in 1:2)
      thres[i] ~ dunif(lb[i], ub[i])
    
    y[1:N,1:D] ~ dbiextmix(theta=theta[1:5], thres=thres[1:D], mu=mu[1:D], 
                           cholesky=U[1:D,1:D],
                           a.ind=a.ind, lam.ind=lam.ind, lamfix=lamfix, 
                           sig.ind=sig.ind[1:D], gamma.ind=gamma.ind[1:D])
  })
  
  
  # why D is unused?
  BivExtMixmodel <- nimbleModel(BivExtMixcode, constants = list(N = 2500, 
                                                                D = 2,
                                                                a.ind = 1,
                                                                lam.ind = 2,
                                                                sig.ind = c(2,3),
                                                                gamma.ind = c(4,5),
                                                                lamfix=TRUE,
                                                                lb = c(4.197212,4.979511),
                                                                ub = c(6.774084,7.069288)), check = FALSE)
  
  
  
  
  BivExtMixmodel$setData(list(y = dat))  ## Set those values as data in the model
  cBivExtMixmodel <- compileNimble(BivExtMixmodel)
  
  
  BivExtMixconf <- configureMCMC(BivExtMixmodel,
                                 enableWAIC = TRUE, time=TRUE)
  BivExtMixMCMC <- buildMCMC(BivExtMixconf)

  cBivExtMixMCMC <- compileNimble(BivExtMixMCMC, project = BivExtMixmodel)

  results <- runMCMC(cBivExtMixMCMC, niter = 25000,nburnin=5000,thin=10,
                     summary = TRUE, WAIC = TRUE,setSeed = seed)
  return(results)
}

# 
# run_MCMC_parallel(1,X.p09)
# 
# t1 <- Sys.time()
# chain_output <- lapply(X=1:3,FUN=run_MCMC_parallel, dat=X.p09)
# t2 <- Sys.time()
# 
# print(t2-t1)
# 

t1 <- Sys.time()
chain_output <- parLapply(cl = this_cluster, X = 1:3, 
                          fun = run_MCMC_parallel, 
                          dat = X.p09)
t2 <- Sys.time()

print(t2-t1)


stopCluster(this_cluster)



library(coda)
coda.options(user.layout=FALSE)
samples <- mcmc.list(mcmc(chain_output[[1]]$samples),
                     mcmc(chain_output[[2]]$samples),
                     mcmc(chain_output[[3]]$samples))

effectiveSize(samples)

autocorr.diag(samples)

autocorr.plot(samples)

par(mfrow = c(2, 2))
traceplot(samples,col=1:3)
densplot(samples)

gelman.diag(samples, confidence = 0.95, transform=FALSE, autoburnin=FALSE,
            multivariate=FALSE)



# --------------------------------------calcualte waic-----------------------------

R_pmnorm_chol <- function(lower, upper, mean, cholesky){
  sigma <- t(cholesky) %*% cholesky
  return(pmvnorm(lower=lower, upper=as.vector(upper), mean=as.vector(mean), sigma=sigma, keepAttr = F))
}
assign('R_pmnorm_chol', R_pmnorm_chol, envir = .GlobalEnv)

R_dmvnorm_chol <- function(x, mean, cholesky, log){
  sigma <- t(cholesky) %*% cholesky
  dvect <- dmvnorm(x, mean=mean, sigma=sigma, log=log)
  if (log) {
    return(sum(dvect))
  }else{
    return(prod(dvect))
  }
  
}
assign('R_dmvnorm_chol', R_dmvnorm_chol, envir = .GlobalEnv)

pmnorm_chol <- nimbleRcall(function(lower = double(1), upper = double(1), 
                                    mean=double(1),cholesky=double(2)){}, 
                           Rfun = 'R_pmnorm_chol',
                           returnType = double(0))

assign('pmnorm_chol', pmnorm_chol, envir = .GlobalEnv)

dmvnorm_chol <- nimbleRcall(function(x = double(2), mean = double(1), 
                                     cholesky=double(2),log=logical(0, default = 0)){}, 
                            Rfun = 'R_dmvnorm_chol',
                            returnType = double(0))

assign('dmvnorm_chol', dmvnorm_chol, envir = .GlobalEnv)

#why switching theta and x will cause error
assign('nll.powunif.GPD', nll.powunif.GPD, envir = .GlobalEnv)
nim_nll_powunif_GPD <- nimbleRcall(function(theta=double(1), x=double(2), u=double(0), a.ind=double(0),
                                            lam.ind=double(0), sig.ind=double(1), gamma.ind=double(1), 
                                            lamfix=logical(0, default = 0), balthresh=logical(0, default = 0), 
                                            marg.scale.ind=double(1), marg.shape.ind=double(1)){}, 
                                   Rfun = 'nll.powunif.GPD',
                                   returnType = double(0))
assign('nim_nll_powunif_GPD', nim_nll_powunif_GPD, envir = .GlobalEnv)

dbiextmix <- nimbleFunction(
  run = function(x=double(2), theta=double(1), thres=double(1), mu=double(1), 
                 cholesky=double(2), D=integer(0, default=2),
                 a.ind=double(0), lam.ind=double(0), lamfix=logical(0, default = 0), 
                 sig.ind=double(1), gamma.ind=double(1),
                 log = logical(0, default = 0)) {
    returnType(double(0))
    
    cond <- (x[,1]>thres[1]) | (x[,2]>thres[2])
    n.tail <- sum(cond)
    n.bulk <- sum(!cond)
    y.tail <- matrix(c(x[cond,1] - thres[1], x[cond,2] - thres[2]), ncol=D)
    y.bulk <- x[!cond,]
    
    pi <- pmnorm_chol(lower=rep(0,D), upper=thres, mean=mu, cholesky = cholesky)
    
    sig <- theta[sig.ind]
    gamma <- theta[gamma.ind]
    eta <- -sig/gamma
    
    dtail <- 0
    dbulk <- 0
    
    if (n.tail>0){
      y.min <- eta
      for (i in 1:D){
        y.min[i] <- min(y.tail[,i])
      }
      if (all(y.min>eta)){
        llt <- -nim_nll_powunif_GPD(x=y.tail, theta=theta, u=min(y.tail)-0.01, a.ind=a.ind,
                                    lam.ind=lam.ind, sig.ind=sig.ind, gamma.ind=gamma.ind, 
                                    lamfix=lamfix, balthresh=FALSE, 
                                    marg.scale.ind=1:2, marg.shape.ind=1:2)
        if (log){
          dtail <- llt
        }else{
          dtail <- exp(llt)
        }
      }else{
        if (log) dtail <- -10^100
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
assign('dbiextmix', dbiextmix, envir = .GlobalEnv)

rbiextmix <- nimbleFunction(
  run = function(n=integer(0), theta=double(1), thres=double(1), mu=double(1), 
                 cholesky=double(2), D=integer(0, default=2),
                 a.ind=double(0), lam.ind=double(0), lamfix=logical(0, default = 0), 
                 sig.ind=double(1), gamma.ind=double(1),
                 log = logical(0, default = 0)) {
    returnType(double(2))
    
    totalProb <- matrix(1,nrow=1000, ncol=D)
    return(totalProb)
  })
assign('rbiextmix', rbiextmix, envir = .GlobalEnv)

registerDistributions(list(
  dbiextmix = list(
    BUGSdist = "dbiextmix(theta, thres, mu, cholesky, a.ind, lam.ind, lamfix, sig.ind, gamma.ind)",
    types = c('value = double(2)', 'theta = double(1)', 'thres = double(1)', 
              'mu = double(1)', 'cholesky = double(2)', 'D = integer(0)', 'a.ind = double(0)', 
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
assign('uppertri_mult_diag', uppertri_mult_diag, envir = .GlobalEnv)

BivExtMixcode <- nimbleCode({
  for (i in 1:D)
    sds[i] ~ dunif(0, 100)
  Ustar[1:D,1:D] ~ dlkj_corr_cholesky(1.3, D)
  U[1:D,1:D] <- uppertri_mult_diag(Ustar[1:D, 1:D], sds[1:D])
  # S[1,1] <- 100
  # S[1,2] <- 0
  # S[2,1] <- 0
  # S[2,2] <- 100
  # U[1:D,1:D] ~ dinvwish(S=S[1:2,1:2], df=3)
  
  for (i in 1:3)
    theta[i] ~ dunif(0,50)
  for (i in 4:5)
    theta[i] ~ dunif(0,1)
  
  for (i in 1:2)
    mu[i] ~ T(dnorm(0, sd=100),0, 4.6)
  
  for (i in 1:2)
    thres[i] ~ dunif(lb[i], ub[i])
  
  y[1:N,1:D] ~ dbiextmix(theta=theta[1:5], thres=thres[1:D], mu=mu[1:D], 
                         cholesky=U[1:D,1:D],
                         a.ind=a.ind, lam.ind=lam.ind, lamfix=lamfix, 
                         sig.ind=sig.ind[1:D], gamma.ind=gamma.ind[1:D])
})


# why D is unused?
BivExtMixmodel <- nimbleModel(BivExtMixcode, constants = list(N = 2500, 
                                                              D = 2,
                                                              a.ind = 1,
                                                              lam.ind = 2,
                                                              sig.ind = c(2,3),
                                                              gamma.ind = c(4,5),
                                                              lamfix=TRUE,
                                                              lb = c(4.197212,4.979511),
                                                              ub = c(6.774084,7.069288)), check = FALSE)



BivExtMixcode <- nimbleCode({
  for (i in 1:D)
    sds[i] ~ dunif(0, 100)
  Ustar[1:D,1:D] ~ dlkj_corr_cholesky(1.3, D)
  U[1:D,1:D] <- uppertri_mult_diag(Ustar[1:D, 1:D], sds[1:D])
  # S[1,1] <- 100
  # S[1,2] <- 0
  # S[2,1] <- 0
  # S[2,2] <- 100
  # U[1:D,1:D] ~ dinvwish(S=S[1:2,1:2], df=3)
  
  for (i in 1:3)
    theta[i] ~ dunif(0,50)
  for (i in 4:5)
    theta[i] ~ dunif(0,1)
  
  for (i in 1:2)
    mu[i] ~ T(dnorm(0, sd=100),0, 4.6)
  
  for (i in 1:2)
    thres[i] ~ dunif(lb[i], ub[i])
  
  y[1:N,1:D] ~ dbiextmix(theta=theta[1:5], thres=thres[1:D], mu=mu[1:D], 
                         cholesky=U[1:D,1:D],
                         a.ind=a.ind, lam.ind=lam.ind, lamfix=lamfix, 
                         sig.ind=sig.ind[1:D], gamma.ind=gamma.ind[1:D])
})


# why D is unused?
BivExtMixmodel <- nimbleModel(BivExtMixcode, constants = list(N = 2500, 
                                                              D = 2,
                                                              a.ind = 1,
                                                              lam.ind = 2,
                                                              sig.ind = c(2,3),
                                                              gamma.ind = c(4,5),
                                                              lamfix=TRUE,
                                                              lb = c(4.197212,4.979511),
                                                              ub = c(6.774084,7.069288)), check = FALSE)

BivExtMixmodel$setData(list(y = X))  ## Set those values as data in the model
cBivExtMixmodel <- compileNimble(BivExtMixmodel)

#different
calculateWAIC(chain_output[[1]]$samples, BivExtMixmodel)
chain_output[[1]]$WAIC


samples.all <- rbind(chain_output[[1]]$samples,
                     chain_output[[2]]$samples,
                     chain_output[[3]]$samples)
calculateWAIC(samples.all, BivExtMixmodel)

