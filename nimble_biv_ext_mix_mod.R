source("KRSW/RevExp_U_Functions.r")
source("KRSW/CommonFunctions.r")
library(nimble, warn.conflicts = F)
library(mvtnorm)


dir.in <- '/home/pgrad2/2448355h/My_PhD_Project/01_Output/Biv_Ext_Mix_Mod/biv_ext_mix_mod_simdat'

load(file.path(dir.in,'simulation_data_pos_shape_fixed.RData'))


X <- X.p09
#X <- X.p05.10t

# total 5%, not marginal 5%
u.x <- u.x.p09

R_pmnorm_chol <- function(lower, upper, mean, cholesky){
  sigma <- t(cholesky) %*% cholesky
  return(pmvnorm(lower=lower, upper=as.vector(upper), mean=as.vector(mean), sigma=sigma, keepAttr = F))
}



pmnorm_chol <- nimbleRcall(function(lower = double(1), upper = double(1), 
                                        mean=double(1),cholesky=double(2)){}, 
                               Rfun = 'R_pmnorm_chol',
                               returnType = double(0))

R_dmvnorm_chol <- function(x, mean, cholesky, log){
  sigma <- t(cholesky) %*% cholesky
  dvect <- dmvnorm(x, mean=mean, sigma=sigma, log=log)
  if (log) {
    return(sum(dvect))
  }else{
    return(prod(dvect))
  }

}


dmvnorm_chol <- nimbleRcall(function(x = double(2), mean = double(1), 
                                     cholesky=double(2),log=logical(0, default = 0)){}, 
                            Rfun = 'R_dmvnorm_chol',
                            returnType = double(0))

# pmnorm_chol(c(0,0), thres, mean=c(3,4), cholesky = diag(2))

#why switching theta and x will cause error
nim_nll_powunif_GPD <- nimbleRcall(function(theta=double(1), x=double(2), u=double(0), a.ind=double(0),
                                            lam.ind=double(0), sig.ind=double(1), gamma.ind=double(1), 
                                            lamfix=logical(0, default = 0), balthresh=logical(0, default = 0), 
                                            marg.scale.ind=double(1), marg.shape.ind=double(1)){}, 
                                   Rfun = 'nll.powunif.GPD',
                                   returnType = double(0))
# 
# thres <- c(5.55,6.213)
# theta <- c(1.656, 0.571, 0.451, 0.253, 0.035)
# a.ind <- c(1)
# lam.ind <- c(1)
# sig.ind <- c(2,3)
# gamma.ind <- c(4,5)
# marg.scale.ind <- c(1,2)
# marg.shape.ind <- c(1,2)
# 
# cond <- (X[,1]>thres[1]) | (X[,2]>thres[2])
# y.tail <- cbind(X[cond,1] - thres[1],
#                 X[cond,2] - thres[2])
# y.single <- y.tail[1,]
# 
# y.bulk <- X[!cond,]

# rbind(y.single,y.single)
# nll.powunif.GPD(theta=theta, x=rbind(y.single,y.single), u=min(y.tail)-0.01, a.ind=a.ind, lam.ind=lam.ind, sig.ind=sig.ind,
#                 gamma.ind=gamma.ind, marg.scale.ind=marg.scale.ind, 
#                 marg.shape.ind = marg.shape.ind)

# nim_nll_powunif_GPD( theta=theta,x=rbind(y.single,y.single), u=min(y.tail)-0.01, a.ind=a.ind, lam.ind=lam.ind, sig.ind=sig.ind,
#                 gamma.ind=gamma.ind, marg.scale.ind=marg.scale.ind, 
#                 marg.shape.ind = marg.shape.ind)


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

# 
# dat <- X[,]
# 
# dbiextmix1(dat,theta=theta, thres=thres, mu=c(3,4),
#           cholesky=diag(2), D=2,
#           a.ind=a.ind, lam.ind=lam.ind, lamfix=FALSE,
#           sig.ind=sig.ind, gamma.ind=gamma.ind,log=TRUE)
# 
# res <- c()
# 
# for (i in 1:nrow(dat)){
#    res <- c(res,dbiextmix(dat[i,],theta=theta, thres=thres, mu=c(3,4),
#              cholesky=diag(2), D=2,
#              a.ind=a.ind, lam.ind=lam.ind, lamfix=FALSE,
#              sig.ind=sig.ind, gamma.ind=gamma.ind,log=TRUE))
#    # if (i>=2) {
#    #   res2 <- dbiextmix1(y.tail[1:i,],theta=theta, thres=thres, mu=c(3,4),
#    #                      cholesky=diag(2), D=2,
#    #                      a.ind=a.ind, lam.ind=lam.ind, lamfix=FALSE,
#    #                      sig.ind=sig.ind, gamma.ind=gamma.ind,log=TRUE)
#    #   if (sum(res)!=res2) break
#    # }
# }
# print(sum(res))


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


BivExtMixconf <- configureMCMC(BivExtMixmodel,
                               enableWAIC = TRUE, time=TRUE)
BivExtMixMCMC <- buildMCMC(BivExtMixconf)
# run it to get R error report
options(error=NULL)

# t1 <- Sys.time()
# BivExtMixMCMC$run(1)
# t2 <- Sys.time()
# 
# print(t2-t1)


cBivExtMixMCMC <- compileNimble(BivExtMixMCMC, project = BivExtMixmodel)

t1 <- Sys.time()
samples <- runMCMC(cBivExtMixMCMC, niter = 30000,nburnin=20000,thin=5,summary = TRUE, WAIC = TRUE)
t2 <- Sys.time()

print(t2-t1)
plot(samples[,'thres[2]'], type = 'l', main = 'thres[2] trace plot')

calculateWAIC(cBivExtMixMCMC)
calculateWAIC(samples, BivExtMixmodel)
