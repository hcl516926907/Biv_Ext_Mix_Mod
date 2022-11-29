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
  return(pmvnorm(lower=lower, upper=upper, mean=mean, sigma=sigma, keepAttr = F))
}



pmnorm_chol <- nimbleRcall(function(lower = double(1), upper = double(1), 
                                        mean=double(1),cholesky=double(2)){}, 
                               Rfun = 'R_pmnorm_chol',
                               returnType = double(0))

# pmnorm_chol(c(0,0), thres, mean=c(3,4), cholesky = diag(2))

#why switching theta and x will cause error
nim_nll_powunif_GPD <- nimbleRcall(function(theta=double(1), x=double(2), u=double(0), a.ind=double(0),
                                            lam.ind=double(0), sig.ind=double(1), gamma.ind=double(1), 
                                            lamfix=logical(0, default = 0), balthresh=logical(0, default = 0), 
                                            marg.scale.ind=double(1), marg.shape.ind=double(1)){}, 
                                   Rfun = 'nll.powunif.GPD',
                                   returnType = double(0))

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
# rbind(y.single,y.single)
# nll.powunif.GPD(theta=theta, x=rbind(y.single,y.single), u=min(y.tail)-0.01, a.ind=a.ind, lam.ind=lam.ind, sig.ind=sig.ind,
#                 gamma.ind=gamma.ind, marg.scale.ind=marg.scale.ind, 
#                 marg.shape.ind = marg.shape.ind)

# nim_nll_powunif_GPD( theta=theta,x=rbind(y.single,y.single), u=min(y.tail)-0.01, a.ind=a.ind, lam.ind=lam.ind, sig.ind=sig.ind,
#                 gamma.ind=gamma.ind, marg.scale.ind=marg.scale.ind, 
#                 marg.shape.ind = marg.shape.ind)


dbiextmix <- nimbleFunction(
  run = function(x=double(1), theta=double(1), thres=double(1), mu=double(1), 
                 chol=double(2), d=integer(0, default=2),
                 a.ind=double(0), lam.ind=double(0), lamfix=logical(0, default = 0), 
                 sig.ind=double(1), gamma.ind=double(1),
                 log = logical(0, default = 0)) {
    returnType(double())
    
    tail.ind <- any(x>thres)
    
    pi <- pmnorm_chol(lower=rep(0,d), upper=thres, mean=mu, cholesky = chol)
    
    sig <- theta[sig.ind]
    gamma <- theta[gamma.ind]
    eta <- -sig/gamma
    
    dtail <- 0
    dbulk <- 0
    
    if(tail.ind){
      if (all((x-thres)>eta)){
        y.tail <- t(matrix(c(x-thres, x-thres),nrow=d))
        twollt <- -nim_nll_powunif_GPD(x=y.tail, theta=theta, u=min(y.tail)-0.01, a.ind=a.ind,
                                      lam.ind=lam.ind, sig.ind=sig.ind, gamma.ind=gamma.ind, 
                                      lamfix=lamfix, balthresh=FALSE, 
                                      marg.scale.ind=1:2, marg.shape.ind=1:2)
        dtail <- exp(twollt/2)
      }
    }else{
    
    dbulk <- dmnorm_chol(x, mean=mu, cholesky = chol, prec_param = FALSE)
    }
    
    totalProb <- (1-pi)*dtail + dbulk
    if (log) return(log(totalProb))
    return(totalProb)
  })

# dbiextmix(y.single+thres,theta=theta, thres=thres, mu=c(3,4),
#           chol=diag(2), d=2,
#           a.ind=a.ind, lam.ind=lam.ind, lamfix=FALSE,
#           sig.ind=sig.ind, gamma.ind=gamma.ind)


registerDistributions(list(
  dbiextmix = list(
    BUGSdist = "dbiextmix(theta, thres, mu, chol, a.ind, lam.ind, lamfix, sig.ind, gamma.ind)",
    types = c('value = double(1)', 'theta = double(1)', 'thres = double(1)', 
              'mu = double(1)', 'chol = double(2)', 'd = integer(0)', 'a.ind = double(0)', 
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

bivextmixcode <- nimbleCode({
  Ustar[1:d,1:d] ~ dlkj_corr_cholesky(1.3, d)
  U[1:d,1:d] <- uppertri_mult_diag(Ustar[1:d, 1:d], sds[1:d])
  
  for (i in 1:3)
    theta[i] ~ dunif(0,50)
  for (i in 4:5)
    theta[i] ~ dunif(0,1)
  
  for (i in 1:2)
    mu[i] ~ dnorm(0,100)
  
  for (i in 1:2)
    thres[i] ~ dunif(lb[i], ub[i])

  for (i in 1:N)
    y[i,1:d] ~ dbiextmix(theta=theta[1:5], thres=thres[1:d], mu=mu[1:d], 
                         chol=U[1:d,1:d],
                         a.ind=a.ind, lam.ind=lam.ind, lamfix=lamfix, 
                         sig.ind=sig.ind[1:d], gamma.ind=gamma.ind[1:d])
})


# why d is unused?
biextmixmodel <- nimbleModel(bivextmixcode, constants = list(N = 2500, 
                                                             d = 2,
                                                             a.ind = 1,
                                                             lam.ind = 2,
                                                             sig.ind = c(2,3),
                                                             gamma.ind = c(4,5),
                                                             lamfix=TRUE,
                                                             lb = c(4.197212,4.979511),
                                                             ub = c(6.774084,7.069288)), check = FALSE)




biextmixmodel$setData(list(y = X))  ## Set those values as data in the model
cbiextmixmodel <- compileNimble(biextmixmodel)


dSS <- nimbleFunction(
  run = function(x = double(1), mu = double(1), cov=double(2), theta=double(0), log = logical(0, default = 0)) {
    returnType(double())
    cov.cho <- chol(cov)
    comp1 <- theta*dmnorm_chol(x, mean=mu[1:2], cholesky=cov.cho, prec_param = FALSE)
    comp2 <- (1-theta)*dmnorm_chol(x, mean=mu[3:4], cholesky=cov.cho, prec_param = FALSE)
    
    totalProb <- comp1 + comp2
    if (log) return(log(totalProb))
    return(totalProb)
  })


test <- nimbleFunction(
  run = function(lower = double(1), upper = double(1), mean=double(1),cholesky=double(2)) {
    returnType(double(0))
    return(pmnorm_chol(lower=lower, upper=upper, mean=mean, cholesky = cholesky))
  })

pmvnorm(lower=rep(0,2), upper=c(5,5), mean=c(3,3), sigma=diag(2), keepAttr = F)
test(lower=rep(0,2), upper=c(5,5), mean=c(3,3), cholesky=diag(2))
## Not run: 
## Say we want an R function that adds 2 to every value in a vector
add2 <- function(x) {
  x + 2 
}
Radd2 <- nimbleRcall(function(x = double(1)){}, Rfun = 'add2',
                     returnType = double(1))
demoCode <- nimbleCode({
  for(i in 1:4) {x[i] ~ dnorm(0,1)} 
  z[1:4] <- Radd2(x[1:4])
})
demoModel <- nimbleModel(demoCode, inits = list(x = rnorm(4)),
                         check = FALSE, calculate = FALSE)
CdemoModel <- compileNimble(demoModel)