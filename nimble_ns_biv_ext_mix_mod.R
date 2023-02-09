source("KRSW/RevExp_U_Functions.r")
source("KRSW/CommonFunctions.r")
library(nimble, warn.conflicts = F)
library(mvtnorm)
library(parallel)



dir.in <- '/home/pgrad2/2448355h/My_PhD_Project/01_Output/Biv_Ext_Mix_Mod/biv_ext_mix_mod_simdat'

load(file.path(dir.in,'simulation_data_non_stationary.RData'))


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

R_map_thres <- function(lower, upper, eta){
  return(lower + sweep(sigmoid(eta),2, upper-lower, "*"))
}

R_mix_prob_sg <- function(upper, mu, cholesky){
  D <- length(mu)
  return(pmnorm_chol(lower=rep(0,D), upper=upper, mean=mu, cholesky = cholesky))
}

R_mix_prob <- function(thres, mu, cholesky){
  prob.all <- apply(thres,1,R_mix_prob_sg, mu=mu, cholesky=cholesky)
  return(sum(log(1-prob.all)))
}

R_mix_prob_1 <- function(thres,mu,cholesky){
  prob.all <- sapply(seq_len(nrow(thres)), function(i) R_mix_prob_sg(thres[i,],mu,cholesky))
  return(sum(log(1-prob.all)))
}

R_mix_prob_2 <- function(thres, mu, cholesky){
  prob.all <- mclapply(seq_len(nrow(thres)), function(i) R_mix_prob_sg(thres[i,],mu,cholesky))
  return(sum(log(1-unlist(prob.all))))
}


system.time(R_mix_prob(thres, mu, cholesky))
t1 <- Sys.time()
system.time(R_mix_prob_1(thres, mu, cholesky))
t2 <- Sys.time()
system.time(R_mix_prob_2(thres, mu, cholesky))
t3 <- Sys.time()
print(c(t3-t2,t2-t1))
# system.time(R_mix_prob_sg(thres[1,], mu, cholesky))

R_pmax <- function(mat, eta){
  return(sweep(mat, 2, eta, 'pmax'))
}

mix_prob_sum <- nimbleRcall(function(thres = double(2), mu = double(1),
                                 cholesky = double(2)){}, 
                        Rfun = 'R_mix_prob',
                        returnType = double(0))


map_thres <- nimbleRcall(function(lower = double(1), upper = double(1),
                                  eta = double(2)){}, 
                         Rfun = 'R_map_thres',
                         returnType = double(2))


pmnorm_chol <- nimbleRcall(function(lower = double(1), upper = double(1), 
                                    mean=double(1),cholesky=double(2)){}, 
                           Rfun = 'R_pmnorm_chol',
                           returnType = double(0))



dmvnorm_chol <- nimbleRcall(function(x = double(2), mean = double(1), 
                                     cholesky=double(2),log=logical(0, default = 0)){}, 
                            Rfun = 'R_dmvnorm_chol',
                            returnType = double(0))


nim_nll_powunif_GPD <- nimbleRcall(function(theta=double(1), x=double(2), u=double(0), a.ind=double(0),
                                            lam.ind=double(0), sig.ind=double(1), gamma.ind=double(1), 
                                            lamfix=logical(0, default = 0), balthresh=logical(0, default = 0), 
                                            marg.scale.ind=double(1), marg.shape.ind=double(1)){}, 
                                   Rfun = 'nll.powunif.GPD',
                                   returnType = double(0))

nim_nrow <- nimbleRcall(function(x = double(2)){}, 
                        Rfun = 'nrow',
                        returnType = double(0))

nim_pmax <- nimbleRcall(function(mat = double(2), eta = double(1)){}, 
                        Rfun = 'R_pmax',
                        returnType = double(2))

dbiextmix <- nimbleFunction(
  run = function(x=double(2), theta=double(1), beta=double(2), X=double(2),
                 lower= double(1), upper= double(1),mu=double(1), 
                 cholesky=double(2), D=integer(0, default=2),
                 a.ind=double(0), lam.ind=double(0), lamfix=logical(0, default = 0), 
                 sig.ind=double(1), gamma.ind=double(1),
                 log = logical(0, default = 0)) {
    returnType(double(0))
    
    lp <-  X %*% beta
    thres <- map_thres(lower,upper,lp)
    cond <- (x[,1]>thres[,1]) | (x[,2]>thres[,2])
    n.tail <- sum(cond)
    n.bulk <- sum(!cond)
    y.tail <- matrix(c(x[cond,1] - thres[cond,1], x[cond,2] - thres[cond,2]), ncol=D)
    y.bulk <- x[!cond,]
    
    log.sum <- mix_prob_sum(thres,mu,cholesky)

    sig <- theta[sig.ind]
    gamma <- theta[gamma.ind]
    eta <- -sig/gamma
    
    dtail <- 0
    dbulk <- 0
    
    if (n.tail>0){
      # max of y.tail and endpoint
      y.tail.max <- nim_pmax(y.tail,eta+10^-10)
      llt <- -nim_nll_powunif_GPD(x=y.tail.max, theta=theta, u=min(y.tail)-0.01, a.ind=a.ind,
                                  lam.ind=lam.ind, sig.ind=sig.ind, gamma.ind=gamma.ind, 
                                  lamfix=lamfix, balthresh=FALSE, 
                                  marg.scale.ind=1:2, marg.shape.ind=1:2)
      if (log){
        dtail <- llt
      }else{
        dtail <- exp(llt)
      }
      
    }
    
    if (n.bulk>0){
      dbulk <- dmvnorm_chol(y.bulk, mean=mu, cholesky = cholesky, log = log )
    }
    
    if (log) {
      totalProb <- log.sum + dtail + dbulk
    }else{
      totalProb <- exp(log.sum) *dtail*dbulk
    }
    
    return(totalProb)
  })






theta <- c(1.656, 0.571, 0.451, 0.253, 0.035)
a.ind <- c(1)
lam.ind <- c(1)
sig.ind <- c(2,3)
gamma.ind <- c(4,5)
marg.scale.ind <- c(1,2)
marg.shape.ind <- c(1,2)
beta <- 0.1*cbind(c(1,2,3,4),c(-2,-3,4,5))
lower <- 5
upper <- 8
mu <- c(5,5.41)
rho=0.5
sigma <- 1.5* matrix(c(1,rho,rho,0.8),ncol=2)
cholesky <- chol(sigma)
x <- Y
D <- 2
log <- TRUE
lamfix=FALSE
dbiextmix(x=Y, theta=theta, beta=beta, X=X,
          lower= lower, upper= upper,mu=mu, 
          cholesky=cholesky,
          a.ind=a.ind, lam.ind=lam.ind, lamfix=0, 
          sig.ind=sig.ind, gamma.ind=gamma.ind,
          log =1)



rbiextmix <- nimbleFunction(
  run = function(n=integer(0), theta=double(1), beta=double(2), X=double(2), 
                 lower= double(1), upper= double(1), mu=double(1), 
                 cholesky=double(2), D=integer(0, default=2),
                 a.ind=double(0), lam.ind=double(0), lamfix=logical(0, default = 0), 
                 sig.ind=double(1), gamma.ind=double(1)) {
    returnType(double(2))
    
    totalProb <- matrix(1,nrow=n, ncol=D)
    return(totalProb)
  })
         
registerDistributions(list(
  dbiextmix = list(
    BUGSdist = "dbiextmix(theta, beta, X, lower, upper, mu, cholesky, D, 
                a.ind, lam.ind, lamfix, sig.ind, gamma.ind)",
    types = c('value = double(2)', 'theta = double(1)', 'beta = double(2)', 'X = double(2)',
              'lower = double(1)', 'upper = double(1)',
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
  
  for (i in 1:3)
    theta[i] ~ dunif(0,50)
  for (i in 4:5)
    theta[i] ~ dunif(0,1)
  for (i in 1:2)
    mu[i] ~ T(dnorm(0, sd=100),0, lower[i])
  
  for (i in 1:D)
    mu_beta[i] <- 0
  for (i in 1:D.pred)
    beta[i,1:D] ~ dmnorm(mu_beta[1:D], cov=cov_beta[1:D,1:D])
  
  y[1:N,1:D] ~ dbiextmix(theta=theta[1:5], beta=beta[1:D.pred,1:D],X=X[1:N,1:D.pred],
                         lower=lower[1:D], upper=upper[1:D], mu=mu[1:D], 
                         cholesky=U[1:D,1:D], D=D,
                         a.ind=a.ind, lam.ind=lam.ind, lamfix=lamfix, 
                         sig.ind=sig.ind[1:D], gamma.ind=gamma.ind[1:D])
})

BivExtMixmodel <- nimbleModel(BivExtMixcode, constants = list(N = 2500, 
                                                              D = 2,
                                                              D.pred = 4,
                                                              cov_beta = 1000*diag(2),
                                                              X = X,
                                                              lower = c(5,5),
                                                              upper = c(8,8),
                                                              a.ind = 1,
                                                              lam.ind = 2,
                                                              sig.ind = c(2,3),
                                                              gamma.ind = c(4,5),
                                                              lamfix=1),
                                                              check = FALSE)


BivExtMixmodel$setData(list(y = Y))  ## Set those values as data in the model
cBivExtMixmodel <- compileNimble(BivExtMixmodel)

BivExtMixconf <- configureMCMC(BivExtMixmodel,
                               enableWAIC = TRUE, time=TRUE)
BivExtMixMCMC <- buildMCMC(BivExtMixconf)

cBivExtMixMCMC <- compileNimble(BivExtMixMCMC, project = BivExtMixmodel)

t1 <- Sys.time()
results <- runMCMC(cBivExtMixMCMC, niter = 5000,nburnin=2500,thin=1,
                   summary = TRUE, WAIC = TRUE,setSeed = 1234)
t2 <- Sys.time()
print(t2-t1)



