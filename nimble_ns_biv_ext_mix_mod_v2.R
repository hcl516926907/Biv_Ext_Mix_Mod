source("KRSW/RevExp_U_Functions.r")
source("KRSW/CommonFunctions.r")
library(nimble, warn.conflicts = F)
library(mvtnorm)
library(parallel)


dir.in <- '/home/pgrad2/2448355h/My_PhD_Project/01_Output/Biv_Ext_Mix_Mod/biv_ext_mix_mod_simdat'

load(file.path(dir.in,'simulation_data_non_stationary_extr_dep_seed_1235_N_2500.RData'))


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
  if (is.null(dim(x))){
    d <- length(x)
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

x <- c(-10.58, 2.07)
u <- min(x)-0.01
a.par <- c(3,4)
lam <- 5
sig<-c(0.571,0.451)
gamma<-c(0.253, 0.135)
theta <- c(a.par, lam, sig, gamma)
a.ind <- c(1,2)
lam.ind <- c(3)
sig.ind <- c(4,5)
gamma.ind <- c(6,7)
marg.scale.ind <- c(1,2)
marg.shape.ind <- c(1,2)
lamfix=FALSE
balthresh=FALSE

nll.powunif.GPD.1(theta,x,u,a.ind,lam.ind,sig.ind,gamma.ind, lamfix=FALSE, balthresh=FALSE, marg.scale.ind,marg.shape.ind)



nim_nll_powunif_GPD.SG <- nimbleRcall(function(theta=double(1), x=double(1), u=double(0), a.ind=double(1),
                                            lam.ind=double(0), sig.ind=double(1), gamma.ind=double(1), 
                                            lamfix=logical(0, default = 0), balthresh=logical(0, default = 0), 
                                            marg.scale.ind=double(1), marg.shape.ind=double(1)){}, 
                                   Rfun = 'nll.powunif.GPD.1',
                                   returnType = double(0))

nim_nll_powunif_GPD.MAT <- nimbleRcall(function(theta=double(1), x=double(2), u=double(0), a.ind=double(1),
                                               lam.ind=double(0), sig.ind=double(1), gamma.ind=double(1), 
                                               lamfix=logical(0, default = 0), balthresh=logical(0, default = 0), 
                                               marg.scale.ind=double(1), marg.shape.ind=double(1)){}, 
                                      Rfun = 'nll.powunif.GPD.1',
                                      returnType = double(0))

x<- rbind(c(1,2),3,4)
theta <- c(2, 3, 1, 0.571, 0.451, 0.253 ,0.035)
u <- min(x)-0.01
a.ind <- c(1,2)
lam.ind <- c(3)
sig.ind <- c(4,5)
gamma.ind <- c(6,7)
marg.scale.ind <- c(1,2)
marg.shape.ind <- c(1,2)
lamfix=FALSE
balthresh=FALSE
nim_nll_powunif_GPD.SG(x=x, theta=theta, u=u, a.ind=a.ind,
                    lam.ind=lam.ind, sig.ind=sig.ind, gamma.ind=gamma.ind, 
                    lamfix=lamfix, balthresh=FALSE, 
                    marg.scale.ind=1:2, marg.shape.ind=1:2)


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
    
    pi <- pmnorm_chol(lower=rep(0,D), upper=thres, mean=mu, cholesky = cholesky)
    
    sig <- para.mg[1:D]
    gamma <- para.mg[(D+1):(2*D)]
    eta <- -sig/gamma
    
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
          llt.sg <- -nim_nll_powunif_GPD.SG(x=y.tail[i,], theta=theta, u=min(y.tail)-0.01, a.ind=a.ind,
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
        llt <- -nim_nll_powunif_GPD.MAT(x=y.tail, theta=theta, u=min(y.tail)-0.01, a.ind=a.ind,
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

test1 <- function(){
  -nim_nll_powunif_GPD(x=y.tail[i,], theta=theta, u=min(y.tail)-0.01, a.ind=a.ind,
                       lam.ind=lam.ind, sig.ind=sig.ind, gamma.ind=gamma.ind, 
                       lamfix=lamfix, balthresh=FALSE, 
                       marg.scale.ind=1:2, marg.shape.ind=1:2)
}

X <- cbind(X1,X2)
para.mg <- c(0.571, 0.451, 0.253, 0.135)
beta.a <- cbind(c(0.5, 0.2),c(0.5, -0.4))
beta.b <- cbind(c(0.25,-0.25),c(0,0))
# beta.a <- cbind(c(0.5, 0),c(-0.4, 0))
# beta.b <- cbind(c(0.25,0),c(0,0))
a.ind <- c(1,2)
lam.ind <- c(3)
sig.ind <- c(4,5)
gamma.ind <- c(6,7)
marg.scale.ind <- c(1,2)
marg.shape.ind <- c(1,2)
thres <- c(12,12)
mu <- c(5,5.41)
rho=0.5
sigma <- 1.5* matrix(c(1,rho,rho,0.8),ncol=2)
cholesky <- chol(sigma)
x <- Y
D <- 2
log <- TRUE
lamfix=FALSE
var.mg <- diag(c(1/sqrt(1.5),1/sqrt(1.2)))
chol.corr <- chol(var.mg %*% sigma %*% var.mg)

t1 <- Sys.time()

dbiextmix_ns(x=Y, para.mg=para.mg, beta.a=beta.a, beta.b=beta.b, X=X,
          thres=thres, mu=mu,
          cholesky=cholesky,
          a.ind=a.ind, lam.ind=lam.ind, lamfix=0,
          sig.ind=sig.ind, gamma.ind=gamma.ind,
          log =1)

# -7339.51

# -7223.138
t2 <- Sys.time()
print(t2-t1)


######################poster sampling debug###################
index <- 1:22
para.mg <- colMeans(results$samples[index,c('para.mg[1]', 'para.mg[2]',
                                       'para.mg[3]', 'para.mg[4]')])
beta.a <- matrix(colMeans(results$samples[index,c('beta.a[1, 1]', 'beta.a[2, 1]',
                                             'beta.a[1, 2]', 'beta.a[2, 2]')]), ncol=2)
# beta.a <- matrix(rnorm(4),ncol=2)
beta.b <- matrix(colMeans(results$samples[index,c('beta.b[1, 1]', 'beta.b[2, 1]',
                                             'beta.b[1, 2]', 'beta.b[2, 2]')]), ncol=2)
beta.b[,2] <- rep(0,2)
thres <- colMeans(results$samples[index,c('thres[1]','thres[2]')])
dbiextmix_ns(x=Y, para.mg=para.mg, beta.a=beta.a, beta.b=beta.b, X=X,
          thres=thres, mu=mu,
          cholesky=cholesky,
          a.ind=a.ind, lam.ind=lam.ind, lamfix=0,
          sig.ind=sig.ind, gamma.ind=gamma.ind,
          log =1)
############################################################

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

BivExtMixcode <- nimbleCode({
  for (i in 1:D)
    sds[i] ~ dunif(0, 100)
  Ustar[1:D,1:D] ~ dlkj_corr_cholesky(1.3, D)
  U[1:D,1:D] <- uppertri_mult_diag(Ustar[1:D, 1:D], sds[1:D])

  for (i in 1:D)
    mu[i] ~ dnorm(0, sd=50)
  
  for (i in 1:D)
    thres[i] ~ T(dnorm(mu.thres[i], sd=sd.thres[i]),0,max.thres[i])
  
  if (station.ind){
    for (i in 1:5)
      theta[i] ~ dunif(0,50)
    for (i in 6:7)
      theta[i] ~ dunif(-1,1)
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

X1.c <- sweep(X1, 2, c(0,mean(X1[,2])), '-')
X2.c <- sweep(X2, 2, c(0,mean(X2[,2])), '-')
BivExtMixmodel <- nimbleModel(BivExtMixcode, constants = list(N = nrow(Y), 
                                                              D = 2,
                                                              station.ind = FALSE,
                                                              D.pred = 2,
                                                              X = cbind(X1.c,X2.c),
                                                              mu.thres = c(6.6,6.9),#90 quantile
                                                              sd.thres = c(10,10), 
                                                              max.thres = apply(Y,2,quantile,0.999),
                                                              mu_beta = rep(0,2),
                                                              cov_beta = 100*diag(2),
                                                              a.ind = 1:2,
                                                              lam.ind = 3,
                                                              sig.ind = c(4,5),
                                                              gamma.ind = c(6,7),
                                                              lamfix=FALSE),check = FALSE)


BivExtMixmodel$setData(list(y = Y))  ## Set those values as data in the model
cBivExtMixmodel <- compileNimble(BivExtMixmodel)

BivExtMixconf <- configureMCMC(BivExtMixmodel,
                               enableWAIC = TRUE, time=TRUE)
# BivExtMixconf$removeSamplers(c('beta[1:4, 1]', 'beta[1:4, 2]'))
# BivExtMixconf$addSampler(target = c('beta[1:4, 1]', 'beta[1:4, 2]'), type = 'AF_slice')

BivExtMixMCMC <- buildMCMC(BivExtMixconf)
# BivExtMixMCMC$run(1)
cBivExtMixMCMC <- compileNimble(BivExtMixMCMC, project = BivExtMixmodel)

t1 <- Sys.time()
results <- runMCMC(cBivExtMixMCMC, niter = 10000,nburnin=0,thin=1,
                   summary = TRUE, WAIC = TRUE,setSeed = 1235)
t2 <- Sys.time()
print(t2-t1)

plot(results$samples[8000:10000,'thres[2]'],type='l')

var.name <- 'thres[2]'
plot(results$samples[8000:10000, var.name],type='l', main=paste('Traceplot of',var.name))
abline(h=u.x[2],col='red')

pairs(results$samples[20000:30000,c('beta.a[1, 1]','beta.a[2, 1]',
                         'beta.a[1, 2]','beta.a[2, 2]')])

pairs(results$samples[,c('beta.b[1, 1]','thres[2]')])

pairs(results$samples[,c('beta.b[1, 2]','beta.b[2, 2]')])

dir.out <- '/home/pgrad2/2448355h/My_PhD_Project/01_Output/Biv_Ext_Mix_Mod/nimble_ns_biv_ext_mix_mod_v2'


save(results, file=file.path(dir.out,'results_non_stationary_seed1235_n2500_itr80k.RData'))



nll.powunif.GPD.1(theta = c(3, 4, 5, 0.95, 0.65, 0.12, 
                            0.05), x = c(-3.58, 1.07
                            ), u = -3.592, a.ind = c(1, 2), lam.ind = 3, sig.ind = c(4,5), 
                             gamma.ind = c(6, 7), lamfix = FALSE, balthresh = FALSE, 
                            marg.scale.ind = c(1, 2), marg.shape.ind = c(1, 2))