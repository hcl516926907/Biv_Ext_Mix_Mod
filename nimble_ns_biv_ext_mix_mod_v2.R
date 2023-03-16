source("KRSW/RevExp_U_Functions.r")
source("KRSW/CommonFunctions.r")
library(nimble, warn.conflicts = F)
library(mvtnorm)
library(parallel)


dir.in <- '/home/pgrad2/2448355h/My_PhD_Project/01_Output/Biv_Ext_Mix_Mod/biv_ext_mix_mod_simdat'

load(file.path(dir.in,'simulation_data_non_stationary_extr_dep.RData'))


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
  for(j in 1:d)
  {
    rej[j]<-gamma[j]<0 && any(x[,j]>-sig[j]/gamma[j])
  }
  
  if(any(lam<0.01)||any(a<=0.01)||any(sig<=0.001)||any(rej)){return(10e10)}
  
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
  if (is.nan(nll.uc)|is.nan(nll.pc)){
    return(10e10)
  }
  nll<-nll.uc+nll.pc
  
  return(nll)
}

# x <- Y.tail
# u <- min(Y.tail)-0.01
# 
# a.par <- a[1,]
# lam <- exp(b[1,1])
# sig<-c(0.571,0.451)
# gamma<-c(0.253,0.135)
# theta <- c(a.par, lam, sig, gamma)
# a.ind <- c(1,2)
# lam.ind <- c(3)
# sig.ind <- c(4,5)
# gamma.ind <- c(6,7)
# marg.scale.ind <- c(1,2)
# marg.shape.ind <- c(1,2)
# lamfix=FALSE
# balthresh=FALSE

# nll.powunif.GPD.1(theta,x,u,a.ind,lam.ind,sig.ind,gamma.ind, lamfix=FALSE, balthresh=FALSE, marg.scale.ind,marg.shape.ind)

nim_nll_powunif_GPD <- nimbleRcall(function(theta=double(1), x=double(1), u=double(0), a.ind=double(1),
                                            lam.ind=double(0), sig.ind=double(1), gamma.ind=double(1), 
                                            lamfix=logical(0, default = 0), balthresh=logical(0, default = 0), 
                                            marg.scale.ind=double(1), marg.shape.ind=double(1)){}, 
                                   Rfun = 'nll.powunif.GPD.1',
                                   returnType = double(0))

x<- c(14.17519,3.667644)
theta <- c(0.01137369, 0.01461362, 1, 0.571, 0.451, 0.253 ,0.035)
u <- min(x)-0.01
a.ind <- c(1,2)
lam.ind <- c(3)
sig.ind <- c(4,5)
gamma.ind <- c(6,7)
marg.scale.ind <- c(1,2)
marg.shape.ind <- c(1,2)
lamfix=FALSE
balthresh=FALSE
nim_nll_powunif_GPD(x=x, theta=theta, u=u, a.ind=a.ind,
                    lam.ind=lam.ind, sig.ind=sig.ind, gamma.ind=gamma.ind, 
                    lamfix=lamfix, balthresh=FALSE, 
                    marg.scale.ind=1:2, marg.shape.ind=1:2)


dbiextmix <- nimbleFunction(
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
            alpha.inv[j] <- exp(inprod(X.tail[i,(1+(j-1)*K):(j*K)], beta.a[,D]))
            if (j<D){
              lam[j] <- exp(inprod(X.tail[i,(1+(j-1)*K):(j*K)] ,beta.b[,D]))
            }else{
              lam[j] <- 1
            }
            
          }
          theta <- c(alpha.inv,lam[1:(D-1)],para.mg)
          llt.sg <- -nim_nll_powunif_GPD(x=y.tail[i,], theta=theta, u=min(y.tail)-0.01, a.ind=a.ind,
                                      lam.ind=lam.ind, sig.ind=sig.ind, gamma.ind=gamma.ind, 
                                      lamfix=lamfix, balthresh=FALSE, 
                                      marg.scale.ind=1:2, marg.shape.ind=1:2)
          if (is.nan(llt.sg)){
            print(c(y.tail[i,],theta))
          }
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

test1 <- function(){
  -nim_nll_powunif_GPD(x=y.tail[i,], theta=theta, u=min(y.tail)-0.01, a.ind=a.ind,
                       lam.ind=lam.ind, sig.ind=sig.ind, gamma.ind=gamma.ind, 
                       lamfix=lamfix, balthresh=FALSE, 
                       marg.scale.ind=1:2, marg.shape.ind=1:2)
}

X <- cbind(X1,X2)
# X <- matrix(1,nrow=1000,ncol=4)
# para.mg <- c(0.571, 0.451, 0.253, 0.035)
para.mg <- c(9.3556,2.02373,0.820941,0.456565)
# beta.a <- cbind(c(0.5, 0.2),c(0.5, -0.4))
beta.a <- cbind(c(15.1993,-12.5238),c(-0.578615,0.323246))
# beta.b <- cbind(c(0.25,-0.25),c(0,0))
beta.b <- cbind(c(6.05148,-1.62846),c(1,1))
a.ind <- c(1,2)
lam.ind <- c(3)
sig.ind <- c(4,5)
gamma.ind <- c(6,7)
marg.scale.ind <- c(1,2)
marg.shape.ind <- c(1,2)
# thres <- u.x
thres <- c(9.44768,7.22659)
mu <- c(5,5.41)
rho=0.5
sigma <- 1.5* matrix(c(1,rho,rho,0.8),ncol=2)
cholesky <- chol(sigma)
x <- Y
D <- 2
log <- TRUE
lamfix=FALSE

t1 <- Sys.time()
# dbiextmix(x=Y, theta=theta, thres = c(7.049958,6.851115),mu=mu,
#           cholesky=cholesky,
#           a.ind=a.ind, lam.ind=lam.ind, lamfix=0,
#           sig.ind=sig.ind, gamma.ind=gamma.ind,
#           log =1)

beta.a <- rmvnorm(2, mean=rep(0,2), sigma=10*diag(2))
# beta.a.issue <- matrix(c(1.902316,-4.395543,2.021196,-3.572355),nrow=2)
# beta.a.issue <- matrix(c(3.39632,-3.91076,-12.0432,-12.3432),nrow=2)
beta.a.issue <- matrix(c(3.39632,-12.0432,-3.91076,-12.3432),,nrow=2)
dbiextmix(x=Y, para.mg=para.mg, beta.a=beta.a.issue, beta.b=beta.b, X=X,
          thres=thres, mu=mu,
          cholesky=cholesky,
          a.ind=a.ind, lam.ind=lam.ind, lamfix=0,
          sig.ind=sig.ind, gamma.ind=gamma.ind,
          log =1)

# -7392.573
t2 <- Sys.time()
print(t2-t1)



rbiextmix <- nimbleFunction(
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
  dbiextmix = list(
    BUGSdist = "dbiextmix(para.mg, beta.a, beta.b, thres, X, mu, cholesky, D, 
                a.ind, lam.ind, lamfix, sig.ind, gamma.ind)",
    types = c('value = double(2)', 'para.mg = double(1)', 'beta.a = double(2)', 
              'beta.b = double(2)', 'thres = double(1)', 'X = double(2)',
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
  
  for (i in 1:2)
    para.mg[i] ~ dunif(0,50)
  for (i in 3:4)
    para.mg[i] ~ dunif(0,1)
  for (i in 1:2)
    mu[i] ~ dnorm(0, sd=50)
  
  for (i in 1:D.pred)
    mu_beta[i] <- 0
  for (i in 1:2)
    thres[i] ~ dunif(lb[i], ub[i])
  for (i in 1:D){
    beta.a[1:D.pred,i] ~ dmnorm(mu_beta[1:D.pred], cov=cov_beta[1:D.pred,1:D.pred])
  }
  beta.b[1:D.pred,1] ~ dmnorm(mu_beta[1:D.pred], cov=cov_beta[1:D.pred,1:D.pred])
  beta.b[1:D.pred,2] <-  rep(1, D.pred)
  
  
  
  y[1:N,1:D] ~ dbiextmix(para.mg=para.mg[1:4], beta.a=beta.a[1:D.pred,1:D],
                         beta.b=beta.b[1:D.pred,1:D], thres=thres[1:2],
                         X=X[1:N,1:(D*D.pred)], mu=mu[1:D], 
                         cholesky=U[1:D,1:D], D=D,
                         a.ind=a.ind[1:D], lam.ind=lam.ind, lamfix=lamfix, 
                         sig.ind=sig.ind[1:D], gamma.ind=gamma.ind[1:D])
})

X1.c <- sweep(X1, 2, c(0,mean(X1[,2])), '-')
X2.c <- sweep(X2, 2, c(0,mean(X2[,2])), '-')
BivExtMixmodel <- nimbleModel(BivExtMixcode, constants = list(N = 2500, 
                                                              D = 2,
                                                              D.pred = 2,
                                                              cov_beta = 100*diag(2),
                                                              X = cbind(X1.c,X2.c),
                                                              lb = c(5.85, 6.15),
                                                              ub = c(21.17, 11.30),
                                                              a.ind = 1:2,
                                                              lam.ind = 3,
                                                              sig.ind = c(4,5),
                                                              gamma.ind = c(6,7),
                                                              lamfix=1),check = FALSE)


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
results <- runMCMC(cBivExtMixMCMC, niter = 250,nburnin=0,thin=1,
                   summary = TRUE, WAIC = TRUE,setSeed = 1234)
t2 <- Sys.time()
print(t2-t1)

plot(results$samples[,'beta.b[1, 1]'],type='l')
plot(results$samples[,'thres[2]'],type='l')

# pairs(results$samples[,c('beta[1, 1]','beta[2, 1]','beta[3, 1]','beta[4, 1]',
#                          'beta[1, 2]','beta[2, 2]','beta[3, 2]','beta[4, 2]')])
pairs(results$samples[,c('beta[1, 1]','beta[2, 1]',
                         'beta[1, 2]','beta[2, 2]')])
plot(results$samples[,c('beta[2, 1]','beta[3, 1]')])

dir.out <- '/home/pgrad2/2448355h/My_PhD_Project/01_Output/Biv_Ext_Mix_Mod/nimble_ns_biv_ext_mix_mod'

save(results, thres.samples, file=file.path(dir.out,'sp_rw_l5dot5_l5dot7_u8_u8.RData'))

