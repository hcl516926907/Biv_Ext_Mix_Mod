source("KRSW/RevExp_U_Functions.r")
source("KRSW/RevExp_T_Functions.r")
source("KRSW/MVGaussian_T_Functions.r")
source("KRSW/CommonFunctions.r")
source("KRSW/Gumbel_U_Functions.r")
source("KRSW/Gumbel_T_Functions.r")
source("KRSW/ModelDiagnosticsNewNames.r")
library(evd)
library(tmvtnorm)

dir.in <- '/home/pgrad2/2448355h/My_PhD_Project/01_Output/Biv_Ext_Mix_Mod/biv_ext_mix_mod_simdat'

load(file.path(dir.in,'simulation_data_pos_shape.RData'))

X <- X.p05

u.x <- u.x.p05

cond <- (X[,1]>u.x[1]) | (X[,2]>u.x[2])

X.gpd <- cbind(X[cond,1] - u.x[1],
               X[cond,2] - u.x[2])


#------------------------------MLE of bivariate GP model----------------------
#with censoring
fit1<-fit.MGPD.RevExpU(x=X.gpd, u=rep(0,2), std=F,
                         dep.scale.fix=T, dep.loc.fix=T, marg.shape.ind=1:2, marg.scale.ind=1:2,
                         maxit=5000)
for (i in 1:5){
  fit1<-fit.MGPD.RevExpU(x=X.gpd, u=rep(0,2), std=F,
                         dep.scale.fix=T, dep.loc.fix=T, marg.shape.ind=1:2, marg.scale.ind=1:2,
                         marg.scale.start=fit1$mle[2:3], marg.shape.start=fit1$mle[4:5],
                         dep.start=fit1$mle[1],maxit=5000,hessian=TRUE)
}


sd1 <- sqrt(diag(solve(fit1$hess)))

upp1 <- fit1$mle + 1.96*sd
low1 <- fit1$mle - 1.96*sd

#without censoring
fit1.1<-fit.MGPD.RevExpU(x=X.gpd, u=apply(X.gpd,2,min)-0.1, std=F,
                         dep.scale.fix=T, dep.loc.fix=T, marg.shape.ind=1:2, marg.scale.ind=1:2,
                         maxit=5000)

for (i in 1:10){
  fit1.1<-fit.MGPD.RevExpU(x=X.gpd, u=apply(X.gpd,2,min)-0.1, std=F,
                           dep.scale.fix=T, dep.loc.fix=T, marg.shape.ind=1:2, marg.scale.ind=1:2,
                           marg.scale.start=fit1.1$mle[2:3], marg.shape.start=fit1.1$mle[4:5],
                           dep.start=fit1.1$mle[1],maxit=15000,hessian=TRUE)
}


sd1.1 <- sqrt(diag(solve(fit1.1$hess)))

upp1.1 <- fit1.1$mle + 1.96*sd
low1.1 <- fit1.1$mle - 1.96*sd


#-----------------------Bayesian inference of the mixture model----------------

#prior for a
prior.a <- function(theta,a.ind){
  d <- length(a.ind)
  prior <- 1
  for (i in 1:d){
    prior <- prior*dunif(theta[a.ind][i],0,100)
  }
  return(prior)
}

#prior for lambda
prior.lam <- function(theta, lam.ind){
  d <- length(lam.ind)
  prior <- 1
  for (i in 1:d){
    #prior <- prior*dgamma(theta[lam.ind][i], shape=1, rate=1)
    prior <- prior*dunif(theta[lam.ind][i],0,100)
  }
  return(prior)
}

#prior for sigma
prior.sig <- function(theta, sig.ind){
  d <- length(sig.ind)
  prior <- 1
  for (i in 1:d){
    prior <- prior*dunif(theta[sig.ind][i],0,100)
  }
}

#prior for gamma
prior.gamma <- function(theta, gamma.ind){
  d <- length(gamma.ind)
  prior <- 1
  for (i in 1:d){
    prior <- prior*dunif(theta[gamma.ind][i],0,1)
  }
}

#prior for thresholds
prior.u <- function(theta, u.ind, y){
  d <- length(u.ind)
  lb <- apply(y, 2, quantile, prob=0.9)
  ub <- apply(y, 2, max)
  prior <- 1
  for (i in 1:d){
    prior <- prior*dunif(theta[u.ind][i], lb[i], ub[i])
  }
}

#full likelihood
ll.tgt <- function(theta, y, u.ind, a.ind,lam.ind,lamfix=FALSE, sig.ind, gamma.ind,
                   marg.scale.ind, marg.shape.ind, balthresh=FALSE){
  u <- theta[u.ind]
  llt <- -nll.powunif.GPD(theta=theta, y=y, u=u, a.ind=a.ind, lam.ind=lam.ind,
                         marg.scale.ind=marg.scale.ind, marg.shape.ind=marg.shape.ind,
                         lamfix=lamfix, balthresh=balthresh)
  # to be done
  llb <- dtmvnorm(2375, mean=u.x-c(2,1.8), sigma=1.5*sigma, lower=c(0,0),upper=u)
  lp1 <- log(prior.a(theta,a.ind))
  if (lamfix==FALSE){
    lp2 <- log(prior.lam(theta, lam.ind))
    return(ll + lp1 + lp2)
  }
  else{
    return(ll + lp1)
  }
  
}


p.log <- function(x){
  return(ll.tgt(theta=x, y=exp(y), u.ind = 1:2, a.ind=3, sig.ind=4:5,
                gamma.ind=6:7, marg.scale.ind=1:2, marg.shape.ind=1:2, lamfix=TRUE))
}

