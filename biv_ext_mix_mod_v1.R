source("KRSW/RevExp_U_Functions.r")
source("KRSW/RevExp_T_Functions.r")
source("KRSW/MVGaussian_T_Functions.r")
source("KRSW/CommonFunctions.r")
source("KRSW/Gumbel_U_Functions.r")
source("KRSW/Gumbel_T_Functions.r")
source("KRSW/ModelDiagnosticsNewNames.r")
library(evd)
library(tmvtnorm)
library(mvtnorm)

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

prior.jef <- function(theta,cov.ind, cov.dim){
  mat <- matrix(theta[cov.ind], nrow=cov.dim)
  return(det(mat)^(-cov.dim/2-1))
}

#full likelihood
ll.tgt <- function(theta, y, u.ind, mu.ind, cov.ind, cov.dim=2,
                   a.ind, lam.ind, lamfix=F, sig.ind, gamma.ind,
                   marg.scale.ind, marg.shape.ind, balthresh=F){
  u <- theta[u.ind]
  #likelihood of the tail
  llt <- -nll.powunif.GPD(theta=theta, y=y, u=u, a.ind=a.ind, lam.ind=lam.ind,
                         marg.scale.ind=marg.scale.ind, marg.shape.ind=marg.shape.ind,
                         lamfix=lamfix, balthresh=balthresh)
  # to be done
  #likelihood of the bulk
  llb <- dtmvnorm(y, mean=theta[mu.ind], sigma=sigma, lower=c(0,0),upper=u, log=T)
  lp.a <- log(prior.a(theta,a.ind))
  lp.lam <- log(prior.lam(theta, lam.ind))
  lp.u <- log(prior.u(theta, u.ind, y))
  lp.jef <- log(prior.jef(theta,cov.ind, cov.dim))
  if (lamfix==F){
    return(llb + lp.jef + lp.u + llt + lp.a + lp.lam)
  }
  else{
    return(llb + lp.jef + lp.u + llt + lp.a)
  }
}


p.log <- function(x){
  return(ll.tgt(theta=x, y=exp(y), mu.ind = 8:9, cov.ind=10:13, u.ind = 1:2, 
                a.ind=3, lamfix=TRUE, 
                sig.ind=4:5, gamma.ind=6:7, 
                marg.scale.ind=1:2, marg.shape.ind=1:2))
}

#----------------------------MH Algorithm---------------------------
res <- MCMC.parallel(p.log, n=itr,
                     init=c(runif(1,0,5),runif(1,0,5),runif(1,0,2)),n.chain=3,n.cpu=6,
                     scale=c(0.08,0.09,0.09),adapt=FALSE)

mh_mcmc <- function(ll, itr, init=c(runif(9,0,1),c(1,0,0,1)),scale, dim.cov=2 ){
  for (i in 1:itr){
    x.old <- init
    # parameters other than matrix elements
    vec.idx <- 1:(length(init)-dim.cov^2)
    x.old.b1 <- x.old[vec.idx]
    x.old.b2 <- x.old[-vec.idx]
    
    x.new.b1 <- rmvnorm(1,mean=x.old.b1, sigma=diag(scale[vec.idx]))
    # move from current state to proposed state
    trans.out <- ll.tgt(x.old,exp(y),exp(u),a.ind,lam.ind,lamfix=FALSE) + 
      dtmvnorm(x.new,mean=x.old,
               sigma=eta*diag(length(x.old)),
               lower=c(0,0,0),
               log=TRUE)
    # move from proposed state to current state
    trans.in <- ll.tgt(x.new,exp(y),exp(u),a.ind,lam.ind,lamfix=FALSE) + 
      dtmvnorm(x.old,mean=x.new,
               sigma=eta*diag(length(x.new)),
               lower=c(0,0,0),
               log=TRUE)
    log.accept.ratio <- trans.in - trans.out
    
    #acceptance rate on a log scale
    if (log(unif[i]) < min(0,log.accept.ratio)){
      x.old <- x.new
      
    }else{
      x.old <- x.old
    }
    samples = rbind(samples,x.old)
    if (i%%500==0){
      cat(c("Done",i,'itertations \n'))
    }
  }
  
}

