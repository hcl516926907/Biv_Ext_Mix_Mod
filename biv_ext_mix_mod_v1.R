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

y <- X.p05

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

upp1 <- fit1$mle + 1.96*sd1
low1 <- fit1$mle - 1.96*sd1

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

upp1.1 <- fit1.1$mle + 1.96*sd1.1
low1.1 <- fit1.1$mle - 1.96*sd1.1


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
  return(prior)
}

#prior for gamma
prior.gamma <- function(theta, gamma.ind){
  d <- length(gamma.ind)
  prior <- 1
  for (i in 1:d){
    prior <- prior*dunif(theta[gamma.ind][i],0,1)
  }
  return(prior)
}

#prior for thresholds
prior.thres <- function(theta, thres.ind, X){
  d <- length(thres.ind)
  lb <- apply(X, 2, quantile, prob=0.9)
  ub <- apply(X, 2, max)
  prior <- 1
  for (i in 1:d){
    prior <- prior*dunif(theta[thres.ind][i], lb[i], ub[i])
  }
  return(prior)
}

prior.jef <- function(theta,cov.ind, d){
  mat <- matrix(theta[cov.ind], nrow=d)
  return(det(mat)^(-d/2-1))
}

#full likelihood
ll.tgt <- function(theta.all, X, thres.ind, mu.ind, cov.ind, d=2,
                   a.ind, lam.ind, lamfix=F, sig.ind, gamma.ind,
                   marg.scale.ind, marg.shape.ind, balthresh=F){
  thres <- theta.all[thres.ind]
  
  cond <- (X[,1]>thres[1]) | (X[,2]>thres[2])
  
  y.tail <- cbind(X[cond,1] - thres[1],
                 X[cond,2] - thres[2])
  
  y.bulk <- X[!cond,]
  
  #likelihood of the tail
  theta <- theta.all[-c(thres.ind, mu.ind, cov.ind)]
  llt <- -nll.powunif.GPD(theta=theta, x=y.tail, u=rep(0,d), a.ind=a.ind, lam.ind=lam.ind,
                         marg.scale.ind=marg.scale.ind, marg.shape.ind=marg.shape.ind,
                         lamfix=lamfix, balthresh=balthresh)
  #likelihood of the bulk
  Sigma <- matrix(theta.all[cov.ind], nrow=2)
  llb <- sum(dtmvnorm(y.bulk, mean=theta.all[mu.ind], sigma=Sigma, lower=c(0,0),upper=thres, log=T))
  #priors
  lp.a <- log(prior.a(theta.all, a.ind))
  lp.thres <- log(prior.thres(theta.all, thres.ind, X))
  lp.jef <- log(prior.jef(theta.all, cov.ind, d))
  
  if (lamfix==F){
    lp.lam <- log(prior.lam(theta.all, lam.ind))
    return(llb + lp.jef + lp.thres + llt + lp.a + lp.lam)
  }
  else{
    return(llb + lp.jef + lp.thres + llt + lp.a)
  }
}


p.log <- function(x){
  return(ll.tgt(theta.all=x, X=X, mu.ind = 8:9, cov.ind=10:13, d=2, thres.ind = 1:2, 
                a.ind=3, lamfix=TRUE, 
                sig.ind=4:5, gamma.ind=6:7, 
                marg.scale.ind=1:2, marg.shape.ind=1:2))
}

p.log(c(u.x.p05, runif(7),1,0,0,1))
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

