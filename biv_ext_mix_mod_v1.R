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
library(LaplacesDemon)

dir.in <- '/home/pgrad2/2448355h/My_PhD_Project/01_Output/Biv_Ext_Mix_Mod/biv_ext_mix_mod_simdat'

load(file.path(dir.in,'simulation_data_pos_shape.RData'))

X <- X.p05

# total 5%, not marginal 5%
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
  lb <- apply(X, 2, quantile, prob=0.5)
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
  
  #proportation of the bulk
  pi <- sum(!cond)/length(cond)
  n.bulk <- sum(!cond)
  n.tail <- sum(cond)
  
  y.tail <- cbind(X[cond,1] - thres[1],
                 X[cond,2] - thres[2])
  
  y.bulk <- X[!cond,]
  
  #likelihood of the tail
  theta <- theta.all[-c(thres.ind, mu.ind, cov.ind)]
  llt <- -nll.powunif.GPD(theta=theta, x=y.tail, u=rep(0,d), a.ind=a.ind, lam.ind=lam.ind,
                          sig.ind=sig.ind-d, gamma.ind=gamma.ind-d,
                         marg.scale.ind=marg.scale.ind, marg.shape.ind=marg.shape.ind,
                         lamfix=lamfix, balthresh=balthresh)
  #likelihood of the bulk
  Sigma <- matrix(theta.all[cov.ind], nrow=2)
  llb <- sum(dtmvnorm(y.bulk, mean=theta.all[mu.ind], sigma=Sigma, lower=c(0,0),upper=thres, log=T))
  #priors
  lp.a <- log(prior.a(theta.all, a.ind))
  lp.sig <- log(prior.sig(theta.all, sig.ind))
  lp.gamma <- log(prior.gamma(theta.all, gamma.ind))
  
  lp.thres <- log(prior.thres(theta.all, thres.ind, X))
  lp.jef <- log(prior.jef(theta.all, cov.ind, d))
  
  if (lamfix==F){
    lp.lam <- log(prior.lam(theta.all, lam.ind))
    return(llb + llt + n.bulk*log(pi) + n.tail*log(1-pi) + lp.jef + lp.thres + 
           lp.a + lp.lam + lp.sig + lp.gamma)
  }
  else{
    return(llb + llt + n.bulk*log(pi) + n.tail*log(1-pi) + lp.jef + lp.thres +
           lp.a + lp.sig + lp.gamma)
  }
}


p.log <- function(x){
  return(ll.tgt(theta.all=x, X=X, mu.ind = 8:9, cov.ind=10:13, d=2, thres.ind = 1:2, 
                a.ind=3, lamfix=TRUE, 
                sig.ind=4:5, gamma.ind=6:7, 
                marg.scale.ind=1:2, marg.shape.ind=1:2))
}

p.log(c( c(5.4,6), rep(1,7), 1,0,0,1))

p.log(c(5.55,6.213, 1.656, 0.571, 0.451, 0.253, 0.035, 3.550, 4.413, 1.5, 0.975, 0.975, 1.2))

#----------------------------MH Algorithm---------------------------
res <- MCMC.parallel(p.log, n=itr,
                     init=c(runif(1,0,5),runif(1,0,5),runif(1,0,2)),n.chain=3,n.cpu=6,
                     scale=c(0.08,0.09,0.09),adapt=FALSE)

#paramter order: u1,u2,a,sigma1,sigma2,gamma1,gamma2,mu1,m2,COV
mh_mcmc <- function(ll, n.itr, init, scale, dim.cov=2 ){
  samples <- init
  x.old <- init
  unif.b1 <- runif(n.itr)
  unif.b2 <- runif(n.itr)
  for (i in 1:n.itr){
    # parameters other than matrix elements
    vec.idx <- 1:(length(init)-dim.cov^2)
    x.old.b1 <- x.old[vec.idx]
    x.old.b2 <- x.old[-vec.idx]
    #-------------------------------------------------------------------
    # update block 1 with random walk proposal. 
    x.new.b1 <- rmvnorm(1,mean=x.old.b1, sigma=diag(scale[vec.idx]))
    x.new <- x.old
    x.new[vec.idx] <- x.new.b1
    # move from current state to proposed state
    trans.out.b1 <- ll(x.old) + 
      dmvnorm(x.new.b1, mean=x.old.b1,
               sigma=diag(scale[vec.idx]),
               log=TRUE)
    # move from proposed state to current state
    trans.in.b1 <- ll(x.new) + 
      dmvnorm(x.old.b1,mean=x.new.b1,
               sigma=diag(scale[vec.idx]),
               log=TRUE)
    log.accept.ratio.b1 <- trans.in.b1 - trans.out.b1
    #acceptance rate on a log scale
    if (log(unif.b1[i]) < min(0,log.accept.ratio.b1)){
      x.old <- x.new
    }else{
      x.old <- x.old
    }
    #------------------------------------------------------------------
    # update block 2 with whishart proposal. 
    x.new.b2 <- as.vector(rwishart(nu=scale[-vec.idx], matrix(x.old.b2, nrow=dim.cov)/scale[-vec.idx]))
    x.new <- x.old
    x.new[-vec.idx] <- x.new.b2
    # move from current state to proposed state
    trans.out.b2 <- ll(x.old) + 
      dwishart(matrix(x.new.b2, nrow=dim.cov), 
               nu=scale[-vec.idx],
               S=matrix(x.old.b2, nrow=dim.cov)/scale[-vec.idx],
               log=TRUE)
    # move from proposed state to current state
    trans.in.b2 <- ll(x.new) + 
      dwishart(matrix(x.old.b2, nrow=dim.cov),
              nu=scale[-vec.idx],
              S=matrix(x.new.b2, nrow=dim.cov)/scale[-vec.idx],
              log=TRUE)
    log.accept.ratio.b2 <- trans.in.b2 - trans.out.b2
    #acceptance rate on a log scale
    if (log(unif.b2[i]) < min(0,log.accept.ratio.b2)){
      x.old <- x.new
    }else{
      x.old <- x.old
    }
    samples = rbind(samples,x.old)
    if (i%%500==0){
      cat(c("Done",i,'itertations \n'))
    }
  }
  return(samples)
  
}

set.seed(1234)
init <- c( c(6,7), runif(7),1,0,0,1)
scale <- c(rep(0.00025,2), #u1 u2 (1,2)   .549946 6.212749
           0.001, #a (3)  (1.656)
           0.001, #sigma1 (4) (0.571)
           0.001,  #sigma2 (5) (0.451)
           0.01,  #gamma1 (6) (0.253)
           0.001,   #gamma2 (7) (0.035)
           0.005,    #mu1 (8)
           0.005,    #mu2 (9)
           80)       #cov (10-13)
n.itr <- 10000
res <- mh_mcmc(p.log, n.itr=n.itr, init=init,scale=scale)

plot(res[,1],type='l', main='Traceplot of u1')
abline(h=u.x[1],col="red")

plot(res[,2],type='l', main='Traceplot of u2')
abline(h=u.x[2],col="red")

plot(res[,3],type='l', main='Traceplot of a', ylim=c(0,2))
abline(h=1.656,col="red")

plot(res[,4],type='l', main='Traceplot of sigma1')
abline(h=0.571,col="red")

plot(res[,5],type='l', main='Traceplot of sigma2', ylim=c(0,1))
abline(h=0.451,col="red")

plot(res[,6],type='l', main='Traceplot of gamma1', ylim=c(0,1))
abline(h=0.253,col="red")

plot(res[,7],type='l', main='Traceplot of gamma2', ylim=c(0,1))
abline(h=0.035,col="red")

plot(res[,8],type='l', main='Traceplot of mu1')
abline(h=3.550,col="red")

plot(res[,9],type='l', main='Traceplot of mu2')
abline(h=4.413,col="red")

plot(res[,10],type='l', main='Traceplot of a11')
abline(h=1.5,col="red")

plot(res[,11],type='l', main='Traceplot of a12')
abline(h=0.975,col="red")

plot(res[,13],type='l', main='Traceplot of a22')
abline(h=1.2,col="red")



#plot(samples[2000:10000,3],type='l')

#issues in the prior of u 
#issues in p.log
#issues in data 
#issues in the prior of gamma2