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

load(file.path(dir.in,'simulation_data_pos_shape_fixed.RData'))


X <- X.p09
#X <- X.p05.10t

# total 5%, not marginal 5%
u.x <- u.x.p09
# u.x <- u.x.p05.10t

cond <- (X[,1]>u.x[1]) | (X[,2]>u.x[2])

X.gpd <- cbind(X[cond,1] - u.x[1],
               X[cond,2] - u.x[2])

y <- X.p09

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

n.sp <- 10000
plot(density(runif(n.sp,0,100)- runif(n.sp)*runif(n.sp, 3.5, 10 )))
plot(density(runif(n.sp,0,100)- runif(n.sp)*runif(n.sp, 4.3, 8 )))
#prior for sigma* that is not dependent on the threshold
prior.sig1 <- function(theta, sig.ind){
  d <- length(sig.ind)
  prior <- 1
  for (i in 1:d){
    prior <- prior*dunif(theta[sig.ind][i],-15,100)
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
  lb <- apply(X, 2, quantile, prob=0.7)
  ub <- apply(X, 2, quantile, prob=0.99)
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
  Sigma <- matrix(theta.all[cov.ind], nrow=d)
  
  cond <- (X[,1]>thres[1]) | (X[,2]>thres[2])
  
  #proportation of the bulk
  pi <- pmvnorm(lower=rep(0,d), upper=thres, mean=theta.all[mu.ind], sigma=Sigma, keepAttr = F)
  n.bulk <- sum(!cond)
  n.tail <- sum(cond)
  
  y.tail <- cbind(X[cond,1] - thres[1],
                 X[cond,2] - thres[2])
  
  y.bulk <- X[!cond,]
  
  #likelihood of the tail
  theta <- theta.all[-c(thres.ind, mu.ind, cov.ind)]
  llt <- -nll.powunif.GPD(theta=theta, x=y.tail, u=rep(0,d), a.ind=a.ind-d, lam.ind=lam.ind-d,
                          sig.ind=sig.ind-d, gamma.ind=gamma.ind-d,
                         marg.scale.ind=marg.scale.ind, marg.shape.ind=marg.shape.ind,
                         lamfix=lamfix, balthresh=balthresh)
  #likelihood of the bulk
  
  llb <- sum(dmvnorm(y.bulk, mean=theta.all[mu.ind], sigma=Sigma, log=T))
  #priors
  lp.a <- log(prior.a(theta.all, a.ind))
  lp.sig <- log(prior.sig(theta.all, sig.ind))
  lp.gamma <- log(prior.gamma(theta.all, gamma.ind))
  
  lp.thres <- log(prior.thres(theta.all, thres.ind, X))
  lp.jef <- log(prior.jef(theta.all, cov.ind, d))

  if (lamfix==F){
    lp.lam <- log(prior.lam(theta.all, lam.ind))
    return(llb + llt + n.tail*log(1-pi) + lp.jef + lp.thres + 
           lp.a + lp.lam + lp.sig + lp.gamma)
  }
  else{
    return(llb + llt + n.tail*log(1-pi) + lp.jef + lp.thres +
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

#optimize the running time by: 1. avoid calculating log likelihood of rejected points repeatedly.
# 2. combine the two blocks
mh_mcmc_1 <- function(ll, n.itr, init, scale, dim.cov=2 ){
  samples <- matrix(NA, nrow=n.itr, ncol=length(init))
  samples[1,] <- init
  unif <- runif(n.itr)
  vec.idx <- 1:(length(init)-dim.cov^2)
  ll.val <- rep(NA, n.itr)
  ll.val[1] <- ll(samples[1,])
  for (i in 2:n.itr){
    #-------------------------------------------------------------------
    # update block 1 with random walk proposal. 
    x.new.b1 <- rmvnorm(1,mean=samples[i-1,][vec.idx], sigma=diag(scale[vec.idx]))
    x.new.b2 <- as.vector(rwishart(nu=scale[-vec.idx], matrix(samples[i-1,][-vec.idx], nrow=dim.cov)/scale[-vec.idx]))
    x.new <- c(x.new.b1, x.new.b2)
    # move from current state to proposed state
    trans.out <- ll.val[i-1] + 
      dmvnorm(x.new.b1, 
              mean=samples[i-1,][vec.idx],
              sigma=diag(scale[vec.idx]),
              log=TRUE) + 
      dwishart(matrix(x.new.b2, nrow=dim.cov), 
               nu=scale[-vec.idx],
               S=matrix(samples[i-1,][-vec.idx], nrow=dim.cov)/scale[-vec.idx],
               log=TRUE)
    # move from proposed state to current state
    ll.val[i] <- ll(x.new)
    trans.in <- ll.val[i] + 
      dmvnorm(samples[i-1,][vec.idx],
              mean=x.new.b1,
              sigma=diag(scale[vec.idx]),
              log=TRUE) + 
      dwishart(matrix(samples[i-1,][-vec.idx], nrow=dim.cov),
               nu=scale[-vec.idx],
               S=matrix(x.new.b2, nrow=dim.cov)/scale[-vec.idx],
               log=TRUE)
    log.accept.ratio <- trans.in - trans.out
    #acceptance rate on a log scale
    if (log(unif[i]) < min(0,log.accept.ratio)){
      samples[i,] <- x.new
    }else{
      samples[i,] <- samples[i-1,]
      ll.val[i] <- ll.val[i-1]
    }
  }
  return(samples)
}


set.seed(1234)
init <- c( c(6,7), runif(7),1,0,0,1)
scale <- c(rep(0.002,2), #u1 u2 (1,2)   .549946 6.212749
           0.064, #a (3)  (1.656)
           0.0056, #sigma1 (4) (0.571)
           0.0008,  #sigma2 (5) (0.451)
           0.0056,  #gamma1 (6) (0.253)
           0.0004,   #gamma2 (7) (0.035)
           0.0025,    #mu1 (8)
           0.0025,    #mu2 (9)
           1000)       #cov (10-13)
n.itr <- 10000

t1 <- Sys.time()
res <- mh_mcmc(p.log, n.itr=n.itr, init=init,scale=scale)
t2 <- Sys.time()
print(t2-t1)

t1 <- Sys.time()
res1 <- mh_mcmc_1(p.log, n.itr=n.itr, init=init,scale=scale)
t2 <- Sys.time()
print(t2-t1)


plot(res1[,1],type='l', main='Traceplot of u1')
abline(h=u.x[1],col="red")

plot(res1[,2],type='l', main='Traceplot of u2')
abline(h=u.x[2],col="red")

plot(res1[,3],type='l', main='Traceplot of a', ylim=c(0,2))
abline(h=1.656,col="red")

plot(res1[,4],type='l', main='Traceplot of sigma1')
abline(h=0.571,col="red")

plot(res1[,5],type='l', main='Traceplot of sigma2', ylim=c(0,1))
abline(h=0.451,col="red")

plot(res1[,6],type='l', main='Traceplot of gamma1', ylim=c(0,1))
abline(h=0.253,col="red")

plot(res1[,7],type='l', main='Traceplot of gamma2', ylim=c(0,1))
abline(h=0.035,col="red")

plot(res1[,8],type='l', main='Traceplot of mu1')
abline(h=3.550,col="red")

plot(res1[,9],type='l', main='Traceplot of mu2')
abline(h=4.413,col="red")

plot(res1[,10],type='l', main='Traceplot of a11')
abline(h=1.5,col="red")

plot(res1[,11],type='l', main='Traceplot of a12')
abline(h=0.975,col="red")

plot(res1[,13],type='l', main='Traceplot of a22')
abline(h=1.2,col="red")


#fix the threshold
mh_mcmc_1_fix <- function(ll, n.itr, init, scale, dim.cov=2 ){
  samples <- matrix(NA, nrow=n.itr, ncol=length(init))
  samples[1,] <- init
  unif <- runif(n.itr)
  vec.idx <- 1:(length(init)-dim.cov^2)
  ll.val <- rep(NA, n.itr)
  ll.val[1] <- ll(c(u.x,samples[1,]))
  for (i in 2:n.itr){
    #-------------------------------------------------------------------
    # update block 1 with random walk proposal. 
    x.new.b1 <- rmvnorm(1,mean=samples[i-1,][vec.idx], sigma=diag(scale[vec.idx]))
    x.new.b2 <- as.vector(rwishart(nu=scale[-vec.idx], matrix(samples[i-1,][-vec.idx], nrow=dim.cov)/scale[-vec.idx]))
    x.new <- c(x.new.b1, x.new.b2)
    # move from current state to proposed state
    trans.out <- ll.val[i-1] + 
      dmvnorm(x.new.b1, 
              mean=samples[i-1,][vec.idx],
              sigma=diag(scale[vec.idx]),
              log=TRUE) + 
      dwishart(matrix(x.new.b2, nrow=dim.cov), 
               nu=scale[-vec.idx],
               S=matrix(samples[i-1,][-vec.idx], nrow=dim.cov)/scale[-vec.idx],
               log=TRUE)
    # move from proposed state to current state
    ll.val[i] <- ll(c(u.x,x.new))
    trans.in <- ll.val[i] + 
      dmvnorm(samples[i-1,][vec.idx],
              mean=x.new.b1,
              sigma=diag(scale[vec.idx]),
              log=TRUE) + 
      dwishart(matrix(samples[i-1,][-vec.idx], nrow=dim.cov),
               nu=scale[-vec.idx],
               S=matrix(x.new.b2, nrow=dim.cov)/scale[-vec.idx],
               log=TRUE)
    log.accept.ratio <- trans.in - trans.out
    #acceptance rate on a log scale
    if (log(unif[i]) < min(0,log.accept.ratio)){
      samples[i,] <- x.new
    }else{
      samples[i,] <- samples[i-1,]
      ll.val[i] <- ll.val[i-1]
    }
  }
  return(samples)
}


burnin <- 5000
n.itr <-  10000
set.seed(1234)
init <- c(runif(7),1,0,0,1)
scale <-c(#rep(0.00025,2), #u1 u2 (1,2)   .549946 6.212749
  0.08, #a (3)  (1.656)
  0.007, #sigma1 (4) (0.571)
  0.003,  #sigma2 (5) (0.451)
  0.005,  #gamma1 (6) (0.253)
  0.0003,   #gamma2 (7) (0.035)
0.0025,    #mu1 (8) (3.550)
0.0025,    #mu2 (9) (4.413)
1000)       #cov (10-13)
t1 <- Sys.time()
res1.fix <- mh_mcmc_1_fix(ll=p.log, n.itr=n.itr, init=init,scale=scale)
t2 <- Sys.time()
print(t2-t1)
plot(res1.fix[,4],type='l')
length(unique(res1.fix[burnin:n.itr,3]))/(n.itr - burnin)




# decrease the scale after several iterations
mh_mcmc_11_fix <- function(ll, n.itr, init, scale, dim.cov=2 ){
  samples <- matrix(NA, nrow=n.itr, ncol=length(init))
  samples[1,] <- init
  unif <- runif(n.itr)
  vec.idx <- 1:(length(init)-dim.cov^2)
  ll.val <- rep(NA, n.itr)
  ll.val[1] <- ll(c(u.x,samples[1,]))
  for (i in 2:n.itr){
    if (i>5000) {
      c <- 0.4
      scale.ada <- c(c*scale[-length(scale)], scale[length(scale)]/c)
      scale.ada[(length(scale)-2):(length(scale)-1)] <- scale[(length(scale)-2):(length(scale)-1)]*0.1
    }else{
      scale.ada <- scale
    }
    #-------------------------------------------------------------------
    # update block 1 with random walk proposal. 
    x.new.b1 <- rmvnorm(1,mean=samples[i-1,][vec.idx], sigma=diag(scale.ada[vec.idx]))
    x.new.b2 <- as.vector(rwishart(nu=scale.ada[-vec.idx], matrix(samples[i-1,][-vec.idx], nrow=dim.cov)/scale.ada[-vec.idx]))
    x.new <- c(x.new.b1, x.new.b2)
    # move from current state to proposed state
    trans.out <- ll.val[i-1] + 
      dmvnorm(x.new.b1, 
              mean=samples[i-1,][vec.idx],
              sigma=diag(scale.ada[vec.idx]),
              log=TRUE) + 
      dwishart(matrix(x.new.b2, nrow=dim.cov), 
               nu=scale.ada[-vec.idx],
               S=matrix(samples[i-1,][-vec.idx], nrow=dim.cov)/scale.ada[-vec.idx],
               log=TRUE)
    # move from proposed state to current state
    ll.val[i] <- ll(c(u.x,x.new))
    trans.in <- ll.val[i] + 
      dmvnorm(samples[i-1,][vec.idx],
              mean=x.new.b1,
              sigma=diag(scale.ada[vec.idx]),
              log=TRUE) + 
      dwishart(matrix(samples[i-1,][-vec.idx], nrow=dim.cov),
               nu=scale.ada[-vec.idx],
               S=matrix(x.new.b2, nrow=dim.cov)/scale.ada[-vec.idx],
               log=TRUE)
    log.accept.ratio <- trans.in - trans.out
    #acceptance rate on a log scale
    if (log(unif[i]) < min(0,log.accept.ratio)){
      samples[i,] <- x.new
    }else{
      samples[i,] <- samples[i-1,]
      ll.val[i] <- ll.val[i-1]
    }
  }
  return(samples)
}


burnin <- 5000
n.itr <-  10000
set.seed(1234)
init <- c(runif(7),1,0,0,1)
scale <-c(#rep(0.00025,2), #u1 u2 (1,2)   .549946 6.212749
  0.12, #a (3)  (1.656)
  0.010, #sigma1 (4) (0.571)
  0.005,  #sigma2 (5) (0.451)
  0.030,  #gamma1 (6) (0.253)
  0.008,   #gamma2 (7) (0.035)
  0.0025,    #mu1 (8) (3.550)
  0.0025,    #mu2 (9) (4.413)
  2000)       #cov (10-13)
t1 <- Sys.time()
res11.fix <- mh_mcmc_11_fix(ll=p.log, n.itr=n.itr, init=init,scale=scale)
t2 <- Sys.time()
print(t2-t1)
plot(res11.fix[,6],type='l')
length(unique(res11.fix[burnin:n.itr,3]))/(n.itr - burnin)

#---------------MH for mixture model with fixed threshold---------------


p.log.bulk.fix <- function(x){
  return(ll.tgt(theta.all=c(u.x, x , u.x-c(2,1.8), 1,0,0,1), X=X, mu.ind = 8:9, cov.ind=10:13, d=2, thres.ind = 1:2, 
                a.ind=3, lamfix=TRUE, 
                sig.ind=4:5, gamma.ind=6:7, 
                marg.scale.ind=1:2, marg.shape.ind=1:2))
}



t1 <- Sys.time()
system.time(res.bulk.fix <- MCMC.parallel(p.log.bulk.fix, n=13000,
                           init=runif(5),n.chain=3,n.cpu=6,
                           scale=0.5*c(0.08, 0.015, 0.0015, 0.01, 0.002),adapt=FALSE,
                           packages=c('mvtnorm','tmvtnorm'))
)
t2 <- Sys.time()
print(t2-t1)


test.para <- 0.8*c(1.656, 0.571, 0.451, 0.253, 0.035)
test.para <- rep(0.1,5)
print(p.log.bulk.fix(test.para) - p.log.gpd.2(test.para))

p.log.bulk.fix(test.para)
p.log.gpd.2(test.para)
#----------------------only GP distribution-------------------------
library(adaptMCMC)
library(posterior)
p.log.gpd <- function(x){
  return(ll.tgt(theta.all=c(u.x, x ,0.571, 0.451, 0.253, 0.035, u.x-c(2,1.8), 1,0,0,1), X=X, mu.ind = 8:9, cov.ind=10:13, d=2, thres.ind = 1:2, 
                a.ind=3, lamfix=TRUE, 
                sig.ind=4:5, gamma.ind=6:7, 
                marg.scale.ind=1:2, marg.shape.ind=1:2))
}

p.log.gpd(c(1.656))

p.log.gpd.1 <- function(x){
  theta <- c(x,0.571, 0.451, 0.253, 0.035)
  d <- 2
  thres <- u.x
  
  cond <- (X[,1]>thres[1]) | (X[,2]>thres[2])
  
  #proportation of the bulk
  
  y.tail <- cbind(X[cond,1] - thres[1],
                  X[cond,2] - thres[2])
  llt <- -nll.powunif.GPD(theta=theta, x=y.tail, u=rep(0,d), a.ind=1, lam.ind=1,
                          sig.ind=2:3, gamma.ind=4:5,
                          marg.scale.ind=1:2, marg.shape.ind=1:2,
                          lamfix=T, balthresh=F)
  lp.a <- log(prior.a(theta, 1))
  lp.sig <- log(prior.sig(theta, 2:3))
  lp.gamma <- log(prior.gamma(theta, 4:5))
  return(llt+lp.a+lp.sig+lp.gamma)
}

p.log.gpd.1(c(1.656))


t1 <- Sys.time()
itr <-  5000
#number of burn-in samples
burnin <- 2000
#c(0.09,0.095,0.095)
res <- MCMC.parallel(p.log.gpd.1, n=itr,
                     init=runif(1,0,1),n.chain=3,n.cpu=6,
                     scale=c(0.08),adapt=FALSE,
                     packages=c('mvtnorm','tmvtnorm'))
t2 <- Sys.time()
print(t2-t1)

plot(burnin:itr,res[[1]]$samples[burnin:itr,1], type='l', xlab='iterations', ylab='sample',
     main='Traceplot of a1')
lines(burnin:itr,res[[2]]$samples[burnin:itr,1], type='l', col='red')
lines(burnin:itr,res[[3]]$samples[burnin:itr,1], type='l', col='blue')
################ GP distribution with scale and shape parameters ################ 
p.log.gpd.2 <- function(x){
  theta <- x
  d <- 2
  thres <- u.x
  
  cond <- (X[,1]>thres[1]) | (X[,2]>thres[2])
  
  #proportation of the bulk
  
  y.tail <- cbind(X[cond,1] - thres[1],
                  X[cond,2] - thres[2])
  llt <- -nll.powunif.GPD(theta=theta, x=y.tail, u=rep(0,d), a.ind=1, lam.ind=1,
                          sig.ind=2:3, gamma.ind=4:5,
                          marg.scale.ind=1:2, marg.shape.ind=1:2,
                          lamfix=T, balthresh=F)
  lp.a <- log(prior.a(theta, 1))
  lp.sig <- log(prior.sig(theta, 2:3))
  lp.gamma <- log(prior.gamma(theta, 4:5))
  return(llt+lp.a+lp.sig+lp.gamma)
}

p.log.gpd.2(c(1.656, 0.571, 0.451, 0.253, 0.035))

t1 <- Sys.time()
itr <-  13000
#number of burn-in samples
burnin <- 3000
#c(0.09,0.095,0.095)
#c(0.08, 0.006, 0.001, 0.008, 0.002)
res.gpd <- MCMC.parallel(p.log.gpd.2, n=itr,
                     init=runif(5,0,1),n.chain=3,n.cpu=6,
                     scale=0.5*c(0.08, 0.015, 0.0015, 0.01, 0.002),adapt=FALSE,
                     packages=c('mvtnorm','tmvtnorm'))
t2 <- Sys.time()
print(t2-t1)


acc1 <- length(unique(res[[1]]$samples[burnin:itr,2]))/(itr-burnin)
acc2 <- length(unique(res[[2]]$samples[burnin:itr,2]))/(itr-burnin)
acc3 <- length(unique(res[[3]]$samples[burnin:itr,2]))/(itr-burnin)
print(c(acc1,acc2,acc3))


# fix the bulk bulk parameters but set threshold flexible
p.log.bulk.fix1 <- function(x){
  return(ll.tgt(theta.all=c(x , u.x-c(2,1.8), 1.5,0.975,0.975,1.2), X=X, mu.ind = 8:9, cov.ind=10:13, d=2, thres.ind = 1:2, 
                a.ind=3, lamfix=TRUE, 
                sig.ind=4:5, gamma.ind=6:7, 
                marg.scale.ind=1:2, marg.shape.ind=1:2))
}

t1 <- Sys.time()
itr <-  5000
#number of burn-in samples
burnin <- 2000
res.bulk.fix1 <- MCMC.parallel(p.log.bulk.fix1, n=itr,
                         init=c(5.5,6.2,runif(5,0,1)),n.chain=3,n.cpu=6,
                         scale=0.5*c(0.00005,0.00005, 0.08, 0.015, 0.0015, 0.01, 0.002),adapt=FALSE,
                         packages=c('mvtnorm','tmvtnorm'))
t2 <- Sys.time()
print(t2-t1)


#fix the bulk, shape and scale parameters
p.log.bulk.fix2 <- function(x){
  return(ll.tgt(theta.all=c(x ,  0.571, 0.451, 0.253, 0.035, u.x-c(2,1.8), 1.5,0.975,0.975,1.2), X=X, mu.ind = 8:9, cov.ind=10:13, d=2, thres.ind = 1:2, 
                a.ind=3, lamfix=TRUE, 
                sig.ind=4:5, gamma.ind=6:7, 
                marg.scale.ind=1:2, marg.shape.ind=1:2))
}

t1 <- Sys.time()
itr <-  50000
#number of burn-in samples
burnin <- 2000
res.bulk.fix2 <- MCMC.parallel(p.log.bulk.fix2, n=itr,
                               init=c(5.5,6.2,runif(1,0,1)),n.chain=3,n.cpu=6,
                               scale=0.5*c(0.00001,0.00001, 0.08),adapt=FALSE,
                               packages=c('mvtnorm','tmvtnorm'))
t2 <- Sys.time()
print(t2-t1)


#fix the bulk, shape, scale and one of the threshold parameters
p.log.bulk.fix3 <- function(x){
  if (x>quantile(X[,2],0.999) | (x<quantile(X[,2],0.8))){
    return(-Inf)
  }else{
  return(ll.tgt(theta.all=c(u.x[1], x , 0.571, 0.451, 0.253, 0.035, u.x-c(2,1.8), 1.5,0.975,0.975,1.2), X=X, mu.ind = 8:9, cov.ind=10:13, d=2, thres.ind = 1:2, 
                a.ind=3, lamfix=TRUE, 
                sig.ind=4:5, gamma.ind=6:7, 
                marg.scale.ind=1:2, marg.shape.ind=1:2))
  }
}


t1 <- Sys.time()
itr <-  5000
#number of burn-in samples
burnin <- 2000
res.bulk.fix3 <- MCMC.parallel(p.log.bulk.fix3, n=itr,
                               init=c(5.5,runif(1,0,1)),n.chain=3,n.cpu=6,
                               scale=0.5*c(0.01, 0.08),adapt=FALSE,
                               packages=c('mvtnorm','tmvtnorm'))
t2 <- Sys.time()
print(t2-t1)



#only threshold of margin1 is flexible
p.log.bulk.fix4 <- function(x){
  if ((x>quantile(X[,1],0.999)) | ((x<quantile(X[,1],0.8)) )){
    return(-Inf)
  }else{
    return(ll.tgt(theta.all=c(x , u.x[2], 1.656, 0.571, 0.451, 0.253, 0.035, u.x-c(2,1.8), 1.5,0.975,0.975,1.2), X=X, mu.ind = 8:9, cov.ind=10:13, d=2, thres.ind = 1:2, 
                  a.ind=3, lamfix=TRUE, 
                  sig.ind=4:5, gamma.ind=6:7, 
                  marg.scale.ind=1:2, marg.shape.ind=1:2))
  }
}


#only threshold of margin2 is flexible
p.log.bulk.fix4.1 <- function(x){
  if ((x>quantile(X[,2],0.999)) | ((x<quantile(X[,2],0.8)) )){
    return(-Inf)
  }else{
  return(ll.tgt(theta.all=c(u.x[1], x , 1.656, 0.571, 0.451, 0.253, 0.035, u.x-c(2,1.8), 1.5,0.975,0.975,1.2), X=X, mu.ind = 8:9, cov.ind=10:13, d=2, thres.ind = 1:2, 
                a.ind=3, lamfix=TRUE, 
                sig.ind=4:5, gamma.ind=6:7, 
                marg.scale.ind=1:2, marg.shape.ind=1:2))
  }
}



x.seq <- seq(4.4,8.3,0.01)
post.ll <- sapply(x.seq, p.log.bulk.fix4)
c <- max(post.ll)
plot(x.seq,exp(post.ll-c),type='l',xlim=c(5,6))
plot(x.seq,exp(post.ll-c),type='l')

#--------------test if likelihood from other data can cover true threshold-----------
d<-2
a<-c(1.656,1.656)
beta<-c(0,0)
sig<-c(0.571,0.451)
gamma<-c(0.253,0.035)
n <- 25000
mu <- c(3.549946,4.412749)
u.x <- c(5.549946,6.212749)
rho=0.65
sigma <- 1.5* matrix(c(1,rho,rho,0.8),ncol=2)
p <- pmvnorm(lower=rep(0,2), upper=u.x, mean=mu, sigma=sigma, keepAttr = F)

likelihood <- function(x){
  #thres <- c(x, u.x[2])
  thres <- x
  cond <- (X[,1]<thres[1]) & (X[,2]<thres[2])
  dat1 <- X[cond,]
  dat2 <- sweep(X[!cond,],2,thres,"-") 
  n1 <- dim(dat1)[1]
  n2 <- dim(dat2)[1]
  #p <- pmvnorm(lower=rep(0,2), upper=thres, mean=mu, sigma=sigma, keepAttr = F)
  p <- n1/(n1+n2)
  theta <- c( a[1], sig, gamma)
  ll <- -nll.powunif.GPD(theta=theta, x=dat2, u=rep(0,d), a.ind=1, lam.ind=2,
                          sig.ind=2:3, gamma.ind=4:5,
                          marg.scale.ind=1:2, marg.shape.ind=1:2,
                          lamfix=T, balthresh=F) +
  
          sum(dmvnorm(dat1, mean=mu, sigma=sigma, log=T)) +

    dunif(thres[1],2,8,log=T) +
    #dnorm(thres[1],5,0.5,log=T) +
    dunif(thres[2],2,8,log=T) +
    n2*log(1-p)
  return (ll)
}
for (i in 1:1){
  set.seed(i)
  X.tail<-sim.RevExpU.MGPD(n=n-floor(n*p),d=d, a=a, beta=beta, sig=sig, gamma=gamma, MGPD = T,std=T)
  
  # GP scale tail data combined with the bulk data
  set.seed(i)
  X.bulk <- rtmvnorm(floor(n*p), mean=mu, sigma=sigma, lower=c(0,0),upper=u.x)
  
  X <- rbind(X.bulk, sweep(X.tail$X,2,u.x,"+"))
  
  x.seq <- seq(4,6,0.1)
  
  y.seq <- seq(5,7,0.1)
  sp.grid <- expand.grid(x.seq, y.seq)

 # post.ll <- sapply(x.seq, likelihood)
  post.ll <- apply(sp.grid, 1, FUN=likelihood)
    
  c <- max(post.ll)
  
 # plot(x.seq,exp(post.ll-c),type='l',xlim=c(5,6), main=paste('seed',i))
  contour(x.seq, y.seq, matrix(post.ll-c, nrow=length(x.seq)))

}

qqplot(X.tail$X[X.tail$X[,1]>0,1],rgpd(1000, mu=0, sigma=sig[1], xi=gamma[1]))
abline(coef = c(0,1))

#-------------use reject sampling to sample from my defined likelihood------------




#--------------------------------------------------------------------------------

t1 <- Sys.time()
itr <-  10000
#number of burn-in samples
burnin <- 5000
res.bulk.fix4 <- MCMC.parallel(p.log.bulk.fix4, n=itr,
                               init=runif(1,5.5,7),n.chain=3,n.cpu=6,
                               scale=c(0.01),adapt=FALSE,
                               packages=c('mvtnorm','tmvtnorm'))
t2 <- Sys.time()
print(t2-t1)
res <- res.bulk.fix3
plot(res[[1]]$samples[,2], type='l')

plot(density(res[[1]]$samples[burnin:itr,2]))


plot(burnin:itr,res[[1]]$samples[burnin:itr,1], type='l', xlab='iterations', ylab='sample',
     main='Traceplot of a1')
lines(burnin:itr,res[[2]]$samples[burnin:itr,1], type='l', col='red')
lines(burnin:itr,res[[3]]$samples[burnin:itr,1], type='l', col='blue')

plot(burnin:itr,res[[1]]$samples[burnin:itr,2], type='l', xlab='iterations', ylab='sample',
     main='Traceplot of a1')
lines(burnin:itr,res[[2]]$samples[burnin:itr,2], type='l', col='red')
lines(burnin:itr,res[[3]]$samples[burnin:itr,2], type='l', col='blue')

plot(burnin:itr,res[[1]]$samples[burnin:itr,3], type='l', xlab='iterations', ylab='sample',
     main='Traceplot of a1')
lines(burnin:itr,res[[2]]$samples[burnin:itr,3], type='l', col='red')
lines(burnin:itr,res[[3]]$samples[burnin:itr,3], type='l', col='blue')

plot(burnin:itr,res[[1]]$samples[burnin:itr,4], type='l', xlab='iterations', ylab='sample',
     main='Traceplot of a1')
lines(burnin:itr,res[[2]]$samples[burnin:itr,4], type='l', col='red')
lines(burnin:itr,res[[3]]$samples[burnin:itr,4], type='l', col='blue')

plot(burnin:itr,res[[1]]$samples[burnin:itr,5], type='l', xlab='iterations', ylab='sample',
     main='Traceplot of a1')
lines(burnin:itr,res[[2]]$samples[burnin:itr,5], type='l', col='red')
lines(burnin:itr,res[[3]]$samples[burnin:itr,5], type='l', col='blue')




#only threshold of margin1 is flexible; change the scale parameter in GP distribution
# threshold independent.
sigma.ast <- c(0.571, 0.451) - c(0.253, 0.035)*u.x
p.log.bulk.fix4.2 <- function(x){
  if ((x>quantile(X[,1],0.999)) | ((x<quantile(X[,1],0.8)) )){
    return(-Inf)
  }else{
    return(ll.tgt(theta.all=c(x , u.x[2], 1.656, -0.833+0.253*x, 0.451, 0.253, 0.035, u.x-c(2,1.8), 1.5,0.975,0.975,1.2), X=X, mu.ind = 8:9, cov.ind=10:13, d=2, thres.ind = 1:2, 
                  a.ind=3, lamfix=TRUE, 
                  sig.ind=4:5, gamma.ind=6:7, 
                  marg.scale.ind=1:2, marg.shape.ind=1:2))
  }
}

x.seq <- seq(4.4,8.3,0.01)
post.ll.2 <- sapply(x.seq, p.log.bulk.fix4.2)
c <- max(post.ll.2)
plot(x.seq,exp(post.ll.2-c),type='l',xlim=c(5,6))
plot(x.seq,exp(post.ll.2-c),type='l')

#----------------Check the estimation of u in the univariate case-------------
u <- 5
n <- 50000
x.seq <- seq(-2,12,0.01)
curve1 <- dnorm(x.seq[x.seq<u], 3, 1)
curve2 <- evd::dgpd(x.seq[x.seq>=u], loc=u, scale=1, shape=0.1)
p <- pnorm(u,3,1)
plot(x.seq, c(curve1, (1-p)*curve2), type='l')

for (i in 1:10){
  set.seed(i)
  sim.dat1 <- rtmvnorm(floor(n*p), mean=3, sigma=matrix(1,nrow=1), upper=u)
  set.seed(i)
  sim.dat2 <- evd::rgpd(n=n-floor(n*p), loc=0, scale=1, shape=0.1)
  sim.dat <- c(sim.dat1, sim.dat2+u)
  #plot(density(sim.dat))
  
  loglikeli <- function(x){
    dat1 <- sim.dat[sim.dat<x]
    dat2 <- sim.dat[sim.dat>=x]
    n1 <- length(dat1)
    n2 <- length(dat2)
    pi <- pnorm(x, mean=3, sd=1)
    ll <- sum(dnorm(dat1, mean=3, sd=1, log=T)) +
      sum(evd::dgpd(dat2-u, loc=0, scale=1,  shape=0.1, log=T)) +
      dunif(x,0,10,log=T) +
      # dnorm(x,6,0.1,log=T)+
      n2*log(1-pi)
    return(ll)
  }
  
  
  # loglikeli.1 <- function(x){
  #   dat1 <- sim.dat[sim.dat<x]
  #   dat2 <- sim.dat[sim.dat>=x]
  #   n1 <- length(dat1)
  #   n2 <- length(dat2)
  #   pi <- pnorm(x, mean=3, sd=1)
  #   ll <- sum(dnorm(dat1, mean=3, sd=1, log=T)) +
  #     sum(evd::dgpd(dat2, loc=x, scale=0.5+x*0.1,  shape=0.1, log=T)) +
  #     dunif(x,0,10,log=T) +
  #     # dnorm(x,6,0.1,log=T)+
  #     n2*log(1-pi)
  #   return(ll)
  # }
  x.seq <- seq(3,8,0.01)
  system.time(ll <-sapply(x.seq, loglikeli))
#  system.time(ll.1 <- sapply(x.seq, loglikeli.1))
  
  c <- max(ll)
#  c.1 <- max(ll.1)
  plot(x.seq, exp(ll-c),type='l',xlim=c(4.5,5.5),main=paste('seed',i))
#  plot(x.seq, exp(ll.1-c.1),type='l',xlim=c(3.5,4.5),main=paste('seed',i))
}
