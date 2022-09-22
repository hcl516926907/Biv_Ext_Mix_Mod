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
  sig <- theta.all[sig.ind]
  gamma <- theta.all[gamma.ind]
  
  #left endpoints in each margin
  eta <- -sig/gamma
  for (i in 1:d){
    if (gamma[i] <= 0) eta[i] <- -Inf
  }

  #assgin 0 to points small than left endpoints
  if (any(apply(y.tail,2,min) < eta)) {
    llt <- -Inf
  }else{
  llt <- -nll.powunif.GPD(theta=theta, x=y.tail, u=min(y.tail)-0.01, a.ind=a.ind-d, lam.ind=lam.ind-d,
                          sig.ind=sig.ind-d, gamma.ind=gamma.ind-d,
                          marg.scale.ind=marg.scale.ind, marg.shape.ind=marg.shape.ind,
                          lamfix=lamfix, balthresh=balthresh)
  }
  #likelihood of the bulk
  
  llb <- sum(dmvnorm(y.bulk, mean=theta.all[mu.ind], sigma=Sigma, log=T))
  if (min(y.bulk)<0) llb <- -Inf
  
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

p.log(c( c(5.5,6.2), rep(1,3),0,0, rep(1,2), 1,0,0,1))

p.log(c(5.55,6.213, 1.656, 0.571, 0.451, 0.253, 0.035, 3.550, 4.413, 1.5, 0.975, 0.975, 1.2))

theta.all <- c( c(5.5,6.2), rep(1,7), 1,0,0,1)
thres.ind <- 1:2
mu.ind = 8:9
cov.ind=10:13
d=2
a.ind=3
lamfix=TRUE
sig.ind=4:5
gamma.ind=6:7
marg.scale.ind=1:2
marg.shape.ind=1:2
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
init <- c( c(6,7), rep(1,3),0,0, rep(1,2),1,0,0,1)
scale <- c(rep(0.0008,2), #u1 u2 (1,2)   .549946 6.212749
           0.04, #a (3)  (1.656)
           0.002, #sigma1 (4) (0.571)
           0.002,  #sigma2 (5) (0.451)
           0.0004,  #gamma1 (6) (0.253)
           0.0004,   #gamma2 (7) (0.035)
           0.0004,    #mu1 (8)
           0.0004,    #mu2 (9)
           1500)       #cov (10-13)
n.itr <- 50000

# t1 <- Sys.time()
# res <- mh_mcmc(p.log, n.itr=n.itr, init=init,scale=scale)
# t2 <- Sys.time()
# print(t2-t1)

t1 <- Sys.time()
res1 <- mh_mcmc_1(p.log, n.itr=n.itr, init=init,scale=scale)
t2 <- Sys.time()
print(t2-t1)


plot(res1[,1],type='l', main='Traceplot of u1',ylim=c(5.2,5.8))
abline(h=u.x[1],col="red")

plot(res1[,2],type='l', main='Traceplot of u2', ylim=c(6,6.5))
abline(h=u.x[2],col="red")

plot(res1[,3],type='l', main='Traceplot of a', ylim=c(0,2.5))
abline(h=1.656,col="red")

plot(res1[,4],type='l', main='Traceplot of sigma1')
abline(h=0.571,col="red")

plot(res1[,5],type='l', main='Traceplot of sigma2', ylim=c(0,1))
abline(h=0.451,col="red")

plot(res1[,6],type='l', main='Traceplot of gamma1', ylim=c(0,0.4))
abline(h=0.253,col="red")

plot(res1[,7],type='l', main='Traceplot of gamma2', ylim=c(0,0.4))
abline(h=0.035,col="red")

plot(res1[,8],type='l', main='Traceplot of mu1', ylim=c(3.4,3.8))
abline(h=3.550,col="red")

plot(res1[,9],type='l', main='Traceplot of mu2',ylim=c(4.2,4.6))
abline(h=4.413,col="red")

plot(res1[,10],type='l', main='Traceplot of a11', ylim=c(1,2))
abline(h=1.5,col="red")

plot(res1[,11],type='l', main='Traceplot of a12')
abline(h=0.975,col="red")

plot(res1[,13],type='l', main='Traceplot of a22')
abline(h=1.2,col="red")

length(unique(res1[30000:50000,1]))/20000



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



#----------------Check the estimation of u in the univariate case-------------
u <- 5
n <- 25000
x.seq <- seq(-2,12,0.01)
curve1 <- dnorm(x.seq[x.seq<u], 3, 1)
curve2 <- evd::dgpd(x.seq[x.seq>=u], loc=u, scale=1, shape=0.1)
p <- pnorm(u,3,1)
plot(x.seq, c(curve1, (1-p)*curve2), type='l')

for (i in 1:4){
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
  plot(x.seq, exp(ll-c),type='l',xlim=c(4.5,5.5),main=paste('seed ',i,', sample size ',n, sep=''))
  abline(v=5)
  #  plot(x.seq, exp(ll.1-c.1),type='l',xlim=c(3.5,4.5),main=paste('seed',i))
}
