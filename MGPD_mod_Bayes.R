#################################################################################################
# Generate simulation data
#################################################################################################
# sample code
source("KRSW/RevExp_U_Functions.r")
source("KRSW/RevExp_T_Functions.r")
source("KRSW/MVGaussian_T_Functions.r")
source("KRSW/CommonFunctions.r")
source("KRSW/Gumbel_U_Functions.r")
source("KRSW/Gumbel_T_Functions.r")

set.seed(15)

d<-2
a<-c(2,4) # NB this is 1/\alpha, where \alpha is the parameterization in Kiriliouk, Rootzen, Segers and Wadsworth. 
beta<-c(0.5,0)
sig<-c(2,2)
gamma<-c(0.1,-0.1)

# Final beta parameter fixed at zero in estimation, so needs to be relative to this

# Simulate with (conditionally) GP(sig,gamma) and (conditionally) exponential margins 
X<-sim.RevExpU.MGPD(n=1000,d=d, a=a, beta=beta, sig=sig, gamma=gamma, MGPD = T,std=T)

#  GP(sig,gamma)
plot(X$X, pch=20, xlab=expression(X[1]),ylab=expression(X[2]))

# Exponential
plot(X$Z, pch=20, xlab=expression(Z[1]),ylab=expression(Z[2]))

# Fit standard (expl) form; censored likelihood censoring at (0,0)
fitS<-fit.MGPD.RevExpU(x=X$Z,u=c(0,0),std=T, hessian=TRUE)
# Parameter order: (a1,a2,beta1)
#$mle
#[1] 1.9095402 3.8381261 0.5267898

# Uncensored estimation
u1<-min(X$Z)-0.01
fitS2<-fit.MGPD.RevExpU(x=X$Z, u=c(u1,u1), std=T, hessian=TRUE)
# Parameter order: (a1,a2,beta1)
#$mle
#[1] 2.1742281 3.8906472 0.4348918

# Fit GP form; censored likelihood censoring at (0,0)
# marg.scale.ind = 1:2 states two scale parameters should be fitted (for common scale set marg.scale.ind = c(1,1))
# marg.shape.ind = 1:2 states two shape parameters should be fitted ( " " )

fitGP<-fit.MGPD.RevExpU(x=X$X, u=c(0,0), std=F, marg.scale.ind = c(1,1), marg.shape.ind = 1:2, maxit=5000, hessian=TRUE)
# Parameter order: (a1,a2,beta1,sig,gamma1,gamma2)
#$mle
#[1]  1.89320730  3.64727336  0.51113618  2.09892248  0.04513193 -0.14167273






#################################################################################################
# MCMC
#################################################################################################

y <- X$Z
d <- dim(x)[2]
a.ind <- 1:d
lam.ind <- (d+1):(2*d-1)
u <- rep(0,d)
theta <- c(2,4,2)

nll.powunif(theta,exp(y),exp(u),a.ind,lam.ind,lamfix=FALSE, balthresh=FALSE)

plot(density(rgamma(1000,shape=4,rate=1)))

prior.a <- function(theta,a.ind){
  d <- length(a.ind)
  prior <- 1
  for (i in 1:d){
    # marginally independent, follows a gamma distribution
    prior <- prior*dgamma(theta[a.ind][i], shape=4, rate=1)
  }
  return(prior)
}

prior.lam <- function(theta, lam.ind){
  d <- length(lam.ind)
  prior <- 1
  for (i in 1:d){
    prior <- prior*dgamma(theta[lam.ind][i], shape=1, rate=1)
  }
  return(prior)
}

ll.tgt <- function(theta,y,u,a.ind,lam.ind,lamfix=FALSE, balthresh=FALSE){
  ll <- -nll.powunif(theta,y,u,a.ind,lam.ind,lamfix,balthresh)
  lp1 <- log(prior.a(theta,a.ind))
  if (lamfix==FALSE){
    lp2 <- log(prior.lam(theta, lam.ind))
    return(ll + lp1 + lp2)
  }
  else{
    return(ll + lp1)
  }

}



#################################################################################################
# Metropolis-Hastings algorithm
#################################################################################################


library(tmvtnorm)
set.seed(1234)
#initial state
x.old <- c(runif(1,0,5),runif(1,0,5),runif(1,0,2))
#scale of the proposal distribution
eta <-  c(0.08,0.08,0.075)
#number of iteration
itr <-  15000
#number of burn-in samples
burnin <- 5000
samples <- x.old

set.seed(1234)
unif <- runif(itr)
t1 <- Sys.time()
for (i in 1:itr){
  # propose a new sample
  x.new <- c(rtmvnorm(1,mean=x.old,sigma=eta*diag(length(x.old)),lower=c(0,0,0)))
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
t2 <- Sys.time()
print(t2 - t1)

n_sp <- nrow(samples)
print(colMeans(samples[burnin:n_sp,]))

plot(samples[burnin:n_sp,1],type='l',main=paste(c('par:','eta:'),c('a1', eta[1])))
plot(density(samples[burnin:n_sp,1]),main=paste(c('par:','eta:'),c('a1', eta[1])))
plot(samples[burnin:n_sp,2],type='l',main=paste(c('par:','eta:'),c('a2', eta[2])))
plot(density(samples[burnin:n_sp,2]),main=paste(c('par:','eta:'),c('a2', eta[2])))
plot(samples[burnin:n_sp,3],type='l',main=paste(c('par:','eta:'),c('lam1', eta[3])))
plot(density(samples[burnin:n_sp,3]),main=paste(c('par:','eta:'),c('lam1', eta[3])))

fit1 <- samples


#################################################################################################
# Constrained scale parameter, free location parameters
#################################################################################################

library(tmvtnorm)
set.seed(1234)
#initial state
x.old <- c(runif(1,0,5),runif(1,0,5))
#scale of the proposal distribution
eta <-  c(0.08,0.08)
#number of iteration
itr <-  15000
#number of burn-in samples
burnin <- 5000
samples <- x.old

set.seed(1234)
unif <- runif(itr)
t1 <- Sys.time()
for (i in 1:itr){
  # propose a new sample
  x.new <- c(rtmvnorm(1,mean=x.old,sigma=eta*diag(length(x.old)),lower=c(0,0)))
  # move from current state to proposed state
  trans.out <- ll.tgt(x.old,exp(y),exp(u),a.ind,lam.ind,lamfix=TRUE) + 
    dtmvnorm(x.new,mean=x.old,
             sigma=eta*diag(length(x.old)),
             lower=c(0,0),
             log=TRUE)
  # move from proposed state to current state
  trans.in <- ll.tgt(x.new,exp(y),exp(u),a.ind,lam.ind,lamfix=TRUE) + 
    dtmvnorm(x.old,mean=x.new,
             sigma=eta*diag(length(x.new)),
             lower=c(0,0),
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
t2 <- Sys.time()
print(t2 - t1)

n_sp <- nrow(samples)
print(colMeans(samples[burnin:n_sp,]))

plot(samples[burnin:n_sp,1],type='l',main=paste(c('par:','eta:'),c('a1', eta[1])))
plot(density(samples[burnin:n_sp,1]),main=paste(c('par:','eta:'),c('a1', eta[1])))
plot(samples[burnin:n_sp,2],type='l',main=paste(c('par:','eta:'),c('a2', eta[2])))
plot(density(samples[burnin:n_sp,2]),main=paste(c('par:','eta:'),c('a2', eta[2])))

fit1.1 <- samples


#################################################################################################
# Constrained scale parameter, Constrained location parameters
#################################################################################################

library(tmvtnorm)
set.seed(1234)
#initial state
x.old <- c(runif(1,0,5))
#scale of the proposal distribution
eta <-  c(0.08)
#number of iteration
itr <-  15000
#number of burn-in samples
burnin <- 5000
samples <- x.old

a.ind <- 1

set.seed(1234)
unif <- runif(itr)
t1 <- Sys.time()
for (i in 1:itr){
  # propose a new sample
  x.new <- c(rtmvnorm(1,mean=x.old,sigma=eta*diag(length(x.old)),lower=c(0)))
  # move from current state to proposed state
  trans.out <- ll.tgt(x.old,exp(y),exp(u),a.ind,lam.ind,lamfix=TRUE) + 
    dtmvnorm(x.new,mean=x.old,
             sigma=eta*diag(length(x.old)),
             lower=c(0),
             log=TRUE)
  # move from proposed state to current state
  trans.in <- ll.tgt(x.new,exp(y),exp(u),a.ind,lam.ind,lamfix=TRUE) + 
    dtmvnorm(x.old,mean=x.new,
             sigma=eta*diag(length(x.new)),
             lower=c(0),
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
t2 <- Sys.time()
print(t2 - t1)

n_sp <- nrow(samples)
print(mean(samples[burnin:n_sp]))

plot(samples[burnin:n_sp,1],type='l',main=paste(c('par:','eta:'),c('a1', eta[1])))
plot(density(samples[burnin:n_sp,1]),main=paste(c('par:','eta:'),c('a1', eta[1])))

fit1.3 <- samples






#################################################################################################
# DIC
#################################################################################################

DIC <- function(post.sp,y,u,a.ind,lam.ind,lamfix=FALSE, balthresh=FALSE,...){
  theta <- colMeans(post.sp)
  pDIC <- 2*(-nll.powunif(theta,y,u,a.ind,lam.ind,lamfix,balthresh)+
               mean(apply(post.sp,MARGIN=1,FUN=nll.powunif,y=y,u=u,a.ind=a.ind,lam.ind=lam.ind,lamfix=lamfix,balthresh=balthresh)) )
  elpdDIC <- -nll.powunif(theta,y,u,a.ind,lam.ind,lamfix,balthresh) - pDIC
  return(-2*elpdDIC)
}

dic1 <- DIC(fit1[burnin:n_sp,], exp(y),exp(u),a.ind=1:2,lam.ind)
t1 <- Sys.time()
dic1.1 <- DIC(fit1.1[burnin:n_sp,], exp(y),exp(u),a.ind=1:2,lam.ind,lamfix=TRUE)
t2 <- Sys.time()
print(t2 - t1)
dic1.3 <- DIC(matrix(fit1.3[burnin:n_sp],nrow=length(burnin:n_sp)), exp(y),exp(u),a.ind=1,lam.ind,lamfix=TRUE)
