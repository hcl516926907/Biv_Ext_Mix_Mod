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
dir.out <- '/home/pgrad2/2448355h/My_PhD_Project/01_Output/Biv_Ext_Mix_Mod/mGPD_mod_Bayes'

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


# Use simulation to check the unidentifiable of beta
set.seed(1234)
X1<-sim.RevExpU.MGPD(n=1000,d=d, a=a, beta=beta+100, sig=sig, gamma=gamma, MGPD = T,std=T)
plot(X1$Z, pch=20, xlab=expression(Z[1]),ylab=expression(Z[2]), main='a=c(2,4), beta=c(100.5,100)')

plot(X$Z, pch=20, xlab=expression(Z[1]),ylab=expression(Z[2]), main='a=c(2,4), beta=c(0.5,0)')

biv.pp <- function(X, u){
  n <- dim(X)[1]
  # table == vector is operating down columns
  is.below <- t(t(X)<=u)
  cnt.below <- sum(rowSums(is.below)==2)
  return(cnt.below/n)
}

nx <- 50
ny <- 50
grid.mat <- matrix(0,nx,ny)
xrange <- seq(from=min(X$Z[,1], X1$Z[,1]), to=max(X$Z[,1], X1$Z[,1]), length.out=nx)
yrange <- seq(from=min(X$Z[,2], X1$Z[,2]), to=max(X$Z[,2], X1$Z[,2]), length.out=ny)
px <- c()
py <- c()
for (i in 1:nx){
  for (j in 1:ny){
    f <- biv.pp(X$Z, c(xrange[i], yrange[j]))
    f1 <- biv.pp(X1$Z, c(xrange[i], yrange[j]))
    px <- c(px, f)
    py <- c(py, f1)
    grid.mat[nx-i+1,j] <- f - f1
  }
  
}

# for better plot, simulate more points, e.g. 10k
plot(px, py, main="PP plot of the empirical cdf of the simulations",
     xlab="a=c(2,4),beta=c(0.5,0)",
     ylab="a=c(2,4),beta=c(100.5,100)")
abline(a=0, b=1)

library(corrplot)
corrplot(grid.mat, method="color",is.corr=FALSE, tl.cex=0.01)
#################################################################################################
# MCMC
#################################################################################################

y <- X$Z
d <- dim(y)[2]
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
abline(v=a[1],col='blue')
plot(samples[burnin:n_sp,2],type='l',main=paste(c('par:','eta:'),c('a2', eta[2])))
plot(density(samples[burnin:n_sp,2]),main=paste(c('par:','eta:'),c('a2', eta[2])))
abline(v=a[2],col='blue')
plot(samples[burnin:n_sp,3],type='l',main=paste(c('par:','eta:'),c('lam1', eta[3])))
plot(density(samples[burnin:n_sp,3]),main=paste(c('par:','eta:'),c('lam1', eta[3])))
abline(v=exp(beta[1]),col='blue')

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


save(fit1,fit1.1,fit1.3, file=file.path(dir.out,'my_mcmc_samples.RData'))


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

#################################################################################################
# WAIC # fail to run due to censoring integral in the nll.powunif 
#################################################################################################

WAIC.sg <- function(yi,post.sp,u,a.ind,lam.ind,lamfix,balthresh){
  yi <- matrix(yi ,nrow=1,ncol=length(yi),byrow=TRUE)
  ll <- -apply(post.sp,MARGIN=1,FUN=nll.powunif,y=yi,u=u,a.ind=a.ind,lam.ind=lam.ind,lamfix=lamfix,balthresh=balthresh)
  lppd <- log(mean(exp(ll)))
  waic.sg <- -lppd + 2*mean(ll)
return(waic.sg)
}

WAIC <- function(post.sp,y,u,a.ind,lam.ind,lamfix=FALSE, balthresh=FALSE,...){
  waic <- sum(apply(y,MARGIN=1,FUN=WAIC.sg,post.sp=post.sp,u=u,a.ind=a.ind,lam.ind=lam.ind,lamfix=lamfix,balthresh=balthresh))
  return(WAIC)
}

waic1 <- WAIC(fit1[burnin:n_sp,], exp(y),exp(u),a.ind=1:2,lam.ind)

#################################################################################################
# Parallelly run MCMC by package
#################################################################################################

set.seed(15)
d<-2
a<-c(2,4) # NB this is 1/\alpha, where \alpha is the parameterization in Kiriliouk, Rootzen, Segers and Wadsworth. 
beta<-c(0.5,0)
sig<-c(2,2)
gamma<-c(0.1,-0.1)

X<-sim.RevExpU.MGPD(n=1000,d=d, a=a, beta=beta, sig=sig, gamma=gamma, MGPD = T,std=T)
y <- X$Z
d <- dim(y)[2]
u <- rep(0,d)

p.log <- function(x){
  return(ll.tgt(x,exp(y),exp(u),a.ind=1:2,lam.ind=3,lamfix=FALSE))
}

library(adaptMCMC)
library(posterior)
t1 <- Sys.time()
itr <-  25000
#number of burn-in samples
burnin <- 5000
#c(0.09,0.095,0.095)
res <- MCMC.parallel(p.log, n=itr,
              init=c(runif(1,0,5),runif(1,0,5),runif(1,0,2)),n.chain=3,n.cpu=6,
              scale=c(0.08,0.09,0.09),adapt=FALSE)
t2 <- Sys.time()
print(t2-t1)

save(res,itr,burnin, file=file.path(dir.out,'pack_mcmc_samples_inform_prior.RData'))

plot(burnin:itr,res[[1]]$samples[burnin:itr,1], type='l', xlab='iterations', ylab='sample',
     main='Traceplot of a1')
lines(burnin:itr,res[[2]]$samples[burnin:itr,1], type='l', col='red')
lines(burnin:itr,res[[3]]$samples[burnin:itr,1], type='l', col='blue')

plot(burnin:itr,res[[1]]$samples[burnin:itr,2], type='l', xlab='iterations', ylab='sample',
     main='Traceplot of a2')
lines(burnin:itr,res[[2]]$samples[burnin:itr,2], type='l', col='red')
lines(burnin:itr,res[[3]]$samples[burnin:itr,2], type='l', col='blue')

plot(burnin:itr,res[[1]]$samples[burnin:itr,3], type='l', xlab='iterations', ylab='sample',
     main='Traceplot of lam')
lines(burnin:itr,res[[2]]$samples[burnin:itr,3], type='l', col='red')
lines(burnin:itr,res[[3]]$samples[burnin:itr,3], type='l', col='blue')


sims.a1 <- cbind(res[[1]]$samples[burnin:itr,1],res[[2]]$samples[burnin:itr,1],
                 res[[3]]$samples[burnin:itr,1])
sims.a2 <- cbind(res[[1]]$samples[burnin:itr,2],res[[2]]$samples[burnin:itr,2],
                 res[[3]]$samples[burnin:itr,2])
sims.lam <- cbind(res[[1]]$samples[burnin:itr,3],res[[2]]$samples[burnin:itr,3],
                  res[[3]]$samples[burnin:itr,3])

thin <- 10
seq.thin <- seq(1, dim(sims.a1)[1], by=thin)

rhat(sims.a1)
rhat(sims.a1[seq.thin,])
ess_basic(sims.a1)
ess_basic(sims.a1[seq.thin,])

rhat(sims.a2)
rhat(sims.a2[seq.thin,])
ess_basic(sims.a2)
ess_basic(sims.a2[seq.thin,])

rhat(sims.lam)
rhat(sims.lam[seq.thin,])
ess_basic(sims.lam)
ess_basic(sims.lam[seq.thin,])

sims.vec.a1 <- c(sims.a1[seq.thin,])
sims.vec.a2 <- c(sims.a2[seq.thin,])
sims.vec.lam <- c(sims.lam[seq.thin,])

plot(density(sims.vec.a1),main=paste(c('par:'),c('a1')))
abline(v=2,col='blue')
plot(density(sims.vec.a2),main=paste(c('par:'),c('a2')))
abline(v=4,col='blue')
plot(density(sims.vec.lam),main=paste(c('par:'),c('lam')))
abline(v=exp(0.5),col='blue')
#################################################################################################
# Posterior predictive check
#################################################################################################

library(evd)
chiplot(cbind(X$Z[,1],X$Z[,2]))
chi.org <- chiplot(X$Z,qlim=c(0.6,0.99),ask=FALSE)

ncomb<- length(sims.vec.a1)
set.seed(1234)
sp.seq <- sample(1:ncomb, size=ncomb, replace=TRUE)
chi.sim <- c()
chibar.sim <- c()

for (i in sp.seq){
  a <- c(sims.vec.a1[i],sims.vec.a2[i])
  beta <- c(log(sims.vec.lam[i]),0)
  set.seed(1234)
  X1<-sim.RevExpU.MGPD(n=1000,d=d, a=a, beta=beta, sig=sig, gamma=gamma, MGPD = T,std=T)
  dev.new()
  chi <- chiplot(X1$Z,qlim=c(0.6,0.99),ask=FALSE)
  dev.off()
  chi.sim <- rbind(chi.sim, chi$chi[,2])
  chibar.sim <- rbind(chibar.sim, chi$chibar[,2])
}

quants <- c(0.025,0.975)
chi.CI <- apply( chi.sim , 2 , quantile , probs = quants , na.rm = TRUE )
chibar.CI <- apply( chibar.sim , 2 , quantile , probs = quants , na.rm = TRUE )

plot(chi.org$quantile, chi.org$chi[,2], type='l', xlim=c(0.6,1), ylim=c(-1,1))
lines(chi.org$quantile, chi.CI[1,], type='l', lty=2, col='red')
lines(chi.org$quantile, chi.CI[2,], type='l', lty=2, col='red')

#################################################################################################
# try vague prior
#################################################################################################

set.seed(15)
d<-2
a<-c(2,4) # NB this is 1/\alpha, where \alpha is the parameterization in Kiriliouk, Rootzen, Segers and Wadsworth. 
beta<-c(0.5,0)
sig<-c(2,2)
gamma<-c(0.1,-0.1)

X<-sim.RevExpU.MGPD(n=1000,d=d, a=a, beta=beta, sig=sig, gamma=gamma, MGPD = T,std=T)
y <- X$Z
d <- dim(y)[2]
u <- rep(0,d)


prior.a <- function(theta,a.ind){
  d <- length(a.ind)
  prior <- 1
  for (i in 1:d){
    # marginally independent, follows a gamma distribution
    #prior <- prior*dgamma(theta[a.ind][i], shape=4, rate=1)
    
    prior <- prior*dunif(theta[a.ind][i],0,100)
  }
  return(prior)
}

prior.lam <- function(theta, lam.ind){
  d <- length(lam.ind)
  prior <- 1
  for (i in 1:d){
    #prior <- prior*dgamma(theta[lam.ind][i], shape=1, rate=1)
    prior <- prior*dunif(theta[lam.ind][i],0,100)
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


p.log <- function(x){
  return(ll.tgt(x,exp(y),exp(u),a.ind=1:2,lam.ind=3,lamfix=FALSE))
}

library(adaptMCMC)
library(posterior)
t1 <- Sys.time()
itr <-  25000
#number of burn-in samples
burnin <- 5000
#c(0.09,0.095,0.095)
res <- MCMC.parallel(p.log, n=itr,
                     init=c(runif(1,0,5),runif(1,0,5),runif(1,0,2)),n.chain=3,n.cpu=6,
                     scale=c(0.08,0.09,0.09),adapt=FALSE)
t2 <- Sys.time()
print(t2-t1)

plot(burnin:itr,res[[1]]$samples[burnin:itr,1], type='l', xlab='iterations', ylab='sample',
     main='Traceplot of a1')
lines(burnin:itr,res[[2]]$samples[burnin:itr,1], type='l', col='red')
lines(burnin:itr,res[[3]]$samples[burnin:itr,1], type='l', col='blue')

plot(burnin:itr,res[[1]]$samples[burnin:itr,2], type='l', xlab='iterations', ylab='sample',
     main='Traceplot of a2')
lines(burnin:itr,res[[2]]$samples[burnin:itr,2], type='l', col='red')
lines(burnin:itr,res[[3]]$samples[burnin:itr,2], type='l', col='blue')

plot(burnin:itr,res[[1]]$samples[burnin:itr,3], type='l', xlab='iterations', ylab='sample',
     main='Traceplot of lam')
lines(burnin:itr,res[[2]]$samples[burnin:itr,3], type='l', col='red')
lines(burnin:itr,res[[3]]$samples[burnin:itr,3], type='l', col='blue')


sims.a1 <- cbind(res[[1]]$samples[burnin:itr,1],res[[2]]$samples[burnin:itr,1],
                 res[[3]]$samples[burnin:itr,1])
sims.a2 <- cbind(res[[1]]$samples[burnin:itr,2],res[[2]]$samples[burnin:itr,2],
                 res[[3]]$samples[burnin:itr,2])
sims.lam <- cbind(res[[1]]$samples[burnin:itr,3],res[[2]]$samples[burnin:itr,3],
                  res[[3]]$samples[burnin:itr,3])

thin <- 10
seq.thin <- seq(1, dim(sims.a1)[1], by=thin)

rhat(sims.a1)
rhat(sims.a1[seq.thin,])
ess_basic(sims.a1)
ess_basic(sims.a1[seq.thin,])

rhat(sims.a2)
rhat(sims.a2[seq.thin,])
ess_basic(sims.a2)
ess_basic(sims.a2[seq.thin,])

rhat(sims.lam)
rhat(sims.lam[seq.thin,])
ess_basic(sims.lam)
ess_basic(sims.lam[seq.thin,])

sims.vec.a1 <- c(sims.a1[seq.thin,])
sims.vec.a2 <- c(sims.a2[seq.thin,])
sims.vec.lam <- c(sims.lam[seq.thin,])

plot(density(sims.vec.a1),main=paste(c('par:'),c('a1')))
abline(v=2,col='blue')
plot(density(sims.vec.a2),main=paste(c('par:'),c('a2')))
abline(v=4,col='blue')
plot(density(sims.vec.lam),main=paste(c('par:'),c('lam')))
abline(v=exp(0.5),col='blue')


library(evd)
chi.org <- chiplot(X$Z,qlim=c(0.6,0.985),ask=FALSE)

ncomb<- length(sims.vec.a1)
set.seed(1234)
sp.seq <- sample(1:ncomb, size=10*ncomb, replace=TRUE)
chi.sim <- c()
chibar.sim <- c()

for (i in sp.seq){
  a <- c(sims.vec.a1[i],sims.vec.a2[i])
  beta <- c(log(sims.vec.lam[i]),0)
  #set.seed(1234) #half of the posterior samples are the same. 
  #need to consider the randomness introduced by generating data.
  X1<-sim.RevExpU.MGPD(n=1000,d=d, a=a, beta=beta, sig=sig, gamma=gamma, MGPD = T,std=T)
  dev.new()
  chi <- chiplot(X1$Z,qlim=c(0.6,0.985),ask=FALSE)
  dev.off()
  chi.sim <- rbind(chi.sim, chi$chi[,2])
  chibar.sim <- rbind(chibar.sim, chi$chibar[,2])
}

quants <- c(0.025,0.975)
chi.CI <- apply( chi.sim , 2 , quantile , probs = quants , na.rm = TRUE )
chibar.CI <- apply( chibar.sim , 2 , quantile , probs = quants , na.rm = TRUE )

plot(chi.org$quantile, chi.org$chi[,2], type='l', xlim=c(0.6,1), ylim=c(-1,1))
lines(chi.org$quantile, chi.CI[1,], type='l', lty=2, col='red')
lines(chi.org$quantile, chi.CI[2,], type='l', lty=2, col='red')
