library(dmvnorm)
library(adaptMCMC)


mu <- c(3,3)
set.seed(1234)
dat <- rmvnorm(20000, mean=mu)


p.log <- function(mu){
  ll <- dmvnorm(dat,mean=mu,log=T)
  prior <- dunif(mu[1],0,5,log=T) + dunif(mu[2],0,5,log=T) 
  return (sum(ll)+prior)
}

p.log(c(2,3))




t1 <- Sys.time() 
itr <-  15000
#number of burn-in samples
burnin <- 5000
#c(0.09,0.095,0.095)
res <- MCMC.parallel(p.log, n=itr,
                     init=runif(2,0,5),n.chain=3,n.cpu=6,
                     scale=c(0.005,0.005),adapt=FALSE,packages='mvtnorm')
t2 <- Sys.time()
print(t2-t1)

plot(burnin:itr,res[[1]]$samples[burnin:itr,1], type='l', xlab='iterations', ylab='sample',
     main='Traceplot of m1')
lines(burnin:itr,res[[2]]$samples[burnin:itr,1], type='l', col='red')
lines(burnin:itr,res[[3]]$samples[burnin:itr,1], type='l', col='blue')


plot(burnin:itr,res[[1]]$samples[burnin:itr,2], type='l', xlab='iterations', ylab='sample',
     main='Traceplot of m2')
lines(burnin:itr,res[[2]]$samples[burnin:itr,2], type='l', col='red')
lines(burnin:itr,res[[3]]$samples[burnin:itr,2], type='l', col='blue')


p.log.1 <- function(x){
  # x is in the form c(theta, z, mu1, mu2)
  if (x[2]==0){
    #mu1 !=mu2
    dz <- dbinom(x[2],size=1,prob=x[1],log=T)
    dtheta <- dbeta(x[1],2,2,log=T)
    return(p.log(x[3:4]) + dz + dtheta)
    
  }else if (x[2]==1){
    dz <- dbinom(x[2],size=1,prob=x[1],log=T)
    dtheta <- dbeta(x[1],2,2,log=T)
    dy <- sum(dmvnorm(dat,mean=x[c(5,5)],log=T))
    dmu <- dunif(x[5],0,5,log=T)
    return(dy+ dmu + dz + dtheta)
    
  }else
  return (-Inf)
}


sigmoid.t <- function(x){
  x[1] <- exp(x[1])/(1+exp(x[1]))
  return(x)
}
sigmoid.inv.t <- function(x){
  x[1] <- log(x[1]/(1-x[1]))
  return(x)
}
mh_mcmc <- function(ll, n.itr, init, scale){
  samples <- matrix(NA, nrow=n.itr, ncol=length(init))
  samples[1,] <- init
  unif <- runif(n.itr)
  ll.val <- rep(NA, n.itr)
  ll.val[1] <- ll(samples[1,])
  for (i in 2:n.itr){
    #-------------------------------------------------------------------
    # update latent state. 
    x.old <- samples[i-1,]
    x.old[2] <- rbinom(1, size=1, prob=x.old[1])
    
    ll.old <- ll(x.old)
    #update other parameters
    x.new <- rmvnorm(1, mean=sigmoid.inv.t(x.old), sigma=diag(scale))
    
    #keep the latent state same
    x.new[2] <- x.old[2]
    x.new <- sigmoid.t(x.new)
    
    
    # move from current state to proposed state
    trans.out <- ll.old + 
      dmvnorm(sigmoid.inv.t(x.new), 
              mean=sigmoid.inv.t(x.old),
              sigma=diag(scale),
              log=TRUE) 
    # move from proposed state to current state
    ll.val[i] <- ll(x.new)
    trans.in <- ll.val[i] + 
      dmvnorm(sigmoid.inv.t(x.old),
              mean=sigmoid.inv.t(x.new),
              sigma=diag(scale),
              log=TRUE)
    log.accept.ratio <- trans.in - trans.out
    #acceptance rate on a log scale
    if (log(unif[i]) < min(0,log.accept.ratio)){
      samples[i,] <- x.new
    }else{
      samples[i,] <- x.old
      ll.val[i] <- ll.old
    }
  }
  return(samples)
}

# remove the transformation on the theta
mh_mcmc_1 <- function(ll, n.itr, init, scale){
  samples <- matrix(NA, nrow=n.itr, ncol=length(init))
  samples[1,] <- init
  unif <- runif(n.itr)
  ll.val <- rep(NA, n.itr)
  ll.val[1] <- ll(samples[1,])
  for (i in 2:n.itr){
    #-------------------------------------------------------------------
    # update latent state. 
    x.old <- samples[i-1,]
    x.old[2] <- rbinom(1, size=1, prob=x.old[1])
    
    ll.old <- ll(x.old)
    #update other parameters
    x.new <- rtmvnorm(1, mean=x.old, lower=c(0,rep(-Inf,length(x.old)-1)),
                      upper=c(1,rep(Inf,length(x.old)-1)),sigma=diag(scale))
    
    #keep the latent state same
    x.new[2] <- x.old[2]
    
    # move from current state to proposed state
    trans.out <- ll.old + 
      dmvnorm(x.new, 
              mean=x.old,
              sigma=diag(scale),
              log=TRUE) 
    # move from proposed state to current state
    ll.val[i] <- ll(x.new)
    trans.in <- ll.val[i] + 
      dmvnorm(x.old,
              mean=x.new,
              sigma=diag(scale),
              log=TRUE)
    log.accept.ratio <- trans.in - trans.out
    #acceptance rate on a log scale
    if (log(unif[i]) < min(0,log.accept.ratio)){
      samples[i,] <- x.new
    }else{
      samples[i,] <- x.old
      ll.val[i] <- ll.old
    }
  }
  return(samples)
}

n.burnin <- 10000
n.itr <- 30000
res <- mh_mcmc_1(ll=p.log.1, n.itr=n.itr, init=c(0.5,1,2,2,2), scale=c(0.45,1,0.001,0.001,0.001))

plot(res[n.burnin:n.itr,2],type='l')

plot(density(res[n.burnin:n.itr,1]))
plot(density(rbeta(20000,2,2)))

     
# another representation of the mixture model 
p.log.2 <- function(x){
  # x is in the form c(theta, z, mu1, mu2)
  dy1 <- dmvnorm(dat,mean=x[c(2,3)])
  dy2 <- dmvnorm(dat,mean=x[c(4,4)])
  dmu1 <- dunif(x[2],0,5,log=T) + dunif(x[3],0,5,log=T) 
  dmu2 <- dunif(x[4],0,5,log=T)
  dtheta <- dbeta(x[1],2,2,log=T)
  ll <- x[1]*dy1 + (1-x[1])*dy2
  
  return(sum(log(ll)) + dmu1 + dmu2 + dtheta)
}



t1 <- Sys.time() 
itr <-  40000
#number of burn-in samples
burnin <- 20000
#scale=0.001*c(0.05,0.05,0.05,10) when data are generated by MVN with mean c(3,4)

# res <- MCMC.parallel(p.log.2, n=itr,
#                      init=c(0.5,3,3,3),n.chain=3,n.cpu=6,
#                      scale=0.001*c(0.05,0.05,0.05,10),adapt=FALSE,packages='mvtnorm')

res <- MCMC.parallel(p.log.2, n=itr,
                     init=c(0.5,3,3,3),n.chain=3,n.cpu=6,
                     scale=0.02*c(0.2,0.2,0.2,0.1),adapt=FALSE,packages='mvtnorm')
t2 <- Sys.time()
print(t2-t1)

plot(burnin:itr,res[[1]]$samples[burnin:itr,4], type='l', xlab='iterations', ylab='sample',
     main='Traceplot of c ')

plot(density(res[[1]]$samples[burnin:itr,1]), main='density of theta')



K <- 
u <- 5
n <- 2500
x.seq <- seq(-2,12,0.01)
curve1 <- dnorm(x.seq[x.seq<u], 3, 1)
curve2 <- evd::dgpd(x.seq[x.seq>=u], loc=u, scale=1, shape=0.1)
f.unsmooth <- function(x){
  if (x<=u){
    c <- dnorm(x,3,1)
  }else{
    p <- pnorm(5,3,1)
    c <- evd::dgpd(x, loc=5, scale=1, shape=0.1)
    c <- (1-p)*c
  }

  return(c)
}

f.smooth <- function(x){
  K <- 0.1
  smooth.fac <- (1+tanh(K*(x-5)))/2
  if (x<=u){
    c <- (1-smooth.fac)*dnorm(x,3,1)
  }else{
    p <- pnorm(5,3,1)
    c <- evd::dgpd(x, loc=5, scale=1, shape=0.1)
    c <- (1-p)*c*smooth.fac
  }
  
  return(c)
}


f.smooth <- function(x){
  K <- 0.1
  smooth.fac <- (1+tanh(K*(x-5)))/2
  c <- (1-smooth.fac)*dnorm(x,3,1) + smooth.fac*
  if (x<=u){
    c <- (1-smooth.fac)*dnorm(x,3,1)
  }else{
    p <- pnorm(5,3,1)
    c <- evd::dgpd(x, loc=5, scale=1, shape=0.1)
    c <- (1-p)*c*smooth.fac
  }
  
  return(c)
}



p <- pnorm(u,3,1)
plot(x.seq, c(curve1, (1-p)*curve2), type='l')
plot(x.seq, sapply(x.seq, f.unsmooth), type='l')
plot(x.seq, sapply(x.seq, f.smooth), type='l')


library(pracma)
f.smooth <- function(x){
  K <- 0.1
  smooth.fac <- sigmoid(K*(x-5))
  if (x<=u){
    c <- (1-smooth.fac)*(-1)
  }else{
    p <- pnorm(5,3,1)
    c <- 1
    c <- c*smooth.fac
  }
  
  return(c)
}
plot(x.seq, sapply(x.seq, f.smooth), type='l')
 