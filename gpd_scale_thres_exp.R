library(evd)
library(ismev)
library(tmvtnorm)
library(nimble)
N <- 10000
#original parameter
sig.t <- runif(N,0,10)
gamma <- runif(N,0,1)
u <- runif(N,6,7)
#scaled scale parameter
sig.s <- sig.t-gamma*u
plot(density(sig.t-gamma*u))

sig<-c(0.571,0.451)
gamma<-c(0.253,0.135)

u.x <- c(7,7.2)

0.571 - 0.253*7
0.451 - 0.135*7.2

plot(density(rnorm(N,1,3) + gamma*u))

sig.t <- sig[1]
xi <- 0.5
u <- 6

sig.s <- sig.t - xi*u

set.seed(1234)
# n <- 1500  #1500 is not large enough to have stationary estimation. Might show multimodality
n <- 15000
p <- pnorm(u, mean=2, sd=sqrt(6))
y.tail <- rgpd(n-floor(n*p), loc=u, scale=sig.t, shape=xi)

# gpd.fit(xdat=y.tail, threshold=8)
# gpd.fitrange(y.tail, umin=6,umax=8,nint=20)
par(mfrow=c(1,1))
y.bulk <- rtmvnorm(n=floor(n*p), mean=2, sigma=6, upper=u)

y <- c(y.bulk, y.tail)

plot(density(y))

ll.extr.mix <- function(y,u,sig,xi){
  y.tail <- y[y>u]
  y.bulk <- y[y<=u]
  n.tail <- length(y.tail)
  pi <- pnorm(u, mean=2, sd=sqrt(6))
  ll.bulk <- sum(dnorm(y.bulk, mean=2, sd=sqrt(6),log=T))
  ll.tail <- sum(dgpd(y.tail, loc=u, scale=sig, shape=xi, log=T))
  ll <- ll.bulk + ll.tail + n.tail*log(1-pi)
  return(ll)
}

ll.extr.mix(y,u=u,sig=sig.t,xi=xi)
ll.extr.mix(y,u=5.6,sig=sig.t,xi=xi)

y.range <- seq(-10,15,0.01)
sp.prob <- rep(NA, length(y.range))
for(i in 1:length(y.range)){
  sp.prob[i] <- exp(ll.extr.mix(y.range[i],u=u,sig=sig.t,xi=xi))
}

plot(density(sample(y.range, 20000, replace=T, prob=sp.prob/sum(sp.prob))))


nim_ll.extr.mix <- nimbleRcall(function(y=double(1), u=double(0), sig=double(0), xi=double(0)){}, 
                                   Rfun = 'll.extr.mix',
                                   returnType = double(0))

dextrmix <- nimbleFunction(
  run = function(x=double(1), u=double(0), sig=double(0), xi=double(0),
                 log = logical(0, default = 0)) {
    returnType(double(0))
    
    if (log){
      totalProb <- nim_ll.extr.mix(y=x, u=u, sig=sig, xi=xi)
    }else{
      totalProb <- exp(nim_ll.extr.mix(y=x, u=u, sig=sig, xi=xi))
    }
    return(totalProb)
  })

dextrmix(x=y, u=6, sig=0.47, xi=0.3,log=T)

rextrmix <- nimbleFunction(
  run = function(n=integer(0),u=double(0), sig=double(0), xi=double(0)) {
    returnType(double(1))
    
    totalProb <- rep(1,n)
    return(totalProb)
  })

registerDistributions(list(
  dextrmix = list(
    BUGSdist = "dextrmix(u, sig, xi)",
    types = c('value = double(1)', 'u = double(0)', 'sig = double(0)', 
              'xi = double(0)')
  )))



ExtMixcode <- nimbleCode({
  u ~ dunif(5,8)
  sig ~ dunif(0,10)
  xi ~ dunif(0,1)

  y[1:N] ~ dextrmix(u=u, sig=sig, xi=xi)
})

ExtMixmodel <- nimbleModel(ExtMixcode, constants = list(N = n),check = FALSE)


ExtMixmodel$setData(list(y = y))  ## Set those values as data in the model
cExtMixmodel <- compileNimble(ExtMixmodel)


ExtMixconf <- configureMCMC(ExtMixmodel,
                               enableWAIC = TRUE, time=TRUE)

ExtMixMCMC <- buildMCMC(ExtMixconf)
# BivExtMixMCMC$run(1)
cExtMixMCMC <- compileNimble(ExtMixMCMC, project = ExtMixmodel)

t1 <- Sys.time()
results.uni <- runMCMC(cExtMixMCMC, niter = 30000,nburnin=0,thin=1,
                   summary = TRUE, WAIC = TRUE,setSeed = 1235)
t2 <- Sys.time()
print(t2-t1)

plot(results.uni$samples[20000:30000, 'u'],type='l')
plot(results.uni$samples[20000:30000, 'sig'],type='l')
plot(results.uni$samples[20000:30000, 'xi'],type='l')



############################ reparametrize the model #####################

ll.extr.mix.1 <- function(y,u,sig,xi){
  y.tail <- y[y>u]
  y.bulk <- y[y<=u]
  n.tail <- length(y.tail)
  pi <- pnorm(u, mean=2, sd=sqrt(6))
  ll.bulk <- sum(dnorm(y.bulk, mean=2, sd=sqrt(6),log=T))
  
  sig.t <- sig + u*xi
  if (sig.t>0){
    ll.tail <- sum(dgpd(y.tail, loc=u, scale=sig.t, shape=xi, log=T))
    ll <- ll.bulk + ll.tail + n.tail*log(1-pi)
    return(ll)
  }else{
    return(-10e10)
  }

}

ll.extr.mix.1(y,u=u,sig=-1,xi=xi)

y.range <- seq(-10,15,0.01)
sp.prob <- rep(NA, length(y.range))
for(i in 1:length(y.range)){
  sp.prob[i] <- exp(ll.extr.mix.1(y.range[i],u=u,sig=sig.s,xi=xi))
}

plot(density(sample(y.range, 20000, replace=T, prob=sp.prob/sum(sp.prob))))


nim_ll.extr.mix.1 <- nimbleRcall(function(y=double(1), u=double(0), sig=double(0), xi=double(0)){}, 
                               Rfun = 'll.extr.mix.1',
                               returnType = double(0))

dextrmix1 <- nimbleFunction(
  run = function(x=double(1), u=double(0), sig=double(0), xi=double(0),
                 log = logical(0, default = 0)) {
    returnType(double(0))
    
    if (log){
      totalProb <- nim_ll.extr.mix.1(y=x, u=u, sig=sig, xi=xi)
    }else{
      totalProb <- exp(nim_ll.extr.mix.1(y=x, u=u, sig=sig, xi=xi))
    }
    return(totalProb)
  })

dextrmix1(x=y, u=6, sig=-1, xi=0.3,log=T)

rextrmix.1 <- nimbleFunction(
  run = function(n=integer(0),u=double(0), sig=double(0), xi=double(0)) {
    returnType(double(1))
    
    totalProb <- rep(1,n)
    return(totalProb)
  })

registerDistributions(list(
  dextrmix1 = list(
    BUGSdist = "dextrmix1(u, sig, xi)",
    types = c('value = double(1)', 'u = double(0)', 'sig = double(0)', 
              'xi = double(0)')
  )))



ExtMixcode1 <- nimbleCode({
  u ~ dunif(5,8)
  sig.s ~ dnorm(-2, sd=5)
  xi ~ dunif(0,1)
  sig.t <- sig.s + u*xi
  y[1:N] ~ dextrmix1(u=u, sig=sig.s, xi=xi)
})

ExtMixmodel1 <- nimbleModel(ExtMixcode1, constants = list(N = n),check = FALSE)


ExtMixmodel1$setData(list(y = y))  ## Set those values as data in the model
cExtMixmodel1 <- compileNimble(ExtMixmodel1)


ExtMixconf1 <- configureMCMC(ExtMixmodel1,
                            enableWAIC = TRUE, time=TRUE,
                            monitors=c('sig.s', 'sig.t','u', 'xi'))

ExtMixMCMC1 <- buildMCMC(ExtMixconf1)
# BivExtMixMCMC$run(1)
cExtMixMCMC1 <- compileNimble(ExtMixMCMC1, project = ExtMixmodel1)

t1 <- Sys.time()
results.uni.1 <- runMCMC(cExtMixMCMC1, niter = 30000,nburnin=0,thin=1,
                         summary = TRUE, WAIC = TRUE,setSeed = 1235)
t2 <- Sys.time()
print(t2-t1)

plot(results.uni.1$samples[20000:30000, 'u'],type='l')
plot(results.uni.1$samples[20000:30000, 'sig.t'],type='l')
plot(results.uni.1$samples[20000:30000, 'sig.s'],type='l')
plot(results.uni.1$samples[20000:30000, 'xi'],type='l')

