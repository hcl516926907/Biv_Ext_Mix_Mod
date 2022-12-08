library(nimble, warn.conflicts = F)

#-----------------------------------Mixture Model Example---------------------------------

ZIPcode <- nimbleCode({
  p ~ dunif(0,1)
  lambda ~ dunif(0,10)
  for (i in 1:N)
    y[i] ~ dZIP(lambda, zeroProb = p) ## Note NIMBLE allows R-like named-parameter syntax
})



dZIP <- nimbleFunction(
  run = function(x = integer(), lambda = double(), 
                 zeroProb = double(), log = logical(0, default = 0)) {
    returnType(double())
    ## First handle non-zero data
    if (x != 0) {
      ## return the log probability if log = TRUE
      if (log) return(dpois(x, lambda, log = TRUE) + log(1 - zeroProb))
      ## or the probability if log = FALSE
      else return((1 - zeroProb) * dpois(x, lambda, log = FALSE))
    }
    ## From here down we know x is 0
    totalProbZero <- zeroProb + (1 - zeroProb) * dpois(0, lambda, log = FALSE)
    if (log) return(log(totalProbZero))
    return(totalProbZero)
  })


rZIP <- nimbleFunction(
  run = function(n = integer(), lambda = double(), zeroProb = double()) {
    returnType(integer())
    isStructuralZero <- rbinom(1, prob = zeroProb, size = 1)
    if (isStructuralZero) return(0)
    return(rpois(1, lambda))
  })

registerDistributions(list(
  dZIP = list(
    BUGSdist = "dZIP(lambda, zeroProb)",
    discrete = TRUE,
    range = c(0, Inf),
    types = c('value = integer()', 'lambda = double()', 'zeroProb = double()')
  )))


ZIPmodel <- nimbleModel(ZIPcode, constants = list(N = 1000), check = FALSE)


ZIPmodel$p <- .4             ## Choose values of p and lambda
ZIPmodel$lambda <- 1.8
ZIPmodel$simulate('y')       ## Simulate values of y[1]...y[100]
simulatedData <- ZIPmodel$y
simulatedData


ZIPmodel$setData(list(y = simulatedData))  ## Set those values as data in the model
cZIPmodel <- compileNimble(ZIPmodel)  


ZIPmcmc <- buildMCMC(ZIPmodel)
cZIPmcmc <- compileNimble(ZIPmcmc, project = ZIPmodel)
samples <- runMCMC(cZIPmcmc, niter = 10000)
plot(samples[,'lambda'], type = 'l', main = 'lambda trace plot')




#---------------------------------------spike and slab--------------------------------
mu <- c(3,4)
set.seed(1234)
dat <- rmvnorm(2000, mean=mu)

# SScode <- nimbleCode({
#   mycov[1:2,1:2] <-  diag(rep(1,2))
#   for (i in 1:2)
#     mu[i] ~ dunif(0,5)
#   for (i in 1:N)
#     y[i,1:2] ~ dmnorm(mu[1:2], cov = mycov[1:2,1:2]) ## Note NIMBLE allows R-like named-parameter syntax
# })

# 
# SScode <- nimbleCode({
#   mycov[1:2,1:2] <-  diag(rep(1,2))
#   mu[1] ~ dunif(0,5)
#   mu[2] <- mu[1]
#   for (i in 1:N)
#     y[i,1:2] ~ dmnorm(mu[1:2], cov = mycov[1:2,1:2]) ## Note NIMBLE allows R-like named-parameter syntax
# })

SScode <- nimbleCode({
  mycov[1:2,1:2] <-  diag(rep(1,2))
  theta ~ dunif(0,1)
  for (i in 1:3)
    mu[i] ~ dunif(0,5)
  mu[4] <- mu[3]
  
  for (i in 1:N)
    y[i,1:2] ~ dSS(mu=mu[1:4],cov=mycov[1:2,1:2], theta=theta) ## Note NIMBLE allows R-like named-parameter syntax
})

dSS <- nimbleFunction(
  run = function(x = double(1), mu = double(1), cov=double(2), theta=double(0), log = logical(0, default = 0)) {
    returnType(double())
    cov.cho <- chol(cov)
    comp1 <- theta*dmnorm_chol(x, mean=mu[1:2], cholesky=cov.cho, prec_param = FALSE)
    comp2 <- (1-theta)*dmnorm_chol(x, mean=mu[3:4], cholesky=cov.cho, prec_param = FALSE)

    totalProb <- comp1 + comp2
    if (log) return(log(totalProb))
    return(totalProb)
  })

dSS(x=c(3,4), mu=c(3,4,3,3), cov=diag(2),theta=0.5)

registerDistributions(list(
  dSS = list(
    BUGSdist = "dSS(mu, cov, theta)",
    types = c('value = double(1)', 'mu = double(1)', 'cov = double(2)', 'theta = double(0)')
  )))

SSmodel <- nimbleModel(SScode, constants = list(N = 2000), check = FALSE)


SSmodel$setData(list(y = dat))  ## Set those values as data in the model
cSSmodel <- compileNimble(SSmodel)


SSmcmc <- buildMCMC(SSmodel)
cSSmcmc <- compileNimble(SSmcmc, project = SSmodel)
samples <- runMCMC(cSSmcmc, niter = 50000,nburnin=10000,thin=10)
plot(samples[,'mu[3]'], type = 'l', main = 'c trace plot')

plot(density(samples[,'theta']), main = 'density of theta')


#---------------------------------Gaussia mixture model-------------------------------
mu1 <- c(3,4)
mu2 <- c(7,8)
# mu2 <- c(10,15)
cov1 <- matrix(c(1,0.8,0.8,1),nrow=2)
cov2 <- matrix(c(2,0.6,0.6,1),nrow=2)

set.seed(1234)

dat1 <- rmvnorm(2000, mean=mu1,sigma=cov1)
dat2 <- rmvnorm(2000, mean=mu2,sigma=cov2)
mix.pi <- runif(2000)

p <- 0.7
dat <- rbind(dat1[mix.pi>p,], dat2[mix.pi<=p,])
plot(dat)

probs <- c(1/4, 1/10, 1 - 1/4 - 1/10)
x <- rcat(n = 30, probs)
dcat(x, probs)

uppertri_mult_diag <- nimbleFunction(
  run = function(mat = double(2), vec = double(1)) {
    returnType(double(2))
    p <- length(vec)
    out <- matrix(nrow = p, ncol = p, init = FALSE)
    for(i in 1:p)
      out[ , i] <- mat[ , i] * vec[i]
    return(out)
  })


dSS <- nimbleFunction(
  run = function(x = double(1), mu = double(1), cov=double(2), theta=double(1), log = logical(0, default = 0)) {
    returnType(double(0))
    D <- length(x)
    K <- length(theta)
    totalProb <- 0
    for (i in 1:K){
      mu.k <- mu[((i-1)*D+1):(i*D)]
      cov.k <- cov[1:D, ((i-1)*D+1):(i*D)] 
      comp.k <- theta[i]*dmnorm_chol(x, mean=mu.k, cholesky=cov.k, prec_param = FALSE)
      totalProb <- totalProb + comp.k
    }
    if (log) return(log(totalProb))
    return(totalProb)
  })

dSS(x=dat[2,], mu=c(3,4,3,3), cov=cbind(diag(2),2*diag(2)),theta=c(1,0))

registerDistributions(list(
  dSS = list(
    BUGSdist = "dSS(mu, cov, theta)",
    types = c('value = double(1)', 'mu = double(1)', 'cov = double(2)', 'theta = double(1)')
  )))

SScode <- nimbleCode({
  for (i in 1:(D*K))
    sds[i] ~ dunif(0, 100)

  for ( i in 1:K) {
    Ustar[1:D,((i-1)*D+1):(i*D)] ~ dlkj_corr_cholesky(1.3, D)
    U[1:D,((i-1)*D+1):(i*D)] <- uppertri_mult_diag(Ustar[1:D, ((i-1)*D+1):(i*D)], sds[((i-1)*D+1):(i*D)])
  }

  theta[1:K] ~ ddirch(alpha[1:K])
  for (i in 1:(D*K))
    mu[i] ~  dunif(0, 100)

  for (i in 1:N)
    y[i,1:D] ~ dSS(mu=mu[1:(D*K)],cov=U[1:D,1:(D*K)], theta=theta[1:K]) ## Note NIMBLE allows R-like named-parameter syntax
})

# SScode <- nimbleCode({
#   for (i in 1:(D*K))
#     sds[i] ~ dunif(0, 100)
#   
#   for ( i in 1:K) {
#     Ustar[1:D,((i-1)*D+1):(i*D)] ~ dlkj_corr_cholesky(1.3, D)
#     U[1:D,((i-1)*D+1):(i*D)] <- uppertri_mult_diag(Ustar[1:D, ((i-1)*D+1):(i*D)], sds[((i-1)*D+1):(i*D)])
#   }
#   
#   for (i in 1:(D*K))
#     mu[i] ~  dunif(0, 100)
#   
#   for (i in 1:N)
#     y[i,1:D] ~ dmnorm(mean=mu[1:(D*K)],cholesky =U[1:D,1:(D*K)], prec_param=0) 
# })



SSmodel <- nimbleModel(SScode, constants = list(N = 2000,
                                                D = 2,
                                                K = 2,
                                                alpha = rep(0.5,2)), check = FALSE)


SSmodel$setData(list(y = dat))  ## Set those values as data in the model
cSSmodel <- compileNimble(SSmodel)


SSmcmc <- buildMCMC(SSmodel)
cSSmcmc <- compileNimble(SSmcmc, project = SSmodel)
samples <- runMCMC(cSSmcmc, niter = 45000,nburnin=5000,thin=20)
plot(samples[,'mu[1]'], type = 'l', main = 'c trace plot')

plot(density(samples[,'theta[1]']), main = 'density of theta')

calculateWAIC(samples, SSmodel)
# K = 1: 12059.89
# K = 2: 11965.95
# K = 3: 12061.66
library(coda)
samples.mcmc <- mcmc(samples)

effectiveSize(samples.mcmc)

autocorr.plot(samples.mcmc)
