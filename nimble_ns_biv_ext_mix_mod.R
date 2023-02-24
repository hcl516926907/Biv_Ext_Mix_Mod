source("KRSW/RevExp_U_Functions.r")
source("KRSW/CommonFunctions.r")
library(nimble, warn.conflicts = F)
library(mvtnorm)
library(parallel)
library(pracma)
library(Rcpp)
library(rbenchmark)



dir.in <- '/home/pgrad2/2448355h/My_PhD_Project/01_Output/Biv_Ext_Mix_Mod/biv_ext_mix_mod_simdat'

load(file.path(dir.in,'simulation_data_non_stationary.RData'))


R_pmnorm_chol <- function(lower, upper, mean, cholesky){
  sigma <- t(cholesky) %*% cholesky
  return(pmvnorm(lower=lower, upper=as.vector(upper), mean=as.vector(mean), sigma=sigma, keepAttr = F))
}


R_dmvnorm_chol <- function(x, mean, cholesky, log){
  sigma <- t(cholesky) %*% cholesky
  dvect <- dmvnorm(x, mean=mean, sigma=sigma, log=log)
  if (log) {
    return(sum(dvect))
  }else{
    return(prod(dvect))
  }
  
}

R_map_thres <- function(lower, upper, eta){
  return(lower + sweep(sigmoid(eta),2, upper-lower, "*"))
}

R_mix_prob_sg <- function(upper, mu, cholesky){
  D <- length(mu)
  return(pmnorm_chol(lower=rep(0,D), upper=upper, mean=mu, cholesky = cholesky))
}

R_mix_prob <- function(thres, mu, cholesky){
  prob.all <- apply(thres,1,R_mix_prob_sg, mu=mu, cholesky=cholesky)
  return(sum(log(1-prob.all)))
}

R_mix_prob_1 <- function(thres,mu,cholesky){
  prob.all <- sapply(seq_len(nrow(thres)), function(i) R_mix_prob_sg(thres[i,],mu,cholesky))
  return(sum(log(1-prob.all)))
}

R_mix_prob_2 <- function(thres, mu, cholesky){
  prob.all <- mclapply(seq_len(nrow(thres)), function(i) R_mix_prob_sg(thres[i,],mu,cholesky))
  return(sum(log(1-unlist(prob.all))))
}

src <-
  '
NumericVector R_mix_prob_3(NumericMatrix UB, NumericVector mu, NumericMatrix M){
  NumericVector lb (mu.length());
  NumericVector prob_vec (UB.nrow());
  
  Environment pkg = Environment::namespace_env("mvtnorm");
  Function pmvnorm = pkg["pmvnorm"];
  
  for(int i=0; i<UB.nrow(); ++i){
    SEXP val = pmvnorm(Named("lower",lb), Named("upper",UB(i,_)),Named("mean",mu),
                Named("sigma",M), Named("keepAttr", false));
    double res =  Rcpp::as<double>(val);
    prob_vec[i] = res;
  }
  return(prob_vec);
}
'


Rcpp::cppFunction(code = src)

t0 <- Sys.time()
system.time(R_mix_prob(thres, mu, cholesky))
t1 <- Sys.time()
system.time(R_mix_prob_1(thres, mu, cholesky))
t2 <- Sys.time()
system.time(R_mix_prob_2(thres, mu, cholesky))
t3 <- Sys.time()
system.time(R_mix_prob_3(thres, mu, t(cholesky) %*% cholesky))
t4 <- Sys.time()
print(c(t1-t0,t2-t1,t3-t2,t4-t3))
rbenchmark::benchmark(R_mix_prob(thres, mu, cholesky),R_mix_prob_1(thres, mu, cholesky),
                      R_mix_prob_2(thres, mu, cholesky),R_mix_prob_3(thres, mu, t(cholesky) %*% cholesky))

# system.time(R_mix_prob_sg(thres[1,], mu, cholesky))

R_pmax <- function(mat, eta){
  return(sweep(mat, 2, eta, 'pmax'))
}

mix_prob_sum <- nimbleRcall(function(thres = double(2), mu = double(1),
                                 cholesky = double(2)){}, 
                        Rfun = 'R_mix_prob_2',
                        returnType = double(0))


map_thres <- nimbleRcall(function(lower = double(1), upper = double(1),
                                  eta = double(2)){}, 
                         Rfun = 'R_map_thres',
                         returnType = double(2))


pmnorm_chol <- nimbleRcall(function(lower = double(1), upper = double(1), 
                                    mean=double(1),cholesky=double(2)){}, 
                           Rfun = 'R_pmnorm_chol',
                           returnType = double(0))



dmvnorm_chol <- nimbleRcall(function(x = double(2), mean = double(1), 
                                     cholesky=double(2),log=logical(0, default = 0)){}, 
                            Rfun = 'R_dmvnorm_chol',
                            returnType = double(0))


nim_nll_powunif_GPD <- nimbleRcall(function(theta=double(1), x=double(2), u=double(0), a.ind=double(0),
                                            lam.ind=double(0), sig.ind=double(1), gamma.ind=double(1), 
                                            lamfix=logical(0, default = 0), balthresh=logical(0, default = 0), 
                                            marg.scale.ind=double(1), marg.shape.ind=double(1)){}, 
                                   Rfun = 'nll.powunif.GPD',
                                   returnType = double(0))

nim_pmax <- nimbleRcall(function(mat = double(2), eta = double(1)){}, 
                        Rfun = 'R_pmax',
                        returnType = double(2))

dbiextmix <- nimbleFunction(
  run = function(x=double(2), theta=double(1), beta=double(2), X=double(2),
                 lower= double(1), upper= double(1),mu=double(1), 
                 cholesky=double(2), D=integer(0, default=2),
                 a.ind=double(0), lam.ind=double(0), lamfix=logical(0, default = 0), 
                 sig.ind=double(1), gamma.ind=double(1),
                 log = logical(0, default = 0)) {
    returnType(double(0))
    
    K <- length(beta[,1])
    lp <- matrix(0,nrow=dim(X)[1], ncol=D)
    for (i in 1:D){
      lp[,i] <- X[,(1+(i-1)*K):(i*K)] %*% beta[,i]
    }
    thres <- map_thres(lower,upper,lp)
    cond <- (x[,1]>thres[,1]) | (x[,2]>thres[,2])
    n.tail <- sum(cond)
    n.bulk <- sum(!cond)
    y.tail <- matrix(c(x[cond,1] - thres[cond,1], x[cond,2] - thres[cond,2]), ncol=D)
    y.bulk <- x[!cond,]
    
    thres.tail <- thres[cond,]
    log.sum <- mix_prob_sum(thres.tail, mu, cholesky)

    sig <- theta[sig.ind]
    gamma <- theta[gamma.ind]
    eta <- -sig/gamma
    
    dtail <- 0
    dbulk <- 0
    
    if (n.tail>0){
      # max of y.tail and endpoint
      y.tail.max <- nim_pmax(y.tail,eta+10^-10)
      llt <- -nim_nll_powunif_GPD(x=y.tail.max, theta=theta, u=min(y.tail)-0.01, a.ind=a.ind,
                                  lam.ind=lam.ind, sig.ind=sig.ind, gamma.ind=gamma.ind, 
                                  lamfix=lamfix, balthresh=FALSE, 
                                  marg.scale.ind=1:2, marg.shape.ind=1:2)
      if (log){
        dtail <- llt
      }else{
        dtail <- exp(llt)
      }
      
    }
    
    if (n.bulk>0){
      dbulk <- dmvnorm_chol(y.bulk, mean=mu, cholesky = cholesky, log = log )
    }
    
    if (log) {
      totalProb <- log.sum + dtail + dbulk
    }else{
      totalProb <- exp(log.sum) *dtail*dbulk
    }
    
    return(totalProb)
  })






theta <- c(1.656, 0.571, 0.451, 0.253, 0.035)
a.ind <- c(1)
lam.ind <- c(1)
sig.ind <- c(2,3)
gamma.ind <- c(4,5)
marg.scale.ind <- c(1,2)
marg.shape.ind <- c(1,2)
beta <- 0.1*cbind(c(1,2,3,4),c(-2,-3,4,5))
lower <- 6
upper <- 8
mu <- c(5,5.41)
rho=0.5
sigma <- 1.5* matrix(c(1,rho,rho,0.8),ncol=2)
cholesky <- chol(sigma)
x <- Y
D <- 2
log <- TRUE
lamfix=FALSE

t1 <- Sys.time()
# dbiextmix(x=Y, theta=theta, thres = c(5.5,5.5),mu=mu, 
#           cholesky=cholesky,
#           a.ind=a.ind, lam.ind=lam.ind, lamfix=0, 
#           sig.ind=sig.ind, gamma.ind=gamma.ind,
#           log =1)

dbiextmix(x=Y, theta=theta, beta=beta, X=cbind(X1,X2),
          lower= lower, upper= upper,mu=mu,
          cholesky=cholesky,
          a.ind=a.ind, lam.ind=lam.ind, lamfix=0,
          sig.ind=sig.ind, gamma.ind=gamma.ind,
          log =1)
t2 <- Sys.time()
print(t2-t1)

beta <- cbind(colMeans(results$samples[,c('beta[1, 1]','beta[2, 1]','beta[3, 1]','beta[4, 1]')]),
              colMeans(results$samples[,c('beta[1, 2]','beta[2, 2]','beta[3, 2]','beta[4, 2]')]))
dbiextmix(x=Y, theta=theta, beta=beta, X=cbind(X1,X2),
          lower= lower, upper= upper,mu=mu,
          cholesky=cholesky,
          a.ind=a.ind, lam.ind=lam.ind, lamfix=0,
          sig.ind=sig.ind, gamma.ind=gamma.ind,
          log =1)

rbiextmix <- nimbleFunction(
  run = function(n=integer(0), theta=double(1), beta=double(2), X=double(2), 
                 lower= double(1), upper= double(1), mu=double(1), 
                 cholesky=double(2), D=integer(0, default=2),
                 a.ind=double(0), lam.ind=double(0), lamfix=logical(0, default = 0), 
                 sig.ind=double(1), gamma.ind=double(1)) {
    returnType(double(2))
    
    totalProb <- matrix(1,nrow=n, ncol=D)
    return(totalProb)
  })
         
registerDistributions(list(
  dbiextmix = list(
    BUGSdist = "dbiextmix(theta, beta, X, lower, upper, mu, cholesky, D, 
                a.ind, lam.ind, lamfix, sig.ind, gamma.ind)",
    types = c('value = double(2)', 'theta = double(1)', 'beta = double(2)', 'X = double(2)',
              'lower = double(1)', 'upper = double(1)',
              'mu = double(1)', 'cholesky = double(2)', 'D = integer(0)', 'a.ind = double(0)', 
              'lam.ind = double(0)', 'lamfix = logical(0)', 'sig.ind = double(1)',
              'gamma.ind = double(1)')
  )))



uppertri_mult_diag <- nimbleFunction(
  run = function(mat = double(2), vec = double(1)) {
    returnType(double(2))
    p <- length(vec)
    out <- matrix(nrow = p, ncol = p, init = FALSE)
    for(i in 1:p)
      out[ , i] <- mat[ , i] * vec[i]
    return(out)
  })

BivExtMixcode <- nimbleCode({
  for (i in 1:D)
    sds[i] ~ dunif(0, 100)
  Ustar[1:D,1:D] ~ dlkj_corr_cholesky(1.3, D)
  U[1:D,1:D] <- uppertri_mult_diag(Ustar[1:D, 1:D], sds[1:D])
  
  for (i in 1:3)
    theta[i] ~ dunif(0,50)
  for (i in 4:5)
    theta[i] ~ dunif(0,1)
  for (i in 1:2)
    mu[i] ~ T(dnorm(0, sd=100),0, lower[i])
  
  for (i in 1:D.pred)
    mu_beta[i] <- 0
  for (i in 1:D)
    beta[1:D.pred,i] ~ dmnorm(mu_beta[1:D.pred], cov=cov_beta[1:D.pred,1:D.pred])
  # for (i in 1:D.pred){
  #   for (j in 1:D){
  #     beta[i,j] ~ dnorm(0, sd=100)
  #   }
  #}
  
  
  y[1:N,1:D] ~ dbiextmix(theta=theta[1:5], beta=beta[1:D.pred,1:D],X=X[1:N,1:D.pred],
                         lower=lower[1:D], upper=upper[1:D], mu=mu[1:D], 
                         cholesky=U[1:D,1:D], D=D,
                         a.ind=a.ind, lam.ind=lam.ind, lamfix=lamfix, 
                         sig.ind=sig.ind[1:D], gamma.ind=gamma.ind[1:D])
})

X1.c <- sweep(X1, 2, colMeans(X1), '-')
X2.c <- sweep(X2, 2, colMeans(X2), '-')
BivExtMixmodel <- nimbleModel(BivExtMixcode, constants = list(N = 2500, 
                                                              D = 2,
                                                              D.pred = 4,
                                                              cov_beta = 1000*diag(4),
                                                              X = cbind(X1.c,X2.c),
                                                              lower = c(6,6),
                                                              upper = c(8,8),
                                                              a.ind = 1,
                                                              lam.ind = 2,
                                                              sig.ind = c(2,3),
                                                              gamma.ind = c(4,5),
                                                              lamfix=1),
                                                              check = FALSE)


BivExtMixmodel$setData(list(y = Y))  ## Set those values as data in the model
cBivExtMixmodel <- compileNimble(BivExtMixmodel)

BivExtMixconf <- configureMCMC(BivExtMixmodel,
                               enableWAIC = TRUE, time=TRUE)
# BivExtMixconf$removeSamplers(c('beta[1:4, 1]', 'beta[1:4, 2]'))
# BivExtMixconf$addSampler(target = c('beta[1:4, 1]', 'beta[1:4, 2]'), type = 'AF_slice')

BivExtMixMCMC <- buildMCMC(BivExtMixconf)

cBivExtMixMCMC <- compileNimble(BivExtMixMCMC, project = BivExtMixmodel)

t1 <- Sys.time()
initsList <- list( beta=matrix(0,ncol=2,nrow=4))
results <- runMCMC(cBivExtMixMCMC, niter = 10000,nburnin=2500,thin=1,
                   summary = TRUE, WAIC = TRUE,setSeed = 1234, inits = initsList)
t2 <- Sys.time()
print(t2-t1)

plot(results$samples[,'beta[4, 1]'],type='l')
plot(results$samples[,'mu[1]'],type='l')

pairs(results$samples[,c('beta[1, 1]','beta[2, 1]','beta[3, 1]','beta[4, 1]',
                         'beta[1, 2]','beta[2, 2]','beta[3, 2]','beta[4, 2]')])

plot(results$samples[,c('beta[2, 1]','beta[3, 1]')])



b0 <- seq(-54,-44, 0.5)
b1 <- seq(88, 106, 0.5)
ll <- function(b0,b1){
  set.seed(1)
  X <-rnorm(1000,mean=3)
  X <- X 
  y <- rnorm(1000, 3+5*X)
  lp <- b0 + X*b1
  return(sum(dnorm(y, mean=lp,log=TRUE  )))
}
z <- matrix(NA, nrow=length(b0), ncol=length(b1))

# need to double check the xlab and ylab

for (i in 1:length(b0)){
  for (j in 1:length(b1)){
    beta.tmp <- 0.1*cbind(c(1,2,3,4),c(-2,-3,4,5))
    beta.tmp[2,1] <- b0[i]
    beta.tmp[3,1] <- b1[j]
    z[i,j] <- dbiextmix(x=Y, theta=theta, beta=beta.tmp, X=cbind(X1.c,X2.c),
                        lower= lower, upper= upper,mu=mu,
                        cholesky=cholesky,
                        a.ind=a.ind, lam.ind=lam.ind, lamfix=0,
                        sig.ind=sig.ind, gamma.ind=gamma.ind,
                        log =1)
    
  }
}


contour(b0, b1, z, main='contour of beta[1,2] by beta[2,2]',xlab='beta[1,2]',ylab='beta[2,2]')
