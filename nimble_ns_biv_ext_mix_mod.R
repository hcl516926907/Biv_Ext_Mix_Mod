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

# t0 <- Sys.time()
# system.time(R_mix_prob(thres, mu, cholesky))
# t1 <- Sys.time()
# system.time(R_mix_prob_1(thres, mu, cholesky))
# t2 <- Sys.time()
# system.time(R_mix_prob_2(thres, mu, cholesky))
# t3 <- Sys.time()
# system.time(R_mix_prob_3(thres, mu, t(cholesky) %*% cholesky))
# t4 <- Sys.time()
# print(c(t1-t0,t2-t1,t3-t2,t4-t3))
# rbenchmark::benchmark(R_mix_prob(thres, mu, cholesky),R_mix_prob_1(thres, mu, cholesky),
#                       R_mix_prob_2(thres, mu, cholesky),R_mix_prob_3(thres, mu, t(cholesky) %*% cholesky))

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
      lp[,i] <- X[,(1+(i-1)*K):(i*K)] %*% matrix(beta[,i],ncol=1)
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
      y.min <- eta
      for (i in 1:D){
        y.min[i] <- min(y.tail[,i])
      }
      if (all(y.min>eta)){
        llt <- -nim_nll_powunif_GPD(x=y.tail, theta=theta, u=min(y.tail)-0.01, a.ind=a.ind,
                                    lam.ind=lam.ind, sig.ind=sig.ind, gamma.ind=gamma.ind, 
                                    lamfix=lamfix, balthresh=FALSE, 
                                    marg.scale.ind=1:2, marg.shape.ind=1:2)
        if (log){
          dtail <- llt
        }else{
          dtail <- exp(llt)
        }
      }else{
        if (log) dtail <- -10^10
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



X <- cbind(X1,X2)
# X <- matrix(1,nrow=1000,ncol=4)
theta <- c(1.656, 0.571, 0.451, 0.253, 0.035)
a.ind <- c(1)
lam.ind <- c(1)
sig.ind <- c(2,3)
gamma.ind <- c(4,5)
marg.scale.ind <- c(1,2)
marg.shape.ind <- c(1,2)
beta <- cbind(c(0.1,0.2),c(0.3,-0.4))
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
# dbiextmix(x=Y, theta=theta, thres = c(7.049958,6.851115),mu=mu,
#           cholesky=cholesky,
#           a.ind=a.ind, lam.ind=lam.ind, lamfix=0,
#           sig.ind=sig.ind, gamma.ind=gamma.ind,
#           log =1)

dbiextmix(x=Y, theta=theta, beta=beta, X=X,
          lower= lower, upper= upper,mu=mu,
          cholesky=cholesky,
          a.ind=a.ind, lam.ind=lam.ind, lamfix=0,
          sig.ind=sig.ind, gamma.ind=gamma.ind,
          log =1)

# -2959.791
t2 <- Sys.time()
print(t2-t1)

burnin <- 1000
n.itr <- nrow(results$samples)
# beta <- cbind(colMeans(results$samples[,c('beta[1, 1]','beta[2, 1]','beta[3, 1]','beta[4, 1]')]),
#               colMeans(results$samples[,c('beta[1, 2]','beta[2, 2]','beta[3, 2]','beta[4, 2]')]))
beta <- cbind(colMeans(results$samples[burnin:n.itr,c('beta[1, 1]','beta[2, 1]')]),
              colMeans(results$samples[burnin:n.itr,c('beta[1, 2]','beta[2, 2]')]))

beta[1,1] <- 10
beta[2,1] <- 3
# beta[1,2] <- 81.1322
#problematic beta 
# [,1]    [,2]
# beta[1, 1] -86.3079 81.1322
# beta[2, 1]   0.0000  0.0000

theta <- colMeans(results$samples[burnin:n.itr,c('theta[1]','theta[2]','theta[3]',
                                                 'theta[4]','theta[5]')])

# theta <- c(17.237010676, 0.902869762, 0.001457657, 0.021203881, 0.980473508 )
#problematic theta 
# theta[1]     theta[2]     theta[3]     theta[4]     theta[5] 
# 17.237010676  0.902869762  0.001457657  0.021203881  0.980473508 

#problematic mu
mu <- colMeans(results$samples[burnin:n.itr,c('mu[1]','mu[2]')])
# mu <- c(5.131380 ,5.347119 )
# mu[1]    mu[2] 
# 5.131380 5.347119 


dbiextmix(x=Y, theta=theta, beta=beta, X=X,
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
  for (i in 1:D){
    beta[1:D.pred,i] ~ dmnorm(mu_beta[1:D.pred], cov=cov_beta[1:D.pred,1:D.pred])
  }
  
  
  y[1:N,1:D] ~ dbiextmix(theta=theta[1:5], beta=beta[1:D.pred,1:D],X=X[1:N,1:D.pred],
                         lower=lower[1:D], upper=upper[1:D], mu=mu[1:D], 
                         cholesky=U[1:D,1:D], D=D,
                         a.ind=a.ind, lam.ind=lam.ind, lamfix=lamfix, 
                         sig.ind=sig.ind[1:D], gamma.ind=gamma.ind[1:D])
})

X1.c <- sweep(X1, 2, c(0,mean(X1[,2])), '-')
X2.c <- sweep(X2, 2, c(0,mean(X2[,2])), '-')
BivExtMixmodel <- nimbleModel(BivExtMixcode, constants = list(N = 1000, 
                                                              D = 2,
                                                              D.pred = 2,
                                                              cov_beta = 25*diag(2),
                                                              X = cbind(X1.c,X2.c),
                                                              # X = matrix(1,nrow=1000,ncol=4),
                                                              lower = c(5.5,5.7),
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
results <- runMCMC(cBivExtMixMCMC, niter = 10000,nburnin=5000,thin=1,
                   summary = TRUE, WAIC = TRUE,setSeed = 1234)
t2 <- Sys.time()
print(t2-t1)

plot(results$samples[,'beta[2, 2]'],type='l', main='Traceplot of beta[2, 2]')
plot(results$samples[,'beta[1, 2]'],type='l')

# pairs(results$samples[,c('beta[1, 1]','beta[2, 1]','beta[3, 1]','beta[4, 1]',
#                          'beta[1, 2]','beta[2, 2]','beta[3, 2]','beta[4, 2]')])
pairs(results$samples[,c('beta[1, 1]','beta[2, 1]',
                         'beta[1, 2]','beta[2, 2]')])
plot(results$samples[,c('beta[2, 1]','beta[3, 1]')])

X <- cbind(X1.c,X2.c)
K <- length(beta[,1])
lp.samples <- list()
for (i in 1:nrow(results$samples)){
  beta.samples <- results$samples[i,c('beta[1, 1]','beta[2, 1]',
                                   'beta[1, 2]','beta[2, 2]')] 
  lp.tmp <- matrix(NA, nrow=nrow(X), ncol=2)
  for (j in 1:D){
    lp.tmp[,j] <- X[,(1+(j-1)*K):(j*K)] %*% matrix(beta.samples[(1+(j-1)*K):(j*K)],ncol=1)
  }
  lp.samples[[i]] <- lp.tmp
}

thres.samples <-  lapply(lp.samples,map_thres,lower=lower,upper=upper)

dir.out <- '/home/pgrad2/2448355h/My_PhD_Project/01_Output/Biv_Ext_Mix_Mod/nimble_ns_biv_ext_mix_mod'

save(results, thres.samples, file=file.path(dir.out,'sp_rw_l5dot5_l5dot7_u8_u8.RData'))

library(ggplot2)
load(file.path(dir.out,'sp_rw_l5dot5_l5dot7_u8_u8.RData'))
thres.sp.5.5_5.7 <- thres.samples
load(file.path(dir.out,'sp_rw_l6_l6_u8_u8.RData'))
thres.sp.6_6 <- thres.samples


obs.idx <- 650
df1 <- as.data.frame(thres.sp.5.5_5.7[[obs.idx]])
df1$catg <- '5.5_5.7'
df2 <- as.data.frame(thres.sp.6_6[[obs.idx]])
df2$catg <- '6_6'
df <- rbind(df1, df2)
p <- ggplot(df, mapping = aes(x = V1, y = V2)) 
p + geom_density2d(aes(color = catg)) + labs(title=paste('observation ',obs.idx))

p2 <- ggplot(as.data.frame(thres.sp.6_6[[1]]), mapping = aes(x = V1, y = V2)) 
p2 + geom_density2d()

x.axis <- seq(0, 0.6, 0.02)
y.axis <- seq(-0.5,0.1, 0.02)
u <- min(c(x.axis,y.axis))-0.01
z <- matrix(NA,nrow=length(x.axis), ncol=length(y.axis))
for (i in 1:length(x.axis)){
  for (j in 1:length(y.axis)){
      beta.tmp <- beta
      beta.tmp[1,2] <- x.axis[i]
      beta.tmp[2,2] <- y.axis[j]
      
      z[i,j] <- dbiextmix(x=Y, theta=theta, beta=beta.tmp, X=cbind(X1,X2),
                          lower= lower, upper= upper,mu=mu,
                          cholesky=cholesky,
                          a.ind=a.ind, lam.ind=lam.ind, lamfix=0,
                          sig.ind=sig.ind, gamma.ind=gamma.ind,
                          log =1)
    }

}
z.new <- z
z.new[which(z.new < -7400)] <- -7400
contour(x.axis, y.axis, z.new, main='contour of bivariate GP density',xlab='x1',ylab='x2')
