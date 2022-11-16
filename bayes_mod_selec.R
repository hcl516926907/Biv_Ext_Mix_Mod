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
library(reticulate)
library(DEoptim)
library(numDeriv)
library(nimble)


dir.in <- '/home/pgrad2/2448355h/My_PhD_Project/01_Output/Biv_Ext_Mix_Mod/biv_ext_mix_mod_simdat'

load(file.path(dir.in,'simulation_data_pos_shape_fixed.RData'))


X <- X.p09
#X <- X.p05.10t

# total 5%, not marginal 5%
u.x <- u.x.p09
# u.x <- u.x.p05.10t



#prior for a
prior.a <- function(theta,a.ind){
  d <- length(a.ind)
  prior <- 0
  for (i in 1:d){
    prior <- prior + dunif(theta[a.ind][i],0,100,log=T)
  }
  return(prior)
}

#prior for lambda
prior.lam <- function(theta, lam.ind){
  d <- length(lam.ind)
  prior <- 0
  for (i in 1:d){
    #prior <- prior*dgamma(theta[lam.ind][i], shape=1, rate=1)
    prior <- prior + dunif(theta[lam.ind][i],0,100, log=T)
  }
  return(prior)
}

#prior for sigma
prior.sig <- function(theta, sig.ind){
  d <- length(sig.ind)
  prior <- 0
  for (i in 1:d){
    prior <- prior+dunif(theta[sig.ind][i],0,100,log=T)
  }
  return(prior)
}


#prior for sigma* that is not dependent on the threshold
prior.sig1 <- function(theta, sig.ind){
  d <- length(sig.ind)
  prior <- 0
  for (i in 1:d){
    prior <- prior+dunif(theta[sig.ind][i],-15,100,log=T)
  }
  return(prior)
}

#prior for gamma
prior.gamma <- function(theta, gamma.ind){
  d <- length(gamma.ind)
  prior <- 0
  for (i in 1:d){
    prior <- prior+dunif(theta[gamma.ind][i],0,1,log=T)
  }
  return(prior)
}

#prior for thresholds
prior.thres <- function(theta, thres.ind, X){
  d <- length(thres.ind)
  lb <- apply(X, 2, quantile, prob=0.7)
  ub <- apply(X, 2, quantile, prob=0.99)
  prior <- 0
  for (i in 1:d){
    prior <- prior+dunif(theta[thres.ind][i], lb[i], ub[i],log=T)
  }
  return(prior)
}


# decompose covariance matrix as the product of standard deviations and correclations
# Use LKJ distribution as the prior for correlation matrix
# reparamatrize correlation by cholesky decomposition.
cov.cho.ind <- c(12,14,15)

prior.bulk <- function(theta, mu.ind, cov.cho.ind, d){
  mu <- theta[mu.ind]
  p_mu <- dmvnorm(mu, mean=rep(0,d), sigma=100*diag(d),log=T)
  L <- matrix(0, nrow=d, ncol=d)
  L[upper.tri(L, diag = TRUE)] <- theta[cov.cho.ind]
  
  # ensure positive value for the diagonal items
  diag(L) <- exp(diag(L))

  D <- sqrt(diag(t(L) %*% L))
  p_sig <- 0
  for (i in 1:d){
    p_sig_d <- dcauchy(D[i],0,2.5,log=T)
    p_sig <- p_sig + p_sig_d
  }
  cor.cho <- L %*% diag(1/D, nrow=d) 
  p_corr <- dlkj_corr_cholesky(cor.cho, eta=2, p=d, log = T)
  
  return(p_mu+p_sig+p_corr)
}

check.thres <- function(thres, X){
  lb <- apply(X, 2, quantile, prob=0.7)
  ub <- apply(X, 2, quantile, prob=0.99)
  if (any(thres<lb) | any(thres>ub)){
    return('Invalid Thresholds')
  }else{
    return('Pass')
  }
}
#full likelihood
ll.tgt <- function(theta.all, X, thres.ind, mu.ind, cov.ind, d=2,
                   a.ind, lam.ind, lamfix=F, sig.ind, gamma.ind,
                   marg.scale.ind, marg.shape.ind, balthresh=F){
  thres <- theta.all[thres.ind]
  
  #reparamatrize the matrix by Cholesky decomposition
  L <- matrix(0, nrow=d, ncol=d)
  L[upper.tri(L, diag = TRUE)] <- theta.all[cov.ind]
  diag(L) <- exp(diag(L))
  Sigma <- t(L) %*% L
  
  #exclude cases where any of the diagnal elements is 0. 
  if (any(diag(Sigma)==0)) {return(-Inf)}
  

  if (check.thres(thres,X)=='Invalid Thresholds') {return(-Inf)}
  
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
    llt <- -10^100
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
  lp.a <- prior.a(theta.all, a.ind)
  lp.sig <- prior.sig(theta.all, sig.ind)
  lp.gamma <- prior.gamma(theta.all, gamma.ind)
  
  lp.thres <- prior.thres(theta.all, thres.ind, X)
  # lp.jef <- log(prior.jef(theta.all, cov.ind, d))
  lp.bulk <- prior.bulk(theta.all, mu.ind, cov.ind,d)
  if (lamfix==F){
    lp.lam <- prior.lam(theta.all, lam.ind)
    return(llb + llt + n.tail*log(1-pi) + lp.bulk + lp.thres + 
             lp.a + lp.lam + lp.sig + lp.gamma)
  }
  else{
    return(llb + llt + n.tail*log(1-pi) + lp.bulk + lp.thres +
             lp.a + lp.sig + lp.gamma)
  }
}


p.log <- function(x){
  return(ll.tgt(theta.all=x, X=X, mu.ind = 8:9, cov.ind=10:12, d=2, thres.ind = 1:2, 
                a.ind=3, lamfix=TRUE, 
                sig.ind=4:5, gamma.ind=6:7, 
                marg.scale.ind=1:2, marg.shape.ind=1:2))
}

p.log(c( c(5.5,6.2), rep(1,3),0,0, rep(1,2), 1,0,1))

dec_2_bin <- function(x, d){
  if(x<2^d){
    bin <- paste(rev(as.integer(intToBits(x))), collapse="")
    return(substr(bin, nchar(bin)-d+1, nchar(bin)))
  }else{
    print(paste('The integer exceeds the largest decimal in dimension',d))
  }
}

idx_2_ind <- function(idx,d){
  bin <- dec_2_bin(idx,d)
  bin.by.char <- unlist(strsplit(bin,split=''))
  equal.idx <- which(bin.by.char=='1')
  ind <- 1:d
  if (length(equal.idx)>0){
    ind[equal.idx] <- min(equal.idx)
  }
  return(ind)
}


thres.idx <- 1
a.idx <- 2
lam.idx <- 1
sig.idx <- 1
gamma.idx <- 1
d <- 2
n.simu <- 100000

#-----------------------Laplace Approximaiton------------------

a.idx <- 3
lam.idx <- 3
sig.idx <- 1
gamma.idx <- 1

thres.equ.ind <- idx_2_ind(thres.idx,d)
a.equ.ind <- idx_2_ind(a.idx,d) 
lam.equ.ind <- idx_2_ind(lam.idx,d) 
sig.equ.ind <- idx_2_ind(sig.idx,d)  
gamma.equ.ind <- idx_2_ind(gamma.idx,d) 

# mu.ind <- idx_2_ind(mu.idx,d) + (4*d):(5*d-1)
# cov.ind <- idx_2_ind(cov.idx,d) + (5*d):(6*d-1)


par.a <-  runif(d,0,5)
par.lam <- runif(d,0,5)
par.sig <- c(3,3)
par.gamma <- runif(d,0,1)
lb <- apply(X, 2, quantile, prob=0.9)
ub <- apply(X, 2, quantile, prob=0.99)
par.thres <- runif(1, lb[1], ub[1])
for (j in 2:d){
  par.thres <- c(par.thres, runif(1, lb[j], ub[j]))
}

par.mu <- rmvnorm(1, mean=rep(3,d), sigma=1*diag(d))
par.cov <- rep(1,d*(d+1)/2)

# the last component of lam is fixed to be 1 to resolve identification issue
if (which(sum(lam.equ.ind == tail(lam.equ.ind, n=1))>=1)){
  equal.last.idx <- which(lam.equ.ind == tail(lam.equ.ind, n=1))
  par.lam[equal.last.idx] <- 1
}

par.all <- c(par.thres[thres.equ.ind],par.a[a.equ.ind], par.lam[lam.equ.ind[1:(d-1)]],
             par.sig[sig.equ.ind], par.gamma[gamma.equ.ind], par.mu, par.cov)
thres.ind <- thres.equ.ind
a.ind <- a.equ.ind + d
lam.ind <-  lam.equ.ind[1:(d-1)] + 2*d
sig.ind <- sig.equ.ind + 3*d - 1
gamma.ind <- gamma.equ.ind + 4*d - 1
mu.ind <- 1:d + 5*d - 1
# only for dimension 2
cov.ind <- 1:(d*(d+1)/2) + 6*d - 1

lam.is.1 <- which(par.lam[lam.equ.ind[1:(d-1)]]==1)

#fix the lam parameter, change the sign of ll.tgt
ll.tgt.1 <- function(theta.all, X.dat, thres.ind, mu.ind, cov.ind, d=2,
                           a.ind, lam.ind, lamfix=F, sig.ind, gamma.ind,
                           marg.scale.ind, marg.shape.ind, balthresh=F, lam.is.1){
  if (length(lam.is.1)>0){
    theta.all[lam.ind][lam.is.1] <- 1
  }
  return(-ll.tgt(theta.all, X.dat, thres.ind, mu.ind, cov.ind, d,
                  a.ind, lam.ind, lamfix=F, sig.ind, gamma.ind,
                  marg.scale.ind=1:d, marg.shape.ind=1:d, balthresh=F))
}

#Use differential evolution to find the optimal point

# lower <-  c(lb,rep(0,2),0,rep(0,6), -1,rep(-5,2))
# upper<- c(ub,rep(100,2),1,rep(100,2),rep(1,2), rep(50,5))
# 
# t1 <- Sys.time()
# res <- DEoptim(fn=ll.tgt.1, lower=lower, upper=upper, 
#                control = DEoptim.control(itermax = 1000,parallelType='auto',trace=20,
#                                          NP=200,CR=0.8),
#       X.dat=X, thres.ind=thres.ind,mu.ind=mu.ind, cov.ind=cov.ind, d=d,
#       a.ind=a.ind, lam.ind=lam.ind, lamfix=F, sig.ind=sig.ind, gamma.ind=gamma.ind,
#       marg.scale.ind=1:d, marg.shape.ind=1:d, balthresh=F, lam.is.1=lam.is.1)
# t2 <- Sys.time()
# print(t2-t1)

# res$optim$bestmem
# 
# hess.mat <- hessian(ll.tgt.1,  c(5.205481,6.358980, 1.367193,2.809423,0.742046, 0.532623, 0.382785, 0.278435, 0.092925, 3.541122, 4.400020, 0.218236, 0.806774, -0.308059), method="Richardson", method.args=list(), 
#         X.dat=X, thres.ind=thres.ind,mu.ind=mu.ind, cov.ind=cov.ind, d=d,
#         a.ind=a.ind, lam.ind=lam.ind, lamfix=F, sig.ind=sig.ind, gamma.ind=gamma.ind,
#         marg.scale.ind=1:d, marg.shape.ind=1:d, balthresh=F, lam.is.1=lam.is.1)



  

#true parameters
# c(5.55,6.213, 1.656, 0.571, 0.451, 0.253, 0.035, 3.550, 4.413, 1.5, 0.975, 0.975, 1.2)



#-----------------------Nested Sampling-----------------------
np = import("numpy")
un = import("ultranest")  # pip or conda install it if you don't have it
p.log <- function(x){
  if (length(lam.is.1)>0){
    x[lam.ind][lam.is.1] <- 1
  }
  L <- ll.tgt(theta.all=x, X=X, mu.ind = mu.ind, cov.ind=cov.ind, d=d, thres.ind = thres.ind, 
                a.ind=a.ind,lam.ind=lam.ind, lamfix=F, 
                sig.ind=sig.ind, gamma.ind=gamma.ind, 
                marg.scale.ind=1:2, marg.shape.ind=1:2)
  return(np$asarray(L))
}


mylikelihood <- function(params){
  L <- apply(params,MARGIN=1, FUN=p.log)
  return(np$asarray(L))
  }
  

mytransform <- function(params) {
  n.params <- ncol(params)
  lower <- c(lb,rep(0,2),0,rep(0,6), -1,rep(-5,2))
  upper <- c(ub,rep(100,2),1,rep(100,2),rep(1,2), rep(50,5))
  
  inv.trans <- function(x, lower,upper){
    return(lower + (upper-lower)*x)
  }
  
  new.params <- apply(params,MARGIN=1,FUN=inv.trans, lower=lower, upper=upper)
  return(t(new.params))
}

t1=Sys.time()
uns <- import('ultranest.stepsampler')
paramnames <- c('u1','u2','a1','a2','b1','sigma1','sigma2','gamma1','gamma2','mu1','mu2',
              'cov1','cov2','cov3')
sampler = un$ReactiveNestedSampler(paramnames, mylikelihood, transform=mytransform,vectorized = T)
nsteps = 2 * length(paramnames)
sampler$stepsampler <- uns$SliceSampler(
  nsteps=nsteps,
  generate_direction=uns$generate_mixture_random_direction
  # adaptive_nsteps=False,
  # max_nsteps=400
)
results = sampler$run()

t2 <- Sys.time()

print(t2-t1)

dir.out <- '/home/pgrad2/2448355h/My_PhD_Project/01_Output/Biv_Ext_Mix_Mod/bayes_mod_selec'
save(results,  file.path(dir.out,'ultranest_res.RData'))
