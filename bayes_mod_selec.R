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

# proper prior for the bulk
prior.bulk <- function(theta, mu.ind, cov.ind,d){
  mu <- theta[mu.ind]
  mat <- matrix(theta[cov.ind], nrow=d)
  p1 <- dmvnorm(mu, mean=rep(0,d), sigma=100*diag(d))
  p2 <- dwishart(mat, nu=100,S=diag(d))
  return(p1*p2)
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
  # lp.jef <- log(prior.jef(theta.all, cov.ind, d))
  lp.bulk <- log(prior.bulk(theta.all, mu.ind, cov.ind,d))
  
  if (lamfix==F){
    lp.lam <- log(prior.lam(theta.all, lam.ind))
    return(llb + llt + n.tail*log(1-pi) + lp.bulk + lp.thres + 
             lp.a + lp.lam + lp.sig + lp.gamma)
  }
  else{
    return(llb + llt + n.tail*log(1-pi) + lp.bulk + lp.thres +
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



p.y.cond.m <- function(thres.idx,a.idx,lam.idx,sig.idx,gamma.idx,d,n.simu){
  # indicator for which parameters are equal
  thres.equ.ind <- idx_2_ind(thres.idx,d)
  a.equ.ind <- idx_2_ind(a.idx,d) 
  lam.equ.ind <- idx_2_ind(lam.idx,d) 
  sig.equ.ind <- idx_2_ind(sig.idx,d)  
  gamma.equ.ind <- idx_2_ind(gamma.idx,d) 

  # mu.ind <- idx_2_ind(mu.idx,d) + (4*d):(5*d-1)
  # cov.ind <- idx_2_ind(cov.idx,d) + (5*d):(6*d-1)

  log.prob.y.val <- rep(NA, n.simu)
  
  for (i in 1:n.simu){ 
    par.a <-  runif(d,0,5)
    par.lam <- runif(d,0,5)
    par.sig <- runif(d,0,5)
    par.gamma <- runif(d,0,1)
    lb <- apply(X, 2, quantile, prob=0.9)
    ub <- apply(X, 2, quantile, prob=0.99)
    par.thres <- runif(1, lb[1], ub[1])
    for (j in 2:d){
      par.thres <- c(par.thres, runif(1, lb[j], ub[j]))
    }
    
    par.mu <- rmvnorm(1, mean=rep(3,d), sigma=1*diag(d))
    par.cov <- rwishart(nu=2,S=diag(d))
    
    # the last component of lam is fixed to be 1 to resolve identification issue
    if (which(sum(lam.equ.ind == tail(lam.equ.ind, n=1))>=1)){
      equal.last.idx <- which(lam.equ.ind == tail(lam.equ.ind, n=1))
      par.lam[equal.last.idx] <- 1
    }
    
    par.all <- c(par.thres[thres.equ.ind],par.a[a.equ.ind], par.lam[lam.equ.ind[1:(d-1)]],
                 par.sig[sig.equ.ind], par.gamma[gamma.equ.ind], par.mu, par.cov)
    thres.ind <- 1:d
    a.ind <- 1:d + d
    lam.ind <-  1:(d-1) + 2*d
    sig.ind <- 1:d + 3*d - 1
    gamma.ind <- 1:d + 4*d - 1
    mu.ind <- 1:d + 5*d - 1
    # only for dimension 2
    cov.ind <- 1:(d^2) + 6*d - 1
   
   
    log.prob.y <- ll.tgt(par.all, X, thres.ind, mu.ind, cov.ind, d,
                     a.ind, lam.ind, lamfix=F, sig.ind, gamma.ind,
                     marg.scale.ind=1:d, marg.shape.ind=1:d, balthresh=F)
    log.prob.y.val[i] <- log.prob.y
  }

  return(log.prob.y.val) 
}

thres.idx <- 1
a.idx <- 2
lam.idx <- 1
sig.idx <- 1
gamma.idx <- 1
d <- 2
n.simu <- 100000

t1 <- Sys.time()
res <- p.y.cond.m(thres.idx,a.idx,lam.idx,sig.idx,gamma.idx,d,n.simu)
t2<- Sys.time()
print(t2-t1)


t1 <- Sys.time()
res1 <- p.y.cond.m(thres.idx,a.idx,lam.idx,sig.idx,gamma.idx,d,n.simu)
t2<- Sys.time()
print(t2-t1)


t1 <- Sys.time()
res2 <- p.y.cond.m(thres.idx,a.idx,lam.idx,sig.idx,gamma.idx,d,n.simu)
t2<- Sys.time()
print(t2-t1)


c.norm <- 9930
mean(exp(log.prob.y.val+c.norm))
c.norm <- 7470.141
mean(exp(res1+c.norm),na.rm=T)

#-----------------------Nested Sampling-----------------------
np = import("numpy")
un = import("ultranest")  # pip or conda install it if you don't have it

#-----------------------Laplace Approximaiton------------------



log.evidence.la <- function(thres.idx,a.idx,lam.idx,sig.idx,gamma.idx,d){
  # indicator for which parameters are equal
  thres.equ.ind <- idx_2_ind(thres.idx,d)
  a.equ.ind <- idx_2_ind(a.idx,d) 
  lam.equ.ind <- idx_2_ind(lam.idx,d) 
  sig.equ.ind <- idx_2_ind(sig.idx,d)  
  gamma.equ.ind <- idx_2_ind(gamma.idx,d) 
  
  # mu.ind <- idx_2_ind(mu.idx,d) + (4*d):(5*d-1)
  # cov.ind <- idx_2_ind(cov.idx,d) + (5*d):(6*d-1)
  

  par.a <-  runif(d,0,5)
  par.lam <- runif(d,0,5)
  par.sig <- runif(d,0,5)
  par.gamma <- runif(d,0,1)
  lb <- apply(X, 2, quantile, prob=0.9)
  ub <- apply(X, 2, quantile, prob=0.99)
  par.thres <- runif(1, lb[1], ub[1])
  for (j in 2:d){
    par.thres <- c(par.thres, runif(1, lb[j], ub[j]))
  }
  
  par.mu <- rmvnorm(1, mean=rep(3,d), sigma=1*diag(d))
  par.cov <- rwishart(nu=2,S=diag(d))
  
  # the last component of lam is fixed to be 1 to resolve identification issue
  if (which(sum(lam.equ.ind == tail(lam.equ.ind, n=1))>=1)){
    equal.last.idx <- which(lam.equ.ind == tail(lam.equ.ind, n=1))
    par.lam[equal.last.idx] <- 1
  }
  
  par.all <- c(par.thres[thres.equ.ind],par.a[a.equ.ind], par.lam[lam.equ.ind[1:(d-1)]],
               par.sig[sig.equ.ind], par.gamma[gamma.equ.ind], par.mu, par.cov)
  thres.ind <- 1:d
  a.ind <- 1:d + d
  lam.ind <-  1:(d-1) + 2*d
  sig.ind <- 1:d + 3*d - 1
  gamma.ind <- 1:d + 4*d - 1
  mu.ind <- 1:d + 5*d - 1
  # only for dimension 2
  cov.ind <- 1:(d^2) + 6*d - 1
  
  optim(ll.tgt, par=par.all, X=X, thres.ind=thres.ind,mu.ind=mu.ind, cov.ind=cov.ind, d=d,
        a.ind=a.ind, lam.ind=lam.ind, lamfix=F, sig.ind=sig.ind, gamma.ind=gamma.ind,
        marg.scale.ind=1:d, marg.shape.ind=1:d, balthresh=F,control=list(maxit=1000,reltol=1e-6))
  
  ll.tgt <- function(theta.all, X, thres.ind, mu.ind, cov.ind, d=2,
                     a.ind, lam.ind, lamfix=F, sig.ind, gamma.ind,
                     marg.scale.ind, marg.shape.ind, balthresh=F)
    
  log.prob.y <- ll.tgt(par.all, X, thres.ind, mu.ind, cov.ind, d,
                       a.ind, lam.ind, lamfix=F, sig.ind, gamma.ind,
                       marg.scale.ind=1:d, marg.shape.ind=1:d, balthresh=F)

  
  return(log.prob.y) 
}

