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
  
  thres.ind <- idx_2_ind(thres.idx,d)
  a.ind <- idx_2_ind(a.idx,d) 
  lam.ind <- idx_2_ind(lam.idx,d) 
  sig.ind <- idx_2_ind(sig.idx,d)  
  gamma.ind <- idx_2_ind(gamma.idx,d) 

  # mu.ind <- idx_2_ind(mu.idx,d) + (4*d):(5*d-1)
  # cov.ind <- idx_2_ind(cov.idx,d) + (5*d):(6*d-1)
  
  par.a <-  runif(d,0,100)
  par.lam <- runif(d,0,100)
  par.sig <- runif(d,0,100)
  par.gamma <- runif(d,0,1)
  lb <- apply(X, 2, quantile, prob=0.7)
  ub <- apply(X, 2, quantile, prob=0.99)
  for (i in 1:d){
    par.thres <- c(runif(1, lb[i], ub[i]))
  }
  
  par.mu <- rmvnorm(1, mean=rep(0,d), sigma=100*diag(d))
  par.cov <- rwishart(nu=100,S=diag(d))
  
  # the last component of lam is fixed to be 1 to resolve identificaiton issue
  if (which(sum(lam.ind == tail(lam.ind, n=1))==1)){
    par.all <- c(par.thres[thres.ind],par.a[a.ind], par.lam[lam.ind[1:(d-1)]],
                 par.sig[sig.ind], par.gamma[gamma.ind])
    #TBD
    a.ind <- a.ind + d
    lam.ind <- lam.ind + d-1
    sig.ind <- sig.ind
  }else{
    equal.last.idx <- which(lam.ind == tail(lam.ind, n=1))
    par.lam[equal.last.idx] <- 1
    par.all <- c(par.thres[thres.ind],par.a[a.ind], par.lam[lam.ind[1:(d-1)]],
                 par.sig[sig.ind], par.gamma[gamma.ind])
  }
  
}
  

