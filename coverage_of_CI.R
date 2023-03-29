source("KRSW/RevExp_U_Functions.r")
source("KRSW/RevExp_T_Functions.r")
source("KRSW/MVGaussian_T_Functions.r")
source("KRSW/CommonFunctions.r")
source("KRSW/Gumbel_U_Functions.r")
source("KRSW/Gumbel_T_Functions.r")
source("KRSW/ModelDiagnosticsNewNames.r")
library(mvtnorm)
library(numDeriv)
library(nimble)

######### Multivariate normal distribution##########
set.seed(1234)
D <- 7
N <- 200

mu <- 1:D
A <- matrix(runif(D^2),nrow=D)
cov <- t(A)%*%A

y <- rmvnorm(200, mean=mu, sigma=cov)

ll <- function(x){
  return(-sum(dmvnorm(y, mean=x,sigma=cov, log=T)))
}

opt1 <- optim(rep(1,7),ll, control=list(maxit=10000,reltol=1e-6))
opt1.1 <- optim(opt1$par, ll,control=list(maxit=10000,reltol=1e-6))
opt1.2 <- optim(opt1.1$par, ll, control=list(maxit=10000,reltol=1e-6))
opt1.3 <- optim(opt1.2$par, ll, control=list(maxit=10000,reltol=1e-6),hessian=TRUE)
opt1.3$par

hessian(ll,opt1.3$par)
numDeriv::hessian(ll,opt1.3$par)
width <- 1.96*sqrt(diag(solve(opt1.3$hessian)))
ll(1:7)
ll(opt1.3$par)

N.sim <- 100
par.mat <- matrix(NA, nrow=N.sim, ncol=length(mu))
CI.width <- matrix(NA, nrow=N.sim, ncol=length(mu))
for (i in 1:N.sim){
  set.seed(i)
  y <- rmvnorm(200, mean=mu, sigma=cov)
  opt1 <- optim(mu,ll, control=list(maxit=10000,reltol=1e-6))
  opt1.1 <- optim(opt1$par, ll,control=list(maxit=10000,reltol=1e-6))
  opt1.2 <- optim(opt1.1$par, ll, control=list(maxit=10000,reltol=1e-6))
  opt1.3 <- optim(opt1.2$par, ll, control=list(maxit=10000,reltol=1e-6),hessian=TRUE)
  par.mat[i,] <- opt1.3$par
  CI.width[i,] <- 1.96*sqrt(diag(solve(opt1.3$hessian)))
  print(i)
}

upper <- par.mat + CI.width
lower <- par.mat - CI.width
within.ind <- rep(NA, length(theta))
for (i in 1:length(mu)){
  cond <- (mu[i]>=lower[,i])&(mu[i]<=upper[,i])
  cond.nona <- cond[!is.na(cond)]
  within.ind[i] <- sum(cond.nona) /length(cond.nona)
}



################   standardized   RevExp_U_Functions  ###################
theta <- c(1.649, 0.670, 2, 0.571, 0.451, 0.253, 0.135)
a.ind <- 1:2
lam.ind <- 3
sig.ind <- 4:5
gamma.ind <- 6:7
a <- theta[a.ind]
lam <- theta[lam.ind]
sig <- theta[sig.ind]
gamma <- theta[gamma.ind]
y <- sim.RevExpU.MGPD(n=200,d=2, a=a, beta=c(log(lam),0), sig=sig, gamma=gamma, MGPD = T,std=T)$Z
plot(y)


u <- min(y)-0.01

theta <- theta[1:3]
nll.powunif.1<-function(theta,y,u,a.ind,lam.ind,lamfix=FALSE, balthresh=FALSE)
{
  d<-dim(y)[2]
  a<-theta[a.ind]
  
  if(length(a)==1){a<-rep(a,d)}
  
  if(lamfix){lam<-rep(1,d)}else{
    lam<-c(theta[lam.ind],1)
  }
  
  if(balthresh){
    lam<-1/(1+a)
  }
  
  if(any(lam<0.01)||any(a<=0.01)){return(10e10)}
  
  ind<-apply(y,1,comp.gt,u=u)
  y.uc<-y[ind,]
  y.pc<-y[!ind,]
  
  L<-apply(y.uc,1,fY.powunif,a=a,lam=lam)
  nll.uc<--sum(log(L))
  
  if(sum(!ind)>0)
  {
    L2<-apply(y.pc,1,fY.powunif.cens,a=a,lam=lam,u=u)
    nll.pc<--sum(log(L2))
  }else{nll.pc<-0}
  nll<-nll.uc+nll.pc
  return(nll)
}

nll.powunif.1(theta[1:3],exp(y),u,a.ind,lam.ind)

par <- theta
opt1 <- optim(nll.powunif.1, par=theta, u=u, y=exp(y), a.ind=a.ind, lam.ind=lam.ind, lamfix=F, control=list(maxit=10000,reltol=1e-6),balthresh=FALSE)
opt1.1 <- optim(nll.powunif.1, par=opt1$par, u=u, y=exp(y), a.ind=a.ind, lam.ind=lam.ind, lamfix=F, control=list(maxit=10000,reltol=1e-6),balthresh=FALSE)
opt1.2 <- optim(nll.powunif.1, par=opt1.1$par, u=u, y=exp(y), a.ind=a.ind, lam.ind=lam.ind, lamfix=F, control=list(maxit=10000,reltol=1e-6),balthresh=FALSE)
opt1.3 <- optim(nll.powunif.1, par=opt1.2$par, u=u, y=exp(y), a.ind=a.ind, lam.ind=lam.ind, lamfix=F, control=list(maxit=10000,reltol=1e-6),balthresh=FALSE,hessian=T)

nll.powunif.1(theta,exp(y),u,a.ind,lam.ind)
nll.powunif.1(opt1.3$par,exp(y),u,a.ind,lam.ind)


N.sim <- 100
par.mat <- matrix(NA, nrow=N.sim, ncol=length(theta))
CI.width <- matrix(NA, nrow=N.sim, ncol=length(theta))
for (i in 1:N.sim){
  set.seed(i)
  y <- sim.RevExpU.MGPD(n=200,d=2, a=a, beta=c(log(lam),0), sig=sig, gamma=gamma, MGPD = T,std=T)$Z
  u <- min(y)-0.01
  par <- theta
  opt1<-optim(nll.powunif.1, par=par, u=u, y=exp(y), lamfix=F, a.ind=a.ind, lam.ind=lam.ind,
              control=list(maxit=10000,reltol=1e-6),balthresh=FALSE)
  par <- opt1$par
  opt1.1 <- optim(nll.powunif.1, par=par, u=u,y=exp(y), lamfix=F, a.ind=a.ind, lam.ind=lam.ind,
              control=list(maxit=10000,reltol=1e-6),balthresh=FALSE)
  par <- opt1.1$par
  opt1.2 <- optim(nll.powunif.1, par=par, u=u, y=exp(y), lamfix=F, a.ind=a.ind, lam.ind=lam.ind,
              control=list(maxit=10000,reltol=1e-6),balthresh=FALSE)
  par <- opt1.2$par
  opt1.3 <- optim(nll.powunif.1, par=par, u=u, y=exp(y), lamfix=F, a.ind=a.ind, lam.ind=lam.ind,
              control=list(maxit=10000,reltol=1e-6),balthresh=FALSE,hessian=TRUE)
  
  par.mat[i,] <- opt1.3$par
  CI.width[i,] <- 1.96*sqrt(diag(solve(opt1.3$hessian)))
  print(i)
}

upper <- par.mat + CI.width
lower <- par.mat - CI.width
within.ind <- rep(NA, length(theta))
for (i in 1:length(theta)){
  cond <- (theta[i]>=lower[,i])&(theta[i]<=upper[,i])
  cond.nona <- cond[!is.na(cond)]
  within.ind[i] <- sum(cond.nona) /length(cond.nona)
}
within.ind
############### profile of the lam #####################################
set.seed(56)
y <- sim.RevExpU.MGPD(n=200,d=2, a=a, beta=c(log(lam),0), sig=sig, gamma=gamma, MGPD = T,std=T)$Z
u <- min(y)-0.01
par <- theta
opt1<-optim(nll.powunif.1, par=par, u=u, y=exp(y), lamfix=F, a.ind=a.ind, lam.ind=lam.ind,
            control=list(maxit=10000,reltol=1e-6),balthresh=FALSE)
par <- opt1$par
opt1.1 <- optim(nll.powunif.1, par=par, u=u,y=exp(y), lamfix=F, a.ind=a.ind, lam.ind=lam.ind,
                control=list(maxit=10000,reltol=1e-6),balthresh=FALSE)
par <- opt1.1$par
opt1.2 <- optim(nll.powunif.1, par=par, u=u, y=exp(y), lamfix=F, a.ind=a.ind, lam.ind=lam.ind,
                control=list(maxit=10000,reltol=1e-6),balthresh=FALSE)
par <- opt1.2$par
opt1.3 <- optim(nll.powunif.1, par=par, u=u, y=exp(y), lamfix=F, a.ind=a.ind, lam.ind=lam.ind,
                control=list(maxit=10000,reltol=1e-6),balthresh=FALSE,hessian=TRUE)

opt1.3$par
1.96*sqrt(diag(solve(opt1.3$hessian)))


nll.powunif.1(theta, y=exp(y), u=u, a.ind=a.ind, lam.ind=lam.ind)
nll.powunif.1(opt1.3$par, y=exp(y), u=u, a.ind=a.ind, lam.ind=lam.ind)

x.seq <- seq(1,2.5,0.1)
y.seq <- rep(NA, length(x.seq))
for (i in 1:length(x.seq)){
  theta.tmp <- c(theta[1:2],x.seq[i])
  y.seq[i] <- nll.powunif.1(theta.tmp,y=exp(y), u=u, a.ind=a.ind, lam.ind=lam.ind)
}
plot(x.seq, y.seq,type='l')

################# mcmc samples from standardized RevExp_U_Funciton#####3

nim_nll_powunif <- nimbleRcall(function( theta=double(1), y=double(2), u=double(0), a.ind=double(1),
                                            lam.ind=double(0),  lamfix=logical(0, default = 0),balthresh=logical(0, default=0)){}, 
                                   Rfun = 'nll.powunif.1',
                                   returnType = double(0))


dpowunif <- nimbleFunction(
  run = function(x=double(2), theta=double(1), 
                 log = logical(0, default = 0)) {
    returnType(double(0))
    
    ll <- -nim_nll_powunif(theta=theta,y=x,  u=min(x)-0.01, a.ind=1:2,
                             lam.ind=3,
                             lamfix=FALSE, balthresh=FALSE)

    if (log) {
      totalProb <- ll
    }else{
      totalProb <- exp(ll)
    }
    return(totalProb)
  })

dpowunif(x=exp(y), theta=theta, log=1)


rpowunif <- nimbleFunction(
  run = function(n=integer(0), theta=double(1)) {
    returnType(double(2))
    
    totalProb <- matrix(1,nrow=n, ncol=2)
    return(totalProb)
  })

registerDistributions(list(
  dpowunif = list(
    BUGSdist = "dpowunif(theta)",
    types = c('value = double(2)', 'theta = double(1)')  )))



BivGPcode <- nimbleCode({
  
  for (i in 1:3)
    theta[i] ~ dunif(0,50)
  
  y[1:N,1:2] ~ dpowunif(theta=theta[1:3])
})

BivGPmodel <- nimbleModel(BivGPcode, constants = list(N = 200),check = FALSE)
y.exp <- exp(y)
BivGPmodel$setData(list(y = y.exp)) 

cBivGPmodel <- compileNimble(BivGPmodel)

BivGPconf <- configureMCMC(BivGPmodel,
                               enableWAIC = TRUE, time=TRUE)

BivGPMCMC <- buildMCMC(BivGPconf)
BivGPMCMC$run(10)
cBivGPMCMC <- compileNimble(BivGPMCMC, project = BivGPmodel)

t1 <- Sys.time()
results <- runMCMC(cBivGPMCMC, niter = 5000,nburnin=2500,thin=1,
                   summary = TRUE, WAIC = TRUE,setSeed = 1235)
t2 <- Sys.time()
print(t2-t1)

plot(results$samples[,'theta[3]'],type='l')
plot(density(results$samples[,'theta[1]']))
abline(v=theta[1],col='red')

plot(density(results$samples[,'theta[2]']))
abline(v=theta[2],col='red')

plot(density(results$samples[,'theta[3]']))
abline(v=theta[3],col='red')
################   standardized   RevExp_T_Functions  ###################
a.ind <- 1:2
lam.ind <- 3
sig.ind <- 4:5
gamma.ind <- 6:7
a <- theta[a.ind]
lam <- theta[lam.ind]
sig <- theta[sig.ind]
gamma <- theta[gamma.ind]
y <- sim.RevExpT.MGPD(n=200,d=2, a=a, beta=c(log(lam),0), sig=sig, gamma=gamma, MGPD = T,std=T)$Z
plot(y)


u <- min(y)-0.01

theta <- theta[1:3]

nll.powunif.Linf.1<-function(theta,y,u,a.ind,lam.ind,lamfix=FALSE, balthresh=FALSE)
{
  d<-dim(y)[2]
  a<-theta[a.ind]
  
  if(length(a)==1){a<-rep(a,d)}
  
  if(lamfix){lam<-rep(1,d)}
  else{
    lam<-c(theta[lam.ind],1)
  }
  
  if(balthresh){
    lam<-1/(1+a)
  }
  
  if(any(lam<0.01)||any(a<=0.01)){return(10e10)}
  
  ind<-apply(y,1,comp.gt,u=u)
  y.uc<-y[ind,]
  y.pc<-y[!ind,]
  
  L<-apply(y.uc,1,fY.powunif.Linf,a=a,lam=lam)
  nll.uc<--sum(log(L))
  
  if(sum(!ind)>0)
  {
    L2<-apply(y.pc,1,fY.powunif.Linf.cens,a=a,lam=lam,u=u)
    nll.pc<--sum(log(L2))
  }
  else{nll.pc<-0}
  nll<-nll.uc+nll.pc
  return(nll)
}

nll.powunif.Linf.1(theta,exp(y),u,a.ind,lam.ind)

par <- theta
opt1 <- optim(nll.powunif.Linf.1, par=theta, u=u, y=exp(y), a.ind=a.ind, lam.ind=lam.ind, lamfix=F, control=list(maxit=10000,reltol=1e-6),balthresh=FALSE)
opt1.1 <- optim(nll.powunif.Linf.1, par=opt1$par, u=u, y=exp(y), a.ind=a.ind, lam.ind=lam.ind, lamfix=F, control=list(maxit=10000,reltol=1e-6),balthresh=FALSE)
opt1.2 <- optim(nll.powunif.Linf.1, par=opt1.1$par, u=u, y=exp(y), a.ind=a.ind, lam.ind=lam.ind, lamfix=F, control=list(maxit=10000,reltol=1e-6),balthresh=FALSE)
opt1.3 <- optim(nll.powunif.Linf.1, par=opt1.2$par, u=u, y=exp(y), a.ind=a.ind, lam.ind=lam.ind, lamfix=F, control=list(maxit=10000,reltol=1e-6),balthresh=FALSE,hessian=T)

nll.powunif.Linf.1(theta,exp(y),u,a.ind,lam.ind)
nll.powunif.Linf.1(opt1.3$par,exp(y),u,a.ind,lam.ind)


N.sim <- 100
par.mat <- matrix(NA, nrow=N.sim, ncol=length(theta))
CI.width <- matrix(NA, nrow=N.sim, ncol=length(theta))
for (i in 1:N.sim){
  set.seed(i)
  y <- sim.RevExpT.MGPD(n=200,d=2, a=a, beta=c(log(lam),0), sig=sig, gamma=gamma, MGPD = T,std=T)$Z
  u <- min(y)-0.01
  par <- theta
  opt1<-optim(nll.powunif.Linf.1, par=par, u=u, y=exp(y), lamfix=F, a.ind=a.ind, lam.ind=lam.ind,
              control=list(maxit=10000,reltol=1e-6),balthresh=FALSE)
  par <- opt1$par
  opt1.1 <- optim(nll.powunif.Linf.1, par=par, u=u,y=exp(y), lamfix=F, a.ind=a.ind, lam.ind=lam.ind,
                  control=list(maxit=10000,reltol=1e-6),balthresh=FALSE)
  par <- opt1.1$par
  opt1.2 <- optim(nll.powunif.Linf.1, par=par, u=u, y=exp(y), lamfix=F, a.ind=a.ind, lam.ind=lam.ind,
                  control=list(maxit=10000,reltol=1e-6),balthresh=FALSE)
  par <- opt1.2$par
  opt1.3 <- optim(nll.powunif.Linf.1, par=par, u=u, y=exp(y), lamfix=F, a.ind=a.ind, lam.ind=lam.ind,
                  control=list(maxit=10000,reltol=1e-6),balthresh=FALSE,hessian=TRUE)
  
  par.mat[i,] <- opt1.3$par
  CI.width[i,] <- 1.96*sqrt(diag(solve(opt1.3$hessian)))
  print(i)
}

upper <- par.mat + CI.width
lower <- par.mat - CI.width
within.ind <- rep(NA, length(theta))
for (i in 1:length(theta)){
  cond <- (theta[i]>=lower[,i])&(theta[i]<=upper[,i])
  cond.nona <- cond[!is.na(cond)]
  within.ind[i] <- sum(cond.nona) /length(cond.nona)
}
within.ind


################   standardized   RevExp_U_Functions  ###################

theta <- c(1.649, 0.670, 2, 1,  0.571, 0.451, 0.253, 0.135)
a.ind <- 1:2
lam.ind <- 3:4
sig.ind <- 5:6
gamma.ind <- 7:8
a <- theta[a.ind]
lam <- theta[lam.ind]
sig <- theta[sig.ind]
gamma <- theta[gamma.ind]
y <- sim.GumbelT.MGPD(n=2000,d=2, a=a, beta=log(lam), sig=sig, gamma=gamma, MGPD = T,std=T)$Z
plot(y)


u <- min(y)-0.01


theta <- theta[1:4]

nll.Frechet.a.Linf(theta,exp(y),u)

par <- theta
opt1 <- optim(nll.Frechet.a.Linf, par=theta, u=u, y=exp(y),  control=list(maxit=10000,reltol=1e-6))
opt1.1 <- optim(nll.Frechet.a.Linf, par=opt1$par, u=u, y=exp(y), control=list(maxit=10000,reltol=1e-6))
opt1.2 <- optim(nll.Frechet.a.Linf, par=opt1.1$par, u=u, y=exp(y),  control=list(maxit=10000,reltol=1e-6))
opt1.3 <- optim(nll.Frechet.a.Linf, par=opt1.2$par, u=u, y=exp(y), control=list(maxit=10000,reltol=1e-6),hessian=T)

nll.Frechet.a.Linf(theta,exp(y),u)
nll.Frechet.a.Linf(opt1.3$par,exp(y),u)

