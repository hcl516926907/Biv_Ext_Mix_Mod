dir.work <- '/home/pgrad2/2448355h/My_PhD_Project/Biv_Ext_Mix_Mod'
dir.data <- '/home/pgrad2/2448355h/My_PhD_Project/00_Dataset'
source(file.path(dir.work, "KRSW/RevExp_U_Functions.r"))
source(file.path(dir.work, "KRSW/CommonFunctions.r"))

load_install_packages <- function(packages) {
  for(package in packages){
    # If the package is not installed, install it
    if(!require(package, character.only = TRUE)) {
      install.packages(package, dependencies = TRUE,repos='http://cran.us.r-project.org')
      # Load the package after installation
      library(package, character.only = TRUE)
    } else {
      # If the package is already installed, just load it
      library(package, character.only = TRUE)
    }
  }
}

packages <- c("nimble",'copula','posterior')  


load_install_packages(packages)

#######################################load the data
load(file.path(dir.data, 'index.daily.RData'))

dat <- index.daily.dropna[,c('return_dax', 'return_cac')]
row.names(dat) <- NULL

dat.tf <- 100*abs(1-dat)

#return_ftse contains 0 in dat.tf

plot(dat.tf)
#######################################Nimble code for MCMC#################
# params is in order Copula parameter, shape 1, shape 2, scale 1, scale 2
R_pbulk <- function(x, params){
  norm.cop <- ellipCopula("normal", param = params[1], dim = 2)
  gaussian <- mvdc(norm.cop, margins = c('gamma','gamma'), 
                   paramMargins=list(list(shape=params[2],scale=params[4]),
                                     list(shape=params[3],scale=params[5])))
  p <- pMvdc(x, gaussian)
  return(p)
}


R_dbulk <- function(x, params, log){
  norm.cop <- ellipCopula("normal", param = params[1], dim = 2)
  gaussian <- mvdc(norm.cop, margins = c('gamma','gamma'), 
                   paramMargins=list(list(shape=params[2],scale=params[4]),
                                     list(shape=params[3],scale=params[5])))
  dvect <- dMvdc(x, gaussian, log=log)
  if (log) {
    return(sum(dvect))
  }else{
    return(prod(dvect))
  }
}


nim_pbulk<- nimbleRcall(function(x = double(1), params = double(1)){}, 
                        Rfun = 'R_pbulk',
                        returnType = double(0))



nim_dbulk <- nimbleRcall(function(x = double(2), params = double(1),
                                  log=logical(0, default = 0)){}, 
                         Rfun = 'R_dbulk',
                         returnType = double(0))



nll.powunif.GPD.1<-function(theta,x,u,a.ind,lam.ind,sig.ind,gamma.ind, lamfix=FALSE, balthresh=FALSE, marg.scale.ind,marg.shape.ind)
{ 
  x.mat.ind <- 1
  if (is.null(dim(x))){
    d <- length(x)
    x.mat.ind <- 0
  }else{
    d<-dim(x)[2]
  }
  
  a<-theta[a.ind]
  if(length(a)==1)
  {
    a<-rep(a,d)
  }
  
  if(lamfix){
    lam<-rep(1,d)
  }else{
    lam<-c(theta[lam.ind],1)
  }
  
  if(balthresh){
    lam<-1/(1+a)
  }
  
  sig<-theta[sig.ind]
  gamma<-theta[gamma.ind]
  
  sig<-sig[marg.scale.ind]
  gamma<-gamma[marg.shape.ind]
  
  rej<-NULL
  # upper bound when xi is greater than 0
  if(x.mat.ind){
    for(j in 1:d)
    {
      rej[j]<- gamma[j]<0 && any(x[,j]>-sig[j]/gamma[j])
    }
  }else{
    for(j in 1:d)
    {
      rej[j]<- gamma[j]<0 && any(x[j]>-sig[j]/gamma[j])
    }
  }
  
  if(any(lam<0.01)||any(a<=0.01)||any(sig<=0.001)||any(rej)){return(10e7)}
  
  nll.uc <- 0
  nll.pc <- 0
  if (!x.mat.ind){
    uc <- comp.gt(x, u)
    if (uc){
      L <- fX.powunif(x=x, a=a, lam=lam, sig=sig, gamma=gamma)
      nll.uc <- -log(L)
    }else{
      L2 <- fX.powunif.cens(x=x, u=u, lam=lam, a=a, sig=sig, gamma=gamma)
      nll.pc <- -log(L2)
    }
  }else{
    uc<-apply(x,1,comp.gt,u=u)
    
    x.uc<-x[uc,]
    x.pc<-x[!uc,]
    
    L<-apply(x.uc,1,fX.powunif,a=a,lam=lam,sig=sig,gamma=gamma)
    nll.uc<--sum(log(L))
    
    if(sum(!uc)>0)
    {
      L2<-apply(x.pc,1,fX.powunif.cens,a=a,lam=lam,u=u,sig=sig,gamma=gamma)
      nll.pc<--sum(log(L2))
    }
  }
  if (is.nan(nll.uc)|is.nan(nll.pc)|(nll.uc==-Inf)|(nll.uc==Inf)){
    return(10e7)
  }
  nll<-nll.uc+nll.pc
  
  return(nll)
}





nim_nll_powunif_GPD_SG <- nimbleRcall(function(theta=double(1), x=double(1), u=double(0), a.ind=double(1),
                                               lam.ind=double(0), sig.ind=double(1), gamma.ind=double(1), 
                                               lamfix=logical(0, default = 0), balthresh=logical(0, default = 0), 
                                               marg.scale.ind=double(1), marg.shape.ind=double(1)){}, 
                                      Rfun = 'nll.powunif.GPD.1',
                                      returnType = double(0))

nim_nll_powunif_GPD_MAT <- nimbleRcall(function(theta=double(1), x=double(2), u=double(0), a.ind=double(1),
                                                lam.ind=double(0), sig.ind=double(1), gamma.ind=double(1), 
                                                lamfix=logical(0, default = 0), balthresh=logical(0, default = 0), 
                                                marg.scale.ind=double(1), marg.shape.ind=double(1)){}, 
                                       Rfun = 'nll.powunif.GPD.1',
                                       returnType = double(0))

dbiextmix <- nimbleFunction(
  run = function(x=double(2), theta=double(1), thres=double(1), params=double(1), D=integer(0, default=2),
                 a.ind=double(1), lam.ind=double(0), lamfix=logical(0, default = 0), 
                 sig.ind=double(1), gamma.ind=double(1),
                 log = logical(0, default = 0)) {
    returnType(double(0))
    
    cond <- (x[,1]>thres[1]) | (x[,2]>thres[2])
    n.tail <- sum(cond)
    n.bulk <- sum(!cond)
    
    pi <- nim_pbulk(thres, params)
    
    sig <- theta[sig.ind]
    gamma <- theta[gamma.ind]
    eta <- -sig/gamma
    eta[which(gamma<=0)] <- -Inf
    
    dtail <- 0
    dbulk <- 0
    
    if (n.tail>0){
      y.min <- eta
      y.tail <- matrix(c(x[cond,1] - thres[1], x[cond,2] - thres[2]), ncol=D)
      for (i in 1:D){
        y.min[i] <- min(y.tail[,i])
      }
      if (all(y.min>eta)){
        
        llt <- -nim_nll_powunif_GPD_MAT(x=y.tail, theta=theta, u=min(y.tail)-0.01, a.ind=a.ind,
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
      y.bulk <- x[!cond,]
      dbulk <- nim_dbulk(y.bulk, params, log = log )
      # use is.finite will have compilation issue: Error in if (a1$nDim == 0)
      if ((dbulk==-Inf)|(dbulk==Inf)) {
        dbulk <- -10^10
      }
    }
    if (log) {
      totalProb <- n.tail*log(1-pi) + dtail + dbulk
    }else{
      totalProb <- (1-pi)^n.tail *dtail*dbulk
    }
    return(totalProb)
  })

x <- as.matrix(dat.tf)
theta=c(1,2,exp(1),0.5,0.5,0.2,0.2)
thres=c(2,2)
params <- c(0.8,0.2,0.2,1,1)
a.ind=c(1,2)
lam.ind=3
sig.ind=c(4,5)
gamma.ind=c(6,7)
log=TRUE
D <- 2
lamfix=TRUE
dbiextmix(x, theta=theta,thres=thres,params=params, a.ind=a.ind,lam.ind=3,sig.ind=sig.ind,gamma.ind=gamma.ind,log=TRUE)


rbiextmix <- nimbleFunction(
  run = function(n=integer(0), theta=double(1), thres=double(1), params=double(1), D=integer(0, default=2),
                 a.ind=double(1), lam.ind=double(0), lamfix=logical(0, default = 0), 
                 sig.ind=double(1), gamma.ind=double(1)) {
    returnType(double(2))
    
    totalProb <- matrix(1,nrow=n, ncol=D)
    return(totalProb)
  })

registerDistributions(list(
  dbiextmix = list(
    BUGSdist = "dbiextmix(theta, thres, params, D, a.ind, lam.ind, lamfix, sig.ind, gamma.ind)",
    types = c('value = double(2)', 'theta = double(1)', 'thres = double(1)', 
              'params = double(1)', 'D = integer(0)', 'a.ind = double(1)', 
              'lam.ind = double(0)', 'lamfix = logical(0)', 'sig.ind = double(1)',
              'gamma.ind = double(1)')
  )))




BivExtMixcode <- nimbleCode({
  params[1] ~ dunif(-1,1)
  for (i in 2:5)
    params[i] ~ dunif(0,50)
  
  for (i in 1:D)
    thres[i] ~ T(dnorm(mu.thres[i], sd=sd.thres[i]),min.thres[i],max.thres[i])
  
  
  # priors for sig
  for (i in 1:5)
    theta[i] ~ dunif(0,50)
  # priors for gamma 
  for (i in 6:7)
    theta[i] ~ dunif(-1,1)
  y[1:N,1:D] ~ dbiextmix(theta=theta[1:7], thres=thres[1:D], params= params[1:5],
                         D=D, a.ind=a.ind[1:D], lam.ind=lam.ind, lamfix=lamfix, 
                         sig.ind=sig.ind[1:D], gamma.ind=gamma.ind[1:D])
  
  
})

BivExtMixmodel <- nimbleModel(BivExtMixcode, constants = list(N = nrow(dat.tf), 
                                                              D = 2,
                                                              mu.thres = apply(dat.tf,2,quantile,0.9),#90 quantile
                                                              sd.thres = c(10,10), 
                                                              min.thres = apply(dat.tf,2,quantile,0.8),
                                                              max.thres = apply(dat.tf,2,quantile,0.99),
                                                              a.ind = 1:2,
                                                              lam.ind = 3,
                                                              sig.ind = c(4,5),
                                                              gamma.ind = c(6,7),
                                                              lamfix=TRUE),check = FALSE)


BivExtMixmodel$setData(list(y = as.matrix(dat.tf)))  ## Set those values as data in the model
cBivExtMixmodel <- compileNimble(BivExtMixmodel)

BivExtMixconf <- configureMCMC(BivExtMixmodel,
                               enableWAIC = TRUE, time=TRUE)
# BivExtMixconf$removeSamplers(c('theta[1:7]', 'thres[1:2]'))
# BivExtMixconf$addSampler(target = c('theta[1:7]', 'thres[1:2]'), type = 'AF_slice')

BivExtMixMCMC <- buildMCMC(BivExtMixconf)
# BivExtMixMCMC$run(10)
cBivExtMixMCMC <- compileNimble(BivExtMixMCMC, project = BivExtMixmodel)

t1 <- Sys.time()
results <- runMCMC(cBivExtMixMCMC, niter = 30000, nburnin=0,thin=10,
                   summary = TRUE, WAIC = TRUE,setSeed = 1234)
t2 <- Sys.time()
print(t2-t1)

dir.out <- '/home/pgrad2/2448355h/My_PhD_Project/01_Output/Biv_Ext_Mix_Mod/financial_application'
# save(results, file=file.path(dir.out, 'dax_cac_gaussian_gamma_gamma_thres_0.8_0.99.RData'))

load(file=file.path(dir.out, 'dax_cac_gaussian_gamma_gamma_thres_0.8_0.99.RData'))
plot(results$samples[,'thres[2]'],type='l')
