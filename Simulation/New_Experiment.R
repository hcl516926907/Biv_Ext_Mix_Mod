dir.work <- '/home/pgrad2/2448355h/My_PhD_Project/Biv_Ext_Mix_Mod'
dir.out <- "/home/pgrad2/2448355h/My_PhD_Project/01_Output/Biv_Ext_Mix_Mod/Simulation"
source(file.path(dir.work, "Simulation/RevExp_U_Functions.r"))
source(file.path(dir.work, "Simulation/CommonFunctions.r"))


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

# List the packages you want to load
packages <- c("nimble", "foreach","doSNOW","parallel",'copula')  


load_install_packages(packages)

print(detectCores())

######################################################################################
# Parallelled code for running the simulation multiple times
# The reason of dividing the code into two parts instead of one is because large number 
# iterations in the foreach cound cause some unexpected error.
######################################################################################
i <- 1
seed <- i
d <- 2
a <- c(2, 2.5)
beta <- c(0, 0)
sig <- c(1, 0.8)
gamma <- c(-0.1, -0.2)

n <- 2000

par <- c(1.5,  7, 2, 2, 4)

u.x <- c(6, 6.5)

lbound.tail=c(0,0)
ubound.tail=c(Inf,Inf)

gum.cop <- gumbelCopula(param = 1.5, dim = 2)

gum.dist <- mvdc(gum.cop, margins = c('gamma','weibull'), 
                 
                 paramMargins=list(list(shape=par[2], rate=par[3]),
                                   
                                   list(shape=par[4],scale=par[5])))

p <- pMvdc(u.x, gum.dist)


set.seed(1111)
Y.tail<-sim.RevExpU.MGPD(n=n-floor(n*p),d=d, a=a, beta=beta, sig=sig, gamma=gamma, MGPD = T,std=T)$X

# GP scale tail data combined with the bulk data

Y.bulk <- rMvdc(floor(n*p), gum.dist)
Y.bulk <- Y.bulk[Y.bulk[,1]<u.x[1] & Y.bulk[,2]<u.x[2],]

Y <- rbind(Y.bulk, sweep(Y.tail,2,u.x,"+"))
apply(Y.bulk,2,min) + u.x
plot(Y)
bound.cond.1 <- Y.tail[,1] > (lbound.tail-u.x)[1] & Y.tail[,1] < (ubound.tail-u.x)[1]
bound.cond.2 <- Y.tail[,2] > (lbound.tail-u.x)[2] & Y.tail[,2] < (ubound.tail-u.x)[2]
Y.tail <- Y.tail[bound.cond.1&bound.cond.2,]
Y <- rbind(Y.bulk, sweep(Y.tail,2,u.x,"+"))
plot(Y)

sum(dMvdc(Y.bulk, gum.dist,log=T))


R_dbulk <- function(x,params, bulk.dist.name, lbound.bulk, ubound.bulk ){
  cop.name.map <- c("normal","clayton","gumbel","frank","joe","plackett")
  cop <- switch(cop.name.map[bulk.dist.name[1]],
                "normal" = normalCopula(param = params[1], dim = 2),
                "clayton" = claytonCopula(param = params[1], dim = 2),
                "gumbel" = gumbelCopula(param = params[1], dim = 2),
                "frank" = frankCopula(param = params[1], dim = 2),
                "joe" = joeCopula(param = params[1], dim = 2),
                "plackett" = plackettCopula(param = params[1]),
                stop("Unsupported copula name"))
  
  margin.name.map <- c("norm","exp","gamma","lnorm","weibull")
  margins.name <- margin.name.map[bulk.dist.name[2:3]]
  params.margin <- params[2:5]
  transformed_params <- list()
  for (i in 1:length(margins.name)) {
    m <- margins.name[i]
    p <- params.margin[(2*i-1) : (2*i)]
    
    if (m == "norm") {
      transformed_params[[i]] <- list(mean = p[1], sd = p[2]) # ensure sd > 0
    } else if (m == "exp") {
      transformed_params[[i]] <- list(rate = p[1]) # ensure rate > 0
    } else if (m == "gamma") {
      transformed_params[[i]] <- list(shape = p[1], rate = p[2]) # ensure shape, rate > 0
    } else if (m == "lnorm") {
      transformed_params[[i]] <- list(meanlog = p[1], sdlog = p[2]) # ensure sdlog > 0
    } else if (m == "weibull") {
      transformed_params[[i]] <- list(shape = p[1], scale = p[2]) # ensure shape, scale > 0
    } else {
      stop("Unsupported marginal distribution")
    }
  }

  
  joint.dist <- mvdc(cop, margins = margins.name, 
                     
                     paramMargins=list(transformed_params[[1]],
                                       
                                       transformed_params[[2]])
  )
  
  p1 <- pMvdc(ubound.bulk, joint.dist)
  p2 <- pMvdc(lbound.bulk, joint.dist)
  p <- p1 - p2 
  ll.bulk <- sum(dMvdc(x, joint.dist,log=T))
  return(c(p, ll.bulk))

}

# gumbel, gamma, weibull
R_dbulk(Y.bulk,c(1.5,  7, 2, 2, 4),c(3,3,5),c(0,0),u.x)



dbulk <- nimbleRcall(function(x = double(2), params = double(1),bulk.dist.name = double(1),
                              lbound.bulk = double(1), ubound.bulk = double(1)){}, 
                            Rfun = 'R_dbulk',
                            returnType = double(1))

dbulk(Y.bulk,c(1.5,  7, 2, 2, 4),c(3,3,5),c(0,0),u.x)


nll.powunif.GPD.1<-function(theta,x)
{ 
  
  u <- min(x)-0.01
  x.mat.ind <- 1
  if (is.null(dim(x))){
    d <- length(x)
    x.mat.ind <- 0
  }else{
    d<-dim(x)[2]
  }
  
  a<-theta[1:2]
  
  lam<-rep(1,d)
  
  
  sig<-theta[3:4]
  gamma<-theta[5:6]
  
  
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





nim_nll_powunif_GPD_SG <- nimbleRcall(function(theta=double(1), x=double(1)){}, 
                                      Rfun = 'nll.powunif.GPD.1',
                                      returnType = double(0))

nim_nll_powunif_GPD_MAT <- nimbleRcall(function(theta=double(1), x=double(2)){}, 
                                       Rfun = 'nll.powunif.GPD.1',
                                       returnType = double(0))

trunc.norm.const <- function(a,lower,upper){
  lb1 <- lower[1]
  lb2 <- lower[2]
  ub1 <- upper[1]
  ub2 <- upper[2]
  EM <- EM.pu(a=a,lam=c(1,1))
  b <- 1/a
  
  den1 <- EM*(1+sum(b))*(1+b[2])
  num1 <- b[1]*(1-exp(-(1+b[2])*ub1))*(1-exp(b[2]*lb2))
  
  den2 <- EM*(1+sum(b))*(1+b[1])
  num2 <- (1+b[1])*b[2]*(1-exp(-ub2))-b[2]*exp(b[1]*lb1)*(1-exp(-(1+b[1])*ub2))
  
  den3 <- EM*(1+sum(b))*(1+b[2])
  num3 <- prod(b)*(1-exp(-ub2))-b[1]*(exp(-ub1)-exp(-(1+b[2])*ub1))
  
  return(num1/den1+num2/den2 + num3/den3)
  
}

trunc.norm.const.GPD <- function(a, sig, gamma, lower, upper ){
  upp.para1 <- if (gamma[1] <0) -sig[1]/gamma[1] else Inf  
  upp.para2 <- if (gamma[2] <0) -sig[2]/gamma[2] else Inf 
  upper <- pmin(upper, c(upp.para1, upp.para2))
  
  low.para1 <- if( gamma[1]>0) -sig[1]/gamma[1] else -Inf  
  low.para2 <- if( gamma[2]>0) -sig[2]/gamma[2] else -Inf  
  lower <- pmax(lower, c(low.para1,low.para2) )
  # can't do it vectorized because the formula for dim 1 and dim 2 may differ depending on the value of gamma
  for (i in 1:length(sig)){
    if (abs(gamma[i])<10^-6){
      lower[i] <- lower[i]/sig[i]
      upper[i] <- upper[i]/sig[i]
    }else{
      if(is.na(log(gamma[i]/sig[i]*upper[i]+1+ 10^-6))) print(list(gamma=gamma,sig=sig,upper=upper))
      lower[i] <- log(gamma[i]/sig[i]*lower[i]+1 + 10^-6)/gamma[i]
      upper[i] <- log(gamma[i]/sig[i]*upper[i]+1 + 10^-6)/gamma[i]
      
    }
  }
  return(trunc.norm.const(a,lower,upper))
}

nim_trunc.norm.const.GPD <- nimbleRcall(function(a=double(1), sig=double(1), gamma=double(1),
                                                 lower=double(1), upper=double(1)){}, 
                                       Rfun = 'trunc.norm.const.GPD',
                                       returnType = double(0))

dbiextmix <- nimbleFunction(
  run = function(x=double(2), thres=double(1), params.bulk=double(1), bulk.dist.name=double(1),
                 theta=double(1),  
                 lower=double(1),upper=double(1),
                 log = logical(0, default = 0)) {
    returnType(double(0))
    
    cond <- (x[,1]>thres[1]) | (x[,2]>thres[2])
    n.tail <- sum(cond)
    n.bulk <- sum(!cond)
    y.tail <- matrix(c(x[cond,1] - thres[1], x[cond,2] - thres[2]), ncol=2)
    y.bulk <- x[!cond,]
    
    lbound.bulk <- lower[1:2]
    lbound.tail <- lower[3:4]
    
    ubound.bulk <- upper[1:2]
    ubound.tail <- upper[3:4]
    
    dbulk.res <- dbulk(y.bulk, params.bulk, bulk.dist.name ,lbound.bulk,thres)
    pi <- dbulk.res[1]
    
    sig <- theta[3:4]
    gamma <- theta[5:6]
    eta <- -sig/gamma
    eta[which(gamma<=0)] <- -Inf
    
    ll.tail <- -10^10
    ll.bulk <- -10^10
    
    if (n.tail>1){
      y.min <- eta
      for (i in 1:2){
        y.min[i] <- min(y.tail[,i])
      }
      if (all(y.min>eta)){
        ll.tail <- -nim_nll_powunif_GPD_MAT(x=y.tail, theta=theta)
        
        den <- nim_trunc.norm.const.GPD(a=theta[1:2], sig=theta[3:4], gamma=theta[5:6],
                                    lower=lbound.tail-thres,upper=ubound.tail-thres)
        
        ll.tail <- ll.tail - n.tail*log(den)
      }else{
        ll.tail <- -10^10
      }
    }

    
    if (n.bulk>0){
      ll.bulk <- max(dbulk.res[2],ll.bulk)
    }
    ll.all <- n.tail*log(1-pi) + ll.tail + ll.bulk 
    if (log) {
      totalProb <- ll.all
    }else{
      totalProb <- exp(ll.all)
    }
    
    return(totalProb)
  })

# x=matrix(c(1,6,6,-4,6,-3),ncol=2)
x <- Y
theta=c(a,sig,gamma)
thres=u.x
params.bulk=par
bulk.dist.name=c(3,3,5)
lower=rep(0,4)
upper=rep(Inf,4)
log = TRUE
dbiextmix(x=Y, thres=thres, params.bulk=par, bulk.dist.name=c(3,3,5),
                 theta=theta,  
                 lower=rep(0,4),upper=rep(Inf,4),
                 log = TRUE)




rbiextmix <- nimbleFunction(
  run = function(n=integer(0), thres=double(1), 
                 params.bulk = double(1), bulk.dist.name=double(1),
                 theta=double(1), 
                 lower=double(1), upper=double(1)) {
    returnType(double(2))
    
    totalProb <- matrix(1,nrow=n, ncol=2)
    return(totalProb)
  })

registerDistributions(list(
  dbiextmix = list(
    BUGSdist = "dbiextmix(thres, params.bulk, bulk.dist.name,
                 theta,  lower, upper)",
    types = c('value = double(2)', 'thres = double(1)', 'params.bulk = double(1)',
              'bulk.dist.name= double(1)', 'theta=double(1)',
              'lower = double(1)', 'upper = double(1)')
  )))



BivExtMixcode <- nimbleCode({
  
  for (i in 1:2)
    thres[i] ~ T(dnorm(0, sd=20),min.thres[i],max.thres[i])

  params.bulk[1] ~  dunif(1,20)
  params.bulk[2] ~  dgamma(shape=1, scale=5)
  params.bulk[3] ~  dgamma(shape=1, scale=5)
  params.bulk[4] ~  dgamma(shape=1, scale=5)
  params.bulk[5] ~  dgamma(shape=1, scale=5)
  
    # priors for sig
  for (i in 1:4)
    theta[i] ~ dgamma(shape=1, scale=5)
  # priors for gamma 
  for (i in 5:6)
    theta[i] ~ dunif(-1,1)

  y[1:N,1:2] ~ dbiextmix(thres=thres[1:2], params.bulk=params.bulk[1:5], bulk.dist.name=bulk.dist.name[1:3],
                         theta=theta[1:6],  
                         lower=lbound[1:4],upper=ubound[1:4])
  
})

dat <- Y
BivExtMixmodel <- nimbleModel(BivExtMixcode, constants = list(N = nrow(dat), 
                                                              bulk.dist.name=c(3,3,5), 
                                                              min.thres = apply(dat,2,quantile,0.8),
                                                              max.thres = apply(dat,2,quantile,0.99),
                                                              lbound = rep(0, 4),
                                                              ubound = rep(Inf, 4)),
                                                              check = FALSE)


BivExtMixmodel$setData(list(y = dat))  ## Set those values as data in the model
cBivExtMixmodel <- compileNimble(BivExtMixmodel, showCompilerOutput = TRUE)

BivExtMixconf <- configureMCMC(BivExtMixmodel,
                               enableWAIC = TRUE, time=TRUE)
BivExtMixconf$removeSamplers(c('theta[1:6]', 'thres[1:2]','params.bulk[1:5]'))
BivExtMixconf$addSampler(target = c('theta[1:6]', 'thres[1:2]','params.bulk[1:5]'), type = 'AF_slice')

BivExtMixMCMC <- buildMCMC(BivExtMixconf)
# BivExtMixMCMC$run(1)
cBivExtMixMCMC <- compileNimble(BivExtMixMCMC, project = BivExtMixmodel, showCompilerOutput = TRUE)

t1 <- Sys.time()
results <- runMCMC(cBivExtMixMCMC, niter = 10000, nburnin=1,thin=10,
                   summary = TRUE, WAIC = TRUE,setSeed = seed)
t2 <- Sys.time() 
print(t2-t1)
save(results, file=file.path(dir.out, 'copula_tGPD_res.RData') )

load(file=file.path(dir.out, 'copula_tGPD_res.RData') )
plot(results$samples[1:999,'params.bulk[1]'],type='l')

# x <- Y
# theta=c(a,sig,gamma)
# thres=u.x
# params.bulk=par
# bulk.dist.name=c(3,3,5)
# lower=rep(0,4)
# upper=rep(Inf,4)
# dbiextmix(x=Y, thres=c(5.5,6.5), params.bulk=par, bulk.dist.name=c(3,3,5),
#           theta=theta,  
#           lower=rep(0,4),upper=rep(Inf,4),
#           log = TRUE)
# 
# pMvdc(c(7.7,7.6), gum.dist)
# 
# estimated.mean <- as.vector(colMeans(results$samples))
# dbiextmix(x=Y, thres=estimated.mean[12:13], params.bulk=estimated.mean[1:5], bulk.dist.name=c(3,3,5),
#           theta=estimated.mean[6:11],  
#           lower=rep(0,4),upper=rep(Inf,4),
#           log = TRUE)
# 
# R_dbulk(Y,estimated.mean[1:5], bulk.dist.name=c(3,3,5),lbound=c(0,0),ubound=u.x)
