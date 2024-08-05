dir.work <- '/home/pgrad2/2448355h/My_PhD_Project/Biv_Ext_Mix_Mod'
dir.out <- "/home/pgrad2/2448355h/My_PhD_Project/01_Output/Biv_Ext_Mix_Mod/Simulation"
source(file.path(dir.work, "Simulation/RevExp_U_Functions.r"))
source(file.path(dir.work, "Simulation/CommonFunctions.r"))


load_install_packages <- function(packages) {
  for(package in packages){
    # If the package is not installed, install it
    if(!require(package, character.only = TRUE)) {
      install.packages(package, dependencies = TRUE,INSTALL_opts = '--no-lock')
      # Load the package after installation
      library(package, character.only = TRUE)
    } else {
      # If the package is already installed, just load it
      library(package, character.only = TRUE)
    }
  }
}

# List the packages you want to load
packages <- c("nimble", "foreach","doSNOW","parallel",'gsl','copula','extraDistr','rugarch','tidyr')  


load_install_packages(packages)

print(detectCores())


#########################################
rep.size <- 1000
set.seed(1234)

params <- c(2, 1.2, 2, 1.5, 3)
thres <- c(1.5,1.6)


total.dist.name <- c(3,3,3,2)

bulk.dist.name <- total.dist.name[1:3]
cop.name.map <- c("normal","clayton","gumbel","frank","joe","plackett")

cop <- switch(cop.name.map[bulk.dist.name[1]],
              "normal" = normalCopula(param = params[1], dim = 2),
              "clayton" = claytonCopula(param = params[1], dim = 2),
              "gumbel" = gumbelCopula(param = params[1], dim = 2),
              "frank" = frankCopula(param = params[1], dim = 2),
              "joe" = joeCopula(param = params[1], dim = 2),
              "plackett" = plackettCopula(param = params[1]),
              stop("Unsupported copula name"))

margin.name.map <- c("norm","exp","gamma","lnorm","weibull",'lst')
margins.name <- margin.name.map[bulk.dist.name[2:3]]
params.margin <- params[-1]
transformed_params <- list()
start.idx <- 1
for (i in 1:length(margins.name)) {
  m <- margins.name[i]
  
  if (m == "norm") {
    n.param <- 2
    p <- params.margin[ start.idx : (start.idx+1)]
    transformed_params[[i]] <- list(mean = p[1], sd = p[2]) # ensure sd > 0
    start.idx <- start.idx + n.param
    
  } else if (m == "exp") {
    n.param <- 1
    p <- params.margin[ start.idx]
    transformed_params[[i]] <- list(rate = p[1]) # ensure rate > 0
    start.idx <- start.idx + n.param
    
  } else if (m == "gamma") {
    n.param <- 2
    p <- params.margin[ start.idx : (start.idx+1)]
    transformed_params[[i]] <- list(shape = p[1], rate = p[2]) # ensure shape, rate > 0
    start.idx <- start.idx + n.param
    
  } else if (m == "lnorm") {
    n.param <- 2
    p <- params.margin[ start.idx : (start.idx+1)]
    transformed_params[[i]] <- list(meanlog = p[1], sdlog = p[2]) # ensure sdlog > 0
    start.idx <- start.idx + n.param
  } else if (m == "weibull") {
    n.param <- 2
    p <- params.margin[ start.idx : (start.idx+1)]
    transformed_params[[i]] <- list(shape = p[1], scale = p[2]) # ensure shape, scale > 0
    start.idx <- start.idx + n.param
  } else if (m == 'lst'){
    n.param <- 3
    p <- params.margin[ start.idx : (start.idx+2)]
    transformed_params[[i]] <- list(df = p[1], mu = p[2], sigma = p[3]) 
    start.idx <- start.idx + n.param
  } else {
    stop("Unsupported marginal distribution")
  }
}


bulk.dist <- mvdc(cop, margins = margins.name, 
                  
                  paramMargins=list(transformed_params[[1]],
                                    
                                    transformed_params[[2]])
)


Y <- rMvdc(rep.size, bulk.dist)



dim(Y)
plot(Y)
#########################################



######################################################################################
# Parallelled code for running the simulation multiple times
# The reason of dividing the code into two parts instead of one is because large number 
# iterations in the foreach cound cause some unexpected error.
######################################################################################
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
  
  margin.name.map <- c("norm","exp","gamma","lnorm","weibull",'lst')
  margins.name <- margin.name.map[bulk.dist.name[2:3]]
  params.margin <- params[-1]
  transformed_params <- list()
  start.idx <- 1
  for (i in 1:length(margins.name)) {
    m <- margins.name[i]
    
    if (m == "norm") {
      n.param <- 2
      p <- params.margin[ start.idx : (start.idx+1)]
      transformed_params[[i]] <- list(mean = p[1], sd = p[2]) # ensure sd > 0
      start.idx <- start.idx + n.param
      
    } else if (m == "exp") {
      n.param <- 1
      p <- params.margin[ start.idx]
      transformed_params[[i]] <- list(rate = p[1]) # ensure rate > 0
      start.idx <- start.idx + n.param
      
    } else if (m == "gamma") {
      n.param <- 2
      p <- params.margin[ start.idx : (start.idx+1)]
      transformed_params[[i]] <- list(shape = p[1], rate = p[2]) # ensure shape, rate > 0
      start.idx <- start.idx + n.param
      
    } else if (m == "lnorm") {
      n.param <- 2
      p <- params.margin[ start.idx : (start.idx+1)]
      transformed_params[[i]] <- list(meanlog = p[1], sdlog = p[2]) # ensure sdlog > 0
      start.idx <- start.idx + n.param
    } else if (m == "weibull") {
      n.param <- 2
      p <- params.margin[ start.idx : (start.idx+1)]
      transformed_params[[i]] <- list(shape = p[1], scale = p[2]) # ensure shape, scale > 0
      start.idx <- start.idx + n.param
    } else if (m == 'lst'){
      n.param <- 3
      p <- params.margin[ start.idx : (start.idx+2)]
      transformed_params[[i]] <- list(df = p[1], mu = p[2], sigma = p[3]) 
      start.idx <- start.idx + n.param
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

# gumbel, lst, lst
# R_dbulk(Y, c(1.5,  1, 2, 2, 1,2,3),c(3,6,6),c(-Inf,-Inf),c(2,3))



dbulk <- nimbleRcall(function(x = double(2), params = double(1),bulk.dist.name = double(1),
                              lbound.bulk = double(1), ubound.bulk = double(1)){}, 
                     Rfun = 'R_dbulk',
                     returnType = double(1))



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
  
  a<- rep(theta[1],2)
  
  lam<-rep(1,d)
  
  
  sig<-theta[2:3]
  gamma<-theta[4:5]
  
  
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

nll.Frechet.a.GPD.1<-function(theta,x)
{
  u <- min(x)-0.01
  d<-2
  a<-theta[1]
  lam<-rep(1,2)
  
  sig<-theta[2:3]
  gamma<-theta[4:5]
  
  
  # check data respect marginal constraints (i.e. no x[,j]>-sig[j]/gamma[j] if gamma[j]<0)
  chk<-rep(FALSE,d)
  for(j in 1:d)
  {
    if(gamma[j]<0){
      chk[j]<-any(x[,j]> (-sig[j]/gamma[j]))
    }
    if(gamma[j]>0){
      chk[j]<-any(x[,j]< (-sig[j]/gamma[j]))
    } 
    
    
  }
  
  if(any(lam<0.01)|any(sig<0.001)|any(a<=1)|any(chk)){return(10e7)}
  
  
  ind<-apply(x,1,comp.gt,u=u)
  x.uc<-x[ind,]
  x.pc<-x[!ind,]
  
  L<-apply(x.uc,1,fX.Fre.a,a=a,lam=lam,sig=sig,gamma=gamma)
  nll.uc<--sum(log(L))
  
  if(sum(!ind)>0)
  {
    L2<-apply(x.pc,1,fX.Fre.a.cens,a=a,lam=lam,u=u,sig=sig,gamma=gamma)
    nll.pc<--sum(log(L2))
  }else{nll.pc<-0}
  
  if (is.nan(nll.uc)|is.nan(nll.pc)|(nll.uc==-Inf)|(nll.uc==Inf)){
    return(10e7)
  }
  
  nll<-nll.uc+nll.pc
  
  return(nll)
  
}







nim_nll_revexp_GPD_MAT <- nimbleRcall(function(theta=double(1), x=double(2)){}, 
                                      Rfun = 'nll.powunif.GPD.1',
                                      returnType = double(0))

nim_nll_gumbel_GPD_MAT <- nimbleRcall(function(theta=double(1), x=double(2)){}, 
                                      Rfun = 'nll.Frechet.a.GPD.1',
                                      returnType = double(0))

cdf.revexp <- function(a,lower,upper){
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

cdf.revexp.GPD <- function(a, sig, gamma, lower, upper ){
  upp.para1 <- if (gamma[1] <0) -sig[1]/gamma[1] else Inf  
  upp.para2 <- if (gamma[2] <0) -sig[2]/gamma[2] else Inf 
  upper <- pmin(upper, c(upp.para1, upp.para2))
  
  low.para1 <- if( gamma[1]>0) -sig[1]/gamma[1] else -Inf  
  low.para2 <- if( gamma[2]>0) -sig[2]/gamma[2] else -Inf  
  lower <- pmax(lower, c(low.para1,low.para2) )
  # can't do it vectorized because the formula for dim 1 and dim 2 may differ depending on the value of gamma
  for (i in 1:length(sig)){
    if (abs(gamma[i])<.Machine$double.eps){
      lower[i] <- lower[i]/sig[i]
      upper[i] <- upper[i]/sig[i]
    }else{
      if(is.na(log(gamma[i]/sig[i]*upper[i]+1+ .Machine$double.eps*10))) print(list(gamma=gamma,sig=sig,upper=upper))
      lower[i] <- log(gamma[i]/sig[i]*lower[i]+1 + .Machine$double.eps*10)/gamma[i]
      upper[i] <- log(gamma[i]/sig[i]*upper[i]+1 + .Machine$double.eps*10)/gamma[i]
      
    }
  }
  return(cdf.revexp(a,lower,upper))
}

nim_cdf.revexp.GPD <- nimbleRcall(function(a=double(1), sig=double(1), gamma=double(1),
                                           lower=double(1), upper=double(1)){}, 
                                  Rfun = 'cdf.revexp.GPD',
                                  returnType = double(0))

cdf.gumbel <- function(a,lower,upper){
  lb1 <- lower[1]
  lb2 <- lower[2]
  ub1 <- upper[1]
  ub2 <- upper[2]
  
  C <- gamma(2-1/a)/gamma(1-1/a)/2^(1/a)/(1-1/a)
  if (exp(-a*lb1) < 10^9){
    comp1 <- 2^(1/a)-(1+exp(-a*ub2))^(1/a)-(1+exp(-a*lb1))^(1/a)+(exp(-a*ub2)+exp(-a*lb1))^(1/a)
  }else{
    comp1 <- 2^(1/a)-(1+exp(-a*ub2))^(1/a)
  }
  
  if(exp(-a*lb2) < 10^9){
    comp2 <- (exp(-a*lb2)+exp(-a*ub1))^(1/a) - (exp(-a*ub2)+exp(-a*ub1))^(1/a) - (exp(-a*lb2)+1)^(1/a) + (exp(-a*ub2)+1)^(1/a)
  }else{
    comp2 <-  - (exp(-a*ub2)+exp(-a*ub1))^(1/a)  + (exp(-a*ub2)+1)^(1/a)
  }
  
  
  return(C*(comp1 + comp2))
  
}


cdf.gumbel.GPD <- function(a, sig, gamma, lower, upper ){
  upp.para1 <- if (gamma[1] <0) -sig[1]/gamma[1] else Inf  
  upp.para2 <- if (gamma[2] <0) -sig[2]/gamma[2] else Inf 
  upper <- pmin(upper, c(upp.para1, upp.para2))
  
  low.para1 <- if( gamma[1]>0) -sig[1]/gamma[1] else -Inf  
  low.para2 <- if( gamma[2]>0) -sig[2]/gamma[2] else -Inf  
  lower <- pmax(lower, c(low.para1,low.para2) )
  # can't do it vectorized because the formula for dim 1 and dim 2 may differ depending on the value of gamma
  for (i in 1:length(sig)){
    if (abs(gamma[i])<.Machine$double.eps){
      lower[i] <- lower[i]/sig[i]
      upper[i] <- upper[i]/sig[i]
    }else{
      if(is.na(log(gamma[i]/sig[i]*upper[i]+1+ .Machine$double.eps*10))) print(list(gamma=gamma,sig=sig,upper=upper))
      lower[i] <- log(gamma[i]/sig[i]*lower[i]+1 + .Machine$double.eps*10)/gamma[i]
      upper[i] <- log(gamma[i]/sig[i]*upper[i]+1 + .Machine$double.eps*10)/gamma[i]
      
    }
  }
  return(cdf.gumbel(a,lower,upper))
}

nim_cdf.gumbel.GPD <- nimbleRcall(function(a=double(0), sig=double(1), gamma=double(1),
                                           lower=double(1), upper=double(1)){}, 
                                  Rfun = 'cdf.gumbel.GPD',
                                  returnType = double(0))





dbiextmix <- nimbleFunction(
  run = function(x=double(2), thres=double(1), params.bulk=double(1), total.dist.name=double(1),
                 theta=double(1),  
                 lower=double(1),upper=double(1),
                 log = logical(0, default = 0)) {
    returnType(double(0))
    
    bulk.dist.name <- total.dist.name[1:3]
    tail.dist.name <- total.dist.name[4]
    
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
    
    sig <- theta[2:3]
    gamma <- theta[4:5]
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
        bound.cond.1 <- y.tail[,1] > (lbound.tail-thres)[1] & y.tail[,1] < (ubound.tail-thres)[1]
        bound.cond.2 <- y.tail[,2] > (lbound.tail-thres)[2] & y.tail[,2] < (ubound.tail-thres)[2]
        y.tail <- y.tail[bound.cond.1&bound.cond.2,]
        
        if (tail.dist.name==1){
          ll.tail <- -nim_nll_revexp_GPD_MAT(x=y.tail, theta=theta)
          
          den <- nim_cdf.revexp.GPD(a=rep(theta[1],2), sig=theta[2:3], gamma=theta[4:5],
                                    lower=lbound.tail-thres,upper=ubound.tail-thres)
          
          ll.tail <- ll.tail - n.tail*log(den)
        }else{
          ## truncation to be done
          ll.tail <- -nim_nll_gumbel_GPD_MAT(x=y.tail, theta=theta)
          den <- nim_cdf.gumbel.GPD(a=theta[1], sig=theta[2:3], gamma=theta[4:5],
                                    lower=lbound.tail-thres,upper=ubound.tail-thres)
          ll.tail <- ll.tail - n.tail*log(den)
        }
        
      }else{
        ll.tail <- -10^10
      }
    }
    
    
    if (n.bulk>0){
      if (!is.nan(dbulk.res[2])){
        ll.bulk <- max(dbulk.res[2],ll.bulk)
      }
      
    }
    ll.all <- n.tail*log(1-pi) + ll.tail + ll.bulk 
    if (log) {
      totalProb <- ll.all
    }else{
      totalProb <- exp(ll.all)
    }
    
    return(totalProb)
  })




rbiextmix <- nimbleFunction(
  run = function(n=integer(0), thres=double(1), 
                 params.bulk = double(1), total.dist.name=double(1),
                 theta=double(1), 
                 lower=double(1), upper=double(1)) {
    returnType(double(2))
    
    totalProb <- matrix(1,nrow=n, ncol=2)
    return(totalProb)
  })

registerDistributions(list(
  dbiextmix = list(
    BUGSdist = "dbiextmix(thres, params.bulk, total.dist.name,
                 theta,  lower, upper)",
    types = c('value = double(2)', 'thres = double(1)', 'params.bulk = double(1)',
              'total.dist.name= double(1)', 'theta=double(1)',
              'lower = double(1)', 'upper = double(1)')
  )))



BivExtMixcode <- nimbleCode({
  
  for (i in 1:2)
    thres[i] ~ T(dnorm(mean.thres[i], sd=sd), min.thres[i],max.thres[i])
  
  # this is for joe copula!
  params.bulk[1] ~  dunif(1,50)
  
  # mu
  params.bulk[2] ~  dunif(0,20)
  # scale
  params.bulk[3] ~  dunif(0,20)
  
  params.bulk[4] ~  dunif(0,20)
  params.bulk[5] ~  dunif(0,20)
  
  
  if (total.dist.name[4]==1){
    theta[1] ~ dunif(0,20)
  }else{
    theta[1] ~ dunif(1,20)
  }
  
  # priors for sig
  for (i in 2:3)
    theta[i] ~ dunif(0,20)
  # priors for gamma 
  for (i in 4:5)
    theta[i] ~ dunif(-1,1)
  
  y[1:N,1:2] ~ dbiextmix(thres=thres[1:2], params.bulk=params.bulk[1:5], total.dist.name=total.dist.name[1:4],
                         theta=theta[1:5],  
                         lower=lbound[1:4],upper=ubound[1:4])
  
})

dat <- Y
total.dist.name <- c(3,3,3,2)
cop.name.map <- c("normal","clayton","gumbel","frank","joe","plackett")
margin.name.map <- c("norm","exp","gamma","lnorm","weibull",'lst')
tail.name.map <- c('revexp','gumbel')
print(c(cop.name.map[total.dist.name[1]], margin.name.map[total.dist.name[2]],margin.name.map[total.dist.name[3]],tail.name.map[total.dist.name[4]]))
BivExtMixmodel <- nimbleModel(BivExtMixcode, constants = list(N = nrow(dat), 
                                                              total.dist.name=total.dist.name, 
                                                              min.thres = apply(dat,2,quantile,0.01),
                                                              mean.thres = apply(dat,2,quantile,0.90),
                                                              max.thres = apply(dat,2,quantile,0.99),
                                                              sd = 1,
                                                              lbound = rep(0, 4),
                                                              ubound = rep(Inf, 4)),
                              check = FALSE)


BivExtMixmodel$setData(list(y = dat))  ## Set those values as data in the model
cBivExtMixmodel <- compileNimble(BivExtMixmodel, showCompilerOutput = TRUE)

BivExtMixconf <- configureMCMC(BivExtMixmodel,
                               enableWAIC = TRUE, time=TRUE)
BivExtMixconf$removeSamplers(c('theta[1:5]', 'thres[1:2]','params.bulk[1:5]'))
BivExtMixconf$addSampler(target = c('theta[1:5]', 'thres[1:2]','params.bulk[1:5]'), type = 'AF_slice')


BivExtMixMCMC <- buildMCMC(BivExtMixconf)
# BivExtMixMCMC$run(1)
cBivExtMixMCMC <- compileNimble(BivExtMixMCMC, project = BivExtMixmodel, showCompilerOutput = TRUE)

inits <- list(thres=apply(Y,2, quantile,0.9))

t1 <- Sys.time()
results <- runMCMC(cBivExtMixMCMC, niter = 5000, nburnin=1,thin=1,
                   summary = TRUE, WAIC = TRUE,setSeed = 1234,inits = inits)
t2 <- Sys.time()
print(t2-t1)
Y.fit <- Y

save(Y.fit ,total.dist.name, results, file=file.path(dir.out, 'Simulation_4_3332_thres_0.01_0.90_0.99_sd_1.RData') )



# 
# 
# post.burnin <- 3000:4999
# x <- Y
# theta=colMeans(results$samples[post.burnin,c('theta[1]','theta[2]','theta[3]','theta[4]'
#                                              ,'theta[5]')])
# thres=colMeans(results$samples[post.burnin,c('thres[1]','thres[2]')])
# params.bulk=colMeans(results$samples[post.burnin,c('params.bulk[1]','params.bulk[2]','params.bulk[3]',
#                                                    'params.bulk[4]','params.bulk[5]')])
# # 
# total.dist.name=c(3,3,3,2)
# lower=rep(0,4)
# upper=rep(Inf,4)
# # thres <- c(2.3,1.7)
# 
# 
# dbiextmix(x=Y, thres=thres, params.bulk=params.bulk, total.dist.name=total.dist.name,
#           theta=theta,
#           lower=rep(-Inf,4),upper=rep(Inf,4),
#           log = TRUE)
# 
