# Functions common to most models
##################################

# comp.gt
#########
# logical function to see if all elements of x are > u
# used via "apply" on matrices of data to determine which rows correspont to completely uncensored observations

comp.gt<-function(x,u){all(x>u)}

# BCi
#####

# Inverse Box--Cox transformation, using limit form if gamma close to 0

BCi<-function(x,gamma,sig)
{
  y<-NULL
  d<-length(x)
  for(j in 1:d)
  {  
    if(abs(gamma[j])>1e-6)
    {
      y[j]<-(1+gamma[j]*x[j]/sig[j])^(1/gamma[j])
    }
    else{
      y[j]<-exp(x[j]/sig[j])
    }
  }
  return(y)
}

# Jac
#####

# Jacobian for transformation from multivartiate Pareto to MGP, using limit for gamma if near 0

Jac<-function(x,gamma,sig)
{
  J<-NULL
  d<-length(x)
  for(j in 1:d)
  {  
    if(abs(gamma[j])>1e-6)
    {
      J[j]<-(1/sig[j])*(1+gamma[j]*x[j]/sig[j])^(1/gamma[j]-1)
    }
    else{
      J[j]<-(1/sig[j])*exp(x[j]/sig[j])
    }
  }
  return(prod(J))
}

##########################################################################################
# MCMC related functions
##########################################################################################

# Bulk density
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


# negative log likelihood of mGPD with reverse exponential generator
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

# negative log likelihood of mGPD with Gumbel generator
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




# CDF of mGPD on standardized scale with reverse exponential generator
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

# CDF of mGPD on GPD scale with reverse exponential generator
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


# CDF of mGPD on standardized scale with Gumbel generator
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

# CDF of mGPD on GPD scale with Gumbel generator
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

##########################################################################################
# Marginal densities related functions
##########################################################################################


#Joint denisty of mGDP with reverse exponential generator.
jointY1.revexp <- function(x,y1,a, sig, gamma){
  if (length(x)==1){
    y <- c(y1,x)
    res <- fX.powunif(y, lam=c(1,1), a=a, sig=sig, gamma=gamma)
    if (is.na(res)) res <- 0
  }else{
    y <- cbind(rep(y1,length(x)),x)
    res <- apply(y, MARGIN = 1, fX.powunif, lam=c(1,1),a=a, sig=sig, gamma=gamma)
    res[is.na(res)] <- 0
  }
  return(res)
}

#Marginal denisty 1 of mGDP with reverse exponential generator. (By numerical integration)
marginY1.revexp <- function(y1, y2.low=-Inf, y2.upp=Inf, a, sig, gamma){
  upp.para <- if (gamma[2] <0) -sig[2]/gamma[2] else Inf  
  upper <- min(y2.upp, upp.para)
  
  low.para <- if( gamma[2]>0) -sig[2]/gamma[2] else -Inf  
  lower <- max(y2.low, low.para )
  
  if(y1 < 0){
    res <- integrate(jointY1.revexp, lower = 0, upper = upper, y1=y1, a=a,sig=sig, gamma=gamma)$value
  }else{
    res <- integrate(jointY1.revexp, lower = lower, upper = upper, y1=y1, a=a,sig=sig, gamma=gamma)$value
  }
  
  return(res)
}


#Joint denisty of mGDP with reverse exponential generator.
jointY2.revexp <- function(x,y2,a, sig, gamma){
  if (length(x)==1){
    y <- c(x,y2)
    res <- fX.powunif(y, lam=c(1,1), a=a, sig=sig, gamma=gamma)
    if (is.na(res)) res <- 0
  }else{
    y <- cbind(x, rep(y2,length(x)))
    res <- apply(y, MARGIN = 1, fX.powunif, lam=c(1,1),a=a, sig=sig, gamma=gamma)
    res[is.na(res)] <- 0
  }
  return(res)
}

#Marginal denisty 2 of mGDP with reverse exponential generator.(By numerical integration)
marginY2.revexp <- function(y2, y1.low=-Inf, y1.upp=Inf, a, sig, gamma){
  upp.para <- if (gamma[1] <0) -sig[1]/gamma[1] else Inf  
  upper <- min(y1.upp, upp.para)
  
  low.para <- if( gamma[1]>0) -sig[1]/gamma[1] else -Inf  
  lower <- max(y1.low, low.para )
  
  if(y2 < 0){
    res <- integrate(jointY2.revexp, lower = 0, upper = upper, y2=y2, a=a,sig=sig, gamma=gamma)$value
  }else{
    res <- integrate(jointY2.revexp, lower = lower, upper = upper, y2=y2, a=a,sig=sig, gamma=gamma)$value
  }
  
  return(res)
}



#Joint denisty of mGDP with Gumbel generator.
jointY1.gumbel <- function(x,y1,a, sig, gamma){
  if (length(x)==1){
    y <- c(y1,x)
    res <- fX.Fre.a(y, lam=c(1,1), a=a, sig=sig, gamma=gamma)
    if (is.na(res)) res <- 0
  }else{
    y <- cbind(rep(y1,length(x)),x)
    res <- apply(y, MARGIN = 1, fX.Fre.a, lam=c(1,1),a=a, sig=sig, gamma=gamma)
    res[is.na(res)] <- 0
  }
  return(res)
}

#Marginal denisty 1 of mGDP with Gumbel generator. (By numerical integration)
marginY1.gumbel <- function(y1, y2.low=-Inf, y2.upp=Inf, a, sig, gamma){
  upp.para <- if (gamma[2] <0) -sig[2]/gamma[2] else Inf  
  upper <- min(y2.upp, upp.para)
  
  low.para <- if( gamma[2]>0) -sig[2]/gamma[2] else -Inf  
  lower <- max(y2.low, low.para )
  
  if(y1 < 0){
    res <- integrate(jointY1.gumbel, lower = 0, upper = upper, y1=y1, a=a,sig=sig, gamma=gamma)$value
  }else{
    res <- integrate(jointY1.gumbel, lower = lower, upper = upper, y1=y1, a=a,sig=sig, gamma=gamma)$value
  }
  
  return(res)
}




# standarized density (Y1<0) of mGPD with reverse exponential generator. (mathematically formula)
fY1_revexp_lt0.margin <- function(y,a,ub){
  EM <- EM.pu(a=a,lam=c(1,1))
  b <- 1/a
  den <- EM*(1+sum(b))
  num <- prod(b)*exp(b[1]*y)*(1-exp(-(1+b[1])*ub))/(1+b[1])
  return(num/den)
}

# standarized density (Y1>0) of mGPD with reverse exponential generator. (mathematically formula)
fY1_revexp_gt0.margin <- function(y,a,lb, ub){
  EM <- EM.pu(a=a,lam=c(1,1))
  b <- 1/a
  num1 <- b[1]*(exp(-y)-exp(b[2]*lb-(1+b[2])*y))
  den1 <- EM*(1+sum(b))
  
  num2 <- prod(b)*exp(-y)-prod(b)*exp(b[1]*y-(1+b[1])*ub)
  den2 <- EM*(1+sum(b))*(1+b[1])
  return(num1/den1+num2/den2)
}


# marginal density of Y1 of mGPD with reverse exponential generator. (mathematically formula)
tail_revexp.margin1 <- function(y, a,sig, gamma, tail.low, tail.upp){
  upp.para1 <- if (gamma[1] <0) -sig[1]/gamma[1] else Inf  
  upp.para2 <- if (gamma[2] <0) -sig[2]/gamma[2] else Inf 
  upper <- pmin(tail.upp, c(upp.para1, upp.para2))
  low.para1 <- if( gamma[1]>0) -sig[1]/gamma[1] else -Inf  
  low.para2 <- if( gamma[2]>0) -sig[2]/gamma[2] else -Inf  
  lower <- pmax(tail.low, c(low.para1,low.para2) )
  
  dtail <- 0
  if ( y> lower[1] & y< upper[1]){
    
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
    
    if (gamma[1]==0){
      z <- y/sig[1]
      jacob <- 1/sig[1]
    }else{
      z <- log(gamma[1]*y/sig[1]+1)/gamma[1]
      jacob <- 1/(y*gamma[1]+sig[1])
    }
    
    if (z<0){
      dtail <- fY1_revexp_lt0.margin(z,a=a,ub=upper[2])*abs(jacob)
    }else{
      dtail <- fY1_revexp_gt0.margin(z,a=a,lb=lower[2], ub=upper[2])*abs(jacob)
    }
  }
  return(dtail)
}




# standarized density (Y2<0) of mGPD with reverse exponential generator. (mathematically formula)
fY2_revexp_lt0.margin <- function(y,a,ub){
  EM <- EM.pu(a=a,lam=c(1,1))
  b <- 1/a
  den <- EM*(1+sum(b))
  num <- prod(b)*exp(b[2]*y)*(1-exp(-(1+b[2])*ub))/(1+b[2])
  return(num/den)
}

# standarized density (Y2>0) of mGPD with reverse exponential generator. (mathematically formula)
fY2_revexp_gt0.margin <- function(y, a, lb, ub){
  EM <- EM.pu(a=a,lam=c(1,1))
  b <- 1/a
  num1 <- b[2]*(exp(-y)-exp(b[1]*lb-(1+b[1])*y))
  den1 <- EM*(1+sum(b))
  
  num2 <- prod(b)*exp(-y)-prod(b)*exp(b[2]*y-(1+b[2])*ub)
  den2 <- EM*(1+sum(b))*(1+b[2])
  
  return(num1/den1+num2/den2)
}



# standarized density (Y1<0) of mGPD with gumbel generator. (mathematically formula)
fY1_gumbel_lt0.margin <- function(y,a,ub){
  C <- gamma(2-1/a)/gamma(1-1/a)/2^(1/a)/(1-1/a)
  
  comp1 <- exp(-a*y)/(exp(-a*ub) + exp(-a*y))^(1-1/a)
  comp2 <- exp(-a*y)/(1 + exp(-a*y))^(1-1/a)
  return(C*(comp1-comp2))
}

# standarized density (Y1>0) of mGPD with gumbel generator. (mathematically formula)
fY1_gumbel_gt0.margin <- function(y,a,lb, ub){
  C <- gamma(2-1/a)/gamma(1-1/a)/2^(1/a)/(1-1/a)
  
  comp1 <- exp(-a*y)/(exp(-a*ub) + exp(-a*y))^(1-1/a)
  comp2 <- exp(-a*y)/(exp(-a*lb) + exp(-a*y))^(1-1/a)
  return(C*(comp1-comp2))
}

# marginal density of Y1 of mGPD with gumbel generator. (mathematically formula)
tail_gumbel.margin1 <- function(y, a,sig, gamma, tail.low, tail.upp){
  upp.para1 <- if (gamma[1] <0) -sig[1]/gamma[1] else Inf  
  upp.para2 <- if (gamma[2] <0) -sig[2]/gamma[2] else Inf 
  upper <- pmin(tail.upp, c(upp.para1, upp.para2))
  low.para1 <- if( gamma[1]>0) -sig[1]/gamma[1] else -Inf  
  low.para2 <- if( gamma[2]>0) -sig[2]/gamma[2] else -Inf  
  lower <- pmax(tail.low, c(low.para1,low.para2) )
  
  dtail <- 0
  if ( y> lower[1] & y< upper[1]){
    
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
    
    if (gamma[1]==0){
      z <- y/sig[1]
      jacob <- 1/sig[1]
    }else{
      z <- log(gamma[1]*y/sig[1]+1)/gamma[1]
      jacob <- 1/(y*gamma[1]+sig[1])
    }
    
    if (z<0){
      dtail <- fY1_gumbel_lt0.margin(z,a=a,ub=upper[2])*abs(jacob)
    }else{
      dtail <- fY1_gumbel_gt0.margin(z,a=a,lb=lower[2], ub=upper[2])*abs(jacob)
    }
  }
  return(dtail)
}


