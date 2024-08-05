dir.work <- '/home/pgrad2/2448355h/My_PhD_Project/Biv_Ext_Mix_Mod'
dir.out <- "/home/pgrad2/2448355h/My_PhD_Project/01_Output/Biv_Ext_Mix_Mod/Simulation"
source(file.path(dir.work, "Simulation/RevExp_U_Functions.r"))
source(file.path(dir.work, "Simulation/Gumbel_U_Functions.r"))
source(file.path(dir.work, "Simulation/CommonFunctions.r"))

install.packages("gsl", dependencies = TRUE, INSTALL_opts = '--no-lock')
# install.packages('copula',dependencies = TRUE, INSTALL_opts = '--no-lock')
library(copula)
library(extraDistr)


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



trunc.prob <- function(lower=rep(-Inf,2),upper=rep(Inf,2),a,sig,gamma){
  upp.para1 <- if (gamma[1] <0) -sig[1]/gamma[1] else Inf  
  upp.para2 <- if (gamma[2] <0) -sig[2]/gamma[2] else Inf 
  upper <- pmin(upper, c(upp.para1, upp.para2))
  
  low.para1 <- if( gamma[1]>0) -sig[1]/gamma[1] else -Inf  
  low.para2 <- if( gamma[2]>0) -sig[2]/gamma[2] else -Inf  
  lower <- pmax(lower, c(low.para1,low.para2) )
  
  double.int.Y1 <- function(x, y2.low=lower[2], y2.upp=upper[2], a=a, sig=sig, gamma=gamma){
    sapply(x, marginY1.revexp, y2.low = y2.low, y2.upp=y2.upp, a=a, sig=sig, gamma=gamma)
  }
  
  p <- integrate(double.int.Y1, lower =lower[1], upper = upper[1] , a=a,sig=sig, gamma=gamma)$value
  return(p)
}



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

trunc.norm.const(a=c(4,5), c(-0.1,-0.1), c(Inf,Inf))
trunc.prob(lower=c(-0.1,-0.1),a=c(4,5),sig=c(1,1),gamma=c(0,0))


trunc.norm.const.GPD <- function(a, sig, gamma, lower, upper ){
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
      if(is.na(log(gamma[i]/sig[i]*upper[i]+1 + .Machine$double.eps*10))) print(list(gamma=gamma,sig=sig,upper=upper))
      lower[i] <- log(gamma[i]/sig[i]*lower[i]+1+.Machine$double.eps*10 )/gamma[i]
      upper[i] <- log(gamma[i]/sig[i]*upper[i]+1+.Machine$double.eps*10 )/gamma[i]
      
    }
  }
  return(trunc.norm.const(a,lower,upper))
}

# 
trunc.norm.const.GPD(a=c(4,1), sig=c(0.1,2), gamma=c(0.2,0.2),lower=c(-0.2,-0.2), upper=c(4,5))
trunc.prob(lower=c(-0.2,-0.2),upper=c(4,5),a=c(4,1),sig=c(0.1,2),gamma=c(0.2,0.2))

trunc.norm.const.GPD(a=c(23,23), sig=c(0.6,0.45),gamma=c(0.3,0),lower=c(-0.5,-0.5), upper=c(Inf,Inf))
trunc.prob(lower=c(-0.5,-0.5),upper=c(Inf,Inf),a=c(23,23),sig=c(0.6,0.45),gamma=c(0.3,0))


nll.powunif.GPD.1<-function(theta,x,u)
{ 
  d <- 2
  a<- rep(theta[1],2)
  
  lam<-rep(1,d)
  
  
  sig<-theta[2:3]
  gamma<-theta[4:5]
  
  rej<-NULL
  # upper bound when xi is greater than 0

  for(j in 1:d)
  {
    rej[j]<- gamma[j]<0 && any(x[j]>-sig[j]/gamma[j])
  }

  
  if(any(lam<0.01)||any(a<=0.01)||any(sig<=0.001)||any(rej)){return(10e7)}
  
  nll.uc <- 0
  nll.pc <- 0

  uc<-apply(x,1,comp.gt,u=u)
  
  x.uc<-x[uc,]
  x.pc<-x[!uc,]
  
  L<-apply(x.uc,1,fX.powunif,a=a,lam=lam,sig=sig,gamma=gamma)
  nll.uc<- -sum(log(L)) 
  
  if(sum(!uc)>0)
  {
    L2<-apply(x.pc,1,fX.powunif.cens,a=a,lam=lam,u=u,sig=sig,gamma=gamma)
    nll.pc<--sum(log(L2))
  }
  
  if (is.nan(nll.uc)|is.nan(nll.pc)|(nll.uc==-Inf)|(nll.uc==Inf)){
    return(10e7)
  }
  nll<-nll.uc+nll.pc
  
  return(nll)
}

nll.Frechet.a.GPD.1<-function(theta,x,u)
{
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

fit.tMGPD.RevExpU <- function(x,lower,upper, init=NULL){
  u <- min(x)-0.01
  if (is.null(init)){
    par <- c(1,1, 1, 0, 0)
  }else{
    par <- init
  }

  a.ind <- 1
  sig.ind <- 2:3
  gamma.ind <- 4:5
  
  
  tgt <- function(par, u=u, x=x, a.ind=a.ind,
                  sig.ind=sig.ind, gamma.ind=gamma.ind, low=lower, upp=upper){
    sig <- par[sig.ind]
    gamma <- par[gamma.ind]
    res <- 0
    cond1 <- any(sig<0)
    cond2 <- rep(FALSE,2)
    if (!cond1){
      for (i in 1:length(sig)){
        if (gamma[i]>0 & -sig[i]/gamma[i]> min(x[,i])){
          cond2[i] <- TRUE
        }
        if (gamma[i]<0 & -sig[i]/gamma[i]< max(x[,i]))
          cond2[i] <- TRUE
      }
    }
    
    if (any(c(cond1,cond2))){
      res <- 10^10
    }
    
    if(res==0){
      nll <-   nll.powunif.GPD.1(theta=par, x=x,u=u)
      den <- trunc.norm.const.GPD(a=rep(par[a.ind],2), sig=par[sig.ind], gamma=par[gamma.ind],lower=low,upper=upp)
      if (is.nan(log(den))) print(par)
      n.obs <- nrow(x)
      res <- nll+ n.obs*log(den)
    }
    return(res)
  }
  
  
  
  
  
  opt<-optim(tgt, par=par, u=u, x=x, a.ind=a.ind,
             sig.ind=sig.ind, gamma.ind=gamma.ind, control=list(maxit=10000,reltol=1e-6), low=lower, upp=upper)
  if(is.null(opt$min))
  {
    mle<-opt$par
    
    nll<-opt$value
    conv<-opt$conv
    hess<-opt$hess
  }
  else{
    mle<-opt$minimum
    nll<-opt$objective
    conv<-NULL
    hess<-NULL
  }
  return(list(mle=mle,nll=nll,conv=conv,hess=hess))
}

log(trunc.norm.const.GPD(a=c(2,2), sig=c(1,2), gamma=c(0.1,0.6),lower=c(-Inf,-Inf),upper=c(Inf,Inf)))


trunc.norm.const.gumbel <- function(a,lower,upper){
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
trunc.norm.const.gumbel(4, c(-4,-3),c(4,5))



trunc.norm.const.gumbel(4, c(-Inf,-Inf),c(Inf,Inf))

trunc.norm.const.gumbel.GPD <- function(a, sig, gamma, lower, upper ){
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
      if(is.na(log(gamma[i]/sig[i]*upper[i]+1 + .Machine$double.eps*10))) print(list(i=i, gamma=gamma,sig=sig,upper=upper,
                                                            val=gamma[i]/sig[i]*upper[i]+1))
      lower[i] <- log(gamma[i]/sig[i]*lower[i]+1 + .Machine$double.eps*10 )/gamma[i]
      upper[i] <- log(gamma[i]/sig[i]*upper[i]+1 + .Machine$double.eps*10 )/gamma[i]
      
    }
  }
  return(trunc.norm.const.gumbel(a,lower,upper))
}



fit.tMGPD.GumbelU <- function(x,lower,upper,init=NULL){
  u <- min(x)-0.01
  if (is.null(init)){
    par <- c(2,1, 1, 0, 0)
  }else{
    par <- init
  }
  
  a.ind <- 1
  sig.ind <- 2:3
  gamma.ind <- 4:5
  
  
  tgt <- function(par, u=u, x=x, a.ind=a.ind,
                  sig.ind=sig.ind, gamma.ind=gamma.ind, low=lower, upp=upper){
    sig <- par[sig.ind]
    gamma <- par[gamma.ind]
    res <- 0
    cond1 <- any(sig<0)
    cond2 <- rep(FALSE,2)
    if (!cond1){
      for (i in 1:length(sig)){
        if (gamma[i]>0 & -sig[i]/gamma[i]> min(x[,i])){
          cond2[i] <- TRUE
        }
        if (gamma[i]<0 & -sig[i]/gamma[i]< max(x[,i]))
          cond2[i] <- TRUE
      }
    }
    
    if (any(c(cond1,cond2))){
      res <- 10^10
    }
    
    if(res==0){
      nll <- nll.Frechet.a.GPD.1(theta=par, x=x, u=u)
      den <- trunc.norm.const.gumbel.GPD(a=par[a.ind], sig=par[sig.ind], gamma=par[gamma.ind],lower=low,upper=upp)
      n.obs <- nrow(x)
      if (is.nan(log(den))) print(par)
      res <- nll+ n.obs*log(den)
    }
    return(res)
  }
  
  
  
  
  
  opt<-optim(tgt, par=par, u=u, x=x, a.ind=a.ind,
             sig.ind=sig.ind, gamma.ind=gamma.ind, control=list(maxit=10000,reltol=1e-6), low=lower, upp=upper)
  if(is.null(opt$min))
  {
    mle<-opt$par
    
    nll<-opt$value
    conv<-opt$conv
    hess<-opt$hess
  }
  else{
    mle<-opt$minimum
    nll<-opt$objective
    conv<-NULL
    hess<-NULL
  }
  return(list(mle=mle,nll=nll,conv=conv,hess=hess))
}



# Y.tail<-sim.RevExpU.MGPD(n=1000,d=2, a=2, beta=c(0,0), sig=c(0.7,0.6), gamma=c(0.1,0.2), MGPD = T,std=T)$X
# lbound.tail <- c(0,0)
# ubound.tail <- c(Inf,Inf)
# u.x <- c(0.5,0.5)
# 
# plot(Y.tail)
# bound.cond.1 <- Y.tail[,1] > (lbound.tail-u.x)[1] & Y.tail[,1] < (ubound.tail-u.x)[1]
# bound.cond.2 <- Y.tail[,2] > (lbound.tail-u.x)[2] & Y.tail[,2] < (ubound.tail-u.x)[2]
# Y.tail.t <- Y.tail[bound.cond.1&bound.cond.2,]
# plot(Y.tail.t)
# 
# fitGP.RevExpU<-fit.MGPD.RevExpU(x=Y.tail, u=c(0,0), std=F, dep.scale.fix=TRUE,dep.loc.fix=TRUE, marg.scale.ind = c(1,2), marg.shape.ind = 1:2, maxit=5000)
# 
# fitTGP.RevExpU <- fit.tMGPD.RevExpU(Y.tail.t,u=min(Y.tail)-0.01, lbound.tail-u.x, c(Inf,Inf))
# fitTGP.RevExpU1 <- fit.tMGPD.RevExpU(Y.tail.t,u=min(Y.tail)-0.01,lbound.tail-u.x, c(Inf,Inf),init=fitTGP.RevExpU$mle)
# fitTGP.RevExpU2 <- fit.tMGPD.RevExpU(Y.tail.t,u=min(Y.tail)-0.01,lbound.tail-u.x, c(Inf,Inf),init=fitTGP.RevExpU1$mle)
# 
# fitTGP.RevExpU2 <- fit.tMGPD.RevExpU(Y.tail.t,u=c(0,0), c(-Inf,-Inf), c(Inf,Inf))
# fitTGP.RevExpU3 <- fit.tMGPD.RevExpU(Y.tail.t,u=c(0,0), c(-Inf,-Inf), c(Inf,Inf),init=fitTGP.RevExpU2$mle)
# 

# 
# trunc.norm.const.GPD(a=rep(par[a.ind],2), sig=par[sig.ind], gamma=par[gamma.ind],lower=c(-Inf,-Inf),upper=c(Inf,Inf))
# 
# tgt(par=fitTGP.RevExpU$mle, u=c(0,0), x=Y.tail, a.ind=1,
#       sig.ind=2:3, gamma.ind=4:5, low= lbound.tail-u.x, upp=c(Inf,Inf))
# 
# opt<-optim(tgt, par=c(1,1,1,0,0), u=c(0,0), x=Y.tail, a.ind=a.ind,
#            sig.ind=sig.ind, gamma.ind=gamma.ind, control=list(maxit=10000,reltol=1e-6), low=lbound.tail-u.x, upp=c(Inf,Inf))
# opt2 <- optim(tgt, par=opt$par, u=c(0,0), x=x, a.ind=a.ind,
#               sig.ind=sig.ind, gamma.ind=gamma.ind, control=list(maxit=10000,reltol=1e-6), low=lbound.tail-u.x, upp=c(Inf,Inf))
# fitTGP.RevExpU1 <- fit.tMGPD.RevExpU(Y.tail, u=c(0,0), lbound.tail-u.x, c(Inf,Inf),init=fitTGP.RevExpU$mle)
# 
# 
# 
# Y.tail<-sim.GumbelU.MGPD(n=1000,d=2, a=3, beta=c(0,0), sig=c(10,20), gamma=c(-0.1,-0.2), MGPD = T,std=T)$X
# lbound.tail <- c(0,0)
# ubound.tail <- c(Inf,Inf)
# u.x <- c(1,2)
# 
# plot(Y.tail)
# bound.cond.1 <- Y.tail[,1] > (lbound.tail-u.x)[1] & Y.tail[,1] < (ubound.tail-u.x)[1]
# bound.cond.2 <- Y.tail[,2] > (lbound.tail-u.x)[2] & Y.tail[,2] < (ubound.tail-u.x)[2]
# Y.tail <- Y.tail[bound.cond.1&bound.cond.2,]
# plot(Y.tail)
# 
# fitGP.GumbelU<-fit.MGPD.GumbelU(x=Y.tail, u=rep(min(Y.tail)-0.01,2), std=F, dep.scale.fix=TRUE,dep.loc.fix=TRUE, marg.scale.ind = c(1,2), marg.shape.ind = 1:2, maxit=5000)
# 
# fitTGP.GumbelU <- fit.tMGPD.GumbelU(Y.tail, lbound.tail-u.x, c(Inf,Inf))
# fitTGP.GumbelU1 <- fit.tMGPD.GumbelU(Y.tail, lbound.tail-u.x, c(Inf,Inf),init=fitTGP.GumbelU$mle)
# 
# 



paralist <- list('gamma'= list(shape=1,rate=1), 
                 'norm'= list(mu=1,sd=1),
                 'exp' = list(rate = 1),
                 'lnorm' = list(meanlog = 0, sdlog = 1),
                 'weibull' = list(shape = 2, scale = 1),
                 'lst' = list(df=5, mu=0, sigma=1)
)

create_copula <- function(copula_name, param, dim = 2) {
  transformed_param <- switch(copula_name,
                              "normal" = tanh(param), # limit to [-1, 1]
                              "clayton" = exp(param), # transform to (0, Inf)
                              "gumbel" = exp(param)+1, # transform to [1, Inf)
                              "frank" = ifelse(abs(param) < 10^-3, 10^-3, param), # no transformation needed, but param != 0
                              "joe" = exp(param)+1, # transform to [1, Inf)
                              "plackett" = exp(param), # transform to (0, Inf)
                              stop("Unsupported copula name"))
  
  cop <- switch(copula_name,
                "normal" = normalCopula(param = transformed_param, dim = dim),
                "clayton" = claytonCopula(param = transformed_param, dim = dim),
                "gumbel" = gumbelCopula(param = transformed_param, dim = dim),
                "frank" = frankCopula(param = transformed_param, dim = dim),
                "joe" = joeCopula(param = transformed_param, dim = dim),
                "plackett" = plackettCopula(param = transformed_param),
                stop("Unsupported copula name"))
  
  return(cop)
}

create_copula('normal',-1.2)

transform_marginal_params <- function(margins, params.margin) {
  transformed_params <- list()
  start.idx <- 1
  for (i in 1:length(margins)) {
    m <- margins[i]
    
    if (m == "norm") {
      n.param <- 2
      p <- params.margin[ start.idx : (start.idx+1)]
      transformed_params[[i]] <- list(mean = p[1], sd = exp(p[2])) # ensure sd > 0
      start.idx <- start.idx + n.param
      
    } else if (m == "exp") {
      n.param <- 1
      p <- params.margin[ start.idx]
      transformed_params[[i]] <- list(rate = exp(p[1])) # ensure rate > 0
      start.idx <- start.idx + n.param
      
    } else if (m == "gamma") {
      n.param <- 2
      p <- params.margin[ start.idx : (start.idx+1)]
      transformed_params[[i]] <- list(shape = exp(p[1]), rate = exp(p[2])) # ensure shape, rate > 0
      start.idx <- start.idx + n.param
      
    } else if (m == "lnorm") {
      n.param <- 2
      p <- params.margin[ start.idx : (start.idx+1)]
      transformed_params[[i]] <- list(meanlog = p[1], sdlog = exp(p[2])) # ensure sdlog > 0
      start.idx <- start.idx + n.param
    } else if (m == "weibull") {
      n.param <- 2
      p <- params.margin[ start.idx : (start.idx+1)]
      transformed_params[[i]] <- list(shape = exp(p[1]), scale = exp(p[2])) # ensure shape, scale > 0
      start.idx <- start.idx + n.param
    } else if (m == 'lst'){
      n.param <- 3
      p <- params.margin[ start.idx : (start.idx+2)]
      transformed_params[[i]] <- list(df = exp(p[1]), mu = p[2], sigma = exp(p[3])) 
      start.idx <- start.idx + n.param
    } else {
      stop("Unsupported marginal distribution")
    }
  }
  
  return(transformed_params)
}

transform_marginal_params(c('norm','lst'),c(1,-22,-3,4,5))

nll.bulk <- function(par, cop.name, margin.name, Y, lbound, ubound){
  cop <- create_copula(cop.name, par[1])
  margin.para <- transform_marginal_params(margin.name, par[-1])
  
  joint.dist <- mvdc(cop, margins = margin.name, 
                     
                     paramMargins=list(margin.para[[1]],
                                       
                                       margin.para[[2]])
  )
  
  p1 <- pMvdc(ubound, joint.dist)
  p2 <- pMvdc(lbound, joint.dist)
  p <- p1 - p2 + .Machine$double.eps*10
  dbulk <- dMvdc(Y, joint.dist,log=T)
  dbulk[which(dbulk==-Inf | is.nan(dbulk))] <- -10^10
  
  return(-sum(dbulk) + nrow(Y)*log(p))
  
}




#################################Preicipitation applicaiton###########################
library(doParallel)
library(foreach)
library(tidyr)
num_cores <- detectCores() - 1
cl <- makeCluster(num_cores,outfile="")
registerDoParallel(cl)

load(file=file.path('/home/pgrad2/2448355h/My_PhD_Project/00_Dataset','datamat5days1218.Rdata'))


df <- drop_na(as.data.frame(datamat[,c(78,79)]))
# df <- drop_na(as.data.frame(datamat[,c(219,217)]))
# df <- drop_na(as.data.frame(datamat[,c(327,328)]))

x1 <- df[,1]/sd(df[,1])
x2 <- df[,2]/sd(df[,2])

x1 <- df[,1]/sd(df[,1])
x2 <- df[,2]/sd(df[,2])


Y <- cbind(x1, x2)
Y <- Y[Y[,1]>0 & Y[,2]>0,]

plot(Y)



cop.family <- c('normal','clayton','gumbel','frank','joe','plackett')
margin.family <- c('norm','exp','gamma','lnorm','weibull','lst')


initial.val <- function(margin.name, data=Y.fit){
  initial <- 1
  for (i in 1:2){
    name <- margin.name[i]
    if (name=='norm'){
      initial <- c(initial, mean(data[,i]), log(sd(data[,i])))
    }else if( name== 'exp'){
      if (mean(data[,i]) > 0 ){
        initial <- c(initial, log(1/mean(data[,i])))
      }else{
        initial <- c(initial, 1)
      }
      
    }else if( name == 'gamma'){
      if (mean(data[,i]) > 0 ){
        initial <- c(initial,  (log(mean(data[,i]) / sd(data[,i]))^2), log(mean(data[,i]) / var(data[,i])))
      }else{
        initial <- c(initial, 1, 1)
      }
    }else if( name=='lnorm'){
      initial <- c(initial,  mean(data[,i]), log(sd(data[,i])))
    }else if( name=='lst'){
      initial <- c(initial,1, mean(data[,i]), log(sd(data[,i])))
    }else{
      shape <- 1.2
      if (mean(data[,i]) > 0 ){
        initial <- c(initial,  log(shape) , log(mean(data[,i]) / gamma(1 + 1/shape) ))
      }else{
        initial <- c(initial,  log(shape) , 1 / gamma(1 + 1/shape) )
      }
    }
  }
  return(initial)
}


thres <- apply(Y, 2, quantile, 0.95)
Y.fit <- Y[Y[,1]<thres[1]&Y[,2]<thres[2],]
plot(Y.fit)



res.bulk <- foreach(margin1.name = margin.family , .combine='rbind',.packages=c('copula','extraDistr')) %:%
  foreach(margin2.name = margin.family, .combine='rbind',.packages=c('copula','extraDistr'))  %:%
    foreach(cop.name = cop.family, .combine='rbind',.packages=c('copula','extraDistr')) %dopar% {
    margin.name <- c(margin1.name,margin2.name)
    print(c(cop.name,margin.name))
    init <- initial.val(margin.name,Y.fit)
    res0 <- optim(par=init, nll.bulk,cop.name=cop.name,margin.name=margin.name, Y=Y.fit, lbound=rep(0,2),ubound=thres )
    res1 <- optim(par=res0$par, nll.bulk,cop.name=cop.name,margin.name=margin.name, Y=Y.fit, lbound=rep(0,2),ubound=thres )
    res2 <- optim(par=res1$par, nll.bulk,cop.name=cop.name,margin.name=margin.name, Y=Y.fit, lbound=rep(0,2),ubound=thres )
    k1 <- if (margin.name[1]=='exp') 1 else if (margin.name[1]=='lst') 3 else 2
    k2 <- if (margin.name[2]=='exp') 1 else if (margin.name[2]=='lst') 3 else 2
    aic <- 2*(k1+k2) + 2*res2$value
    data.frame('copula'=cop.name,'margin1'=margin1.name, 'margin2'=margin2.name, 'nll'= res2$value, 'aic'=aic)
  }

as.character(res.bulk[which.min(res.bulk$aic),c('copula','margin1','margin2')])
res.bulk[res.bulk$nll< 1040,]

cop.name <- 'gumbel'
margin.name <- c('gamma','gamma')
init <- initial.val(margin.name,Y.fit)
res0 <- optim(par=init, nll.bulk,cop.name=cop.name,margin.name=margin.name, Y=Y.fit, lbound=rep(0,2),ubound=thres )
res1 <- optim(par=res0$par, nll.bulk,cop.name=cop.name,margin.name=margin.name, Y=Y.fit, lbound=rep(0,2),ubound=thres )
res2 <- optim(par=res1$par, nll.bulk,cop.name=cop.name,margin.name=margin.name, Y=Y.fit, lbound=rep(0,2),ubound=thres )



res.bulk[res.bulk$copula=='plackett',]




thres <- apply(Y, MARGIN = 2, quantile, 0.9)

tail.cond <- (Y[,1] > thres[1]) | (Y[,2] > thres[2])
Y.tail <- matrix(c(Y[tail.cond,1] - thres[1], Y[tail.cond,2] - thres[2]), ncol=2)
plot(Y.tail)
lbound.tail <- c(0,0)
ubound.tail <- c(Inf,Inf)

bound.cond.1 <- Y.tail[,1] > (lbound.tail-thres)[1] & Y.tail[,1] < (ubound.tail-thres)[1]
bound.cond.2 <- Y.tail[,2] > (lbound.tail-thres)[2] & Y.tail[,2] < (ubound.tail-thres)[2]
Y.tail.trunc <- Y.tail[bound.cond.1&bound.cond.2,]

plot(Y.tail.trunc)
fitGP.RevExpU<-fit.MGPD.RevExpU(x=Y.tail, u=c(0,0), std=F, dep.scale.fix=TRUE,dep.loc.fix=TRUE, marg.scale.ind = c(1,2), marg.shape.ind = 1:2, maxit=5000)

fitTGP.RevExpU <- fit.tMGPD.RevExpU(Y.tail.trunc, lbound.tail-thres, c(Inf,Inf))
fitTGP.RevExpU1 <- fit.tMGPD.RevExpU(Y.tail.trunc,lbound.tail-thres, c(Inf,Inf),init=fitTGP.RevExpU$mle)
fitTGP.RevExpU2 <- fit.tMGPD.RevExpU(Y.tail.trunc,lbound.tail-thres, c(Inf,Inf),init=fitTGP.RevExpU1$mle)


fitGP.GumbelU<-fit.MGPD.GumbelU(x=Y.tail, u=min(Y.tail)-0.01, std=F, dep.scale.fix=TRUE,dep.loc.fix=TRUE, marg.scale.ind = c(1,2), marg.shape.ind = 1:2, maxit=5000)
fitTGP.GumbelU <- fit.tMGPD.GumbelU(Y.tail.trunc, lbound.tail-thres, c(Inf,Inf))
fitTGP.GumbelU1 <- fit.tMGPD.GumbelU(Y.tail.trunc, lbound.tail-thres, c(Inf,Inf),init=fitTGP.GumbelU$mle)
fitTGP.GumbelU2 <- fit.tMGPD.GumbelU(Y.tail.trunc, lbound.tail-thres, c(Inf,Inf),init=fitTGP.GumbelU1$mle)

# 
mle <- fitTGP.GumbelU2$mle
Y.sim <- sim.GumbelU.MGPD(n=1000,d=2, a=mle[1], beta=c(0,0), sig=mle[2:3],
                                           gamma=mle[4:5], MGPD = T,std=F)
# mle <- fitTGP.RevExpU2$mle
# Y.sim <- sim.RevExpU.MGPD(n=1000,d=2, a=mle[1], beta=c(0,0), sig=mle[2:3],
#                  gamma=mle[4:5], MGPD = T,std=F)
plot(Y.sim)
plot(Y.tail.trunc)
plot(Y.sim[Y.sim[,1]>(lbound.tail-thres)[1] & Y.sim[,2] > (lbound.tail-thres)[2] ,])

library(evd)
chiplot(Y.tail.trunc)
chiplot(Y.sim[Y.sim[,1]>(lbound.tail-thres)[1] & Y.sim[,2] > (lbound.tail-thres)[2] ,])

