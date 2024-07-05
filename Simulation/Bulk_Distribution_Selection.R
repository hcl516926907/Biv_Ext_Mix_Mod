dir.work <- '/home/pgrad2/2448355h/My_PhD_Project/Biv_Ext_Mix_Mod'
dir.out <- "/home/pgrad2/2448355h/My_PhD_Project/01_Output/Biv_Ext_Mix_Mod/Simulation"
source(file.path(dir.work, "Simulation/RevExp_U_Functions.r"))
source(file.path(dir.work, "Simulation/Gumbel_U_Functions.r"))
source(file.path(dir.work, "Simulation/CommonFunctions.r"))

# install.packages("gsl", dependencies = TRUE, INSTALL_opts = '--no-lock')
# install.packages('copula')
library(copula)
library(extraDistr)
i <- 1
seed <- i
d <- 2
a <- c(2, 2.5)
beta <- c(0, 0)
sig <- c(1, 0.8)
gamma <- c(-0.1, -0.2)
# a <- c(3, 2)
# beta <- c(0, 0)
# sig <- c(2, 3)
# gamma <- c(-0.2, -0.6)
n <- 5000

mu <- c(7, 8)
sd <- c(2, 2.5)
rho <- 0.4

u.x <- c(5, 5.5)

lbound.tail=c(1,1)
ubound.tail=c(Inf,Inf)

norm.cop <- gumbelCopula(param = 1.5, dim = 2)
# 
# gaussian <- mvdc(norm.cop, margins = c('gamma','weibull'),
# 
#                  paramMargins=list(list(shape=mu[1], rate=sd[1]),
# 
#                                    list(shape=2,scale=4)))
# 

u.x <- c(3,2)
gaussian <- mvdc(norm.cop, margins = c('lst','lst'),

                 paramMargins=list(list(df=4, mu=2,sigma=0.5),

                                   list(df=5, mu=1,sigma=0.6)))

p <- pMvdc(u.x, gaussian)


set.seed(1111)

a <- c(2, 2)
beta <- c(0, 0)
sig <- c(0.3, 0.3)
gamma <- c(0.22, 0.2)

##### tail bound is defined on the observational scale
lbound.tail=c(-Inf,-Inf)
ubound.tail=c(Inf,Inf)


Y.tail<-sim.RevExpU.MGPD(n=n-floor(n*p),d=d, a=a, beta=beta, sig=sig, gamma=gamma, MGPD = T,std=T)$X

# GP scale tail data combined with the bulk data

Y.bulk <- rMvdc(floor(n*p), gaussian)
Y.bulk <- Y.bulk[Y.bulk[,1]<u.x[1] & Y.bulk[,2]<u.x[2],]

Y <- rbind(Y.bulk, sweep(Y.tail,2,u.x,"+"))
apply(Y.bulk,2,min) + u.x
plot(Y)
bound.cond.1 <- Y.tail[,1] > (lbound.tail-u.x)[1] & Y.tail[,1] < (ubound.tail-u.x)[1]
bound.cond.2 <- Y.tail[,2] > (lbound.tail-u.x)[2] & Y.tail[,2] < (ubound.tail-u.x)[2]
Y.tail <- Y.tail[bound.cond.1&bound.cond.2,]
Y <- rbind(Y.bulk, sweep(Y.tail,2,u.x,"+"))
plot(Y)



a <- c(3,2)
y1 <- c(1,2)
y2 <- c(3,4)
fY.powunif(y1,c(1,1),a)
fY.powunif(y2,c(1,1),a)

EM <- EM.pu(a=a, lam=c(1,1))





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


marginY1.revexp(-10, y2.upp=Inf, a=a, sig=sig, gamma=gamma)

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

marginY2.revexp(-0.1,y1.upp=Inf, a=a, sig=sig, gamma=gamma)

double.int.Y1 <- function(x, y2.low=-Inf, y2.upp=Inf, a=c(1,1), sig=c(1,1), gamma=c(0,0)){
  sapply(x, marginY1.revexp, y2.low = y2.low, y2.upp=y2.upp, a=a, sig=sig, gamma=gamma)
}

integrate(double.int.Y1, lower =-Inf, upper = 10 , a=a,sig=c(1,2), gamma=c(-0.1,0))

double.int.Y2 <- function(x, y1.low=-Inf, y1.upp=Inf, a=c(1,1), sig=c(1,1), gamma=c(0,0)){
  sapply(x, marginY2.revexp, y1.low = y1.low, y1.upp=y1.upp, a=a, sig=sig, gamma=gamma)
}

integrate(double.int.Y2, lower =-Inf, upper = Inf , a=a,sig=c(1,2), gamma=c(-0.1,0))


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

norm.const <- trunc.prob(lower=lbound.tail-u.x,a=a,sig=sig,gamma=gamma)

# standarized Y1,Y1<0 with sig=1 and gamma=0
fY1_lt0.margin <- function(y,a,ub){
  EM <- EM.pu(a=a,lam=c(1,1))
  b <- 1/a
  den <- EM*(1+sum(b))
  num <- prod(b)*exp(b[1]*y)*(1-exp(-(1+b[1])*ub))/(1+b[1])
  return(num/den)
}

fY1_lt0.margin(-5,a=a,ub=Inf)

marginY1.revexp(-5, y2.upp=Inf, a=a, sig=c(1,1), gamma=c(0,0))

# standarized Y1,Y1<0 with sig=1 and gamma=0
fY2_lt0.margin <- function(y,a,ub){
  EM <- EM.pu(a=a,lam=c(1,1))
  b <- 1/a
  den <- EM*(1+sum(b))
  num <- prod(b)*exp(b[2]*y)*(1-exp(-(1+b[2])*ub))/(1+b[2])
  return(num/den)
}

fY2_lt0.margin(-5,a=a,ub=5)
marginY2.revexp(-5, y1.upp=5, a=a, sig=c(1,1), gamma=c(0,0))

fY2_gt0.margin <- function(y, a, lb, ub){
  EM <- EM.pu(a=a,lam=c(1,1))
  b <- 1/a
  num1 <- b[2]*(exp(-y)-exp(b[1]*lb-(1+b[1])*y))
  den1 <- EM*(1+sum(b))
  
  num2 <- prod(b)*exp(-y)-prod(b)*exp(b[2]*y-(1+b[2])*ub)
  den2 <- EM*(1+sum(b))*(1+b[2])
  
  return(num1/den1+num2/den2)
}
fY2_gt0.margin(0.1,a=a, lb=-1,ub=5)
marginY2.revexp(0.1, y1.low=-1, y1.upp=5, a=a, sig=c(1,1), gamma=c(0,0))

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


trunc.norm.const.GPD(a=c(4,1), sig=c(0.1,2), gamma=c(0.2,0.2),lower=c(-0.2,-0.2), upper=c(4,5))
trunc.prob(lower=c(-0.2,-0.2),upper=c(4,5),a=c(4,1),sig=c(0.1,2),gamma=c(0.2,0.2))



fitGP<-fit.MGPD.RevExpU(x=Y.tail, u=rep(min(Y.tail)-0.01,2), std=F, dep.loc.fix=TRUE, marg.scale.ind = c(1,2), marg.shape.ind = 1:2, maxit=5000)
# Parameter order: (a1,a2,beta1,sig,gamma1,gamma2)


fit.tMGPD.RevExpU <- function(x,lower,upper){
  par <- c(1,1,1, 1, 0, 0)
  a.ind <- 1:2
  sig.ind <- 3:4
  gamma.ind <- 5:6
  
  
  tgt <- function(par, u=min(x)-0.01, x=x, a.ind=a.ind,
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
      nll <-   nll.powunif.GPD(theta=par, u=min(x)-0.01, x=x, lamfix=T, a.ind=a.ind,
                               sig.ind=sig.ind, gamma.ind=gamma.ind, marg.scale.ind=1:2, marg.shape.ind=1:2)
      den <- trunc.norm.const.GPD(a=par[a.ind], sig=par[sig.ind], gamma=par[gamma.ind],lower=low,upper=upp)
      
      n.obs <- nrow(x)
      res <- nll+ n.obs*log(den)
    }
    return(res)
  }
  


  
  
  opt<-optim(tgt, par=par, u=min(x)-0.01, x=x, a.ind=a.ind,
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

fitTGP <- fit.tMGPD.RevExpU(Y.tail, lbound.tail-u.x, c(Inf,Inf))



mix.p.hat <- nrow(Y.bulk)/nrow(Y)
tail.norm.const <- trunc.norm.const.GPD(a=c(4,1), sig=c(0.1,2), gamma=c(0.2,0.2),lower=c(-Inf,-Inf), upper=c(Inf,Inf))


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
  p <- p1 - p2 + 10^-6
  dbulk <- dMvdc(Y, joint.dist,log=T)
  dbulk[which(dbulk==-Inf)] <- -10^10
  return(-sum(dbulk) + nrow(Y)*log(p))
  
}

nll.bulk(rep(0,5),cop.name='normal',margin.name=c('gamma','gamma'), Y=Y.bulk, lbound=c(0,0),ubound=u.x)
nll.bulk(c(-0.93,log(mu[1]),log(sd[1]),log(2),log(4)),cop.name='gumbel',margin.name=c('gamma','weibull'), Y=Y.bulk, lbound=c(0,0),ubound=u.x)

nll.bulk(c(-0.9,4,2,log(0.6),5,1, log(0.5)),cop.name='gumbel',margin.name=c('lst','lst'), Y=Y.bulk, lbound=rep(-Inf,2),ubound=u.x)


transform_marginal_params(c('norm','norm'),c(2,log(0.6),1, log(0.5)))
nll.bulk(c(log(1.5),2,log(0.6),1, log(0.5)),cop.name='gumbel',margin.name=c('norm','norm'), Y=Y.bulk, lbound=rep(-Inf,2),ubound=u.x)

nll.bulk(c(log(1.5),log(4),2,log(0.6),log(5),1, log(0.5)),cop.name='gumbel',margin.name=c('lst','lst'), Y=Y.bulk, lbound=rep(-Inf,2),ubound=u.x)


test <- optim(par=rep(0,7), nll.bulk,cop.name='gumbel',margin.name=c('lst','lst'), Y=Y.bulk, lbound=rep(-Inf,2),ubound=u.x )

test2 <- optim(par=test$par, nll.bulk,cop.name='gumbel',margin.name=c('lst','lst'), Y=Y.bulk, lbound=rep(-Inf,2),ubound=u.x)
test3 <- optim(par=test2$par, nll.bulk,cop.name='gumbel',margin.name=c('lst','lst'), Y=Y.bulk, lbound=rep(-Inf,2),ubound=u.x )
test4 <- optim(par=test3$par, nll.bulk,cop.name='gumbel',margin.name=c('lst','lst'), Y=Y.bulk, lbound=rep(-Inf,2),ubound=u.x )



nll.bulk(test3$par,cop.name='gumbel',margin.name=c('lst','lst'), Y=Y.bulk, lbound=c(0,0),ubound=u.x)
nll.bulk(c(-0.1872918, 1.2 , 2.0226106, -0.8886660  ,1.5,  1.0140278, -0.5081560),cop.name='gumbel',margin.name=c('lst','lst'), Y=Y.bulk, lbound=c(0,0),ubound=u.x)

# cop.col <- c()
# marg1.col <- c()
# marg2.col <- c()
# nll.col <- c()
# for (cop.name in cop.family){
#   print(cop.name)
#   for (margin1.name in margin.family){
#     for (margin2.name in margin.family){
#       margin.name <- c(margin1.name,margin2.name)
#       res0 <- optim(par=rep(0,5), nll.bulk,cop.name=cop.name,margin.name=margin.name, Y=Y.bulk, lbound=c(0,0),ubound=u.x )
#       res1 <- optim(par=res0$par, nll.bulk,cop.name=cop.name,margin.name=margin.name, Y=Y.bulk, lbound=c(0,0),ubound=u.x )
#       res2 <- optim(par=res1$par, nll.bulk,cop.name=cop.name,margin.name=margin.name, Y=Y.bulk, lbound=c(0,0),ubound=u.x )
#       cop.col <- c(cop.col, cop.name)
#       marg1.col <- c(marg1.col, margin1.name)
#       marg2.col <- c(marg2.col, margin2.name)
#       nll.col <- c(nll.col, res2$value)
#     }
#   }
# }

#################################Simulation###########################
library(doParallel)
library(foreach)
num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)

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
initial.val(c('lst','lnorm'),Y.bulk)

for (margin1.name in margin.family[1]){
  for (margin2.name in margin.family[6]){
    margin.name <- c(margin1.name,margin2.name)
    print(margin.name)
    init <- initial.val(margin.name,Y.fit)
    res0 <- optim(par=init, nll.bulk,cop.name='normal',margin.name=margin.name, Y=Y.fit, lbound=rep(-Inf,2),ubound=u.x )
    res1 <- optim(par=res0$par, nll.bulk,cop.name='normal',margin.name=margin.name, Y=Y.fit, lbound=rep(-Inf,2),ubound=u.x )
    res2 <- optim(par=res1$par, nll.bulk,cop.name='normal',margin.name=margin.name, Y=Y.fit, lbound=rep(-Inf,2),ubound=u.x )
    k1 <- if (margin.name[1]=='exp') 1 else 2
    k2 <- if (margin.name[2]=='exp') 1 else 2
    aic <- 2*(k1+k2) + 2*res2$value
    data.frame('margin1'=margin1.name, 'margin2'=margin2.name, 'nll'= res2$value, 'aic'=aic)
    
  }
  
}




thres <- u.x
Y.fit <- Y.bulk
res.margin <- foreach(margin1.name = margin.family , .combine='rbind',.packages='copula') %:%
  foreach(margin2.name = margin.family, .combine='rbind',.packages='copula') %dopar% {
    margin.name <- c(margin1.name,margin2.name)
    print(margin.name)
    init <- initial.val(margin.name,Y.fit)
    res0 <- optim(par=init, nll.bulk,cop.name='normal',margin.name=margin.name, Y=Y.fit, lbound=rep(-Inf,2),ubound=thres )
    res1 <- optim(par=res0$par, nll.bulk,cop.name='normal',margin.name=margin.name, Y=Y.fit, lbound=rep(-Inf,2),ubound=thres )
    res2 <- optim(par=res1$par, nll.bulk,cop.name='normal',margin.name=margin.name, Y=Y.fit, lbound=rep(-Inf,2),ubound=thres )
    k1 <- if (margin.name[1]=='exp') 1 else 2
    k2 <- if (margin.name[2]=='exp') 1 else 2
    aic <- 2*(k1+k2) + 2*res2$value
    data.frame('margin1'=margin1.name, 'margin2'=margin2.name, 'nll'= res2$value, 'aic'=aic)
}

best.margin <- as.character(res.margin[which.min(res.margin$aic),c('margin1','margin2')])


res.copula <- foreach(cop.name = cop.family, .combine='rbind',.packages='copula') %dopar% {
    margin.name <- best.margin
    print(margin.name)
    init <- initial.val(margin.name,Y.fit)
    res0 <- optim(par=init, nll.bulk,cop.name=cop.name,margin.name=margin.name, Y=Y.fit, lbound=c(0,0),ubound=thres )
    res1 <- optim(par=res0$par, nll.bulk,cop.name=cop.name,margin.name=margin.name, Y=Y.fit, lbound=c(0,0),ubound=thres )
    res2 <- optim(par=res1$par, nll.bulk,cop.name=cop.name,margin.name=margin.name, Y=Y.fit, lbound=c(0,0),ubound=thres )
    k1 <- if (margin.name[1]=='exp') 1 else 2
    k2 <- if (margin.name[2]=='exp') 1 else 2
    aic <- 2*(k1+k2+1) + 2*res2$value
    data.frame('copula'=cop.name,'margin1'=margin.name[1], 'margin2'=margin.name[2], 'nll'= res2$value, 'aic'=aic)
}

best.bulk<- as.character(res.copula[which.min(res.copula$aic),c('copula','margin1','margin2')])


#################################Financial applicaiton###########################
library(doParallel)
library(foreach)
num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)

load(file=file.path('/home/pgrad2/2448355h/My_PhD_Project/00_Dataset','index.daily.RData'))


library("rugarch")
spec <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                   mean.model = list(armaOrder = c(0, 0), include.mean = TRUE),
                   distribution.model = "norm")

# Fit the GARCH model to the daily returns
fit.ftse <- ugarchfit(spec = spec, data = -log(index.daily.dropna$return_ftse))
residuals.ftse <- residuals(fit.ftse, standardize = TRUE)

fit.cac <- ugarchfit(spec = spec, data = -log(index.daily.dropna$return_cac))
residuals.cac <- residuals(fit.cac, standardize = TRUE)

fit.dax <- ugarchfit(spec = spec, data = -log(index.daily.dropna$return_dax))
residuals.dax <- residuals(fit.dax, standardize = TRUE)

plot(as.numeric(residuals.ftse), as.numeric(residuals.dax))

plot(as.numeric(residuals.ftse), as.numeric(residuals.cac))

plot(as.numeric(residuals.dax), as.numeric(residuals.cac))

# Y.daily.return <- as.matrix(100*abs(1-index.daily.dropna[,c('return_ftse','return_dax')]))
Y.daily.return <- cbind('log.ftse.return'=as.numeric(residuals.ftse), 'log.dax.return'=as.numeric(residuals.dax) )
thres <- apply(Y.daily.return, 2, quantile, 0.6)
Y.fit <-  Y.daily.return[Y.daily.return[,1]<=thres[1] & Y.daily.return[,2]<=thres[2],]

plot(Y.fit)


cop.family <- c('normal','clayton','gumbel','frank','joe','plackett')
margin.family <- c('norm','exp','gamma','lnorm','weibull','lst')

# for (margin1.name in margin.family){
#   for (margin2.name in margin.family){
#     margin.name <- c(margin1.name,margin2.name)
#     print(margin.name)
#     init <- initial.val(margin.name,Y.fit)
#     res0 <- optim(par=init, nll.bulk,cop.name='normal',margin.name=margin.name, Y=Y.fit, lbound=rep(-Inf,2),ubound=thres )
#     res1 <- optim(par=res0$par, nll.bulk,cop.name='normal',margin.name=margin.name, Y=Y.fit, lbound=rep(-Inf,2),ubound=thres )
#     res2 <- optim(par=res1$par, nll.bulk,cop.name='normal',margin.name=margin.name, Y=Y.fit, lbound=rep(-Inf,2),ubound=thres )
#     k1 <- if (margin.name[1]=='exp') 1 else if (margin.name[1]=='lst') 3 else 2
#     k2 <- if (margin.name[2]=='exp') 1 else if (margin.name[2]=='lst') 3 else 2
#     aic <- 2*(k1+k2) + 2*res2$value
#     data.frame('margin1'=margin1.name, 'margin2'=margin2.name, 'nll'= res2$value, 'aic'=aic)
# 
#   }
# }


res.margin <- foreach(margin1.name = margin.family , .combine='rbind',.packages='copula') %:%
  foreach(margin2.name = margin.family, .combine='rbind',.packages='copula') %dopar% {
    margin.name <- c(margin1.name,margin2.name)
    init <- initial.val(margin.name,Y.fit)
    res0 <- optim(par=init, nll.bulk,cop.name='normal',margin.name=margin.name, Y=Y.fit, lbound=rep(-Inf,2),ubound=thres )
    res1 <- optim(par=res0$par, nll.bulk,cop.name='normal',margin.name=margin.name, Y=Y.fit, lbound=rep(-Inf,2),ubound=thres )
    res2 <- optim(par=res1$par, nll.bulk,cop.name='normal',margin.name=margin.name, Y=Y.fit, lbound=rep(-Inf,2),ubound=thres )
    k1 <- if (margin.name[1]=='exp') 1 else if (margin.name[1]=='lst') 3 else 2
    k2 <- if (margin.name[2]=='exp') 1 else if (margin.name[2]=='lst') 3 else 2
    aic <- 2*(k1+k2) + 2*res2$value
    data.frame('margin1'=margin1.name, 'margin2'=margin2.name, 'nll'= res2$value, 'aic'=aic)
  }

best.margin <- as.character(res.margin[which.min(res.margin$aic),c('margin1','margin2')])


res.copula <- foreach(cop.name = cop.family, .combine='rbind',.packages='copula') %dopar% {
  margin.name <- best.margin
  init <- initial.val(margin.name,Y.fit)
  res0 <- optim(par=init, nll.bulk,cop.name=cop.name,margin.name=margin.name, Y=Y.fit, lbound=rep(-Inf,2),ubound=thres )
  res1 <- optim(par=res0$par, nll.bulk,cop.name=cop.name,margin.name=margin.name, Y=Y.fit, lbound=rep(-Inf,2),ubound=thres )
  res2 <- optim(par=res1$par, nll.bulk,cop.name=cop.name,margin.name=margin.name, Y=Y.fit, lbound=rep(-Inf,2),ubound=thres )
  k1 <- if (margin.name[1]=='exp') 1 else 2
  k2 <- if (margin.name[2]=='exp') 1 else 2
  aic <- 2*(k1+k2+1) + 2*res2$value
  data.frame('copula'=cop.name,'margin1'=margin.name[1], 'margin2'=margin.name[2], 'nll'= res2$value, 'aic'=aic)
}

best.bulk<- as.character(res.copula[which.min(res.copula$aic),c('copula','margin1','margin2')])

thres <- apply(Y.fit, MARGIN = 2, quantile, 0.8)
tail.cond <- (Y.fit[,1] > thres[1]) | (Y.fit[,2] > thres[2])
Y.tail <- matrix(c(Y.fit[tail.cond,1] - thres[1], Y.fit[tail.cond,2] - thres[2]), ncol=2)
plot(Y.tail)

fit.RevExpU<-fit.MGPD.RevExpU(x=Y.tail, u=rep(min(Y.tail)-0.01,2), std=F, dep.scale.fix=F,dep.loc.fix=TRUE, marg.scale.ind = c(1,2), marg.shape.ind = 1:2, maxit=5000)

fit.RevExpU.1<-fit.MGPD.RevExpU(x=Y.tail, u=rep(min(Y.tail)-0.01,2), std=F, dep.scale.fix=T,dep.loc.fix=TRUE, marg.scale.ind = c(1,2), marg.shape.ind = 1:2, maxit=5000)

fit.RevExpU.2<-fit.MGPD.RevExpU(x=Y.tail, u=rep(0,2), std=F, dep.scale.fix=F,dep.loc.fix=TRUE, marg.scale.ind = c(1,2), marg.shape.ind = 1:2, maxit=5000)

ubound.tail <- rep(8,2)
fit.RevExpU.t <- fit.tMGPD.RevExpU(Y.tail, rep(-Inf,2), ubound.tail)

par <- fit.RevExpU$mle
a <- par[1:2]
sig <- par[3:4]
gamma <- par[5:6]
Y.sim <- sim.RevExpU.MGPD(n=100,d=2, a=a, beta=c(0,0), sig=sig, gamma=gamma, MGPD = T,std=T)$X
plot(Y.sim)
Y.sim.t <- Y.sim[Y.sim[,1]<ubound.tail[1] & Y.sim[,2]<ubound.tail[2],]
plot(Y.sim.t)

optim(nll.powunif.GPD, par=fit.RevExpU$mle, u=rep(min(Y.tail)-0.01,2), x=Y.tail, lamfix=T, a.ind=1:2,
      sig.ind=3:4, gamma.ind=5:6, marg.scale.ind=1:2, marg.shape.ind=1:2, control=list(maxit=5000,reltol=1e-6))

optim(nll.powunif.GPD.1, par=rep(1,6),  x=Y.tail)


-nim_nll_powunif_GPD_MAT(theta,Y.tail)
-nim_nll_powunif_GPD_MAT(fit.RevExpU$mle,Y.tail)

fitGumbelU<-fit.MGPD.GumbelU(x=Y.tail, u=rep(min(Y.tail)-0.01,2), std=F, dep.scale.fix=T, dep.loc.fix = TRUE, marg.scale.ind = c(1,2), marg.shape.ind = c(1,2), maxit=5000)

Y.sim2 <- sim.GumbelU.MGPD(n=1000,d=2,a=fitGumbelU$mle[1],beta=c(0,0),sig=fitGumbelU$mle[2:3],gamma=fitGumbelU$mle[4:5],MGPD=TRUE,std=FALSE)
plot(Y.sim2)


d<-2
a<-c(1.5) # common shape (quicker estimation)
beta<-c(0,0)
sig<-c(2,2)
gamma<-c(0.1,0.1)

# Final beta parameter fixed at zero in estimation, so needs to be relative to this

# Simulate with (conditionally) GP(sig,gamma) and (conditionally) exponential margins 
X<-sim.GumbelU.MGPD(n=1000,d=d, a=a, beta=beta, sig=sig, gamma=gamma, MGPD = T,std=T)
