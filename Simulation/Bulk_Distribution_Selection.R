dir.work <- '/home/pgrad2/2448355h/My_PhD_Project/Biv_Ext_Mix_Mod'
dir.out <- "/home/pgrad2/2448355h/My_PhD_Project/01_Output/Biv_Ext_Mix_Mod/Simulation"
source(file.path(dir.work, "Simulation/RevExp_U_Functions.r"))
source(file.path(dir.work, "Simulation/CommonFunctions.r"))

library(copula)

i <- 1
seed <- i
d <- 2
a <- c(1, 1.2)
beta <- c(0, 0)
sig <- c(1, 1.5)
gamma <- c(-0.1, -0.3)
n <- 5000

mu <- c(7, 8)
sd <- c(2, 2.5)
rho <- 0.4

u.x <- c(5, 5.5)

lbound=c(1,1)
ubound=c(Inf,Inf)

norm.cop <- ellipCopula("normal", param = rho, dim = 2)

gaussian <- mvdc(norm.cop, margins = c('gamma','gamma'), 

                 paramMargins=list(list(shape=mu[1], rate=sd[1]),
                                   
                                   list(shape=mu[2],rate=sd[2])))

p <- pMvdc(u.x, gaussian)

set.seed(1111)
Y.tail<-sim.RevExpU.MGPD(n=n-floor(n*p),d=d, a=a, beta=beta, sig=sig, gamma=gamma, MGPD = T,std=T)$X

# GP scale tail data combined with the bulk data

Y.bulk <- rMvdc(floor(n*p), gaussian)
Y.bulk <- Y.bulk[Y.bulk[,1]<u.x[1] & Y.bulk[,2]<u.x[2],]

Y <- rbind(Y.bulk, sweep(Y.tail,2,u.x,"+"))
apply(Y.bulk,2,min) + u.x
plot(Y)
bound.cond.1 <- Y.tail[,1] > (lbound-u.x)[1] & Y.tail[,1] < (ubound-u.x)[1]
bound.cond.2 <- Y.tail[,2] > (lbound-u.x)[2] & Y.tail[,2] < (ubound-u.x)[2]
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


marginY1.revexp(101, y2.upp=Inf, a=a, sig=sig, gamma=gamma)

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

norm.const <- trunc.prob(lower=lbound-u.x,a=a,sig=sig,gamma=gamma)

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

fitGP<-fit.MGPD.RevExpU(x=Y.tail, u=rep(min(Y.tail)-0.01,2), std=F, dep.loc.fix=TRUE, marg.scale.ind = c(1,2), marg.shape.ind = 1:2, maxit=5000)
# Parameter order: (a1,a2,beta1,sig,gamma1,gamma2)



