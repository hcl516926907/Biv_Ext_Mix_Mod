# ***************************
# ***************************
# Lídia André 
# March 2023
# Lancaster University, UK
# ***************************
# ***************************

rm(list=ls())
library(pracma)
library(copula)
library(parallel)
library(evd)
# library(SpatialADAI)

# =================
# Copula densities
# =================

## Gaussian copula
### rho: correlation parameter
gausdensity<-function(u,v,rho){
  s<-qnorm(u)
  r<-qnorm(v)
  first<-1/(sqrt(1-rho^2))
  second<-exp(-(rho^2*s^2+rho^2*r^2-2*rho*s*r)/(2*(1-rho^2)))
  dens<-first*second
  return(dens)
}

## Student t copula
### eta: degrees of freedom
### rho: correlation parameter
tdensity<-function(u,v,eta,rho){
  s<-qt(u,df=eta)
  r<-qt(v,df=eta)
  numerator<-gamma((eta+2)/2)*gamma(eta/2)*((1+s^2/eta)*(1+r^2/eta))^((eta+1)/2)
  denominator<-sqrt(1-rho^2)*(gamma((eta+1)/2))^2*(1+(s^2+r^2-2*rho*s*r)/(eta*(1-rho^2)))^((eta+2)/2)
  dens<-numerator/denominator
  return(dens)
}


## Frank copula
### alpha: copula parameter
frankdensity<-function(u,v,alpha){
  numerator<-alpha*(1-exp(-alpha))*exp(-alpha*(u+v))
  denominator<-(1-exp(-alpha)-(1-exp(-alpha*u))*(1-exp(-alpha*v)))^2
  dens<-numerator/denominator
  return(dens)
}

## Clayton copula
### alpha: copula parameter
claytondensity<-function(u,v,alpha){
  numerator<-(alpha+1)*(u*v)^alpha
  denominator<-(u^alpha+v^alpha-(u*v)^alpha)^(1/alpha+2)
  dens<-numerator/denominator
  return(dens)
}

## Joe copula
### alpha: copula parameter
joedensity<-function(u,v,alpha){
  w<-1-u
  z<-1-v
  dens<-(w^alpha+z^alpha-(w*z)^alpha)^(1/alpha-2)*(w*z)^(alpha-1)*(alpha-1+w^alpha+z^alpha-(w*z)^alpha)
  return(dens)
}

## Gumbel copula
### alpha: copula parameter
gumbeldensity<-function(u,v,alpha){
  x<--log(u)
  y<--log(v)
  numerator<-exp(-(x^alpha+y^alpha)^(1/alpha))*((x^alpha+y^alpha)^(1/alpha)+alpha-1)*(x^alpha+y^alpha)^(1/alpha-2)*(x*y)^(alpha-1)
  denominator<-u*v
  dens<-numerator/denominator
  return(dens)
}

## Inverted gumbel copula
### alpha: copula parameter
invergumbeldensity<-function(u,v,alpha){
  x<--log(1-u)
  y<--log(1-v)
  numerator<-exp(-(x^alpha+y^alpha)^(1/alpha))*((x^alpha+y^alpha)^(1/alpha)+alpha-1)*(x^alpha+y^alpha)^(1/alpha-2)*(x*y)^(alpha-1)
  denominator<-(1-u)*(1-v)
  dens<-numerator/denominator
  return(dens)
}

## Husler-Reiss copula
### alpha: copula parameter
huslerreiss<-function(u,v,alpha){
  x<-evd::qgev(u,loc=1,scale=1,shape=1)
  y<-evd::qgev(v,loc=1,scale=1,shape=1)
  jacob<-u*v*((log(u))^2)*((log(v))^2)
  w1<-1/x
  w2<-1/y
  hr<-function(x1,x2,model=pnorm){
    x1*model(1/alpha+alpha*log(x1/x2)/2)
  }
  v<-hr(w1,w2)+hr(w2,w1)
  hr2<-hr(w1,w2)*hr(w2,w1)+alpha*hr(w1,w2,model=dnorm)/2
  jac<-log(w1)+log(w2)
  dens<-exp(log(hr2)+jac-v)/jacob
  return(dens)
}

## Galambos copula
### alpha: copula parameter
galambos<-function(u,v,alpha){
  x<-evd::qgev(u,loc=1,scale=1,shape=1)
  y<-evd::qgev(v,loc=1,scale=1,shape=1)
  jacob<-u*v*((log(u))^2)*((log(v))^2)
  w1<-1/x
  w2<-1/y
  ialpha<-1/alpha
  z<-(w1^(-alpha)+w2^(-alpha))^(-ialpha)
  v<-w1+w2-z
  f1<-(-alpha-1)*log(w1)
  f2<-(-alpha-1)*log(w2)
  jac<-2*(log(w1)+log(w2))
  e1<-(1+alpha)*log(z)+log(exp(f1)+exp(f2))
  e2<-f1+f2+(1+2*alpha)*log(z)+log(1+alpha+z)
  dens<-exp(log(1-exp(e1)+exp(e2))-v+jac)/jacob
  return(dens)
}

## Coles-Tawn copula
### alpha, beta: copula parameters
colestawn<-function(u,v,alpha,beta){
  x<-evd::qgev(u,loc=1,scale=1,shape=1)
  y<-evd::qgev(v,loc=1,scale=1,shape=1)
  jacob<-u*v*((log(u))^2)*((log(v))^2)
  mar1<-mar2<-matrix(t(c(1,1,1)),nrow=length(u),ncol=3,byrow=T)
  w1<-1/x
  w2<-1/y
  uw<-(alpha*w2)/(alpha*w2+beta*w1)
  vw<-w2*pbeta(uw,shape1=alpha,shape2=beta+1)+w1*pbeta(uw,shape1=alpha+1,shape2=beta,lower.tail=F)
  jac<-(1+mar1[,3])*log(w1)+(1+mar2[,3])*log(w2)-log(mar1[,2]*mar2[,2])
  exp1<-pbeta(uw,shape1=alpha,shape2=beta+1)*pbeta(uw,shape1=alpha+1,shape2=beta,lower.tail=F)
  exp2<-dbeta(uw,shape1=alpha+1,shape2=beta+1)/(alpha*w2+beta*w1)
  d<-log(exp1+(alpha*beta/(alpha+beta+1))*exp2)-vw+jac
  dens<-exp(d)/jacob
  return(dens)
}

# --------------------------------------
# Not used in the paper but implemented
# --------------------------------------

## Bilogistic copula
### alpha, beta: copula parameters
bilogistic<-function(u,v,alpha,beta){
  x<-evd::qgev(u,loc=1,scale=1,shape=1)
  y<-evd::qgev(v,loc=1,scale=1,shape=1)
  jacob<-u*v*((log(u))^2)*((log(v))^2)
  w1<-1/x
  w2<-1/y
  gma<-c()
  for(i in 1:length(u)){
    gmafn<-function(z){
      (1-alpha)*w1[i]*(1-z)^beta-(1-beta)*w2[i]*z^alpha
    }
    gma[i]<-uniroot(gmafn,lower=0,upper=1,tol=.Machine$double.eps^0.5)$root
  }
  v<-w1*gma^(1-alpha)+w2*(1-gma)^(1-beta)
  jac<-2*(log(w1)+log(w2))
  e1<-exp((1-alpha)*log(gma)+(1-beta)*log(1-gma))
  e2<-exp(log(1-alpha)+log(beta)+(beta-1)*log(1-gma)+log(w1))+
    exp(log(1-beta)+log(alpha)+(alpha-1)*log(gma)+log(w2))
  dens<-exp(log(e1+(1-alpha)*(1-beta)/e2)-v+jac)/jacob
  return(dens)
}

## Negative bilogistic copula
### alpha, beta: copula parameters
negbilogistic<-function(u,v,alpha,beta){
  x<-evd::qgev(u,loc=1,scale=1,shape=1)
  y<-evd::qgev(v,loc=1,scale=1,shape=1)
  jacob<-u*v*((log(u))^2)*((log(v))^2)
  w1<-1/x
  w2<-1/y
  gma<-c()
  for(i in 1:length(u)){
    gmafn<-function(z){
      (1+alpha)*w1[i]*z^alpha-(1+beta)*w2[i]*(1-z)^beta
    }
    gma[i]<-uniroot(gmafn,lower=0,upper=1,tol=.Machine$double.eps^0.5)$root
  }
  v<-w1+w2-w1*gma^(1+alpha)-w2*(1-gma)^(1+beta)
  jac<-2*(log(w1)+log(w2))
  e1<-(1-gma^(1+alpha))*(1-(1-gma)^(1+beta))
  e2<-exp(log(1+alpha)+log(1+beta)+alpha*log(gma)+beta*log(1-gma))
  e3<-exp(log(1+alpha)+log(alpha)+(alpha-1)*log(gma)+log(w1))+
    exp(log(1+beta)+log(beta)+(beta-1)*log(1-gma)+log(w2))
  dens<-exp(log(e1+e2/e3)-v+jac)/jacob
  return(dens)
}

## Asymetric Mixed copula
### alpha, beta: copula parameters
asymix<-function(u,v,alpha,beta){
  x<-evd::qgev(u,loc=1,scale=1,shape=1)
  y<-evd::qgev(v,loc=1,scale=1,shape=1)
  jacob<-u*v*((log(u))^2)*((log(v))^2)
  w1<-1/x
  w2<-1/y
  sumw<-w1+w2
  v<-sumw-(alpha+beta)*w1+alpha*(w1^2)/sumw+beta*(w1^3)/(sumw^2)
  jac<-2*(log(w1)+log(w2))
  z1<-w1/sumw
  z2<-w2/sumw
  v1<-1-alpha*(z2^2)-beta*(3*z2^2-2*z2^3)
  v2<-1-alpha*(z1^2)-2*beta*z1^3
  v12<-(-2*alpha*z1*z2-6*beta*z1^2*z2)/sumw
  dens<-exp(log(v1*v2-v12)-v+jac)/jacob
  return(dens)
}

## Huser and Wadsworth (2019)
### equation (9) in the paper: quantile function
### delta: copula parameter that represents the extremal dependence of the model
### NOTE: this model is implemented but it is very computational expensive and might not give results in time
quanty<-function(q,delta){
  dummyd<-function(y){
    x<-exp(y)
    (1-q)-(delta/(2*delta-1))*x^(-1/delta)+((1-delta)/(2*delta-1))*x^(-1/(1-delta))
  }
  dummy0.5<-function(y){
    x<-exp(y)
    (1-q)-x^(-2)*(2*log(x)+1)
  }
  if(abs(delta-0.5)>1e-05){
    ur<-uniroot(dummyd,interval=c(0,100),extendInt="yes")
  }
  else{
    ur<-uniroot(dummy0.5,interval=c(0,100),extendInt="yes")
  }
  return(ur$root)
}
quanty<-Vectorize(quanty,"q")

### W copula: gaussian
### thetaW: parameter of the gaussian copula
jgaus<-function(s,thetaW){
  z<-qnorm(1-exp(-s))
  if(any(z==Inf)){
    z<- -qnorm(exp(-s))
  }
  inside<--(1/(2*(1-thetaW^2)))*(z[1]^2-2*thetaW*z[1]*z[2]+z[2]^2)
  expinside<-exp(inside)
  dens<-(1/(2*pi*sqrt(1-thetaW^2)))*expinside*prod(exp(-s)/dnorm(z))
  return(dens)
}

### after equation (16) in the paper: joint pdf
jointdist<-function(y,delta,thetaW){
  int<-function(v){
    s<-c(y[1]-delta*v,y[2]-delta*v)/(1-delta)
    out<-jgaus(s=s,thetaW=thetaW)*exp(-v)
    return(out)
  }
  intvect<-Vectorize(int)
  integ<-integrate(intvect,lower=0,upper=min(y)/delta,abs.tol=0)$val
  dens<-integ*(1-delta)^(-2)
  return(dens)
}

### marginal pdf: differentiation of equation (9)
marg<-function(y,delta){
  if(abs(delta-0.5)>1e-05){
    dens<-((1/(2*delta-1))*exp(y)^(-1/delta-1)-(1/(2*delta-1))*exp(y)^(-1/(1-delta)-1))*exp(y)
    return(dens)
  }
  else{
    dens<-((exp(y)^(-3)*(4*y)))*exp(y)
    return(dens)
  }
}

### final density
cop2019<-function(u,v,delta,thetaW){
  y<-matrix(NA,nrow=2,ncol=length(u))
  y[1,]<-apply(as.matrix(u),1,quanty,delta=delta)
  y[2,]<-apply(as.matrix(v),1,quanty,delta=delta)
  numerator<-apply(y,2,jointdist,delta=delta,thetaW=thetaW)
  denominator<-sapply(y[1,],marg,delta=delta)*sapply(y[2,],marg,delta=delta)
  dens<-numerator/denominator
  return(dens)
}

# ===============
# Functions used
# ===============

## Weighting Functions
pifun<-function(u,v,theta){
  return((u*v)^theta)
}

pifun1<-function(u,v,theta){
  return((1/exp((1-u)*(1-v)))^theta)
}

## Identification of the copulas used for each component of the model
copmodel<-function(u,v,par,copula){
  if(copula=="t"){
    rho<-par[1]
    eta<-par[2]
    density<-tdensity(u=u,v=v,rho=rho,eta=eta)
  }
  else if(copula=="gaus"){
    rho<-par[1]
    density<-gausdensity(u=u,v=v,rho=rho)
  }
  else if(copula=="joe"){
    alpha<-par[1]
    density<-joedensity(u=u,v=v,alpha=alpha)
  }
  else if(copula=="frank"){
    alpha<-par[1]
    density<-frankdensity(u=u,v=v,alpha=alpha)
  }
  else if(copula=="gumbel"){
    alpha<-par[1]
    density<-gumbeldensity(u=u,v=v,alpha=alpha)
  }
  else if(copula=="InverGumbel"){
    alpha<-par[1]
    density<-invergumbeldensity(u=u,v=v,alpha=alpha)
  }
  else if(copula=="cop19"){
    delta<-par[1]
    thetaW<-par[2]
    density<-cop2019(u=u,v=v,delta=delta,thetaW=thetaW)
  }
  else if(copula=="clayton"){
    alpha<-par[1]
    density<-claytondensity(u=u,v=v,alpha=alpha)
  }
  else if(copula=="HR"){
    alpha<-par[1]
    density<-huslerreiss(u=u,v=v,alpha=alpha)
  }
  else if(copula=="galambos"){
    alpha<-par[1]
    density<-galambos(u=u,v=v,alpha=alpha)
  }
  else if(copula=="bilogistic"){
    alpha<-par[1]
    beta<-par[2]
    density<-bilogistic(u=u,v=v,alpha=alpha,beta=beta)
  }
  else if(copula=="negbilogistic"){
    alpha<-par[1]
    beta<-par[2]
    density<-negbilogistic(u=u,v=v,alpha=alpha,beta=beta)
  }
  else if(copula=="ct"){
    alpha<-par[1]
    beta<-par[2]
    density<-colestawn(u=u,v=v,alpha=alpha,beta=beta)
  }
  else if(copula=="AsyMix"){
    alpha<-par[1]
    beta<-par[2]
    density<-asymix(u=u,v=v,alpha=alpha,beta=beta)
  }
  else{
    warning("This copula is not implemented")
  }
  return(density)
}

## c^* without the normalising constant K(gamma)
### ct: copula density tailored to the tail
### cb: copula density tailored to the body
### parct: paramter(s) of ct
### parcb: paramter(s) of cb
### weightfun: weighting function (has to be written as a function)
### parweight: parameter of the weighting function
cstarfun<-function(u,v,parct,parcb,ct,cb,weightfun,parweight){
  ctdensity<-copmodel(u=u,v=v,par=parct,copula=ct)
  cbdensity<-copmodel(u=u,v=v,par=parcb,copula=cb)
  dens<-weightfun(u,v,parweight)*ctdensity+(1-weightfun(u,v,parweight))*cbdensity
  return(dens)
}
cstarfun<-Vectorize(cstarfun,vectorize.args=c("u","v"))

## normalising constant K(gamma)
K<-function(parct,parcb,ct,cb,weightfun,parweight,singular=F){
  if(ct=="gumbel" || ct=="t" || ct=="gaus" || ct=="joe" || ct=="HR" || ct=="galambos" || ct=="ct" || cb=="gumbel" || cb=="t" || cb=="gaus" || cb=="joe" || cb=="HR" || cb=="galambos"){
    a<-integral2(cstarfun,0,1,0,1,parct=parct,parcb=parcb,ct=ct,cb=cb,weightfun=weightfun,parweight=parweight,singular=T)
    return(a$Q)
  }
  a<-integral2(cstarfun,0,1,0,1,parct=parct,parcb=parcb,ct=ct,cb=cb,weightfun=weightfun,parweight=parweight)
  return(a$Q)
}

## c^* density
### K: normalising constant
fun<-function(u,v,parct,parcb,ct,cb,weightfun,parweight,K){
  dens<-(1/K)*cstarfun(u,v,parct,parcb,ct,cb,weightfun,parweight)
  return(dens)
}

## marginal cdf F
marginal<-function(x,parct,parcb,ct,cb,weightfun,parweight,K){
  integral<-c()
  for(i in 1:length(x)){
    if(x[i]<0.7){
      integral[i]<-tryCatch(expr=integral2(fun,0.000001,x[i],0,1,parct=parct,parcb=parcb,ct=ct,
                                           cb=cb,weightfun=weightfun,parweight=parweight,K=K,vectorized=T,reltol=0.001,singular=T)$Q,
                            error=function(e){message("Numerical difficulties occurred in the cdf calculations \n NA returned to integration");return(NA)})
      
    }
    else{
      integral[i]<-tryCatch(expr=integral2(fun,0.000001,x[i],0,1,parct=parct,parcb=parcb,
                                           ct=ct,cb=cb,weightfun=weightfun,parweight=parweight,K=K,vectorized=T,singular=T)$Q,
                            error=function(e){message("Numerical difficulties occurred in the cdf calculations \n NA returned to integration");return(NA)})
    }
  }
  return(integral)
}

## marginal pdf f
### the arguments of parct, parcb and weightfun have to be lists
pdfun<-function(u,parct,parcb,ct,cb,weightfun,parweight,K){
  fun1<-Vectorize(fun,"v")
  fun2<-function(v,parct,parcb,ct,cb,weightfun,parweight,K){
    fun1(u,v,parct=parct,ct=ct,parcb=parcb,cb=cb,weightfun=weightfun,parweight=parweight,K=K)
  }
  int<-tryCatch(expr=integrate(fun2,0,1,parct=parct,ct=ct,parcb=parcb,cb=cb,weightfun=weightfun,parweight=parweight,K=K)$value,
                error=function(e){message("Numerical difficulties occurred in the pdf calculations \n NA returned to integration");return(NA)})
  return(int)
}
pdfun<-Vectorize(pdfun)

## inverse of the marginal cdf F^{-1}
quant<-function(u,parct,parcb,ct,cb,weightfun,parweight,K){
  n<-50
  x<-seq(0+.Machine$double.eps^.25,1,len=n)
  y<-marginal(x,parct=parct,parcb=parcb,ct=ct,cb=cb,weightfun=weightfun,parweight=parweight,K=K)
  f<-splinefun(y,x)
  return(f(u))
}

## copula density 
cop<-function(u,v,parct,parcb,ct,cb,weightfun,parweight,K){
  ustar<-quant(u=u,parct=parct,parcb=parcb,ct=ct,cb=cb,weightfun=weightfun,parweight=parweight,K=K)
  vstar<-quant(u=v,parct=parct,parcb=parcb,ct=ct,cb=cb,weightfun=weightfun,parweight=parweight,K=K)
  numerator<-fun(ustar,vstar,parct=parct,parcb=parcb,ct=ct,cb=cb,weightfun=weightfun,parweight=parweight,K=K)
  denominator<-pdfun(ustar,parct=list(parct),parcb=list(parcb),ct=ct,cb=cb,weightfun=list(weightfun),parweight=parweight,K=K)*pdfun(vstar,parct=list(parct),parcb=list(parcb),ct=ct,cb=cb,weightfun=list(weightfun),parweight=parweight,K=K)
  dens<-numerator/denominator
  return(dens)
}

## log-likelihood
### par: vector of model parameters where the parameter(s) for ct are the first to be included, followed
### by the parameter(s) for cb, followed by the parameter for the weighting function
loglik<-function(x,ct,cb,weightfun,par){
  u<-x[,1]
  v<-x[,2]
  if(ct=="t"){
    parct<-c(par[1],par[2])
    if(cb=="t"){
      parcb<-c(par[3],par[4])
      parweight<-par[5]
      if((par[1]<=-1 || par[1] >=1) | (par[2]<=1 || par[2] >=15) | (par[3]<=-1 || par[3] >=1) | (par[4]<=1 || par[4] >=15) | par[5]<=0){return(1000e10)}
    }
    else if(cb=="cop2019"){
      parcb<-c(par[3],par[4])
      parweight<-par[5]
      if((par[1]<=-1 || par[1] >=1) | (par[2]<=1 || par[2] >=15) | (par[3]<=0 || par[3] >=1) | (par[4]<=-1 || par[4] >=1) | par[5]<=0){return(1000e10)}
    }
    else if(cb=="bilogistic"){
      parcb<-c(par[3],par[4])
      parweight<-par[5]
      if((par[1]<=-1 || par[1] >=1) | (par[2]<=1 || par[2] >=15) | (par[3]<=0 || par[3] >=1) | (par[4]<=0 || par[4] >=1) | par[5]<=0){return(1000e10)}
    }
    else if(cb=="negbilogistic" || cb=="ct"){
      parcb<-c(par[3],par[4])
      parweight<-par[5]
      if((par[1]<=-1 || par[1] >=1) | (par[2]<=1 || par[2] >=15) | par[3]<=0 | par[4]<=0 | par[5]<=0){return(1000e10)}
    }
    else if(cb=="AsyMix"){
      parcb<-c(par[3],par[4])
      parweight<-par[5]
      if((par[1]<=-1 || par[1] >=1) | (par[2]<=1 || par[2] >=15) | (par[3]<0 || par[3] >1.5) | (par[4]< -0.5 || par[4] >0.5) | par[5]<=0 | par[3]+3*par[4]<0 | par[3]+2*par[4]>1){return(1000e10)}
    }
    else{
      parcb<-c(par[3])
      if(cb=="gaus"){
        if((par[1]<=-1 || par[1] >=1) | (par[2]<=1 || par[2] >=15) | (par[3]<=-1 || par[3] >=1)){return(1000e10)}
      }
      else if(cb=="joe" | cb=="gumbel" | cb=="InverGumbel"){
        if((par[1]<=-1 || par[1] >=1) | (par[2]<=1 || par[2] >=10) | (par[3]<1 || par[3]>=15)){return(1000e10)}
      }
      else if(cb=="clayton"){
        if((par[1]<=-1 || par[1] >=1) | (par[2]<=1 || par[2] >=15) | (par[3]<0 || par[3]>=15)){return(1000e10)}
      }
      else if(cb=="frank"){
        if((par[1]<=-1 || par[1] >=1) | (par[2]<=1 || par[2] >=15) | par[3]==0){return(1000e10)}
      }
      else if(cb=="HR" || cb=="galambos"){
        if((par[1]<=-1 || par[1] >=1) | (par[2]<=1 || par[2] >=15) | par[3]<=0){return(1000e10)}
      }
      parweight<-par[4]
      if(par[4]<=0){return(1000e10)}
    }
  }
  else if(ct=="cop2019"){
    parct<-c(par[1],par[2])
    if(cb=="t"){
      parcb<-c(par[3],par[4])
      parweight<-par[5]
      if((par[1]<=0 || par[1] >=1) | (par[2]<=-1 || par[2] >=1) | (par[3]<=-1 || par[3] >=1) | (par[4]<=1 || par[4] >=15) | par[5]<=0){return(1000e10)}
    }
    else if(cb=="cop2019"){
      parcb<-c(par[3],par[4])
      parweight<-par[5]
      if((par[1]<=0 || par[1] >=1) | (par[2]<=-1 || par[2] >=1) | (par[3]<=0 || par[3] >=1) | (par[4]<=-1 || par[4] >=1) | par[5]<=0){return(1000e10)}
    }
    else if(cb=="bilogistic"){
      parcb<-c(par[3],par[4])
      parweight<-par[5]
      if((par[1]<=0 || par[1] >=1) | (par[2]<=-1 || par[2] >=1) | (par[3]<=0 || par[3] >=1) | (par[4]<=0 || par[4] >=1) | par[5]<=0){return(1000e10)}
    }
    else if(cb=="negbilogistic" || cb=="ct"){
      parcb<-c(par[3],par[4])
      parweight<-par[5]
      if((par[1]<=0 || par[1] >=1) | (par[2]<=-1 || par[2] >=1) | par[3]<=0 | par[4]<=0 | par[5]<=0){return(1000e10)}
    }
    else if(cb=="AsyMix"){
      parcb<-c(par[3],par[4])
      parweight<-par[5]
      if((par[1]<=0 || par[1] >=1) | (par[2]<=-1 || par[2] >=1) | (par[3]<0 || par[3] >1.5) | (par[4]< -0.5 || par[4] >0.5) | par[5]<=0 | par[3]+3*par[4]<0 | par[3]+2*par[4]>1){return(1000e10)}
    }
    else{
      parcb<-c(par[3])
      if(cb=="gaus"){
        if((par[1]<=0 || par[1] >=1) | (par[2]<=-1 || par[2] >=1) | (par[3]<=-1 || par[3] >=1)){return(1000e10)}
      }
      else if(cb=="joe" | cb=="gumbel" | cb=="InverGumbel"){
        if((par[1]<=0 || par[1] >=1) | (par[2]<=-1 || par[2] >=1) | (par[3]<1 || par[3]>=10)){return(1000e10)}
      }
      else if(cb=="clayton"){
        if((par[1]<=0 || par[1] >=1) | (par[2]<=-1 || par[2] >=1) | (par[3]<0 || par[3]>=15)){return(1000e10)}
      }
      else if(cb=="frank"){
        if((par[1]<=0 || par[1] >=1) | (par[2]<=-1 || par[2] >=1) | par[3]==0){return(1000e10)}
      }
      else if(cb=="HR" || cb=="galambos"){
        if((par[1]<=0 || par[1] >=1) | (par[2]<=-1 || par[2] >=1) | par[3]<=0){return(1000e10)}
      }
      parweight<-par[4]
      if(par[4]<=0){return(1000e10)}
    }
  }
  else if(ct=="bilogistic"){
    parct<-c(par[1],par[2])
    if(cb=="t"){
      parcb<-c(par[3],par[4])
      parweight<-par[5]
      if((par[1]<=0 || par[1] >=1) | (par[2]<=0 || par[2] >=1) | (par[3]<=-1 || par[3] >=1) | (par[4]<=1 || par[4] >=15) | par[5]<=0){return(1000e10)}
    }
    else if(cb=="cop2019"){
      parcb<-c(par[3],par[4])
      parweight<-par[5]
      if((par[1]<=0 || par[1] >=1) | (par[2]<=0 || par[2] >=1) | (par[3]<=0 || par[3] >=1) | (par[4]<=-1 || par[4] >=1) | par[5]<=0){return(1000e10)}
    }
    else if(cb=="bilogistic"){
      parcb<-c(par[3],par[4])
      parweight<-par[5]
      if((par[1]<=0 || par[1] >=1) | (par[2]<=0 || par[2] >=1) | (par[3]<=0 || par[3] >=1) | (par[4]<=0 || par[4] >=1) | par[5]<=0){return(1000e10)}
    }
    else if(cb=="negbilogistic" || cb=="ct"){
      parcb<-c(par[3],par[4])
      parweight<-par[5]
      if((par[1]<=0 || par[1] >=1) | (par[2]<=0 || par[2] >=1) | par[3]<=0 | par[4]<=0 | par[5]<=0){return(1000e10)}
    }
    else if(cb=="AsyMix"){
      parcb<-c(par[3],par[4])
      parweight<-par[5]
      if((par[1]<=0 || par[1] >=1) | (par[2]<=0 || par[2] >=1) | (par[3]<0 || par[3] >1.5) | (par[4]< -0.5 || par[4] >0.5) | par[5]<=0 | par[3]+3*par[4]<0 | par[3]+2*par[4]>1){return(1000e10)}
    }
    else{
      parcb<-c(par[3])
      if(cb=="gaus"){
        if((par[1]<=0 || par[1] >=1) | (par[2]<=0 || par[2] >=1) | (par[3]<=-1 || par[3] >=1)){return(1000e10)}
      }
      else if(cb=="joe" | cb=="gumbel" | cb=="InverGumbel"){
        if((par[1]<=0 || par[1] >=1) | (par[2]<=0 || par[2] >=1) | (par[3]<1 || par[3]>=10)){return(1000e10)}
      }
      else if(cb=="clayton"){
        if((par[1]<=0 || par[1] >=1) | (par[2]<=0 || par[2] >=1) | (par[3]<0 || par[3]>=15)){return(1000e10)}
      }
      else if(cb=="frank"){
        if((par[1]<=0 || par[1] >=1) | (par[2]<=0 || par[2] >=1) | par[3]==0){return(1000e10)}
      }
      else if(cb=="HR" || cb=="galambos"){
        if((par[1]<=0 || par[1] >=1) | (par[2]<=0 || par[2] >=1) | par[3]<=0){return(1000e10)}
      }
      parweight<-par[4]
      if(par[4]<=0){return(1000e10)}
    }
  }
  else if(ct=="negilogistic" || ct=="ct"){
    parct<-c(par[1],par[2])
    if(cb=="t"){
      parcb<-c(par[3],par[4])
      parweight<-par[5]
      if(par[1]<=0 | par[2]<=0 | (par[3]<=-1 || par[3] >=1) | (par[4]<=1 || par[4] >=15) | par[5]<=0){return(1000e10)}
    }
    else if(cb=="cop2019"){
      parcb<-c(par[3],par[4])
      parweight<-par[5]
      if(par[1]<=0 | par[2]<=0 | (par[3]<=0 || par[3] >=1) | (par[4]<=-1 || par[4] >=1) | par[5]<=0){return(1000e10)}
    }
    else if(cb=="bilogistic"){
      parcb<-c(par[3],par[4])
      parweight<-par[5]
      if(par[1]<=0 | par[2]<=0 | (par[3]<=0 || par[3] >=1) | (par[4]<=0 || par[4] >=1) | par[5]<=0){return(1000e10)}
    }
    else if(cb=="negbilogistic" || cb=="ct"){
      parcb<-c(par[3],par[4])
      parweight<-par[5]
      if(par[1]<=0 | par[2]<=0 | par[3]<=0 | par[4]<=0 | par[5]<=0){return(1000e10)}
    }
    else if(cb=="AsyMix"){
      parcb<-c(par[3],par[4])
      parweight<-par[5]
      if(par[1]<=0 | par[2]<=0 | (par[3]<0 || par[3] >1.5) | (par[4]< -0.5 || par[4] >0.5) | par[5]<=0 | par[3]+3*par[4]<0 | par[3]+2*par[4]>1){return(1000e10)}
    }
    else{
      parcb<-c(par[3])
      if(cb=="gaus"){
        if(par[1]<=0 | par[2]<=0 | (par[3]<=-1 || par[3] >=1)){return(1000e10)}
      }
      else if(cb=="joe" | cb=="gumbel" | cb=="InverGumbel"){
        if(par[1]<=0 | par[2]<=0 | (par[3]<1 || par[3]>=10)){return(1000e10)}
      }
      else if(cb=="clayton"){
        if(par[1]<=0 | par[2]<=0 | (par[3]<0 || par[3]>=15)){return(1000e10)}
      }
      else if(cb=="frank"){
        if(par[1]<=0 | par[2]<=0 | par[3]==0){return(1000e10)}
      }
      else if(cb=="HR" || cb=="galambos"){
        if(par[1]<=0 | par[2]<=0 | par[3]<=0){return(1000e10)}
      }
      parweight<-par[4]
      if(par[4]<=0){return(1000e10)}
    }
  }
  else if(ct=="AsyMix"){
    parct<-c(par[1],par[2])
    if(cb=="t"){
      parcb<-c(par[3],par[4])
      parweight<-par[5]
      if((par[1]<0 || par[1] >1.5) | (par[2]< -0.5 || par[2] > 0.5) | (par[3]<=-1 || par[3] >=1) | (par[4]<=1 || par[4] >=15) | par[5]<=0 | par[1]+3*par[2]<0 | par[1]+2*par[2]>1){return(1000e10)}
    }
    else if(cb=="cop2019"){
      parcb<-c(par[3],par[4])
      parweight<-par[5]
      if((par[1]<0 || par[1] >1.5) | (par[2]< -0.5 || par[2] > 0.5) | (par[3]<=0 || par[3] >=1) | (par[4]<=-1 || par[4] >=1) | par[5]<=0 | par[1]+3*par[2]<0 | par[1]+2*par[2]>1){return(1000e10)}
    }
    else if(cb=="bilogistic"){
      parcb<-c(par[3],par[4])
      parweight<-par[5]
      if((par[1]<0 || par[1] >1.5) | (par[2]< -0.5 || par[2] > 0.5) | (par[3]<=0 || par[3] >=1) | (par[4]<=0 || par[4] >=1) | par[5]<=0 | par[1]+3*par[2]<0 | par[1]+2*par[2]>1){return(1000e10)}
    }
    else if(cb=="negbilogistic" || cb=="ct"){
      parcb<-c(par[3],par[4])
      parweight<-par[5]
      if((par[1]<0 || par[1] >1.5) | (par[2]< -0.5 || par[2] > 0.5) | par[3]<=0 | par[4]<=0 | par[5]<=0 | par[1]+3*par[2]<0 | par[1]+2*par[2]>1){return(1000e10)}
    }
    else if(cb=="AsyMix"){
      parcb<-c(par[3],par[4])
      parweight<-par[5]
      if((par[1]<0 || par[1] >1.5) | (par[2]< -0.5 || par[2] > 0.5) | (par[3]<0 || par[3] >1.5) | (par[4]< -0.5 || par[4] >0.5) | par[5]<=0 | par[1]+3*par[2]<0 | par[1]+2*par[2]>1 | par[3]+3*par[4]<0 | par[3]+2*par[4]>1){return(1000e10)}
    }
    else{
      parcb<-c(par[3])
      if(cb=="gaus"){
        if((par[1]<0 || par[1] >1.5) | (par[2]< -0.5 || par[2] > 0.5) | (par[3]<=-1 || par[3] >=1) | par[1]+3*par[2]<0 | par[1]+2*par[2]>1){return(1000e10)}
      }
      else if(cb=="joe" | cb=="gumbel" | cb=="InverGumbel"){
        if((par[1]<0 || par[1] >1.5) | (par[2]< -0.5 || par[2] > 0.5) | (par[3]<1 || par[3]>=10) | par[1]+3*par[2]<0 | par[1]+2*par[2]>1){return(1000e10)}
      }
      else if(cb=="clayton"){
        if((par[1]<0 || par[1] >1.5) | (par[2]< -0.5 || par[2] > 0.5) | (par[3]<0 || par[3]>=15) | par[1]+3*par[2]<0 | par[1]+2*par[2]>1){return(1000e10)}
      }
      else if(cb=="frank"){
        if((par[1]<0 || par[1] >1.5) | (par[2]< -0.5 || par[2] > 0.5) | par[3]==0 | par[1]+3*par[2]<0 | par[1]+2*par[2]>1){return(1000e10)}
      }
      else if(cb=="HR" || cb=="galambos"){
        if((par[1]<0 || par[1] >1.5) | (par[2]< -0.5 || par[2] > 0.5) | par[3]<=0 | par[1]+3*par[2]<0 | par[1]+2*par[2]>1){return(1000e10)}
      }
      parweight<-par[4]
      if(par[4]<=0){return(1000e10)}
    }
  }
  else if(ct=="gaus"){
    parct<-par[1]
    if(cb=="t"){
      parcb<-c(par[2],par[3])
      parweight<-par[4]
      if((par[1]<=-1 || par[1]>=1) | (par[2]<=-1 || par[2] >=1) | (par[3]<=1 || par[3] >=15) | par[4]<=0){return(1000e10)}
    }
    else if(cb=="cop2019"){
      parcb<-c(par[2],par[3])
      parweight<-par[4]
      if((par[1]<=-1 || par[1]>=1) | (par[2]<=0 || par[2] >=1) | (par[3]<=-1 || par[3] >=1) | par[4]<=0){return(1000e10)}
    }
    else if(cb=="bilogistic"){
      parcb<-c(par[2],par[3])
      parweight<-par[4]
      if((par[1]<=-1 || par[1]>=1) | (par[2]<=0 || par[2] >=1) | (par[3]<=0 || par[3] >=1) | par[4]<=0){return(1000e10)}
    }
    else if(cb=="negbilogistic" || cb=="ct"){
      parcb<-c(par[2],par[3])
      parweight<-par[4]
      if((par[1]<=-1 || par[1]>=1) | par[2]<=0 | par[3]<=0 | par[4]<=0){return(1000e10)}
    }
    else if(cb=="AsyMix"){
      parcb<-c(par[2],par[3])
      parweight<-par[4]
      if((par[1]<=-1 || par[1]>=1) | (par[2]<0 || par[2] >1.5) | (par[3]< -0.5 || par[3] >0.5) | par[4]<=0 | par[2]+3*par[3]<0 | par[2]+2*par[3]>1){return(1000e10)}
    }
    else{
      parcb<-c(par[2])
      if(cb=="gaus"){
        if((par[1]<=-1 || par[1]>=1) | (par[2]<=-1 || par[2] >=1)){return(1000e10)}
      }
      else if(cb=="joe" | cb=="gumbel" | cb=="InverGumbel"){
        if((par[1]<=-1 || par[1]>=1) | (par[2]<1 || par[2]>=10)){return(1000e10)}
      }
      else if(cb=="clayton"){
        if((par[1]<=-1 || par[1]>=1) | (par[2]<0 || par[2]>=15)){return(1000e10)}
      }
      else if(cb=="frank"){
        if((par[1]<=-1 || par[1]>=1) | par[2]==0){return(1000e10)}
      }
      else if(cb=="HR" || cb=="galambos"){
        if((par[1]<=-1 || par[1]>=1) | par[2]<=0){return(1000e10)}
      }
      parweight<-par[3]
      if(par[3]<=0){return(1000e10)}
    }
  }
  else if(ct=="joe" || ct=="gumbel" || ct=="InverGumbel"){
    parct<-par[1]
    if(cb=="t"){
      parcb<-c(par[2],par[3])
      parweight<-par[4]
      if((par[1]<=1 || par[1]>=10) | (par[2]<=-1 || par[2] >=1) | (par[3]<=1 || par[3] >=15) | par[4]<=0){return(1000e10)}
    }
    else if(cb=="cop2019"){
      parcb<-c(par[2],par[3])
      parweight<-par[4]
      if((par[1]<=1 || par[1]>=10) | (par[2]<=0 || par[2] >=1) | (par[3]<=-1 || par[3] >=1) | par[4]<=0){return(1000e10)}
    }
    else if(cb=="bilogistic"){
      parcb<-c(par[2],par[3])
      parweight<-par[4]
      if((par[1]<=1 || par[1]>=10) | (par[2]<=0 || par[2] >=1) | (par[3]<=0 || par[3] >=1) | par[4]<=0){return(1000e10)}
    }
    else if(cb=="negbilogistic" || cb=="ct"){
      parcb<-c(par[2],par[3])
      parweight<-par[4]
      if((par[1]<=1 || par[1]>=10) | par[2]<=0 | par[3]<=0 | par[4]<=0){return(1000e10)}
    }
    else if(cb=="AsyMix"){
      parcb<-c(par[2],par[3])
      parweight<-par[4]
      if((par[1]<=1 || par[1]>=10) | (par[2]<0 || par[2] >1.5) | (par[3]< -0.5 || par[3] >0.5) | par[4]<=0 | par[2]+3*par[3]<0 | par[2]+2*par[3]>1){return(1000e10)}
    }
    else{
      parcb<-c(par[2])
      if(cb=="gaus"){
        if((par[1]<=1 || par[1]>=10) | (par[2]<=-1 || par[2] >=1)){return(1000e10)}
      }
      else if(cb=="joe" | cb=="gumbel" | cb=="InverGumbel"){
        if((par[1]<=1 || par[1]>=10) | (par[2]<1 || par[2]>=10)){return(1000e10)}
      }
      else if(cb=="clayton"){
        if((par[1]<=1 || par[1]>=15) | (par[2]<0 || par[2]>=15)){return(1000e10)}
      }
      else if(cb=="frank"){
        if((par[1]<=1 || par[1]>=15) | par[2]==0){return(1000e10)}
      }
      else if(cb=="HR" || cb=="galambos"){
        if((par[1]<=1 || par[1]>=15) | par[2]<=0){return(1000e10)}
      }
      parweight<-par[3]
      if(par[3]<=0){return(1000e10)}
    }
  }
  else if(ct=="frank"){
    parct<-par[1]
    if(cb=="t"){
      parcb<-c(par[2],par[3])
      parweight<-par[4]
      if(par[1]==0 | (par[2]<=-1 || par[2] >=1) | (par[3]<=1 || par[3] >=15) | par[4]<=0){return(1000e10)}
    }
    else if(cb=="cop2019"){
      parcb<-c(par[2],par[3])
      parweight<-par[4]
      if(par[1]==0 | (par[2]<=0 || par[2] >=1) | (par[3]<=-1 || par[3] >=1) | par[4]<=0){return(1000e10)}
    }
    else if(cb=="bilogistic"){
      parcb<-c(par[2],par[3])
      parweight<-par[4]
      if(par[1]==0 | (par[2]<=0 || par[2] >=1) | (par[3]<=0 || par[3] >=1) | par[4]<=0){return(1000e10)}
    }
    else if(cb=="negbilogistic" || cb=="ct"){
      parcb<-c(par[2],par[3])
      parweight<-par[4]
      if(par[1]==0 | par[2]<=0 | par[3]<=0 | par[4]<=0){return(1000e10)}
    }
    else if(cb=="AsyMix"){
      parcb<-c(par[2],par[3])
      parweight<-par[4]
      if(par[1]==0 | (par[2]<0 || par[2] >1.5) | (par[3]< -0.5 || par[3] >0.5) | par[4]<=0 | par[2]+3*par[3]<0 | par[2]+2*par[3]>1){return(1000e10)}
    }
    else{
      parcb<-c(par[2])
      if(cb=="gaus"){
        if(par[1]==0 | (par[2]<=-1 || par[2] >=1)){return(1000e10)}
      }
      else if(cb=="joe" | cb=="gumbel" | cb=="InverGumbel"){
        if(par[1]==0 | (par[2]<1 || par[2]>=10)){return(1000e10)}
      }
      else if(cb=="clayton"){
        if(par[1]==0 | (par[2]<0 || par[2]>=15)){return(1000e10)}
      }
      else if(cb=="frank"){
        if(par[1]==0 | par[2]==0){return(1000e10)}
      }
      else if(cb=="HR" || cb=="galambos"){
        if(par[1]==0 | par[2]<=0){return(1000e10)}
      }
      parweight<-par[3]
      if(par[3]<=0){return(1000e10)}
    }
  }
  else if(ct=="clayton"){
    parct<-par[1]
    if(cb=="t"){
      parcb<-c(par[2],par[3])
      parweight<-par[4]
      if((par[1]<=0 || par[1]>=15) | (par[2]<=-1 || par[2] >=1) | (par[3]<=1 || par[3] >=15) | par[4]<=0){return(1000e10)}
    }
    else if(cb=="cop2019"){
      parcb<-c(par[2],par[3])
      parweight<-par[4]
      if((par[1]<=0 || par[1]>=15) | (par[2]<=0 || par[2] >=1) | (par[3]<=-1 || par[3] >=1) | par[4]<=0){return(1000e10)}
    }
    else if(cb=="bilogistic"){
      parcb<-c(par[2],par[3])
      parweight<-par[4]
      if((par[1]<=0 || par[1]>=15) | (par[2]<=0 || par[2] >=1) | (par[3]<=0 || par[3] >=1) | par[4]<=0){return(1000e10)}
    }
    else if(cb=="negbilogistic" || cb=="ct"){
      parcb<-c(par[2],par[3])
      parweight<-par[4]
      if((par[1]<=0 || par[1]>=15) | par[2]<=0 | par[3]<=0 | par[4]<=0){return(1000e10)}
    }
    else if(cb=="AsyMix"){
      parcb<-c(par[2],par[3])
      parweight<-par[4]
      if((par[1]<=0 || par[1]>=15) | (par[2]<0 || par[2] >1.5) | (par[3]< -0.5 || par[3] >0.5) | par[4]<=0 | par[2]+3*par[3]<0 | par[2]+2*par[3]>1){return(1000e10)}
    }
    else{
      parcb<-c(par[2])
      if(cb=="gaus"){
        if((par[1]<=0 || par[1]>=15) | (par[2]<=-1 || par[2] >=1)){return(1000e10)}
      }
      else if(cb=="joe" | cb=="gumbel" | cb=="InverGumbel"){
        if((par[1]<=0 || par[1]>=15) | (par[2]<1 || par[2]>=10)){return(1000e10)}
      }
      else if(cb=="clayton"){
        if((par[1]<=0 || par[1]>=15) | (par[2]<0 || par[2]>=15)){return(1000e10)}
      }
      else if(cb=="frank"){
        if((par[1]<=0 || par[1]>=15) | par[2]==0){return(1000e10)}
      }
      else if(cb=="HR" || cb=="galambos"){
        if((par[1]<=0 || par[1]>=15) | par[2]<=0){return(1000e10)}
      }
      parweight<-par[3]
      if(par[3]<=0){return(1000e10)}
    }
  }
  else if(ct=="HR" || ct=="galambos"){
    parct<-par[1]
    if(cb=="t"){
      parcb<-c(par[2],par[3])
      parweight<-par[4]
      if(par[1]<=0 | (par[2]<=-1 || par[2] >=1) | (par[3]<=1 || par[3] >=15) | par[4]<=0){return(1000e10)}
    }
    else if(cb=="cop2019"){
      parcb<-c(par[2],par[3])
      parweight<-par[4]
      if(par[1]<=0 | (par[2]<=0 || par[2] >=1) | (par[3]<=-1 || par[3] >=1) | par[4]<=0){return(1000e10)}
    }
    else if(cb=="bilogistic"){
      parcb<-c(par[2],par[3])
      parweight<-par[4]
      if(par[1]<=0 | (par[2]<=0 || par[2] >=1) | (par[3]<=0 || par[3] >=1) | par[4]<=0){return(1000e10)}
    }
    else if(cb=="negbilogistic" || cb=="ct"){
      parcb<-c(par[2],par[3])
      parweight<-par[4]
      if(par[1]<=0 | par[2]<=0 | par[3]<=0 | par[4]<=0){return(1000e10)}
    }
    else if(cb=="AsyMix"){
      parcb<-c(par[2],par[3])
      parweight<-par[4]
      if(par[1]<=0 | (par[2]<0 || par[2] >1.5) | (par[3]< -0.5 || par[3] >0.5) | par[4]<=0 | par[2]+3*par[3]<0 | par[2]+2*par[3]>1){return(1000e10)}
    }
    else{
      parcb<-c(par[2])
      if(cb=="gaus"){
        if(par[1]<=0 | (par[2]<=-1 || par[2] >=1)){return(1000e10)}
      }
      else if(cb=="joe" | cb=="gumbel" | cb=="InverGumbel"){
        if(par[1]<=0 | (par[2]<1 || par[2]>=10)){return(1000e10)}
      }
      else if(cb=="clayton"){
        if(par[1]<=0 | (par[2]<0 || par[2]>=15)){return(1000e10)}
      }
      else if(cb=="frank"){
        if(par[1]<=0 | par[2]==0){return(1000e10)}
      }
      else if(cb=="HR" || cb=="galambos"){
        if(par[1]<=0 | par[2]<=0){return(1000e10)}
      }
      parweight<-par[3]
      if(par[3]<=0){return(1000e10)}
    }
  }
  else{
    warning("The copula is not implemented!")
  }
  K1<-K(parct=parct,parcb=parcb,ct=ct,cb=cb,weightfun=weightfun,parweight=parweight)
  logcop<-log(cop(u,v,parct=parct,parcb=parcb,ct=ct,cb=cb,weightfun=weightfun,parweight=parweight,K=K1))
  return(-sum(logcop))
}

# ==============================
# Functions for the diagnostics
# ==============================

## AIC
AICm<-function(par,logvalue){
  return(2*length(par)-2*logvalue)
}

## Survival function 1-F
survivalf<-function(x,parct,parcb,ct,cb,weightfun,parweight,K){
  l<-quant(x,parct=parct,parcb=parcb,ct=ct,cb=cb,weightfun=weightfun,parweight=parweight,K=K)
  surv<-sapply(l,function(x){integral2(fun,x,1,x,1,parct=parct,parcb=parcb,ct=ct,cb=cb,weightfun=weightfun,parweight=parweight,K=K,vectorize=T,singular=T)$Q})
  return(surv)
}

## Empirical \chi(r)
### thresh: r \in (0,1)
chiu_emp<-function(u,v,thresh){
  chi_emp<-c()
  for(i in 1:length(thresh)){
    chi_emp[i]<-mean(u>thresh[i]&v>thresh[i])/(mean(u>thresh[i]))
  }
  chi_emp[is.na(chi_emp)]<-0
  return(chi_emp)
}

## Model \chi(r)
### survP: survival function
chiu_model<-function(thresh,survP){
  chi_model<-c()
  for(i in 1:length(thresh)){
    chi_model[i]<-survP[i]/(1-thresh[i])
  }
  return(chi_model)
}

## Empirical \eta(r) and \eta(r) with Hill estimator
etau_emp<-function(u,v,thresh){
  uexp<-sapply(u,qexp) # for eta computed with the Hill estimator
  vexp<-sapply(v,qexp) # for eta computed with the Hill estimator
  compmin<-pmin(uexp,vexp) # for eta estimated via the Hill estimator
  eta_emp<-c()
  eta_hill<-c()
  for(i in 1:length(thresh)){
    eta_emp[i]<-(log(mean(u>thresh[i])))/(log(mean(u>thresh[i]&v>thresh[i])))
    eta_hill[i]<-mean(compmin[compmin>thresh[i]]-thresh[i])
  }
  eta_emp[is.na(eta_emp)]<-0
  return(cbind("Empirical"=eta_emp,"Hill"=eta_hill))
}

## Model \eta(r)
etau_model<-function(thresh,survP){
  eta_model<-c()
  for(i in 1:length(thresh)){
    eta_model[i]<-(log(1-thresh[i]))/(log(survP[i]))
  }
  return(eta_model)
}

## Block bootstrap
### k: block length
blockboot<-function(data,k){
  finaldata<-NULL
  for(i in 1:length(unique(data$Year))){
    datasummer<-data[data$Year==unique(data$Year)[i],]
    nyear<-dim(datasummer)[1]
    nblocks<-ceiling(nyear/k)
    nnew<-nblocks*k
    newdata<-matrix(NA,nrow=nnew,ncol=dim(data)[2])
    idx<-1:(nyear-k+1)
    startpoints<-sample(x=idx,size=nblocks,replace=T)
    datasummer<-as.matrix(datasummer)
    for(j in 1:nblocks){
      newdata[((j-1)*k+1):(j*k),]<-datasummer[(startpoints[j]:(startpoints[j]+k-1)),]
    }
    finaldata<-rbind(finaldata,newdata)
  }
  return(finaldata)
}

## Block bootstrap of the original data
bootoriginal<-function(data,nsim,k){
  nboots<-dim(data)[1]
  bootstrap<-lapply(1:nsim,matrix,data=NA,nrow=nboots,ncol=dim(data)[2])
  for(i in 1:nsim){
    boot<-blockboot(data,k)
    u<-as.numeric(boot[,9])
    v<-as.numeric(boot[,10])
    bootstrap[[i]]<-cbind(u,v)
  }
  return(bootstrap)
}

## Confidence intervals
CIs<-function(data,thresh,k){
  nthresh<-length(thresh)
  nboots<-dim(data)[1]
  chibootmodel<-matrix(0,ncol=nthresh,nrow=nboots)
  etabootmodel_e<-matrix(0,ncol=nthresh,nrow=nboots)
  etabootmodel_h<-matrix(0,ncol=nthresh,nrow=nboots)
  for(i in 1:nboots){
    boot<-blockboot(data,k)
    chibootmodel[i,]<-chiu_emp(u=as.numeric(boot[,9]),v=as.numeric(boot[,10]),thresh=thresh)
    etabootmodel_e[i,]<-etau_emp(u=as.numeric(boot[,9]),v=as.numeric(boot[,10]),thresh=thresh)[,1]
    etabootmodel_h[i,]<-etau_emp(u=as.numeric(boot[,9]),v=as.numeric(boot[,10]),thresh=thresh)[,2]
  }
  CIchimodel<-cbind(apply(chibootmodel,2,quantile,0.025),apply(chibootmodel,2,quantile,0.975))
  CIetamodel_e<-cbind(apply(etabootmodel_e,2,quantile,0.025),apply(etabootmodel_e,2,quantile,0.975))
  CIetamodel_h<-cbind(apply(etabootmodel_h,2,quantile,0.025),apply(etabootmodel_h,2,quantile,0.975))
  return(list(CIchimodel,CIetamodel_e,CIetamodel_h))
}

## Simulation from c^*
sim<-function(n=n,nsim=nsim,parct=parct,parcb=parcb,ct=ct,cb=cb,weightfun=weightfun,parweight=parweight){
  if(ct=="t"){
    data_ct<-rCopula(n,tCopula(param=parct[1],df=parct[2]))
  }
  else if(ct=="gaus"){
    data_ct<-rCopula(n,normalCopula(param=parct[1]))
  }
  else if(ct=="joe"){
    data_ct<-rCopula(n,joeCopula(param=parct[1]))
  }
  else if(ct=="frank"){
    data_ct<-rCopula(n,frankCopula(param=parct[1]))
  }
  else if(ct=="gumbel"){
    data_ct<-rCopula(n,gumbelCopula(param=parct[1]))
  }
  else if(ct=="InverGumbel"){
    data_log<-rbvevd(n,dep=1/parct[1],mar1=c(1,1,1))
    data_exp<-1/data_log
    data_ct<-apply(data_exp,2,pexp)
  }
  else if(ct=="HR"){
    data_log<-rbvevd(n,dep=parct[1],model="hr",mar1=c(1,1,1))
    data_ct<-apply(data_log,2,evd::pgev,loc=1,scale=1,shape=1)
  }
  else if(ct=="ct"){
    data_log<-rbvevd(n,alpha=parct[1],beta=parct[2],model="ct",mar1=c(1,1,1))
    data_ct<-apply(data_log,2,evd::pgev,loc=1,scale=1,shape=1)
  }
  else if(ct=="galambos"){
    data_log<-rbvevd(n,dep=parct[1],model="neglog",mar1=c(1,1,1))
    data_ct<-apply(data_log,2,evd::pgev,loc=1,scale=1,shape=1)
  }
  else if(ct=="cop2019"){
    data_ct<-rC2(n,delta=parct[1],theta=parct[2],model="Gauss",scale="unif")
  }
  else if(ct=="clayton"){
    data_ct<-rCopula(n,claytonCopula(param=parct[1]))
  }
  u_ct<-data_ct[,1]
  v_ct<-data_ct[,2]
  pifun_ct<-weightfun(u_ct,v_ct,parweight)
  keep_ct<-matrix(NA,nrow=n,ncol=2)
  for(i in 1:n){
    if(runif(1)<pifun_ct[i]){keep_ct[i,]<-data_ct[i,]}
  }
  keep_ct<-keep_ct[complete.cases(keep_ct),]
  if(cb=="t"){
    data_cb<-rCopula(n,tCopula(param=parcb[1],df=parcb[2]))
  }
  else if(cb=="gaus"){
    data_cb<-rCopula(n,normalCopula(param=parcb[1]))
  }
  else if(cb=="joe"){
    data_cb<-rCopula(n,joeCopula(param=parcb[1]))
  }
  else if(cb=="frank"){
    data_cb<-rCopula(n,frankCopula(param=parcb[1]))
  }
  else if(cb=="gumbel"){
    data_cb<-rCopula(n,gumbelCopula(param=parcb[1]))
  }
  else if(cb=="InverGumbel"){
    data_log<-rbvevd(n,dep=1/parcb[1],mar1=c(1,1,1))
    data_exp<-1/data_log
    data_cb<-apply(data_exp,2,pexp)
  }
  else if(cb=="cop2019"){
    data_cb<-rC2(n,delta=parcb[1],theta=parcb[2],model="Gauss",scale="unif")
  }
  else if(cb=="clayton"){
    data_cb<-rCopula(n,claytonCopula(param=parcb[1]))
  }
  else if(cb=="HR"){
    data_log<-tryCatch(exp=rbvevd(n,dep=parcb[1],model="hr",mar1=c(1,1,1)),error=function(e) e)
    if(is.element("error",class(data_log))){
      data_log<-rbvevd(n,dep=parcb[1],model="hr",mar1=c(1,1,1))
    }
    data_cb<-apply(data_log,2,evd::pgev,loc=1,scale=1,shape=1)
  }
  else if(cb=="ct"){
    data_log<-rbvevd(n,alpha=parcb[1],beta=parcb[2],model="ct",mar1=c(1,1,1))
    data_cb<-apply(data_log,2,evd::pgev,loc=1,scale=1,shape=1)
  }
  else if(cb=="galambos"){
    data_log<-rbvevd(n,dep=parcb[1],model="neglog",mar1=c(1,1,1))
    data_cb<-apply(data_log,2,evd::pgev,loc=1,scale=1,shape=1)
  }
  u_cb<-data_cb[,1]
  v_cb<-data_cb[,2]
  pifun_cb<-weightfun(u_cb,v_cb,parweight)
  keep_cb<-matrix(NA,nrow=n,ncol=2)
  for(i in 1:n){
    if(runif(1)<1-pifun_cb[i]){keep_cb[i,]<-data_cb[i,]}
  }
  keep_cb<-keep_cb[complete.cases(keep_cb),]
  cstar<-rbind(keep_ct,keep_cb)
  cstar<-cstar[sample(nrow(cstar),nsim),]
  return(cstar)
}
