
# T representation where V is a vector of MV Gaussian r.v.s 

#--------------------------------------------------------------------------------------------------------------------------------
# Front end functions (fit.MGPD.*** and Sim.***.MGPD) have syntax to match the paper closely
# However, other back end functions have syntax relating to earlier parameterizations and models constructed on different scales
#--------------------------------------------------------------------------------------------------------------------------------


################################################################################
# Basic structure:

# fY.*** - density in MV Pareto margins (exponential of the standard form margins in KRSW)
# fX.*** - density in MGPD margins

# f*.***.cens - censored density

# nll.*** - negative log-likelihood to be optimised (MV Pareto scale)
# nll.***.GPD - negative log-likelihood to be optimised (MGPD scale)

# fit.MGPD.*** - function to perform the fit

# Sim.***.MGPD - function to simulate from the model
################################################################################

# For identifiability final lambda parameter (scale in pow unif / log(lam) is  in reverse exponential) always = 1

###################################################################################################################
# Functions allowing only a single common shape parameter for the Powunif (scale parameter for the Gumbel)
###################################################################################################################

# MP scale
##########


fY.logGauss.Linf<-function(y,beta,Sig)
{
  d<-length(y)
  one<-rep(1,d)
  SI<-solve(Sig)
  
  q<-sum(SI)
  
  a1<-colSums(SI)
  A<-SI-outer(a1,a1)/c(q)
  
  out<-(2*pi)^((1-d)/2) *sqrt(1/q) *det(SI)^(1/2)*(1/max(y))*(prod(1/y))*exp(-0.5*t(log(y)-beta)%*%(A%*%(log(y)-beta)))
  return(out)
}


##############################################################

# censored density


# Arguments:
# y - data (in "Pareto" margins)
# Sig - covariance matrix for the lognormal
# u - vector (of same length as y) of censoring thresholds

fY.logGauss.Linf.cens<-function(y,Sig,beta,u)
{
  d<-length(y)
  D<-1:d
  C<-which(y<=u)
  
  one<-rep(1,d)
  SI<-solve(Sig)
  q<-sum(SI)
  a1<-colSums(SI)
  A<-SI-outer(a1,a1)/c(q)
  
  oneDC<-rep(1,d-length(C))
  DC<-D[-C]
  
  SIDC<-solve(Sig[DC,DC])
  qDC<-sum(SIDC)
  
  a1DC<-colSums(SIDC)
  ADC<-SIDC-outer(a1DC,a1DC)/c(qDC)
  
  KC<-matrix(0,ncol=length(C),nrow=d)
  KC[-(DC),]<-diag(length(C))
  
  KDC<-matrix(0,ncol=d-length(C),nrow=d)
  KDC[(DC),]<-diag(length(DC))
  
  GI<-A[C,C]
  G<-solve(GI)
  m<--(G%*%(t(KC)%*%A))%*%(KDC%*%(log(y[DC])-beta[DC]))
  
  out<-fY.logGauss.Linf(y=y[DC],beta=beta[DC],Sig=Sig[DC,DC])*pmvnorm(upper=log(u[C])-beta[C],mean=c(m),sigma=G)
  return(out)
}



nll.logGauss.Linf<-function(theta, y, u, structured.cor=FALSE, cor.ind, beta.ind, dist)
{
  if(structured.cor)
  {
    phi<-theta[1]
    if(phi<0.01){return(10e10)}
    Sig<-matern(dist,phi=phi,kappa=1)
  }
  else{
    d<-dim(y)[2]
    rho<-theta[cor.ind]
    Sig<-diag(d)
    Sig[upper.tri(Sig)]<-rho  
    Sig<-forceSymmetric(Sig)
    Sig<-as.matrix(Sig)
    if(any(eigen(Sig)$values<1e-8) || any(rho<=-1) || any(rho>=1)){return(10e10)}
  }
  
  # set final location param value equal to 0 
  if(!is.null(beta.ind)){beta<-c(theta[beta.ind],0)}
  else(beta<-rep(0,d))
  
  # logical vector of completely uncensored values
  uc<-apply(y,1,comp.gt,u=u)
  
  y.uc<-y[uc,]
  y.pc<-y[!uc,]
  
  L<-apply(y.uc,1,fY.logGauss.Linf,Sig=Sig, beta=beta)
  nll.uc<--sum(log(L))
  
  # censored contribution
  
  if(sum(!uc)>0)
  {
    L2<-apply(y.pc,1,fY.logGauss.Linf.cens,u=u,Sig=Sig, beta=beta)
    nll.c<--sum(log(L2))
  }
  else{nll.c<-0}
  
  nll<-nll.uc+nll.c
  return(nll)
}

# MGP scale
###########

# Arguments
# x - data on **observed** margins (but with threshold subtracted, so at least 1 coordinate >0 -- assuming data to be generally +ve)
# Sig - covariance matrix
# beta - vecotr of location parameters for underlying Gaussian r.v.
# eta - d-vector of scale parameter
# gamma - d-vector of shape parameters


fX.logGauss.Linf<-function(x,Sig,beta,eta,gamma)
{
  y<-BCi(x=x,gamma=gamma,sig=eta)
  J<-Jac(x=x,gamma=gamma,sig=eta)
  return(fY.logGauss.Linf(y=y,beta=beta,Sig=Sig)*J)
}
###########################################################

# Censored version

# Arguments
# x - data on **observed** margins (but with threshold subtracted, so at least 1 coordinate >0 -- assuming data to be generally +ve)
# Sig - covariance matrix
# eta - d-vector of scale parameter
# gamma - d-vector of shape parameters
# u - d-vector of censoring threshold

fX.logGauss.Linf.cens<-function(x,Sig,beta,eta,gamma,u)
{
  
  y<-BCi(x=x,gamma=gamma,sig=eta)
  uy<-BCi(x=u,gamma=gamma,sig=eta)
  
  d<-length(y)
  D<-1:d
  C<-which(x<=u)
  
  one<-rep(1,d)
  SI<-solve(Sig)
  q<-sum(SI)
  a1<-colSums(SI)
  A<-SI-outer(a1,a1)/c(q)
  
  
  oneDC<-rep(1,d-length(C))
  DC<-D[-C]
  
  SIDC<-solve(Sig[DC,DC])
  qDC<-sum(SIDC)
  
  a1DC<-colSums(SIDC)
  ADC<-SIDC-outer(a1DC,a1DC)/c(qDC)
  
  KC<-matrix(0,ncol=length(C),nrow=d)
  KC[-(DC),]<-diag(length(C))
  
  KDC<-matrix(0,ncol=d-length(C),nrow=d)
  KDC[(DC),]<-diag(length(DC))
  
  GI<-A[C,C]
  G<-solve(GI)
  m<--G%*%t(KC)%*%A%*%(KDC%*%(log(y[DC])-beta[DC]))
  
  J<-prod((1/eta[DC])*(1+gamma[DC]*x[DC]/eta[DC])^(1/gamma[DC]-1))
  
  Ind<-1+gamma[C]*u[C]/eta[C]<0
  log.u.y<-log(uy[C])
  log.u.y[Ind]<--Inf
  
  out<- J * fY.logGauss.Linf(y=y[DC],beta=beta[DC],Sig=Sig[DC,DC])*pmvnorm(upper=log.u.y-beta[C],mean=c(m),sigma=G)
  return(out)
}


# marg.scale.ind
# marg.shape.ind

nll.logGauss.Linf.GPD<-function(theta, x, u, structured.cor, cor.ind, dist, eta.ind, gamma.ind,
                                marg.scale.ind, marg.shape.ind, beta.ind)
{
  if(structured.cor)
  {
    phi<-theta[1]
    if(phi<0.01){return(10e10)}
    Sig<-matern(dist,phi=phi,kappa=1)
  }
  else{
    d<-dim(x)[2]
    rho<-theta[cor.ind]
    Sig<-diag(d)
    Sig[upper.tri(Sig)]<-rho  
    Sig<-forceSymmetric(Sig)
    Sig<-as.matrix(Sig)
    if(any(eigen(Sig)$values<1e-8) || any(rho<=-1) || any(rho>=1)){return(10e10)}
  }
  
  
  if(!is.null(beta.ind)){beta<-c(theta[beta.ind],0)}
  else{beta<-rep(0,d)}
  
  eta<-theta[eta.ind]
  gamma<-theta[gamma.ind]
  
  # make eta / gamma appropriate vectors of length d
  eta<-eta[marg.scale.ind]
  gamma<-gamma[marg.shape.ind]
  
  d<-dim(x)[2]
  
  if(any(eta<0.001)){return(10e10)}
  
  
  # logical vector of completely uncensored values
  uc<-apply(x,1,comp.gt,u=u)
  
  x.uc<-x[uc,]
  x.pc<-x[!uc,]
  
  # uncensored contribution
  
  L<-apply(x.uc,1,fX.logGauss.Linf,Sig=Sig,beta=beta,eta=eta,gamma=gamma)
  nll.uc<--sum(log(L))
  
  # censored contribution
  
  if(sum(!uc)>0)
  {
    L2<-apply(x.pc,1,fX.logGauss.Linf.cens,u=u,Sig=Sig,beta=beta,eta=eta,gamma=gamma)
    nll.c<--sum(log(L2))
  }
  else{nll.c<-0}
  nll<-nll.uc+nll.c
  return(nll)
}



###################################################################################################################


# fit.MGPD is a wrapper function to fit all coded MGPD models using their individual NLLs

# Arguments:
# x - data (with >=1 exceedance above 0). Matrix of dim n by d.
# u - d-vector of thresholds at which to censor. Max value 0.
# std - logical. If TRUE then fit on standardized (Exp) scale, else fit on GPD scale [actual fitting done after tf to Pareto scale]
# marg.scale.ind, marg.shape.ind - numerical vecotrs of length d denoting which margins should have common scale / shape pars 
#   (e.g. if d=3, say, c(1,1,1) denotes all identical; c(1,2,1), denotes 1st and 3rd the same, 2nd diff; c(1,2,3) denotes all different)
# structured.cor - logical. If TRUE there is a structured (e.g. spatial) correlation structure, and argument
#       dist is required. Currently fits a matern with shape 1 and free scale.
# dist - matrix of "distances" for structured correlation.
# dep.loc.fix - logical should the location parameters of the Gaussian be fixed 

#####################################################################################################################################
#####################################################################################################################################

# Modified to allow input of starting values - needs checking / testing


fit.MGPD.MVGaussT<-function(x, u, std=FALSE, dep.start=NULL, marg.scale.start=NULL,
                   marg.shape.start=NULL, marg.scale.ind, marg.shape.ind, structured.cor, dist, dep.loc.fix=FALSE, maxit=1000,...)
{
  # if std=TRUE then fit multivariate Pareto (i.e. "standardized" scale)
  # tf to Pareto scale for fitting
  if(std){x<-exp(x); u<-exp(u)} 
  
  d<-dim(x)[2]
  
  if(std)
  {
      if(!structured.cor){
        # rho, beta, sigma if applicable
        par<-c(dep.start[!is.null(dep.start)], rep(0,sum(1:(d-1)))[is.null(dep.start)],rep(0,d-1)[is.null(dep.start)&!dep.loc.fix])
        cor.ind<-1:sum(1:(d-1))
        if(!dep.loc.fix){beta.ind<-length(cor.ind)+1:(d-1)}
        else{beta.ind<-NULL}

      }
      else{
        par<-c(dep.start[!is.null(dep.start)], c(1)[is.null(dep.start)],rep(0,d-1)[is.null(dep.start)&!dep.loc.fix])
        if(!dep.loc.fix){beta.ind<-(2):(d)}
        else{beta.ind<-NULL}

      }
      if(length(par)==1)
      {
        if(!structured.cor){intvl<-c(-0.999,0.999)}
        else{intvl<-c(0.001,5*max(dist))}
        opt<-optimize(nll.logGauss.Linf, interval=intvl, y=x, u=u, structured.cor=structured.cor, cor.ind=cor.ind, dist=dist,beta.ind=beta.ind,...)
      }
      else{
        opt<-optim(nll.logGauss.Linf, par=par, y=x, u=u, structured.cor=structured.cor, cor.ind=cor.ind,
                   dist=dist, beta.ind=beta.ind, control=list(maxit=maxit,reltol=1e-6),...)
      }
  }
  
  ##############################
  # GPD scale
  
  else{
    
    n.eta<-length(unique(marg.scale.ind))
    n.gamma<-length(unique(marg.shape.ind))
    
      if(!structured.cor){
        par<-c(dep.start[!is.null(dep.start)], rep(0,sum(1:(d-1)))[is.null(dep.start)],
               rep(0,d-1)[is.null(dep.start)&!dep.loc.fix],
               marg.scale.start[!is.null(marg.scale.start)], rep(1,n.eta)[is.null(marg.scale.start)],
               marg.shape.start[!is.null(marg.shape.start)], rep(0.1,n.gamma)[is.null(marg.shape.start)])
        cor.ind<-1:sum(1:(d-1))
        if(!dep.loc.fix){beta.ind<-length(cor.ind)+1:(d-1)}
        else{beta.ind<-NULL}
        eta.ind<-length(cor.ind)+length(beta.ind)+1:n.eta
        gamma.ind<-length(cor.ind)+length(beta.ind)+length(eta.ind)+1:n.gamma
        
      }
      else{
        par<-c(dep.start[!is.null(dep.start)], c(1)[is.null(dep.start)],rep(0,d-1)[is.null(dep.start)&!dep.loc.fix],
               marg.scale.start[!is.null(marg.scale.start)], rep(1,n.eta)[is.null(marg.scale.start)],
               marg.shape.start[!is.null(marg.shape.start)], rep(0.1,n.gamma)[is.null(marg.shape.start)])
        if(!dep.loc.fix){beta.ind<-(2):(d)}
        else{beta.ind<-NULL}
        eta.ind<-1+length(beta.ind)+1:n.eta
        gamma.ind<-1 +length(beta.ind)+ length(eta.ind)+1:n.gamma

      }
      
      #marg.scale.ind # c(1,1,1) etc means all same; c(1,2,3) means all diff; c(1,2,1) means 1 and 3 equal etc
      #marg.shape.ind
      
      
      opt<-optim(nll.logGauss.Linf.GPD, par=par, x=x, u=u, structured.cor=structured.cor, cor.ind=cor.ind, dist=dist,
                 eta.ind=eta.ind, gamma.ind=gamma.ind, marg.scale.ind=marg.scale.ind, marg.shape.ind=marg.shape.ind,
                 beta.ind=beta.ind, control=list(maxit=maxit,reltol=1e-6),...)
    
  }

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
  }
  return(list(mle=mle,nll=nll,conv=conv,hess=hess))
}



#####################################################################################################################################
#####################################################################################################################################


# sim.logGaussLinf.MGPD - Simulate from log Gaussian (Linf) model

# n - number of vectors to simulate
# d - dimension
# a - dependence parameter a 
# Sig - covariance matrix in MV Gauss
# beta - location parameters in MV Gauss
# sig - marginal scale param (NB elsewhere in the code in this file, the marginal scale param is called eta)
# gamma - marginal shape param
# std / MGPD  - logical: return standardized (Exp) / MGPD / both scales

sim.MVGaussT.MGPD<-function(n,d,Sig,beta,sig,gamma,MGPD=TRUE,std=FALSE)
{
  Z<-rmvnorm(n,sigma=Sig,mean=beta)
  W<-exp(Z)
  mW<-apply(W,1,max)
  R<-runif(n,0,mW)
  Y<-W/R
  
  if(std&&!MGPD){return(log(Y))}
  
  X<-NULL
  for(j in 1:d)
  {
    if(gamma[j]!=0)
    {
      X<-cbind(X,sig[j]*(Y[,j]^gamma[j]-1)/gamma[j])
    }
    else
    {
      X<-cbind(X,sig[j]*log(Y[,j]))
    }
  }
  if(MGPD&&!std){return(X)}
  if(std&&MGPD){return(list(X=X,Z=log(Y)))}
}
