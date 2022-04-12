
# T representation where V is a vector of independendent Gumbel r.v.s 

#--------------------------------------------------------------------------------------------------------------------------------
# Front end functions (fit.MGPD.*** and Sim.***.MGPD) have syntax to match the paper closely
# However, other back end functions have syntax relating to earlier parameterizations and models constructed on different scales
#--------------------------------------------------------------------------------------------------------------------------------

################################################################################
# Basic structure:

# fY.*** - density in MV Pareto margins  (exponential of the standard form margins in KRSW)
# fX.*** - density in MGPD margins

# f*.***.cens - censored density

# nll.*** - negative log-likelihood to be optimised (MV Pareto scale)
# nll.***.GPD - negative log-likelihood to be optimised (MGPD scale)

# fit.MGPD.*** - function to perform the fit

# Sim.***.MGPD - function to simulate from the model
################################################################################

# For identifiability final lambda parameter (scale in Frechet / location in Gumbel) always = 1

###################################################################################################################
# Functions allowing only a single common shape parameter for the Frechet (scale parameter for the Gumbel)
###################################################################################################################

# MP scale
##########



fY.Fre.a.Linf<-function(y,a,lam)
{
  d<-length(y)
  c1<-gamma(d)*a^(d-1) 
  num<-prod(y^(-a-1)*lam^a)
  den<-max(y)*sum((y/lam)^(-a))^(d)
  
  return(c1*num/den)
}


fY.Fre.a.Linf.cens<-function(y,a,lam,u)
{
  
  d<-length(y)
  D<-1:d
  C<-which(y<=u)
  DC<-D[-C]
  
  c1<-gamma(d-length(C))*a^(d-length(C)-1)
  num<-prod(y[DC]^(-a-1)*lam[DC]^a)
  den<-max(y)*(sum((y[DC]/lam[DC])^(-a)) + sum((u[C]/lam[C])^(-a)))^(d-length(C))
  return(c1*num/den)
}



nll.Frechet.a.Linf<-function(theta,y,u,lamfix=FALSE)
{
  d<-dim(y)[2]
  a<-theta[1]
  if(lamfix){lam<-rep(1,d)}
  else{
    lam<-c(theta[2:d],1)
  }
  
  if(any(lam<0.01)){return(10e10)}
  
  ind<-apply(y,1,comp.gt,u=u)
  y.uc<-y[ind,]
  y.pc<-y[!ind,]
  
  L<-apply(y.uc,1,fY.Fre.a.Linf,a=a,lam=lam)
  nll.uc<--sum(log(L))
  
  if(sum(!ind)>0)
  {
    L2<-apply(y.pc,1,fY.Fre.a.Linf.cens,a=a,lam=lam,u=u)
    nll.pc<--sum(log(L2))
  }
  else{nll.pc<-0}
  nll<-nll.uc+nll.pc
  return(nll)
}

###################################################################################################################
# MGPD scale
#############


fX.Fre.a.Linf<-function(x,a,lam,sig,gamma)
{
  y<-BCi(x=x,gamma=gamma,sig=sig)
  J<-Jac(x=x,gamma=gamma,sig=sig)
  return(J*fY.Fre.a.Linf(y=y,a=a,lam=lam))
}


fX.Fre.a.Linf.cens<-function(x,a,lam,u,sig,gamma)
{
  y<-BCi(x=x,gamma=gamma,sig=sig)
  uy<-BCi(x=u,gamma=gamma,sig=sig)
  
  d<-length(y)
  D<-1:d
  C<-which(x<=u)
  DC<-D[-C]
  J<-Jac(x=x[DC],sig=sig[DC],gamma=gamma[DC])
    
  c1<-gamma(d-length(C))*a^(d-length(C)-1)
  num<-prod(y[DC]^(-a-1)*lam[DC]^a)
  den<-max(y,na.rm=T)*(sum((y[DC]/lam[DC])^(-a)) + sum((uy[C]/lam[C])^(-a)))^(d-length(C))
  return(J*c1*num/den)
}



nll.Frechet.a.Linf.GPD<-function(theta,x,u,sig.ind,gamma.ind,lamfix=F,marg.scale.ind,marg.shape.ind)
{
  d<-dim(x)[2]
  a<-theta[1]
  if(lamfix){lam<-rep(1,d)}
  else{
    lam<-c(theta[(2):(d)],1)
  }
  sig<-theta[sig.ind]
  gamma<-theta[gamma.ind]
  
  sig<-sig[marg.scale.ind]
  gamma<-gamma[marg.shape.ind]
  
  # check data respect marginal constraints (i.e. no x[,j]>-sig[j]/gamma[j] if gamma[j]<0)
  chk<-rep(TRUE,d)
  for(j in 1:d)
  {
    if(gamma[j]<0){chk[j]<-all(x[,j]< (-sig[j]/gamma[j]))}
  }
  if(any(lam<0.01)||any(sig<0.001)||any(!chk)){return(10e10)}
  
  ind<-apply(x,1,comp.gt,u=u)
  x.uc<-x[ind,]
  x.pc<-x[!ind,]
  
  L<-apply(x.uc,1,fX.Fre.a.Linf,a=a,lam=lam,sig=sig,gamma=gamma)
  nll.uc<--sum(log(L))
  
  if(sum(!ind)>0)
  {
    L2<-apply(x.pc,1,fX.Fre.a.Linf.cens,a=a,lam=lam,u=u,sig=sig,gamma=gamma)
    nll.pc<--sum(log(L2))
  }
  else{nll.pc<-0}
  nll<-nll.uc+nll.pc
  return(nll)
}



###################################################################################################################
###################################################################################################################

###################################################################################################################
# Functions allowing different shape parameters for the Frechet (scale parameters for the Gumbel)
###################################################################################################################

# MP scale
##########


fY.Fre.Linf<-function(y,a,lam)
{
  to.int<-function(r)
  {
    p<-r*y
    return(prod(a*lam^a *p^(-a)*y^(-1)*exp(-(p/lam)^(-a)))/r)
  }
  vti<-Vectorize(to.int)
  int<-tryCatch(integrate(vti,0,Inf,abs.tol=0), error=function(e) e)
  if(is.element("error",class(int)))
  {
    warning(paste("\n Numerical difficulties (probably an integration error) occured at parameter values:", round(a,4),",",round(lam,4),
                  "\n NA returned to optimization"))
    
    # try once more with different tolerance
    int<-tryCatch(integrate(vti,0,Inf,rel.tol=1e-10), error=function(e) e)
    if(is.element("error",class(int)))
    {
      warning(paste("\n Numerical difficulties at second attempt (probably an integration error) occured at parameter values:", round(a,4),",",round(lam,4),
                    "\n NA returned to optimization"))
    return(NA)
    }
  }
  return(int$value/max(y))
}


fY.Fre.Linf.cens<-function(y,u,a,lam)
{
  d<-length(y)
  D<-1:d
  C<-which(y<=u)
  DC<-D[-C]
  
  to.int<-function(r)
  {
    p<-r*y
    uncens<-prod(a[DC]*lam[DC]^a[DC] *p[DC]^(-a[DC])*y[DC]^(-1)*exp(-(p[DC]/lam[DC])^(-a[DC])))/r
    cens<-prod(exp(-((u[C]*r)/lam[C])^(-a[C])))
    return(uncens*cens)
  }
  vti<-Vectorize(to.int)
  int<-tryCatch(integrate(vti,0,Inf,abs.tol=0), error=function(e) e)
  if(is.element("error",class(int)))
  {
    warning(paste("Numerical difficulties (probably an integration error) occured at parameter values:", round(a,4),",",round(lam,4),
                  "\n NA returned to optimization"))
    return(NA)
  }
  return(int$value/max(y))
}



nll.Frechet.Linf<-function(theta,y,u,lamfix=F)
{
  d<-dim(y)[2]
  a<-theta[1:d]
  # single constraint - make lambdas relative to 1 in dth margin
  if(lamfix){lam<-rep(1,d)}
  else{
    lam<-c(theta[(d+1):(2*d-1)],1)
  }
  
  # can have a<1 here as variables normalized
  if(any(lam<0.01)){return(10e10)}
  
  uc<-apply(y,1,comp.gt,u=u)
  y.uc<-y[uc,]
  y.pc<-y[!uc,]
  
  L<-apply(y.uc,1,fY.Fre.Linf,a=a,lam=lam)
  nll.uc<--sum(log(L))
  
  if(sum(!uc)>1)
  {
    L2<-apply(y.pc,1,fY.Fre.Linf.cens,a=a,lam=lam,u=u)
    nll.pc<--sum(log(L2))
  }
  else{nll.pc<-0}
  nll<-nll.uc+nll.pc
  return(nll)
}



# MGPD scale
############


fX.Fre.Linf<-function(x,a,lam,sig,gamma)
{
  y<-BCi(x=x,gamma=gamma,sig=sig)
  J<-Jac(x=x,gamma=gamma,sig=sig)
  
  return(fY.Fre.Linf(y=y,a=a,lam=lam)*J)
}


fX.Fre.Linf.cens<-function(x,u,a,lam,sig,gamma)
{
  y<-BCi(x=x,gamma=gamma,sig=sig)
  uy<-BCi(x=u,gamma=gamma,sig=sig)
  
  d<-length(y)
  D<-1:d
  C<-which(x<=u)
  DC<-D[-C]
  
  J<-prod((1/sig[DC])*(1+gamma[DC]*x[DC]/sig[DC])^(1/gamma[DC]-1))
  
  to.int<-function(r)
  {
    p<-r*y
    uncens<-prod(a[DC]*lam[DC]^a[DC] *p[DC]^(-a[DC])*y[DC]^(-1)*exp(-(p[DC]/lam[DC])^(-a[DC])))/r
    cens<-prod(exp(-((uy[C]*r)/lam[C])^(-a[C])))
    return(uncens*cens)
  }
  vti<-Vectorize(to.int)
  
  int<-tryCatch(integrate(vti,0,Inf,abs.tol=0), error=function(e) e)
  if(is.element("error",class(int)))
  {
    warning(paste("\n Numerical difficulties (probably an integration error) occured at parameter values:", round(a,4),",",round(lam,4),",",round(sig,4),",",round(gamma,4),
                  "\n NA returned to optimization"))
    
    # try once more with altered precision
    int<-tryCatch(integrate(vti,0,Inf,rel.tol=1e-10), error=function(e) e)
    if(is.element("error",class(int)))
    {
      warning(paste("\n Numerical difficulties at second attempt (probably an integration error) occured at parameter values:", round(a,4),",",round(lam,4),",",round(sig,4),",",round(gamma,4),
                    "\n NA returned to optimization"))
    return(NA)
    }
  }
  return(J*int$value/max(y,na.rm=T))
}


nll.Frechet.Linf.GPD<-function(theta,x,u,lamfix=F,sig.ind,gamma.ind,marg.scale.ind,marg.shape.ind)
{
  d<-dim(x)[2]
  a<-theta[1:d]
  if(lamfix){lam<-rep(1,d)}
  else{
    lam<-c(theta[(d+1):(2*d-1)],1)
  }
  sig<-theta[sig.ind]
  gamma<-theta[gamma.ind]
  
  sig<-sig[marg.scale.ind]
  gamma<-gamma[marg.shape.ind]
  
  # check data respect marginal constraints (i.e. no x[,j]>-sig[j]/gamma[j] if gamma[j]<0)
  chk<-rep(TRUE,d)
  for(j in 1:d)
  {
    if(gamma[j]<0){chk[j]<-all(x[,j]< (-sig[j]/gamma[j]))}
  }
  if(any(lam<0.01)||any(sig<0.001)||any(!chk)){return(10e10)}
  
  uc<-apply(x,1,comp.gt,u=u)
  x.uc<-x[uc,]
  x.pc<-x[!uc,]
  
  L<-apply(x.uc,1,fX.Fre.Linf,a=a,lam=lam,sig=sig,gamma=gamma)
  nll.uc<--sum(log(L))
  
  if(sum(!uc)>0)
  {
    L2<-apply(x.pc,1,fX.Fre.Linf.cens,a=a,lam=lam,u=u,sig=sig,gamma=gamma)
    nll.pc<--sum(log(L2))
  }
  else{nll.pc<-0}
  nll<-nll.uc+nll.pc
  
  return(nll)
  
}


#####################################################################################################################################

# fit.MGPD.GumbelT

# Arguments:
# x - data (with >=1 exceedance above 0). Matrix of dim n by d.
# u - d-vector of thresholds at which to censor. Max value 0.
# std - logical. If TRUE then fit on standardized (Exp) scale, else fit on GPD scale [actual fitting done after tf to Pareto scale]
# dep.scale.fix - logical. If TRUE then where appropriate the variables in the dependence model are driven by a single scale parameter.
# dep.loc.fix - logical. If TRUE then where appropriate the variables in the dependence model are driven by a single location parameter.
# marg.scale.ind, marg.shape.ind - numerical vecotrs of length d denoting which margins should have common scale / shape pars 
#   (e.g. if d=3, say, c(1,1,1) denotes all identical; c(1,2,1), denotes 1st and 3rd the same, 2nd diff; c(1,2,3) denotes all different)


fit.MGPD.GumbelT<-function(x, u, std=FALSE, dep.scale.fix=FALSE, dep.loc.fix=FALSE, dep.start=NULL, marg.scale.start=NULL,
                   marg.shape.start=NULL, marg.scale.ind, marg.shape.ind, maxit=1000,...)
{
  # if std=TRUE then fit multivariate Pareto (i.e. "standardized" scale)
  # tf to Pareto scale for fitting
  if(std){x<-exp(x); u<-exp(u)} 
  
  d<-dim(x)[2]

  # Move from beta scale to lambda scale
  if(!is.null(dep.start)&&!dep.loc.fix){dep.start[(length(dep.start)-d+2):length(dep.start)]<-exp(dep.start[(length(dep.start)-d+2):length(dep.start)])}
  
    
  if(std)
  {
      if(dep.scale.fix)
      {
        if(dep.loc.fix)
        {
          par<-c(dep.start[!is.null(dep.start)], c(2)[is.null(dep.start)])
          a.ind<-1
          opt<-optimize(nll.Frechet.a.Linf, interval=c(par/10,10*par), u=u, y=x, lamfix=T,...)
          
        }
        else{
          par<-c(dep.start[!is.null(dep.start)], c(2,rep(1,d-1))[is.null(dep.start)])
          a.ind<-1
          lam.ind<-2:d
          opt<-optim(nll.Frechet.a.Linf, par=par, u=u, y=x, lamfix=F, control=list(maxit=maxit,reltol=1e-6),...)
        }
      }
      else{
        if(dep.loc.fix)
        {
          par<-c(dep.start[!is.null(dep.start)], c(rep(2,d))[is.null(dep.start)])
          a.ind<-1:d
          opt<-optim(nll.Frechet.Linf, par=par, u=u, y=x, lamfix=T, control=list(maxit=maxit,reltol=1e-6),...)
        }
        else
        {
          par<-c(dep.start[!is.null(dep.start)], c(rep(2,d),rep(1,d-1))[is.null(dep.start)])
          a.ind<-1:d
          lam.ind<-(d+1):(2*d-1)
          opt<-optim(nll.Frechet.Linf, par=par, u=u, y=x, lamfix=F, control=list(maxit=maxit,reltol=1e-6),...)
        }
      }
    }
 
  ##############################
  # GPD scale
  
  else{
    n.sig<-length(unique(marg.scale.ind))
    n.gamma<-length(unique(marg.shape.ind))
    

      if(dep.scale.fix)
      {
        if(dep.loc.fix)
        {
          par<-c(dep.start[!is.null(dep.start)], c(2)[is.null(dep.start)],
                 marg.scale.start[!is.null(marg.scale.start)], rep(1,n.sig)[is.null(marg.scale.start)],
                 marg.shape.start[!is.null(marg.shape.start)], rep(0.1,n.gamma)[is.null(marg.shape.start)])
          a.ind<-1
          sig.ind<-1+1:n.sig
          gamma.ind<-1+length(sig.ind)+1:n.gamma
          
          opt<-optim(nll.Frechet.a.Linf.GPD, par=par, u=u, x=x, lamfix=T, 
                     sig.ind=sig.ind, gamma.ind=gamma.ind, marg.scale.ind=marg.scale.ind, marg.shape.ind=marg.shape.ind, control=list(maxit=maxit,reltol=1e-6),...)
        }
        else{
          par<-c(dep.start[!is.null(dep.start)], c(2,rep(1,d-1))[is.null(dep.start)],
                 marg.scale.start[!is.null(marg.scale.start)], rep(1,n.sig)[is.null(marg.scale.start)],
                 marg.shape.start[!is.null(marg.shape.start)], rep(0.1,n.gamma)[is.null(marg.shape.start)])
          a.ind<-1
          lam.ind<-2:d
          sig.ind<-d+1:n.sig
          gamma.ind<-d+length(sig.ind)+1:n.gamma
          
          opt<-optim(nll.Frechet.a.Linf.GPD, par=par, u=u, x=x, lamfix=F, 
                     sig.ind=sig.ind, gamma.ind=gamma.ind, marg.scale.ind=marg.scale.ind, marg.shape.ind=marg.shape.ind, control=list(maxit=maxit,reltol=1e-6),...)
        }
      }
      else{
        if(dep.loc.fix)
        {
          par<-c(dep.start[!is.null(dep.start)], c(rep(2,d))[is.null(dep.start)],
                 marg.scale.start[!is.null(marg.scale.start)], rep(1,n.sig)[is.null(marg.scale.start)],
                 marg.shape.start[!is.null(marg.shape.start)], rep(0.1,n.gamma)[is.null(marg.shape.start)])
          a.ind<-1:d
          sig.ind<-d+1:n.sig
          gamma.ind<-d+length(sig.ind)+1:n.gamma
          
          opt<-optim(nll.Frechet.Linf.GPD, par=par, u=u, x=x, lamfix=T, 
                     sig.ind=sig.ind, gamma.ind=gamma.ind, marg.scale.ind=marg.scale.ind, marg.shape.ind=marg.shape.ind, control=list(maxit=maxit,reltol=1e-6),...)
        }
        else
        {
          par<-c(dep.start[!is.null(dep.start)], c(rep(2,d),rep(1,d-1))[is.null(dep.start)],
                 marg.scale.start[!is.null(marg.scale.start)], rep(1,n.sig)[is.null(marg.scale.start)],
                 marg.shape.start[!is.null(marg.shape.start)], rep(0.1,n.gamma)[is.null(marg.shape.start)])
          a.ind<-1:d
          lam.ind<-(d+1):(2*d-1)
          sig.ind<-2*d-1 +1:n.sig
          gamma.ind<-2*d-1+length(sig.ind)+1:n.gamma
          
          opt<-optim(nll.Frechet.Linf.GPD, par=par, u=u, x=x, lamfix=F, 
                     sig.ind=sig.ind, gamma.ind=gamma.ind, marg.scale.ind=marg.scale.ind, marg.shape.ind=marg.shape.ind, control=list(maxit=maxit,reltol=1e-6),...)
          
        }
      }  
  }
  
  if(is.null(opt$min))
  {
    mle<-opt$par
    # transform from lambda scale to beta scale
    if(!dep.loc.fix){mle[lam.ind]<-log(mle[lam.ind])}
    nll<-opt$value
    conv<-opt$conv
    hess<-opt$hess
    if(!dep.loc.fix && !is.null(hess)){warn<-"WARNING: Hessian for beta parameters is on lambda=exp(beta) scale. Use delta method."}
    else{warn<-NULL}
  }
  else{
    mle<-opt$minimum
    nll<-opt$objective
    conv<-NULL
    hess<-NULL
    warn<-NULL
  }
  return(list(mle=mle,nll=nll,conv=conv,hess=hess,warn=warn))
}



#############################################################################################################################
#############################################################################################################################

rfrechet<-function(n,a,lam)
{
  lam*(-log(runif(n)))^(-1/a)
}

# sim.GumbelT.MGPD - Simulate from Frechet model

# n - number of vectors to simulate
# d - dimension
# a - dependence parameter a 
# beta - dependence parameter beta
# sig - marginal scale param
# gamma - marginal shape param
# std / MGPD  - logical: return standardized (Exp) / MGPD / both scales

sim.GumbelT.MGPD<-function(n,d,a,beta,sig,gamma,MGPD=TRUE,std=FALSE)
{
  lam<-exp(beta)
  W<-matrix(0,ncol=d,nrow=n)  
  
  for(i in 1:n)
  {
  W[i,]<-rfrechet(d,a=a,lam=lam)
  }
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

