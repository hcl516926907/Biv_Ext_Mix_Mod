
# U representation where V is a vector of independent reverse exponential r.v.s

# !!
# Parameterization note: a is 1/\alpha, where \alpha is the parameterization in Kiriliouk, Rootzen, Segers and Wadsworth. 
# !!

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

fY.powunif<-function(y,lam,a)
{
  EM<-EM.pu(a=a,lam=lam)
  b<-1/a
  num<- prod(lam*b*(lam*y)^(b-1))
  den<-EM*((1+sum(b))*max(lam*y)^(sum(b)+1))
  return(num/den)
}



fY.powunif.cens<-function(y,u,lam,a)
{
  EM<-EM.pu(a=a,lam=lam)
  b<-1/a
  d<-length(y)
  C<-which(y<=u)
  D<-1:d
  DC<-D[-C]
  
  m1<-max(u[C]*lam[C])
  m2<-max(y[DC]*lam[DC])
  if(m1<m2)
  {
    num<-prod((u[C]*lam[C])^b[C])*prod((lam[DC]*b[DC])*(y[DC]*lam[DC])^(b[DC]-1))
    den<-EM*(sum(b)+1)*(max(y[DC]*lam[DC]))^(sum(b)+1)
    return(num/den)
  }
  else
  {
    v<-u[C]*lam[C]
    I<-which(v>max(y[DC]*lam[DC]))
    
    k<-length(I)
    tmp<-sort(v,index.return=T)
    v2<-rev(tmp$x)
    I2<-rev(tmp$ix)
    
    QDC<-prod(lam[DC]*b[DC]*(y[DC]*lam[DC])^(b[DC]-1))
    bit1<-prod((u[C]*lam[C])^b[C])*QDC / ((sum(b)+1)*(v2[1])^(sum(b)+1))
    
    if(k>1)
    {
      num<-den<-int<-NULL
      for(i in 1:(k-1))
      {
        num[i]<-prod((v[-I2[1:i]])^(b[C][-I2[1:i]])) * QDC
        den[i]<-sum(b[C][-I2[1:i]]) + sum(b[DC]) +1
        
        int[i]<-(v[I2[i+1]])^(-sum(b[C][-I2[1:i]]) - sum(b[DC])-1) -(v[I2[i]])^(-sum(b[C][-I2[1:i]]) - sum(b[DC])-1) 
      }
      bit2<-sum(num*int/den)
    }
    else{bit2<-0}
    bit3<- prod((v[-I2[1:k]])^(b[C][-I2[1:k]])) * QDC * ((m2)^(-sum(b[C][-I2[1:k]]) - sum(b[DC])-1) -(v[I2[k]])^(-sum(b[C][-I2[1:k]]) - sum(b[DC])-1)) / (sum(b[C][-I2[1:k]]) + sum(b[DC])+1)
    
    return((bit1+bit2+bit3)/EM)
  }
  
}




nll.powunif<-function(theta,y,u,a.ind,lam.ind,lamfix=FALSE, balthresh=FALSE)
{
  d<-dim(y)[2]
  a<-theta[a.ind]
  
  if(length(a)==1){a<-rep(a,d)}
  
  if(lamfix){lam<-rep(1,d)}
  else{
    lam<-c(theta[lam.ind],1)
  }
  
  if(balthresh){
    lam<-1/(1+a)
  }
  
  if(any(lam<0.01)||a<=0.01){return(10e10)}
  
  ind<-apply(y,1,comp.gt,u=u)
  y.uc<-y[ind,]
  y.pc<-y[!ind,]
  
  L<-apply(y.uc,1,fY.powunif,a=a,lam=lam)
  nll.uc<--sum(log(L))
  
  if(sum(!ind)>0)
  {
    L2<-apply(y.pc,1,fY.powunif.cens,a=a,lam=lam,u=u)
    nll.pc<--sum(log(L2))
  }
  else{nll.pc<-0}
  nll<-nll.uc+nll.pc
  return(nll)
}


# MGP scale
############


fX.powunif<-function(x,lam,a,sig,gamma)
{
  y<-BCi(x=x,gamma=gamma,sig=sig)
  J<-Jac(x=x,gamma=gamma,sig=sig)
  return(fY.powunif(y=y,a=a,lam=lam)*J)
}



fX.powunif.cens<-function(x,u,lam,a,sig,gamma)
{
  EM<-EM.pu(a=a,lam=lam)
  b<-1/a
  d<-length(x)
  C<-which(x<=u)
  D<-1:d
  DC<-D[-C]
  
  y<-BCi(x=x,gamma=gamma,sig=sig)
  uy<-BCi(x=u,gamma=gamma,sig=sig)
  
  J<-prod((1/sig[DC])*(1+gamma[DC]*x[DC]/sig[DC])^(1/gamma[DC]-1))
  
  
  m1<-max(uy[C]*lam[C])
  m2<-max(y[DC]*lam[DC])
  if(m1<m2)
  {
    num<-prod((uy[C]*lam[C])^b[C])*prod((lam[DC]*b[DC])*(y[DC]*lam[DC])^(b[DC]-1))
    # remove NA in max, as any X obs below -sig/gamma will cause NAs
    den<-EM*(sum(b)+1)*(max(y[DC]*lam[DC]))^(sum(b)+1)
    return(J*num/den)
  }
  else
  {
    v<-uy[C]*lam[C]
    I<-which(v>max(y[DC]*lam[DC]))
    
    k<-length(I)
    tmp<-sort(v,index.return=T)
    v2<-rev(tmp$x)
    I2<-rev(tmp$ix)
    
    QDC<-prod(lam[DC]*b[DC]*(y[DC]*lam[DC])^(b[DC]-1))
    
    bit1<-prod((uy[C]*lam[C])^b[C])*QDC / ((sum(b)+1)*(v2[1])^(sum(b)+1))
    
    if(k>1)
    {
      num<-den<-int<-NULL
      for(i in 1:(k-1))
      {
        num[i]<-prod((v[-I2[1:i]])^(b[C][-I2[1:i]])) * QDC
        den[i]<-sum(b[C][-I2[1:i]]) + sum(b[DC]) +1
        
        int[i]<-(v[I2[i+1]])^(-sum(b[C][-I2[1:i]]) - sum(b[DC]) -1) -(v[I2[i]])^(-sum(b[C][-I2[1:i]]) - sum(b[DC])-1) 
      }
      bit2<-sum(num*int/den)
    }
    else{bit2<-0}
    bit3<- prod((v[-I2[1:k]])^(b[C][-I2[1:k]])) * QDC * ((m2)^(-sum(b[C][-I2[1:k]]) - sum(b[DC])-1) -(v[I2[k]])^(-sum(b[C][-I2[1:k]]) - sum(b[DC]) -1)) / (sum(b[C][-I2[1:k]]) + sum(b[DC])+1)    
    
    return(J*(bit1+bit2+bit3)/EM)
  }
  
}



nll.powunif.GPD<-function(theta,x,u,a.ind,lam.ind,sig.ind,gamma.ind, lamfix=FALSE, balthresh=FALSE, marg.scale.ind,marg.shape.ind)
{
  d<-dim(x)[2]
  a<-theta[a.ind]
  if(length(a)==1)
  {
    a<-rep(a,d)
  }
  
  if(lamfix){lam<-rep(1,d)}
  else{
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
  for(j in 1:d)
  {
    rej[j]<-gamma[j]<0 && any(x[,j]>-sig[j]/gamma[j])
  }
  
  if(any(lam<0.01)||any(a<=0.01)||any(sig<=0.001)||any(rej)){return(10e10)}
  
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
  else{nll.pc<-0}
  nll<-nll.uc+nll.pc
  
  return(nll)
}

########################################################################################################
# Expected value of max(W_1,...W_d)

EM.pu<-function(a,lam)
{
  d<-length(a)
  tmp<-sort(lam,index.return=T)
  ls<-rev(tmp$x)
  ix<-rev(tmp$ix)
  as<-a[ix]
  
  bs<-1/as
  
  bit1<-1/ls[d] - ls[1]^(-sum(bs)-1) *prod(ls^bs)/(sum(bs)+1)
  
  bit2<-NULL
  for(i in 2:d)
  {
    bit2[i-1]<- - (ls[i]^(-sum(bs[i:d])-1) -ls[i-1]^(-sum(bs[i:d])-1))*prod(ls[i:d]^bs[i:d])/(sum(bs[i:d])+1)
  }
  
  return(bit1+sum(bit2))
}


#######################################################################################################


# Arguments:
# x - data (with >=1 exceedance above 0). Matrix of dim n by d.
# u - d-vector of thresholds at which to censor. Max value 0.
# std - logical. If TRUE then fit on standardized (Exp) scale, else fit on GPD scale [actual fitting done after tf to Pareto scale]
# dep.scale.fix - logical. If TRUE then where appropriate the variables in the dependence model are driven by a single scale parameter.
# dep.loc.fix - logical. If TRUE then where appropriate the variables in the dependence model are driven by a single location parameter.
# marg.scale.ind, marg.shape.ind - numerical vecotrs of length d denoting which margins should have common scale / shape pars 
#   (e.g. if d=3, say, c(1,1,1) denotes all identical; c(1,2,1), denotes 1st and 3rd the same, 2nd diff; c(1,2,3) denotes all different)
# balthresh - logical; if TRUE then scale parameter set so that E(e^(S_j)) = = const. for all j, thus only
#             appropriate when exceedance rates are the same in each margin.
# maxit - maximuma number of iterations allowed in optimization


#####################################################################################################################################



fit.MGPD.RevExpU<-function(x, u, std=FALSE, dep.scale.fix=FALSE, dep.loc.fix=FALSE, dep.start=NULL, marg.scale.start=NULL,
                   marg.shape.start=NULL, marg.scale.ind, marg.shape.ind, balthresh=FALSE, maxit=1000,...)
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
          opt<-optimize(nll.powunif, interval=c(par/10,10*par), u=u, y=x, a.ind=a.ind, lamfix=T, balthresh=balthresh,...)
        }
        else{
          a.ind<-1
          lam.ind<-2:d
          par<-c(dep.start[!is.null(dep.start)], c(2,rep(1,d-1))[is.null(dep.start)])
          opt<-optim(nll.powunif, par=par, u=u, y=x, a.ind=a.ind, lam.ind=lam.ind, lamfix=F, control=list(maxit=maxit,reltol=1e-6),balthresh=balthresh,...)
        }
      }
      else{
        if(dep.loc.fix)
        {
          par<-c(dep.start[!is.null(dep.start)], c(rep(2,d))[is.null(dep.start)])
          a.ind<-1:d
          opt<-optim(nll.powunif, par=par, u=u, y=x, a.ind=a.ind, lamfix=T, control=list(maxit=maxit,reltol=1e-6),balthresh=balthresh,...)
        }
        else
        {
          par<-c(dep.start[!is.null(dep.start)], c(rep(2,d),rep(1,d-1))[is.null(dep.start)])
          a.ind<-1:d
          lam.ind<-(d+1):(2*d-1)
          opt<-optim(nll.powunif, par=par, u=u, y=x, a.ind=a.ind, lam.ind=lam.ind, lamfix=F, control=list(maxit=maxit,reltol=1e-6),balthresh=balthresh,...)
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
          
          opt<-optim(nll.powunif.GPD, par=par, u=u, x=x, lamfix=T, a.ind=a.ind,
                     sig.ind=sig.ind, gamma.ind=gamma.ind, marg.scale.ind=marg.scale.ind, marg.shape.ind=marg.shape.ind, control=list(maxit=maxit,reltol=1e-6),balthresh=balthresh,...)
        }
        else{
          par<-c(dep.start[!is.null(dep.start)], c(2,rep(1,d-1))[is.null(dep.start)],
                 marg.scale.start[!is.null(marg.scale.start)], rep(1,n.sig)[is.null(marg.scale.start)],
                 marg.shape.start[!is.null(marg.shape.start)], rep(0.1,n.gamma)[is.null(marg.shape.start)])
          a.ind<-1
          lam.ind<-2:d
          sig.ind<-d+1:n.sig
          gamma.ind<-d+length(sig.ind)+1:n.gamma
          
          opt<-optim(nll.powunif.GPD, par=par, u=u, x=x, lamfix=F, a.ind=a.ind, lam.ind=lam.ind,
                     sig.ind=sig.ind, gamma.ind=gamma.ind, marg.scale.ind=marg.scale.ind, marg.shape.ind=marg.shape.ind, control=list(maxit=maxit,reltol=1e-6),balthresh=balthresh,...)
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
          
          opt<-optim(nll.powunif.GPD, par=par, u=u, x=x, lamfix=T, a.ind=a.ind,
                     sig.ind=sig.ind, gamma.ind=gamma.ind, marg.scale.ind=marg.scale.ind, marg.shape.ind=marg.shape.ind, control=list(maxit=maxit,reltol=1e-6),balthresh=balthresh,...)
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
          
          opt<-optim(nll.powunif.GPD, par=par, u=u, x=x, lamfix=F, a.ind=a.ind, lam.ind=lam.ind,
                     sig.ind=sig.ind, gamma.ind=gamma.ind, marg.scale.ind=marg.scale.ind, marg.shape.ind=marg.shape.ind, control=list(maxit=maxit,reltol=1e-6),balthresh=balthresh,...)
          
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



#####################################################################################################################################
#####################################################################################################################################


# sim.RevExpU.MGPD - Simulate from Frechet model

# n - number of vectors to simulate
# d - dimension
# a - dependence parameter a 
# beta - dependence parameter beta
# sig - marginal scale param
# gamma - marginal shape param
# std / MGPD  - logical: return standardized (Exp) / MGPD / both scales

sim.RevExpU.MGPD<-function(n,d,a,beta,sig,gamma,MGPD=TRUE,std=FALSE)
{
  lam<-exp(beta)
  Y<-matrix(0,ncol=d,nrow=n)
  C<-1/min(lam)
  k<-0
  while(k<n)
  {
    U<-runif(d,0,1)
    W<-(U^a)/lam
    R<-runif(1,0,C)
    Ind<-R<max(W)
    if(Ind)
    {
      k<-k+1
      Y[k,]<-W/R
    }
  }
  
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

