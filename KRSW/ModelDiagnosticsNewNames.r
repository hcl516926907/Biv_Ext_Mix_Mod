# Functions to produce diagnostics for fitted MGPD models
#########################################################

# q.gpd - gpd quantile function

q.gpd<-function(q,sig,gamma)
{
  if(abs(gamma)>1e-6){return(sig*((1-q)^(-gamma)-1)/gamma)}
  else{return(-sig*log(1-q))}
}

# GPD.diag: marginal QQ plots for a matrix of MGPD data ("origin" should be 0, i.e. x[,j] | x[,j]>0 assumed GPD)

# GPD.diag.sim - uses simulation for CIs
# GPD.diag.th - uses theory for CIs


GPD.diag.sim<-function(x,marg.scale,marg.shape,nsim=1000)
{
  d<-dim(x)[2]
  n<-NULL
  sig<-marg.scale
  gam<-marg.shape
  par(mfrow=c(1,d))
  for(j in 1:d)
  {
    n[j]<-sum(x[,j]>0)
    tmp<-matrix(0,ncol=n[j],nrow=nsim)
    for(i in 1:nsim)
    {
      tmp[i,]<-sort(q.gpd(runif(n[j]),sig=sig[j],gamma=gam[j]))
    }
    up<-apply(tmp,2,quantile,0.975)
    low<-apply(tmp,2,quantile,0.025)
    
    plot(q.gpd((1:length(x[x[,j]>0,j]))/(length(x[x[,j]>0,j])+1),sig=sig[j],gamma=gam[j]),sort(x[x[,j]>0,j]),ylab="Empirical",xlab="Model",cex.lab=1.5,cex.axis=1.5)
    abline(a=0,b=1)
    lines(q.gpd((1:length(x[x[,j]>0,j]))/(length(x[x[,j]>0,j])+1),sig=sig[j],gamma=gam[j]),up,lty=2,col=4)
    lines(q.gpd((1:length(x[x[,j]>0,j]))/(length(x[x[,j]>0,j])+1),sig=sig[j],gamma=gam[j]),low,lty=2,col=4)
  }
}


GPD.diag.th<- function(x,sig,gam){
  d <- dim(x)[2]
  n <- dim(x)[1]
  par(mfrow=c(1,d))
  for(j in 1:d){
    n[j]<-sum(x[,j]>0)
    gpdQ <- q.gpd((1:n[j])/(n[j]+1),sig=sig[j],gamma=gam[j])
    par(cex.lab=1.5,cex.axis=1.5,cex.main=1.5,mar=c(5,4.4,4,2)+0.6)
    plot(gpdQ,sort(x[x[,j]>0,j]),ylab="Empirical",xlab="Model",cex.lab=2,cex.axis=2)
    abline(a=0,b=1,lwd=2)
    lowUnif <- sapply(c(1:n[j]), function(i) qbeta(0.025,i,n[j]+1-i))
    highUnif <- sapply(c(1:n[j]), function(i) qbeta(0.975,i,n[j]+1-i))
    lines(gpdQ,q.gpd(lowUnif, sig=sig[j],gamma=gam[j]),lty=2,lwd=2,col="gray")
    lines(gpdQ,q.gpd(highUnif, sig=sig[j],gamma=gam[j]),lty=2,lwd=2,col="gray")
  }
}

# MGPD.diag.***

# Create diagnostic plots for dependence fit and combined marginal / dependence features for model ***

# Arguments:
#   x - matrix of MGP data (at least one column contains an exceedance of 0)
#   a, beta, Sig - dependence parameters (or matrix for MVGauss case) from the fit
#   marg.scale / marg.shape - d-vectors of marginal scale / shape parameters; only required if want to use "g" feature
#   nsim - number of simulations to use to establish model features by Monte Carlo
#   chi - logical; compute and plot chi or not?
#   cols - which columns of x do you want chi for? defaults to all of them
#   q1 - minumum marginal quantile (between 0 and 1) to calculate chi above
#   k - number of points at which to calculate chi
#   g - user-defined function, which should take a d-vector as an argument. If supplied, will also produce a plot of P(g(X)>g.thresh) versus g.thresh;
#         model probability is black line, empirical is red dots. NB: does not undo conditioning, so probabilities are on the scale P(... | max(X_j)>0)
#   g.thresh - threshold sequence for P(g(X)>g.thresh) 
#   plot.g - logical; plot P(g(X)>g.thresh)?

MGPD.diag.GumbelU<-function(x,a,beta,marg.scale,marg.shape,nsim=10000,chi=TRUE,cols=NULL,q1=0.5,k=25,chiylabel,
                            g=NULL,g.thresh,plot.g,chiBS=TRUE,nbs=1000)
{
  d<-dim(x)[2]
  lam<-exp(beta)
  
  if(is.null(cols)){cols<-1:d}
  out<-NULL
  if(chi)
  {
    # use simulation to get model fitted value
    
    W<-W2<-matrix(0,ncol=d,nrow=nsim)
    for(i in 1:nsim)
    {
      W[i,]<-rfrechet(d,a=a,lam=lam) 
    }
    EW<-apply(W,2,mean)
    for(j in 1:d)
    {
      W2[,j]<-W[,j]/EW[j]
    }
    chi.est<-mean(apply(W2[,cols],1,min))
    
    # empirical + model
    chiPlot(data=x[,cols], ylabel=chiylabel, chimod=chi.est, nsim=nbs)
    
    out<-c(out,chi.est)
  }
  
  if(is.function(g))
  {
    Xmod<-sim.GumbelU.MGPD(nsim,d=d,a=a,beta=beta,sig=marg.scale,gamma=marg.shape)
    mod<-apply(Xmod,1,g)
    emp<-apply(x,1,g)
    
    mod.prob<-emp.prob<-NULL
    
    for(i in 1:length(g.thresh))
    {
      mod.prob[i]<-mean(mod>g.thresh[i])
      emp.prob[i]<-mean(emp>g.thresh[i])
    }
    
    if(plot.g)
    {
      plot(g.thresh,mod.prob,typ="l",xlab="z",ylab="P(g(X)>z)")
      points(g.thresh,emp.prob,col=2)
    }
    out<-list(out, g.thresh=g.thresh,mod.prob=mod.prob,emp.prob=emp.prob)
  }
  return(out)
}



MGPD.diag.GumbelT<-function(x,a,beta,marg.scale,marg.shape,nsim=10000,chi=TRUE,cols=NULL,q1=0.5,k=25,chiylabel,
                                 g=NULL,g.thresh,plot.g,chiBS=TRUE,nbs=1000)
{
  d<-dim(x)[2]
  lam<-exp(beta)
  
  if(is.null(cols)){cols<-1:d}
  out<-NULL
  if(chi)
  {
    # use simulation to get model fitted value
    
    W<-W2<-matrix(0,ncol=d,nrow=nsim)
    for(i in 1:nsim)
    {
      W[i,]<-rfrechet(d,a=a,lam=lam) 
      W[i,]<-W[i,]/max(W[i,])
    }
    EW<-apply(W,2,mean)
    for(j in 1:d)
    {
      W2[,j]<-W[,j]/EW[j]
    }
    chi.est<-mean(apply(W2[,cols],1,min))
    
    # empirical + model
    chiPlot(data=x[,cols], ylabel=chiylabel, chimod=chi.est, nsim=nbs)
    out<-c(out,chi.est)
  }
  
  if(is.function(g))
  {
    Xmod<-sim.GumbelT.MGPD(nsim,d=d,a=a,beta=beta,sig=marg.scale,gamma=marg.shape)
    mod<-apply(Xmod,1,g)
    emp<-apply(x,1,g)
    
    mod.prob<-emp.prob<-NULL
    
    for(i in 1:length(g.thresh))
    {
      mod.prob[i]<-mean(mod>g.thresh[i])
      emp.prob[i]<-mean(emp>g.thresh[i])
    }
    
    if(plot.g)
    {
      plot(g.thresh,mod.prob,typ="l",xlab="z",ylab="P(g(X)>z)")
      points(g.thresh,emp.prob,col=2)
    }
    out<-list(out, g.thresh=g.thresh,mod.prob=mod.prob,emp.prob=emp.prob)
  }
  return(out)
}



MGPD.diag.RevExpU<-function(x,a,beta,marg.scale,marg.shape,nsim=10000,chi=TRUE,cols=NULL,q1=0.5,k=25,chiylabel,
                            g=NULL,g.thresh,plot.g,chiBS=TRUE,nbs=1000)
{
  d<-dim(x)[2]
  lam<-exp(beta)
  
  if(is.null(cols)){cols<-1:d}
  out<-NULL
  if(chi)
  {

    # use simulation to get model fitted value
    
    W<-W2<-matrix(0,ncol=d,nrow=nsim)
    for(i in 1:nsim)
    {
      W[i,]<-runif(d)^a / lam 
    }
    EW<-apply(W,2,mean)
    for(j in 1:d)
    {
      W2[,j]<-W[,j]/EW[j]
    }
    chi.est<-mean(apply(W2[,cols],1,min))
    
    # empirical + model
    chiPlot(data=x[,cols], ylabel=chiylabel, chimod=chi.est, nsim=nbs)
    
        out<-c(out,chi.est)
  }
  
  if(is.function(g))
  {
    Xmod<-sim.RevExpU.MGPD(nsim,d=d,a=a,beta=beta,sig=marg.scale,gamma=marg.shape)
    mod<-apply(Xmod,1,g)
    emp<-apply(x,1,g)
    
    mod.prob<-emp.prob<-NULL
    
    for(i in 1:length(g.thresh))
    {
      mod.prob[i]<-mean(mod>g.thresh[i])
      emp.prob[i]<-mean(emp>g.thresh[i])
    }
    
    if(plot.g)
    {
      plot(g.thresh,mod.prob,typ="l",xlab="z",ylab="P(g(X)>z)")
      points(g.thresh,emp.prob,col=2)
    }
    out<-list(out, g.thresh=g.thresh,mod.prob=mod.prob,emp.prob=emp.prob)
  }
  return(out)
}



MGPD.diag.RevExpT<-function(x,a,beta,marg.scale,marg.shape,nsim=10000,chi=TRUE,cols=NULL,q1=0.5,k=25,chiylabel,
                            g=NULL,g.thresh,plot.g,chiBS=TRUE,nbs=1000)
{
  d<-dim(x)[2]
  lam<-exp(beta)
  
  if(is.null(cols)){cols<-1:d}
  out<-NULL
  if(chi)
  {
    # use simulation to get model fitted value
    
    W<-W2<-matrix(0,ncol=d,nrow=nsim)
    for(i in 1:nsim)
    {
      W[i,]<-runif(d)^a / lam 
      W[i,]<-W[i,]/max(W[i,])
    }
    EW<-apply(W,2,mean)
    for(j in 1:d)
    {
      W2[,j]<-W[,j]/EW[j]
    }
    chi.est<-mean(apply(W2[,cols],1,min))
    
    # empirical + model
    chiPlot(data=x[,cols], ylabel=chiylabel, chimod=chi.est, nsim=nbs)
    
    out<-c(out,chi.est)
  }
  
  if(is.function(g))
  {
    Xmod<-sim.RevExpT.MGPD(nsim,d=d,a=a,beta=beta,sig=marg.scale,gamma=marg.shape)
    mod<-apply(Xmod,1,g)
    emp<-apply(x,1,g)
    
    mod.prob<-emp.prob<-NULL
    
    for(i in 1:length(g.thresh))
    {
      mod.prob[i]<-mean(mod>g.thresh[i])
      emp.prob[i]<-mean(emp>g.thresh[i])
    }
    
    if(plot.g)
    {
      plot(g.thresh,mod.prob,typ="l",xlab="z",ylab="P(g(X)>z)")
      points(g.thresh,emp.prob,col=2)
    }
    out<-list(out, g.thresh=g.thresh,mod.prob=mod.prob,emp.prob=emp.prob)
  }
  return(out)
}

MGPD.diag.MVGaussT<-function(x,Sig,beta,marg.scale,marg.shape,nsim=10000,chi=TRUE,cols=NULL,q1=0.5,k=25,chiylabel,
                             g=NULL,g.thresh,plot.g,chiBS=TRUE,nbs=1000)
{
  d<-dim(x)[2]
  if(is.null(cols)){cols<-1:d}
  out<-NULL
  if(chi)
  {
    # use simulation to get model fitted value
    
    W<-W2<-matrix(0,ncol=d,nrow=nsim)
    for(i in 1:nsim)
    {
      W[i,]<-exp(rmvnorm(1,sigma=Sig,mean=beta))
      W[i,]<-W[i,]/max(W[i,])
    }
    EW<-apply(W,2,mean)
    for(j in 1:d)
    {
      W2[,j]<-W[,j]/EW[j]
    }
    chi.est<-mean(apply(W2[,cols],1,min))
    
    # empirical + model
    chiPlot(data=x[,cols], ylabel=chiylabel, chimod=chi.est, nsim=nbs)
    
    out<-c(out,chi.est)
  }
  
  if(is.function(g))
  {
    Xmod<-sim.MVGauss.MGPD(nsim,d=d,Sig=Sig,beta=beta,sig=marg.scale,gamma=marg.shape)
    mod<-apply(Xmod,1,g)
    emp<-apply(x,1,g)
    
    mod.prob<-emp.prob<-NULL
    
    for(i in 1:length(g.thresh))
    {
      mod.prob[i]<-mean(mod>g.thresh[i])
      emp.prob[i]<-mean(emp>g.thresh[i])
    }
    
    if(plot.g)
    {
      plot(g.thresh,mod.prob,typ="l",xlab="z",ylab="P(g(X)>z)")
      points(g.thresh,emp.prob,col=2)
    }
    out<-list(out, g.thresh=g.thresh,mod.prob=mod.prob,emp.prob=emp.prob)
  }
  return(out)
}

########################################################################


chiPlot <- function(data, ylabel, chimod, nsim, nq = 35, qmin = 0.5, qmax = 0.99){ #chimod is the model-based value
  tmp <- matrix(,nrow=nsim,ncol=nq)
  n <- nrow(data)
  for(j in 1:nsim){
    nsample <- sample(1:n,size=n,replace=T)
    newdata <- data[nsample,]
    tmp[j,] <- chiEmp(newdata,nq=nq,qmin=qmin,qmax=qmax)[,2] 
  }
  CIlow <- apply(tmp,2,quantile,0.025)
  CIhigh <- apply(tmp,2,quantile,0.975)
  chi <-chiEmp(data,nq=nq,qmin=qmin,qmax=qmax)
  
  par(cex.lab=2,cex.axis=2,cex.main=2,mar=c(5,4.4,4,2)+0.9)
  plot(chi[,1],chi[,2],ylim=c(0,1),xlab="q",
       ylab=ylabel,lwd=2)
  lines(chi[,1],CIlow,lty=3,lwd=2)
  lines(chi[,1],CIhigh,lty=3,lwd=2)
  abline(h = chimod, lwd = 2)
}

chiEmp <-function(data,nq=25,qmin, qmax){
  n<-nrow(data)
  datatr<-(apply(data,2,rank) - 0.5)/n   
  qlim<-c(qmin,qmax)
  u <- seq(qlim[1], qlim[2], length = nq)
  cu<-sapply(c(1:nq), function(i) mean(apply(datatr,1,min) >= u[i]))
  return(cbind(u,cu/(1-u)))
}



# unif - use ranks to transform to uniform scale

unif<-function(x){rank(x)/(length(x)+1)}
