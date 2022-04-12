# Example use of functions for Gumbel T and U models, Reverse exponential T and U models and Gaussian T model
#---------------------------------------------------------------------------------------------------------------

source("Functions/Gumbel_T_Functions.r")
source("Functions/Gumbel_U_Functions.r")
source("Functions/RevExp_T_Functions.r")
source("Functions/RevExp_U_Functions.r")
source("Functions/MVGaussian_T_Functions.r")
source("Functions/CommonFunctions.r")

library(mvtnorm)
library(Matrix)
library(ismev)


# Gumbel T
###########

set.seed(15)

d<-2
a<-c(1.5) # common shape (quicker estimation)
beta<-c(0,0)
sig<-c(2,2)
gamma<-c(0.1,0.1)

# Final beta parameter fixed at zero in estimation, so needs to be relative to this

# Simulate with (conditionally) GP(sig,gamma) and (conditionally) exponential margins 
X<-sim.GumbelT.MGPD(n=1000,d=d, a=a, beta=beta, sig=sig, gamma=gamma, MGPD = T,std=T)

#  GP(sig,gamma)
plot(X$X, pch=20, xlab=expression(X[1]),ylab=expression(X[2]))

# Exponential
plot(X$Z, pch=20, xlab=expression(Z[1]),ylab=expression(Z[2]))

# Fit standard (expl) form; censored likelihood censoring at (0,0); single a parameter
fitS<-fit.MGPD.GumbelT(x=X$Z,u=c(0,0),std=T, dep.scale.fix = T)
# Parameter order: (a,beta1)

# Uncensored estimation
u1<-min(X$Z)-0.01
fitS2<-fit.MGPD.GumbelT(x=X$Z, u=c(u1,u1), std=T, dep.scale.fix = T)
# Parameter order: (a,beta1)

# Fit GP form; censored likelihood censoring at (0,0). Single a and fixed beta parameters
# marg.scale.ind = c(1,1) states single common scale parameter should be fitted (for different scale parameters set marg.scale.ind = c(1,2))
# marg.shape.ind = c(1,1) states single common shape parameter should be fitted ( " " )


fitGP<-fit.MGPD.GumbelT(x=X$X, u=c(0,0), std=F, dep.scale.fix = T, dep.loc.fix = T, marg.scale.ind = c(1,1), marg.shape.ind = c(1,1))

# Parameter order: (a,sig,gamma)

#-----------------------------------------------------------------------------------------------------------------------

set.seed(15)

d<-2
a<-c(2,4) # Different a parameters means numerical integration in the likelihood (slower)
beta<-c(0.5,0)


# Final beta parameter fixed at zero in estimation, so needs to be relative to this

# Simulate with (conditionally) GP(sig,gamma) and (conditionally) exponential margins 
Z<-sim.GumbelT.MGPD(n=1000,d=d, a=a, beta=beta, MGPD = F,std=T)

# Exponential
plot(Z, pch=20, xlab=expression(Z[1]),ylab=expression(Z[2]))

# Fit standard (expl) form; censored likelihood censoring at (0,0)
fitS<-fit.MGPD.GumbelT(x=X$Z,u=c(0,0),std=T)
# Parameter order: (a1,a2,beta1)

# Uncensored estimation
u1<-min(X$Z)-0.01
fitS2<-fit.MGPD.GumbelT(x=X$Z, u=c(u1,u1), std=T)
# Parameter order: (a1,a2,beta1)


#-----------------------------------------------------------------------------------------------------------------------

# Gumbel U
###########

set.seed(15)

d<-2
a<-c(1.5) # common shape (quicker estimation)
beta<-c(0,0)
sig<-c(2,2)
gamma<-c(0.1,0.1)

# Final beta parameter fixed at zero in estimation, so needs to be relative to this

# Simulate with (conditionally) GP(sig,gamma) and (conditionally) exponential margins 
X<-sim.GumbelU.MGPD(n=1000,d=d, a=a, beta=beta, sig=sig, gamma=gamma, MGPD = T,std=T)

#  GP(sig,gamma)
plot(X$X, pch=20, xlab=expression(X[1]),ylab=expression(X[2]))

# Exponential
plot(X$Z, pch=20, xlab=expression(Z[1]),ylab=expression(Z[2]))

# Fit standard (expl) form; censored likelihood censoring at (0,0); single a parameter
fitS<-fit.MGPD.GumbelU(x=X$Z,u=c(0,0),std=T, dep.scale.fix = T)
# Parameter order: (a,beta1)

# Uncensored estimation
u1<-min(X$Z)-0.01
fitS2<-fit.MGPD.GumbelU(x=X$Z, u=c(u1,u1), std=T, dep.scale.fix = T)
# Parameter order: (a,beta1)

# Fit GP form; censored likelihood censoring at (0,0). Single a and fixed beta parameters

fitGP<-fit.MGPD.GumbelU(x=X$X, u=c(0,0), std=F, dep.scale.fix = T, dep.loc.fix = T, marg.scale.ind = c(1,1), marg.shape.ind = c(1,1), maxit=5000)
# Parameter order: (a,sig,gamma)

#-----------------------------------------------------------------------------------------------------------------------

# Reverse Exponential T

set.seed(15)

d<-2
a<-c(2,4) # NB this is 1/\alpha, where \alpha is the parameterization in Kiriliouk, Rootzen, Segers and Wadsworth. 
beta<-c(0.5,0)
sig<-c(2,2)
gamma<-c(0.1,-0.1)

# Final beta parameter fixed at zero in estimation, so needs to be relative to this

# Simulate with (conditionally) GP(sig,gamma) and (conditionally) exponential margins 
X<-sim.RevExpT.MGPD(n=1000,d=d, a=a, beta=beta, sig=sig, gamma=gamma, MGPD = T,std=T)

#  GP(sig,gamma)
plot(X$X, pch=20, xlab=expression(X[1]),ylab=expression(X[2]))

# Exponential
plot(X$Z, pch=20, xlab=expression(Z[1]),ylab=expression(Z[2]))

# Fit standard (expl) form; censored likelihood censoring at (0,0)
fitS<-fit.MGPD.RevExpT(x=X$Z,u=c(0,0),std=T)
# Parameter order: (a1,a2,beta1)

# Uncensored estimation
u1<-min(X$Z)-0.01
fitS2<-fit.MGPD.RevExpT(x=X$Z, u=c(u1,u1), std=T)
# Parameter order: (a1,a2,beta1)

# Fit GP form; censored likelihood censoring at (0,0)
# marg.scale.ind = 1:2 states two scale parameters should be fitted (for common scale set marg.scale.ind = c(1,1))
# marg.shape.ind = 1:2 states two shape parameters should be fitted ( " " )

fitGP<-fit.MGPD.RevExpT(x=X$X, u=c(0,0), std=F, marg.scale.ind = 1:2, marg.shape.ind = 1:2, maxit=5000)

# Parameter order: (a1,a2,beta1,sig1,sig2,gamma1,gamma2)

# Simpler parameterization with Hessian
# Delta method needed for Var-Cov matrix as optimization (and hence Hessian calculation) done on scale lambda=exp(beta)
# for the dependence location parameters

fitGP<-fit.MGPD.RevExpT(x=X$X, u=c(0,0), std=F, marg.scale.ind = c(1,1), marg.shape.ind = 1:2, maxit=5000, hessian=T)

# Parameter order: (a1,a2,beta1,sig,gamma1,gamma2)


# Gradient matrix
# Transformation is the idendity except for component 3 (beta)
# beta = t(lambda) = log(lambda), so derivative is 1/lambda = exp(-beta)
GM<-diag(6)
GM[3,3]<-exp(-fitGP$mle[3])

# Approx Var-Cov matrix
GM%*% solve(fitGP$hess)%*% t(GM)

#-----------------------------------------------------------------------------------------------------------------------

# Reverse Exponential U

set.seed(15)

d<-2
a<-c(2,4) # NB this is 1/\alpha, where \alpha is the parameterization in Kiriliouk, Rootzen, Segers and Wadsworth. 
beta<-c(0.5,0)
sig<-c(2,2)
gamma<-c(0.1,-0.1)

# Final beta parameter fixed at zero in estimation, so needs to be relative to this

# Simulate with (conditionally) GP(sig,gamma) and (conditionally) exponential margins 
X<-sim.RevExpU.MGPD(n=1000,d=d, a=a, beta=beta, sig=sig, gamma=gamma, MGPD = T,std=T)

#  GP(sig,gamma)
plot(X$X, pch=20, xlab=expression(X[1]),ylab=expression(X[2]))

# Exponential
plot(X$Z, pch=20, xlab=expression(Z[1]),ylab=expression(Z[2]))

# Fit standard (expl) form; censored likelihood censoring at (0,0)
fitS<-fit.MGPD.RevExpU(x=X$Z,u=c(0,0),std=T)
# Parameter order: (a1,a2,beta1)

# Uncensored estimation
u1<-min(X$Z)-0.01
fitS2<-fit.MGPD.RevExpU(x=X$Z, u=c(u1,u1), std=T)
# Parameter order: (a1,a2,beta1)

# Fit GP form; censored likelihood censoring at (0,0)
# marg.scale.ind = 1:2 states two scale parameters should be fitted (for common scale set marg.scale.ind = c(1,1))
# marg.shape.ind = 1:2 states two shape parameters should be fitted ( " " )

fitGP<-fit.MGPD.RevExpU(x=X$X, u=c(0,0), std=F, marg.scale.ind = c(1,1), marg.shape.ind = 1:2, maxit=5000)
# Parameter order: (a1,a2,beta1,sig,gamma1,gamma2)

#-----------------------------------------------------------------------------------------------------------------------

# MV Gaussian T

set.seed(15)

d<-2
rho<-c(0.6)
Sig<-matrix(c(1,rho,rho,1),2,2)
beta<-c(0.5,0)
sig<-c(2,2)
gamma<-c(0.1,-0.1)

# Final beta parameter fixed at zero in estimation, so needs to be relative to this

# Simulate with (conditionally) GP(sig,gamma) and (conditionally) exponential margins 
X<-sim.MVGaussT.MGPD(n=1000,d=d, Sig= Sig, beta=beta, sig=sig, gamma=gamma, MGPD = T,std=T)

#  GP(sig,gamma)
plot(X$X, pch=20, xlab=expression(X[1]),ylab=expression(X[2]))

# Exponential
plot(X$Z, pch=20, xlab=expression(Z[1]),ylab=expression(Z[2]))

# Fit standard (expl) form; censored likelihood censoring at (0,0)
# structured.cor=F means the correlation matrix is unstructured (for higher dimensions possibility of structured matrix)
fitS<-fit.MGPD.MVGaussT(x=X$Z,u=c(0,0),std=T, structured.cor=F)
# Parameter order: (rho,beta1)

# Uncensored estimation
u1<-min(X$Z)-0.01
fitS2<-fit.MGPD.MVGaussT(x=X$Z, u=c(u1,u1), std=T, structured.cor = F)
# Parameter order: (rho,beta1)

# Fit GP form; censored likelihood censoring at (0,0)
# marg.scale.ind = 1:2 states two scale parameters should be fitted (for common scale set marg.scale.ind = c(1,1))
# marg.shape.ind = 1:2 states two shape parameters should be fitted ( " " )


fitGP<-fit.MGPD.MVGaussT(x=X$X, u=c(0,0), std=F, marg.scale.ind = c(1,1), marg.shape.ind = 1:2, structured.cor = F, maxit=5000)
# Parameter order: (rho,beta1,sig,gamma1,gamma2)

#-----------------------------------------------------------------------------------------------------------------------
