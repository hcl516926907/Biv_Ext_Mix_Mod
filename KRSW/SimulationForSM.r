rm(list=ls(all=T))

source("Functions/CommonFunctions.r")
source("Functions/Gumbel_T_Functions.r")
source("Functions/Gumbel_U_Functions.r")
source("Functions/MVGaussian_T_Functions.r")
source("Functions/RevExp_T_Functions.r")
source("Functions/RevExp_U_Functions.r")


library(mvtnorm)
library(ismev)
library(Matrix)
library(fields)


# 
# ################################################################################################
# set.seed(123)
# ################################################################################################
# 
 n<-200
 d<-3
 a<-c(2,1,2)
 beta<-log(c(1,1,1)) # lam=exp(beta)
 
 nrep<-100
 
 parsGum<-matrix(0,nrow=nrep,ncol=length(a)+length(beta)-1)
 parsGauss<-matrix(0,nrow=nrep,ncol=length(a)+length(beta)-1)
 parsRevExp<-matrix(0,nrow=nrep,ncol=length(a)+length(beta)-1)
 nllGum<-nllGauss<-nllRevExp<-NULL
 conv<-matrix(0,nrow=nrep,ncol=3)
 
# Gumbel T data
 ###############################################

 for(s in 1:nrep)
 {
   Z<-sim.GumbelT.MGPD(n=n,d=d,a=a,beta=beta,std=T,MGPD=F)

   fitGum<-fit.MGPD.GumbelT(x=Z, u=rep(0,d), std=T, dep.scale.fix=F,maxit=10000)
   fitGauss<-fit.MGPD.MVGaussT(x=Z, u=rep(0,d), std=T, structured.cor=F,maxit=10000)
   fitRevExp<-fit.MGPD.RevExpT(x=Z, u=rep(0,d), std=T, dep.scale.fix=F,maxit=10000)

   parsGum[s,]<-fitGum$mle
   parsGauss[s,]<-fitGauss$mle
   parsRevExp[s,]<-fitRevExp$mle

   nllGum[s]<-fitGum$nll
   nllGauss[s]<-fitGauss$nll
   nllRevExp[s]<-fitRevExp$nll

   conv[s,]<-c(fitGum$conv,fitGauss$conv,fitRevExp$conv)

   print(s)
   print(parsGum[s,])
   print(parsGauss[s,])
   print(parsRevExp[s,])

   save.image("GumbelT.Rdata")
 }


 ################################################################################################
 set.seed(456)
 ################################################################################################


 # Rev Exp T data
 ###############################################

 for(s in 1:nrep)
 {
   Z<-sim.RevExpT.MGPD(n=n,d=d,a=a,beta=beta,std=T,MGPD=F)

   fitGum<-fit.MGPD.GumbelT(x=Z, u=rep(0,d), std=T, dep.scale.fix=F,maxit=10000)
   fitGauss<-fit.MGPD.MVGaussT(x=Z, u=rep(0,d), std=T, structured.cor=F,maxit=10000)
   fitRevExp<-fit.MGPD.RevExpT(x=Z, u=rep(0,d), std=T, dep.scale.fix=F,maxit=10000)

   parsGum[s,]<-fitGum$mle
   parsGauss[s,]<-fitGauss$mle
   parsRevExp[s,]<-fitRevExp$mle

   nllGum[s]<-fitGum$nll
   nllGauss[s]<-fitGauss$nll
   nllRevExp[s]<-fitRevExp$nll

   conv[s,]<-c(fitGum$conv,fitGauss$conv,fitRevExp$conv)

   print(s)
   print(parsGum[s,])
   print(parsGauss[s,])
   print(parsRevExp[s,])

   save.image("RevExpT.Rdata")
 }


 ################################################################################################
 set.seed(789)
 ################################################################################################


 rho<-c(-0.1,0.5,0.7)
 Sig<-matrix(c(1,rho[1:2],rho[1],1,rho[3],rho[2:3],1),3,3)

 # MV Gaussian T data
 ###############################################

 for(s in 1:nrep)
 {
   Z<-sim.MVGaussT.MGPD(n=n,d=d,Sig = Sig, beta=beta,std=T,MGPD=F)

   fitGum<-fit.MGPD.GumbelT(x=Z, u=rep(0,d), std=T, dep.scale.fix=F,maxit=10000)
   fitGauss<-fit.MGPD.MVGaussT(x=Z, u=rep(0,d), std=T, structured.cor=F,maxit=10000)
   fitRevExp<-fit.MGPD.RevExpT(x=Z, u=rep(0,d), std=T, dep.scale.fix=F,maxit=10000)

   parsGum[s,]<-fitGum$mle
   parsGauss[s,]<-fitGauss$mle
   parsRevExp[s,]<-fitRevExp$mle

   nllGum[s]<-fitGum$nll
   nllGauss[s]<-fitGauss$nll
   nllRevExp[s]<-fitRevExp$nll

   conv[s,]<-c(fitGum$conv,fitGauss$conv,fitRevExp$conv)

   print(s)
   print(parsGum[s,])
   print(parsGauss[s,])
   print(parsRevExp[s,])

   save.image("GaussT.Rdata")
 }



################################################################################################
set.seed(101112)
################################################################################################

 # With all positive correlations

rho<-c(0.1,0.5,0.7)
Sig<-matrix(c(1,rho[1:2],rho[1],1,rho[3],rho[2:3],1),3,3)

# log-Gaussian / Gaussian-max data
###############################################

for(s in 1:nrep)
{
  Z<-sim.MVGaussT.MGPD(n=n,d=d,Sig = Sig, beta=beta,std=T,MGPD=F)
  
  fitGum<-fit.MGPD.GumbelT(x=Z, u=rep(0,d), std=T, dep.scale.fix=F,maxit=10000)
  fitGauss<-fit.MGPD.MVGaussT(x=Z, u=rep(0,d), std=T, structured.cor=F,maxit=10000)
  fitRevExp<-fit.MGPD.RevExpT(x=Z, u=rep(0,d), std=T, dep.scale.fix=F,maxit=10000)
  
  parsGum[s,]<-fitGum$mle
  parsGauss[s,]<-fitGauss$mle
  parsRevExp[s,]<-fitRevExp$mle
  
  nllGum[s]<-fitGum$nll
  nllGauss[s]<-fitGauss$nll
  nllRevExp[s]<-fitRevExp$nll
  
  conv[s,]<-c(fitGum$conv,fitGauss$conv,fitRevExp$conv)
  
  print(s)
  print(parsGum[s,])
  print(parsGauss[s,])
  print(parsRevExp[s,])
  
  save.image("GaussT2.Rdata")
}

# Comment:
##############################################################################################################################

# In final simulation, convergence code = 10 for reps s=78 and s=88 for the Gumbel T model. These optimizations have been re-run 
# with starting values equal to the previous finishing values; results below. No change in reported results since:
# (i) Gauss T parameters unaffected
# (ii) Ordering of likelihoods remains the same

# s=78
#----------------------------------------
# 
#> fitGum<-fit.MGPD.GumbelT(x=Z, u=rep(0,d), std=T, dep.scale.fix=F,maxit=10000)
# fitGum
# $mle
# [1]  0.98660375  1.63904982  2.89446022 -0.37600337  0.05814264
# 
# $nll
# [1] 1000.722
# 
# $conv
# [1] 10
# 
# $hess
# NULL
# 
# $warn
# NULL

# fitGum1<-fit.MGPD.GumbelT(x=Z, u=rep(0,d), std=T, dep.scale.fix=F,maxit=10000,dep.start = fitGum$mle)
# > fitGum1
# $mle
# [1]  0.98484756  1.69249318  2.76191714 -0.38358126  0.07017427
# 
# $nll
# [1] 1000.685
# 
# $conv
# [1] 0
# 
# $hess
# NULL
# 
# $warn
# NULL

#----------------------------------------
# s=88
#----------------------------------------

# fitGum<-fit.MGPD.GumbelT(x=Z, u=rep(0,d), std=T, dep.scale.fix=F,maxit=10000)
# > fitGum
# $mle
# [1]  1.1582025  1.3805937  3.7978428 -0.2333193 -0.2034222
# 
# $nll
# [1] 960.1294
# 
# $conv
# [1] 10
# 
# $hess
# NULL
# 
# $warn
# NULL

# fitGum1<-fit.MGPD.GumbelT(x=Z, u=rep(0,d), std=T, dep.scale.fix=F,maxit=10000,dep.start = fitGum$mle)
# > fitGum1
# $mle
# [1]  1.1535859  1.3837609  3.8694092 -0.2364688 -0.2018865
# 
# $nll
# [1] 960.1289
# 
# $conv
# [1] 0
# 
# $hess
# NULL
# 
# $warn
# NULL
