# ***************************
# ***************************
# Lídia André 
# March 2023
# Lancaster University, UK
# ***************************
# ***************************

rm(list=ls())
source("Simulation/Functions.R")
library(latex2exp)

# ====================
# Extremal Dependence
# ====================

theta<-seq(0.2,15,len=10)
thresh<-seq(0.7,1-.Machine$double.eps^0.5,len=20)
palette<-colorRampPalette(c("darkblue","red"))(length(theta))
legendtext<-TeX(paste("$\\theta=$",round(theta,2)))

## ct: asymptotically independent - Gumbel
## cb: asymptotically dependent - Gaussian
alphag<-1.2
rho<-0.5

survct<-1-2*thresh+pCopula(cbind(thresh,thresh),normalCopula(rho))
chict<-survct/(1-thresh)
etact<-log(1-thresh)/log(survct)

survcb<-1-2*thresh+pCopula(cbind(thresh,thresh),gumbelCopula(alphag))
chicb<-survcb/(1-thresh)
etacb<-log(1-thresh)/log(survcb)

chicstar<-etacstar<-NULL
for(i in 1:length(theta)){
  KM<-K(parct=rho,parcb=alphag,ct="gaus",cb="gumbel",
        weightfun=function(u,v,theta){pifun(u,v,theta)},parweight=theta[i])
  r<-quant(thresh,parct=rho,parcb=alphag,ct="gaus",cb="gumbel",
           weightfun=function(u,v,theta){pifun(u,v,theta)},parweight=theta[i])
  surv<-sapply(r,function(x){
    integral2(fun,x,1,x,1,parct=rho,parcb=alphag,ct="gaus",cb="gumbel",
              weightfun=function(u,v,theta){pifun(u,v,theta)},parweight=theta[i],K=KM)$Q
  })
  chicstar[[i]]<-surv/(1-thresh)
  etacstar[[i]]<-log(1-thresh)/log(surv)
}

par(mfrow=c(1,2))
plot(thresh,chict,"l",ylim=c(0,0.6),lwd=3,ylab=expression(chi),xlab="Thresholds")
for(i in 1:length(theta)){
  lines(thresh,chicstar[[i]],col=palette[i])
  lines(thresh,chicb,lty=2)
}
abline(h=0,lty=3) # \chi_t
text(x=0.75,y=0.02,expression(chi[t]),cex=1.5)
abline(h=2-2^(1/alphag),lty=3) # \chi_b
text(x=0.75,y=0.24,expression(chi[b]),cex=1.5)
legend("topright",c("ct","cb",legendtext),col=c(1,1,palette),lty=c(1,2,rep(1,10)),
       lwd=c(3,rep(1,11)),cex=0.5)
plot(thresh,etact,"l",ylim=c(0.45,1),lwd=3,ylab=expression(eta),xlab="Thresholds")
for(i in 1:length(theta)){
  lines(thresh,etacstar[[i]],col=palette[i])
  lines(thresh,etacb,lty=2)
}
legend("bottomright",c("ct","cb",legendtext),col=c(1,1,palette),lty=c(1,2,rep(1,20)),
       lwd=c(3,rep(1,21)),cex=0.45)
abline(h=0.75,lty=3) # \eta_t
text(x=0.75,y=0.77,expression(eta[t]),cex=1.5)
abline(h=1,lty=3) # \eta_b
text(x=0.75,y=0.98,expression(eta[b]),cex=1.5)
mtext("Case 3: Body Gumbel and Tail Gaussian", side = 3, line = -2.5, outer = TRUE,font=2)
par(mfrow=c(1,1))

# ======================
# Simulation 1
# Parameter Estimation
# ======================

## data simulation
n1<-n2<-50000000
n<-100000
rho<-0.6
alpha<-2
theta<-1
nsim<-100
len<-1000
thresh<-seq(0.01,0.99,len=100)

set.seed(123)
data_ct<-rCopula(n1,normalCopula(rho))
u_ct<-data_ct[,1]
v_ct<-data_ct[,2]
pi1_ct<-(u_ct*v_ct)^theta

set.seed(321)
data_cb<-rCopula(n2,joeCopula(alpha))
u_cb<-data_cb[,1]
v_cb<-data_cb[,2]
pi2_cb<-(u_cb*v_cb)^theta

keep_ct<-matrix(NA,nrow=n1,ncol=2)
for(i in 1:n1){
  if(runif(1)<pi1_ct[i]){
    keep_ct[i,]<-data_ct[i,]
  }
}
keep_ct<-keep_ct[complete.cases(keep_ct),]

keep_cb<-matrix(NA,nrow=n2,ncol=2)
for(i in 1:n2){
  if(runif(1)<(1-pi2_cb[i])){
    keep_cb[i,]<-data_cb[i,]
  }
}
keep_cb<-keep_cb[complete.cases(keep_cb),]

cstar<-rbind(keep_ct,keep_cb)
set.seed(123)
cstar<-cstar[sample(nrow(cstar),n),]

Kdata<-K(parct=rho,parcb=alpha,ct="gaus",cb="joe",
         weightfun=function(u,v,theta){pifun(u,v,theta)},parweight=theta)
datalist<-list()
for(i in 1:nsim){
  aux<-cstar[(1+(i-1)*len):(i*len),]
  datalist[[i]]<-apply(aux,2,marginal,parct=rho,parcb=alpha,ct="gaus",cb="joe",
                       weightfun=function(u,v,theta){pifun(u,v,theta)},parweight=theta,K=Kdata)
}

## done in parallel 
## example of one parallelisation
library(doParallel)
cl <- makeCluster(3)
registerDoParallel(cl)
outputM1 <- foreach(j=datalist[1:3],.packages=c("copula","pracma","evd")) %dopar% {
  source("Functions.R")
  rho<-0.6
  alpha<-2
  theta<-1
  start<-Sys.time()
  opt<-optim(loglik,par=c(rho,alpha,theta),x=j,ct="gaus",cb="joe",
             weightfun=function(u,v,theta){pifun(u,v,theta)})
  end<-Sys.time()-start
  return(list("Estimates"=opt$par,"Estimated Loglikelihood"=opt$value,"Time"=end))
}
stopCluster(cl)

## results with all parallelisations
estlog<-time<-NULL
est<-matrix(NA,ncol=3,nrow=length(outputM1))
for(i in 1:length(outputM1)){
  est[i,]<-outputM1[[i]][[1]]
  estlog[i]<-outputM1[[i]][[2]]
  time[i]<-outputM1[[i]][[3]]
}

timemin<-NULL
for(i in 1:length(time)){
  if(time[i]<3){
    timemin[i]<-round(time[i]*60,2)
  }
  else{
    timemin[i]<-round(time[i],2)
  }
}

# ====================
# Simulation 2
# Misspecification 1
# ====================

## data simulation
n<-100000
alpha<-2
nsim<-100
len<-1000

set.seed(123)
datajoe<-rCopula(n,joeCopula(alpha))

datalist<-list()
for(i in 1:nsim){
  datalist[[i]]<-datajoe[(1+(i-1)*len):(i*len),]
}

## done in parallel 
## example of one parallelisation
cl <- makeCluster(10)
registerDoParallel(cl)
outputsim1 <- foreach(j=datalist[1:10],.packages=c("copula","pracma","evd")) %dopar% {
  source("Functions.R")
  alpha<-2
  optTrue<-fitCopula(joeCopula(),j,method="ml",start=alpha,lower=1,upper=50,optim.method="Brent")
  aicTrue<-AICm(optTrue@estimate,optTrue@loglik)
  start<-Sys.time()
  optMixed<-optim(loglik1,par=c(alpha,0.5,1),x=j,ct="joe",cb="gaus",
                  weightfun=function(u,v,theta){pifun(u,v,theta)})
  end<-Sys.time()-start
  aicMixed<-AICm(optMixed$par,-optMixed$value)
  return(list("Estimates True Model"=optTrue@estimate,"Estimated Loglikelihood True Model"=-optTrue@loglik,"AIC True Model"=aicTrue,"Estimates Mixed Model"=optMixed$par,"Estimated Loglikelihood Mixed Model"=optMixed$value,"AIC Mixed Model"=aicMixed,"Time"=end))
}
stopCluster(cl)

## results with all parallelisations
esttrue<-estlogtrue<-aictrue<-estlogmixed<-aicmixed<-time<-NULL
estmixed<-matrix(NA,ncol=3,nrow=length(outputM1))
for(i in 1:length(outputM1)){
  esttrue[i]<-outputM1[[i]][[1]]
  estlogtrue[i]<-outputM1[[i]][[2]]
  aictrue[i]<-outputM1[[i]][[3]]
  estmixed[i,]<-outputM1[[i]][[4]]
  estlogmixed[i]<-outputM1[[i]][[5]]
  aicmixed[i]<-outputM1[[i]][[6]]
  time[i]<-outputM1[[i]][[7]]
}

# ====================
# Simulation 3
# Misspecification 2
# ====================


## data simulation
n<-100000
rho<-0.65
alpha<-2
nsim<-50
len<-1000
thresh<-seq(0.01,0.99,len=100)

set.seed(123)
datagaus<-rCopula(n,normalCopula(rho))

datalist<-list()
for(i in 1:nsim){
  datalist[[i]]<-datagaus[(1+(i-1)*len):(i*len),]
}

## done in parallel 
## example of one parallelisation of one model
cl <- makeCluster(10)
registerDoParallel(cl)
outputM1 <- foreach(j=datalist1[1:10],.packages=c("copula","pracma","evd")) %dopar% {
  source("Functions.R")
  rho<-0.65
  thresh<-seq(0.01,0.99,len=100)
  start<-Sys.time()
  optM1<-optim(loglik1,par=c(2,2,1),x=j,ct="joe",cb="frank",
               weightfun=function(u,v,theta){pifun(u,v,theta)})
  end<-Sys.time()-start
  KM1<-K(parct=optM1$par[1],parcb=optM1$par[2],
         ct="joe",cb="frank",weightfun=function(u,v,theta){pifun(u,v,theta)},parweight=optM1$par[3])
  survM1<-survivalf(thresh=thresh,parct=optM1$par[1],parcb=optM1$par[2],
                    ct="joe",cb="frank",weightfun=function(u,v,theta){pifun(u,v,theta)},parweight=optM1$par[3],K=KM1)
  aicM1<-AICm(optM1$par,-optM1$value)
  return(list("Estimates M1"=optM1$par,"Estimated Loglikelihood M1"=optM1$value,"AIC M1"=aicM1,"Time"=end,"KM1"=KM1,"Survival M1"=survM1))
}
stopCluster(cl)

## results with all parallelisations for one model
Kconst<-NULL
est<-matrix(NA,ncol=4,nrow=length(ooutputM1))
surv<-list()
for(i in 1:length(ooutputM1)){
  est[i,]<-ooutputM1[[i]][[1]]
  Kconst[i]<-ooutputM1[[i]][[5]]
  surv[[i]]<-ooutputM1[[i]][[6]]
}

### model \chi(r) and \eta(r)
chi0.65<-chi0.7<-chi0.75<-chi0.8<-chi0.85<-chi0.9<-chi0.95<-chi0.99<-eta0.65<-eta0.7<-eta0.75<-eta0.8<-eta0.85<-eta0.9<-eta0.95<-eta0.99<-c()
for(i in 1:length(ooutputM1)){
  chi0.65[i]<-chiu_model(thresh=thresh[66],survP=surv[[i]][66])
  chi0.7[i]<-chiu_model(thresh=thresh[71],survP=surv[[i]][71])
  chi0.75[i]<-chiu_model(thresh=thresh[76],survP=surv[[i]][76])
  chi0.8[i]<-chiu_model(thresh=thresh[81],survP=surv[[i]][81])
  chi0.85[i]<-chiu_model(thresh=thresh[86],survP=surv[[i]][86])
  chi0.9[i]<-chiu_model(thresh=thresh[91],survP=surv[[i]][91])
  chi0.95[i]<-chiu_model(thresh=thresh[96],survP=surv[[i]][96])
  chi0.99[i]<-chiu_model(thresh=thresh[100],survP=surv[[i]][100])
  eta0.65[i]<-etau_model(thresh=thresh[66],survP=surv[[i]][66])
  eta0.7[i]<-etau_model(thresh=thresh[71],survP=surv[[i]][71])
  eta0.75[i]<-etau_model(thresh=thresh[76],survP=surv[[i]][76])
  eta0.8[i]<-etau_model(thresh=thresh[81],survP=surv[[i]][81])
  eta0.85[i]<-etau_model(thresh=thresh[86],survP=surv[[i]][86])
  eta0.9[i]<-etau_model(thresh=thresh[91],survP=surv[[i]][91])
  eta0.95[i]<-etau_model(thresh=thresh[96],survP=surv[[i]][96])
  eta0.99[i]<-etau_model(thresh=thresh[100],survP=surv[[i]][100])
}

### \chi(r) for the gaussian with \rho=0.65
chi0.65geral<-(1-2*thresh[66]+pCopula(cbind(thresh[66],thresh[66]),normalCopula(rho)))/(1-thresh[66])
chi0.7geral<-(1-2*thresh[71]+pCopula(cbind(thresh[71],thresh[71]),normalCopula(rho)))/(1-thresh[71])
chi0.75geral<-(1-2*thresh[76]+pCopula(cbind(thresh[76],thresh[76]),normalCopula(rho)))/(1-thresh[76])
chi0.8geral<-(1-2*thresh[81]+pCopula(cbind(thresh[81],thresh[81]),normalCopula(rho)))/(1-thresh[81])
chi0.85geral<-(1-2*thresh[86]+pCopula(cbind(thresh[86],thresh[86]),normalCopula(rho)))/(1-thresh[86])
chi0.9geral<-(1-2*thresh[91]+pCopula(cbind(thresh[91],thresh[91]),normalCopula(rho)))/(1-thresh[91])
chi0.95geral<-(1-2*thresh[96]+pCopula(cbind(thresh[96],thresh[96]),normalCopula(rho)))/(1-thresh[96])
chi0.99geral<-(1-2*thresh[100]+pCopula(cbind(thresh[100],thresh[100]),normalCopula(rho)))/(1-thresh[100])

eta0.65geral<-(log(1-thresh[66]))/(log(1-2*thresh[66]+pCopula(cbind(thresh[66],thresh[66]),normalCopula(rho))))
eta0.7geral<-(log(1-thresh[71]))/(log(1-2*thresh[71]+pCopula(cbind(thresh[71],thresh[71]),normalCopula(rho))))
eta0.75geral<-(log(1-thresh[76]))/(log(1-2*thresh[76]+pCopula(cbind(thresh[76],thresh[76]),normalCopula(rho))))
eta0.8geral<-(log(1-thresh[81]))/(log(1-2*thresh[81]+pCopula(cbind(thresh[81],thresh[81]),normalCopula(rho))))
eta0.85geral<-(log(1-thresh[86]))/(log(1-2*thresh[86]+pCopula(cbind(thresh[86],thresh[86]),normalCopula(rho))))
eta0.9geral<-(log(1-thresh[91]))/(log(1-2*thresh[91]+pCopula(cbind(thresh[91],thresh[91]),normalCopula(rho))))
eta0.95geral<-(log(1-thresh[96]))/(log(1-2*thresh[96]+pCopula(cbind(thresh[96],thresh[96]),normalCopula(rho))))
eta0.99geral<-(log(1-thresh[100]))/(log(1-2*thresh[100]+pCopula(cbind(thresh[100],thresh[100]),normalCopula(rho))))

# model kendall's \tau
taumodel<-c()
nsim<-50000
n<-100000
for(i in 1:length(ooutputM1)){
  cstar<-sim(n=n,nsim=nsim,parct=est[i,1],parcb=c(est[i,2],est[i,3]),ct="InverGumbel",cb="t",
             weightfun=function(u,v,theta){pifun(u,v,theta)},parweight=est[i,4])
  taumodel[i]<-cor(cstar[,1],cstar[,2],method="kendall")
}
# kendall's \tau for the gaussian with \rho=0.65
kendalgeral<-2/pi*asin(rho)

# ====================
# Case Study: Example
# ====================

thresh<-seq(0.01,0.99,len=100) # thresholds at which \chi and \eta are computed
# initial parameters for loglikelihood optimisation for the HR and Clayton and weighting function
alpha<-theta<-2 
rho<-0.5 # initial parameter for loglikelihood optimisation for the Gaussian

# ---------------------------
# Single copula fit: Clayton
# ---------------------------

# data$u; data$v: data in uniform margins
claytonfit<-fitCopula(claytonCopula(),as.matrix(cbind(data$u,data$v)),method="ml",
                      start=alpha,lower=0,upper=50,optim.method="Brent") # mle
AICclayton<-AICm(claytonfit@estimate,claytonfit@loglik) 
survclayton<-1-2*thresh+pCopula(cbind(thresh,thresh),
                                claytonCopula(claytonfit@estimate)) # 1-F
chiclayton<-chiu_model(thresh,survclayton) # model \chi(r)
etaclayton<-etau_model(thresh,survclayton) # model \eta(r)

# -----------------------------------------------------------
# Weighted copula model fit: cb Gaussian and ct Husler-Reiss
# -----------------------------------------------------------

optimM1<-optim(loglik,par=c(alpha,rho,theta),x=as.data.frame(cbind(data$u,data$v)),
               ct="HR",cb="gaus",weightfun=function(u,v,theta){pifun(u,v,theta)},
               control=list(reltol=1e-4)) # mle
# normalising constant
KM1<-K(parct=optimM1$par[1],parcb=optimM1$par[2],ct="HR",cb="gaus",
       weightfun=function(u,v,theta){pifun1(u,v,theta)},parweight=optimM1$par[3])
survPM1<-survivalf(thresh=thresh,parct=optimM1$par[1],parcb=optimM1$par[2],
                   ct="HR",cb="gaus",weightfun=function(u,v,theta){pifun1(u,v,theta)},
                   parweight=optimM1$par[3],K=KM1) # 1-F
AICM1<-AICm(optimM1$par,-optimM1$value)
chiu_modelM1<-chiu_model(thresh=thresh,survP=survPM1) # model \chi(r)
etau_modelM1<-etau_model(thresh=thresh,survP=survPM1) # model \eta(r)

# ------------
# Diagnostics
# ------------

CompCI<-CIs(data=data,thresh=thresh,k=14) # confidence bands
CompChi<-CompCI[[1]] # for \chi(r)
CompEta<-CompCI[[2]] # for \eta(r)

chiemp<-chiu_emp(data$u,data$v,thresh) # empirical \chi(r)
etaemp<-etau_emp(data$u,data$v,thresh) # empirical \eta(r)

# rozone: threshold above which the GPD is fitted; here quantile(data$ozone,0.85)
# qozone: P[Ozone<=rozone]
# xiozone: mle for the shape parameter of the gpd
# sigmaozone: mle for the scale parameter of the gpd
# rtemperature: threshold above which the GPD is fitted; here quantile(data$temperature,0.9)
# qtemp: P[Temperature<=rtemperature]
# xitemp: mle for the shape parameter of the gpd
# sigmatemp: mle for the scale parameter of the gpd 
ozone100<-1-(1-qozone)*(max(0,(1+(xiozone*(100-rozone))/sigmaozone)))^(-1/xiozone)

temp16<-1/(dim(data)[1]+1)*length(data$temperature[data$temperature<=16])
temp23<-1-(1-qtemp)*(max(0,(1+(xitemp*(23-rtemperature))/sigmatemp)))^(-1/xitemp)

prob1M1<-integral2(cop,0,temp16,ozone100,1,parct=optimM1$par[1],parcb=optimM1$par[2],
                   ct="HR",cb="gaus",weightfun=function(u,v,theta){pifun1(u,v,theta)},
                   parweight=optimM1$par[3],K=KM1)
prob2M1<-integral2(cop,temp23,1,ozone100,1,parct=optimM1$par[1],parcb=optimM1$par[2],
                   ct="HR",cb="gaus",weightfun=function(u,v,theta){pifun1(u,v,theta)},
                   parweight=optimM1$par[3],K=KM1)
