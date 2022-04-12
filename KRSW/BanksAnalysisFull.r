###############################################################################################
# R script to reproduce analyses of UK Bank data
###############################################################################################


HSBC<-scan("Data/HSBCW_Adj_Oct16.txt")
HSBCDates<-read.table("Data/HSBCWDates_Adj_Oct16.txt")
HSBCDates<-as.character(HSBCDates[,1])
HSBCDates<-as.Date(c(HSBCDates),format="%Y-%m-%d")
plot(HSBCDates,HSBC)
# Most recent observations are first, so this is -(Z_{t}-Z_{t-1})/Z_{t-1} = 1-Z_{t}/Z_{t-1}
HSBCR<--(HSBC[1:(length(HSBC)-1)]-HSBC[2:(length(HSBC))])/HSBC[2:(length(HSBC))]
HSBCDates<-HSBCDates[-length(HSBCDates)]
plot(HSBCDates,HSBCR)


Lloyds<-scan("Data/LloydsW_Adj_Oct16.txt")
LloydsDates<-read.table("Data/LloydsWDates_Adj_Oct16.txt")
LloydsDates<-as.character(LloydsDates[,1])
LloydsDates<-as.Date(c(LloydsDates),format="%Y-%m-%d")
plot(LloydsDates,Lloyds)
LloydsR<--(Lloyds[1:(length(Lloyds)-1)]-Lloyds[2:(length(Lloyds))])/Lloyds[2:(length(Lloyds))]
LloydsDates<-LloydsDates[-length(LloydsDates)]
plot(LloydsDates,LloydsR)


RBS<-scan("Data/RBSW_Adj_Oct16.txt")
RBSDates<-read.table("Data/RBSWDates_Adj_Oct16.txt")
RBSDates<-as.character(RBSDates[,1])
RBSDates<-as.Date(c(RBSDates),format="%Y-%m-%d")
plot(RBSDates,RBS)
RBSR<--(RBS[1:(length(RBS)-1)]-RBS[2:(length(RBS))])/RBS[2:(length(RBS))]
RBSDates<-RBSDates[-length(RBSDates)]
plot(RBSDates,RBSR)


Barclays<-scan("Data/BarclaysW_Adj_Oct16.txt")
BarclaysDates<-read.table("Data/BarclaysWDates_Adj_Oct16.txt")
BarclaysDates<-as.character(BarclaysDates[,1])
BarclaysDates<-as.Date(c(BarclaysDates),format="%Y-%m-%d")
plot(BarclaysDates,Barclays)
BarclaysR<--(Barclays[1:(length(Barclays)-1)]-Barclays[2:(length(Barclays))])/Barclays[2:(length(Barclays))]
BarclaysDates<-BarclaysDates[-length(BarclaysDates)]
plot(BarclaysDates,BarclaysR)


# To get contemporaneous series:

AllDates<-intersect(HSBCDates,LloydsDates)
AllDates<-intersect(AllDates,RBSDates)
AllDates<-intersect(AllDates,BarclaysDates)

AllDates<-as.POSIXct(AllDates*24*60*60,origin = "1970-01-01")
AllDates<-as.Date(AllDates,format="%Y-%m-%d %H:%M:%S")
is.element(HSBCDates[1],AllDates)

date.logic<-function(date1,dates2)
{
  is.element(date1,dates2)
}

# check
sum(sapply(HSBCDates,date.logic,dates2=AllDates))

Banks<-cbind(HSBCR[sapply(HSBCDates,date.logic,dates2=AllDates)],LloydsR[sapply(LloydsDates,date.logic,dates2=AllDates)],
             RBSR[sapply(RBSDates,date.logic,dates2=AllDates)], BarclaysR[sapply(BarclaysDates,date.logic,dates2=AllDates)])



dim(Banks)
plot(Banks[,1:2])
plot(Banks[,c(1,3)])
plot(Banks[,c(1,4)])
plot(Banks[,2:3])
plot(Banks[,c(2,4)])
plot(Banks[,3:4])

###################################################
# Pairwise plots
###################################################

par(mar=c(4,6,2,2))

pdf("BanksHL.pdf")
par(mar=c(5,6,2,0.5))
plot(Banks[,1:2], xlab=expression(Y[H]),ylab=expression(Y[L]),cex.axis=2,cex.lab=2,pch=20)
dev.off()

pdf("BanksHR.pdf")
par(mar=c(5,6,2,0.5))
plot(Banks[,c(1,3)], xlab=expression(Y[H]),ylab=expression(Y[R]),cex.axis=2,cex.lab=2,pch=20)
dev.off()

pdf("BanksHB.pdf")
par(mar=c(5,6,2,0.5))
plot(Banks[,c(1,4)], xlab=expression(Y[H]),ylab=expression(Y[B]),cex.axis=2,cex.lab=2,pch=20)
dev.off()

pdf("BanksLR.pdf")
par(mar=c(5,6,2,0.5))
plot(Banks[,2:3], xlab=expression(Y[L]),ylab=expression(Y[R]),cex.axis=2,cex.lab=2,pch=20)
dev.off()

pdf("BanksLB.pdf")
par(mar=c(5,6,2,0.5))
plot(Banks[,c(2,4)], xlab=expression(Y[L]),ylab=expression(Y[B]),cex.axis=2,cex.lab=2,pch=20)
dev.off()

pdf("BanksRB.pdf")
par(mar=c(5,6,2,0.5))
plot(Banks[,3:4], xlab=expression(Y[R]),ylab=expression(Y[B]),cex.axis=2,cex.lab=2,pch=20)
dev.off()



#############################################################


source("Functions/ModelDiagnosticsNewNames.r")

# 4-dim chi
chiPlot(data=Banks, ylabel=expression(chi[HLRB]~(q)), chimod=NULL, nsim=1000, nq = 35, qmin = 0.5, qmax = 0.99)

# 3-dim chis
chiPlot(data=Banks[,1:3], ylabel=expression(chi[HLR]~(q)), chimod=NULL, nsim=1000, nq = 35, qmin = 0.5, qmax = 0.99)
chiPlot(data=Banks[,2:4], ylabel=expression(chi[LRB]~(q)), chimod=NULL, nsim=1000, nq = 35, qmin = 0.5, qmax = 0.99)
chiPlot(data=Banks[,c(1,3,4)], ylabel=expression(chi[HRB]~(q)), chimod=NULL, nsim=1000, nq = 35, qmin = 0.5, qmax = 0.99)
chiPlot(data=Banks[c(1,2,4)], ylabel=expression(chi[HLB]~(q)), chimod=NULL, nsim=1000, nq = 35, qmin = 0.5, qmax = 0.99)

# 2-dim chis
chiPlot(data=Banks[,1:2], ylabel=expression(chi[HL]~(q)), chimod=NULL, nsim=1000, nq = 35, qmin = 0.5, qmax = 0.99)
chiPlot(data=Banks[,c(1,3)], ylabel=expression(chi[HR]~(q)), chimod=NULL, nsim=1000, nq = 35, qmin = 0.5, qmax = 0.99)
chiPlot(data=Banks[,c(1,4)], ylabel=expression(chi[HB]~(q)), chimod=NULL, nsim=1000, nq = 35, qmin = 0.5, qmax = 0.99)
chiPlot(data=Banks[,c(2,3)], ylabel=expression(chi[LR]~(q)), chimod=NULL, nsim=1000, nq = 35, qmin = 0.5, qmax = 0.99)
chiPlot(data=Banks[,c(2,4)], ylabel=expression(chi[LB]~(q)), chimod=NULL, nsim=1000, nq = 35, qmin = 0.5, qmax = 0.99)
chiPlot(data=Banks[,c(3,4)], ylabel=expression(chi[RB]~(q)), chimod=NULL, nsim=1000, nq = 35, qmin = 0.5, qmax = 0.99)

# produce 4-dim chi figure
pdf("ChiAll.pdf")
chiPlot(data=Banks, ylabel=expression(hat(chi)[HLRB]~(q)), chimod=NULL, nsim=1000, nq = 35, qmin = 0.5, qmax = 0.99)
dev.off()

#########################################################################
#########################################################################

# Initially transform and fit on standardized scale. NB code for standardized data is written on standard Pareto scale (exponential
# of standard exponential scale), but works exactly the same as if it were on standard exponential scale.

x<-Banks

# Use empirical probability integral transform to put on uniform scale

u1<-unif(x=x[,1])
hist(u1)

u2<-unif(x=x[,2])
hist(u2)

u3<-unif(x=x[,3])
hist(u3)

u4<-unif(x=x[,4])
hist(u4)

# Transform to standard Pareto scale

xpb<-cbind(1/(1-u1),1/(1-u2),1/(1-u3),1/(1-u4))


# Create matrix of data on MVP scale (at least one threshold exc)
# (exponential of Y_E-u_E|Y_E \not\leq u_E)

logic<-xpb[,1]>quantile(xpb[,1],0.83)|xpb[,2]>quantile(xpb[,2],0.83)|xpb[,3]>quantile(xpb[,3],0.83)|xpb[,4]>quantile(xpb[,4],0.83)
xpb2<-cbind(xpb[logic,1]/quantile(xpb[,1],0.83),xpb[logic,2]/quantile(xpb[,2],0.83),xpb[logic,3]/quantile(xpb[,3],0.83),xpb[logic,4]/quantile(xpb[,4],0.83))
plot(log(xpb2))
dim(xpb2)

# exponential margins
xeb2<-log(xpb2)


# Same but on MGPD scale (i.e. scale of the observations)

Banks2<-cbind(Banks[logic,1]-quantile(Banks[,1],0.83),Banks[logic,2]-quantile(Banks[,2],0.83),Banks[logic,3]-quantile(Banks[,3],0.83),Banks[logic,4]-quantile(Banks[,4],0.83))
plot(Banks2)
dim(Banks2)


# Initially examine "most complex" dependence models to home in on best family of models
########################################################################################
source("Functions/CommonFunctions.r")
source("Functions/Gumbel_T_Functions.r")
source("Functions/MVGaussian_T_Functions.r")
source("Functions/RevExp_T_Functions.r")
source("Functions/RevExp_U_Functions.r")
source("Functions/Gumbel_U_Functions.r")

#######################################################################################################
# NB: all fits below are written twice, the second version with starting values obtained from the first 
# fit. For simple models, the optimum will be obtained with a single optimization, but for others, a few
# iterations is necessary.
#######################################################################################################

# h_T based on the Gumbel

fit1<-fit.MGPD.GumbelT(x=xeb2, u=rep(0,4), std=T, dep.scale.fix=F, maxit=2000)
fit1
fit1<-fit.MGPD.GumbelT(x=xeb2, u=rep(0,4), std=T, dep.scale.fix=F, maxit=2000, dep.start=fit1$mle)
fit1

fit1.1<-fit.MGPD.GumbelT(x=xeb2, u=rep(0,4), std=T, dep.scale.fix=T, maxit=2000)
fit1.1
fit1.1<-fit.MGPD.GumbelT(x=xeb2, u=rep(0,4), std=T, dep.scale.fix=T, maxit=2000,dep.start=fit1.1$mle)
fit1.1


#####################################################################################################
#library(geoR)
library(mvtnorm)
library(Matrix)

# h_T based on the Gaussian

fit2<-fit.MGPD.MVGaussT(x=xeb2, u=rep(0,4), std=T, structured.cor=F)
# Run some repeats of this:
fit2<-fit.MGPD.MVGaussT(x=xeb2, u=rep(0,4), std=T, structured.cor=F, dep.start=fit2$mle,maxit=2000)


# Issues with "degeneracy of Nelder-Mead simplex"; BFGS gives apparent convergence
fit2<-fit.MGPD.MVGaussT(x=xeb2, u=rep(0,4), std=T, structured.cor=F, dep.start=fit2$mle, method="BFGS")

# Without free Gaussian location parameters
fit2.1<-fit.MGPD.MVGaussT(x=xeb2, u=rep(0,4), std=T, structured.cor=F, dep.start=fit2$mle[1:6],dep.loc.fix=T,maxit=2000)
fit2.1<-fit.MGPD.MVGaussT(x=xeb2, u=rep(0,4), std=T, structured.cor=F, dep.start=fit2.1$mle,
                               dep.loc.fix=T,maxit=2000,method="BFGS")


#####################################################################################################
# h_T based on the reverse exponential

fit3<-fit.MGPD.RevExpT(x=xeb2, u=rep(0,4), std=T, dep.scale.fix=F)
fit3
fit3<-fit.MGPD.RevExpT(x=xeb2, u=rep(0,4), std=T, dep.scale.fix=F, dep.start = fit3$mle)
fit3

fit3.1<-fit.MGPD.RevExpT(x=xeb2, u=rep(0,4), std=T, dep.scale.fix=T)
fit3.1
fit3.1<-fit.MGPD.RevExpT(x=xeb2, u=rep(0,4), std=T, dep.scale.fix=T, dep.start = fit3.1$mle)
fit3.1


#####################################################################################################
# h_U based on the reverse exponential

fit4<-fit.MGPD.RevExpU(x=xeb2, u=rep(0,4), std=T, dep.scale.fix=F)
fit4
fit4<-fit.MGPD.RevExpU(x=xeb2, u=rep(0,4), std=T, dep.scale.fix=F,dep.start = fit4$mle)
fit4

fit4.1<-fit.MGPD.RevExpU(x=xeb2, u=rep(0,4), std=T, dep.scale.fix=T)
fit4.1
fit4.1<-fit.MGPD.RevExpU(x=xeb2, u=rep(0,4), std=T, dep.scale.fix=T,dep.start = fit4.1$mle)
fit4.1


#####################################################################################################
# h_U based on the Gumbel

fit5<-fit.MGPD.GumbelU(x=xeb2, u=rep(0,4), std=T, dep.scale.fix=F)
fit5
fit5<-fit.MGPD.GumbelU(x=xeb2, u=rep(0,4), std=T, dep.scale.fix=F, dep.start=fit5$mle, maxit=2000)
fit5

fit5.1<-fit.MGPD.GumbelU(x=xeb2, u=rep(0,4), std=T, dep.scale.fix=T)
fit5.1
fit5.1<-fit.MGPD.GumbelU(x=xeb2, u=rep(0,4), std=T, dep.scale.fix=T,dep.start = fit5.1$mle)
fit5.1

########################################################################################

# Minimized log-likelihoods and AICs for most complicated models

fit1$nll
fit2$nll
fit3$nll
fit4$nll
fit5$nll

fit1$nll+2*length(fit1$mle)
fit2$nll+2*length(fit2$mle)
fit3$nll+2*length(fit3$mle)
fit4$nll+2*length(fit4$mle)
fit5$nll+2*length(fit5$mle)

########################################################################################

# Gumbel model based on h_T is preferred

fit1
fit1.1

# Further parameterizations of this model

# Free scale parameter, constrained location parameters
fit1.2<-fit.MGPD.GumbelT(x=xeb2, u=rep(0,4), std=T, dep.scale.fix=F, dep.loc.fix=T)
fit1.2
fit1.2<-fit.MGPD.GumbelT(x=xeb2, u=rep(0,4), std=T, dep.scale.fix=F, dep.loc.fix=T, dep.start = fit1.2$mle)
fit1.2

# Single scale parameter, fixed location
fit1.3<-fit.MGPD.GumbelT(x=xeb2, u=rep(0,4), std=T, dep.scale.fix=T, dep.loc.fix=T)
fit1.3
fit1.3<-fit.MGPD.GumbelT(x=xeb2, u=rep(0,4), std=T, dep.scale.fix=T, dep.loc.fix=T, dep.start = fit1.3$mle)
fit1.3


# NLLs to 1d.p.

# \alpha1,2,3,4; \beta1,2,3
round(fit1$nll,1)
# \alpha; \beta1,2,3
round(fit1.1$nll,1)
# \alpha1,2,3,4
round(fit1.2$nll,1)
# \alpha
round(fit1.3$nll,1)

# Likelihood ratio tests

# Sequence: 1->1.2->1.3
1-pchisq(2*(fit1.2$nll-fit1$nll), df=3)
1-pchisq(2*(fit1.3$nll-fit1.2$nll), df=3)

# Sequence 1->1.1->1.3
1-pchisq(2*(fit1.1$nll-fit1$nll), df=3)
1-pchisq(2*(fit1.3$nll-fit1.1$nll), df=3)

1-pchisq(2*(fit1.3$nll-fit1$nll), df=6)

#################################################################################################
# Marginal parameter stability plots to understand if 0.83 quantile might be OK
#################################################################################################

library(ismev)

plot(Banks[,1])
quantile(Banks[,1],0.83)
gpd.fitrange(Banks[,1],umin=0,umax=0.1,nint=20)
gpd.fit(Banks[,1],thresh=0.06)
m1<-gpd.fit(Banks[,1],thresh=quantile(Banks[,1],0.83))
gpd.diag(m1)

plot(Banks[,2])
quantile(Banks[,2],0.83)
gpd.fitrange(Banks[,2],umin=0,umax=0.2,nint=20)
gpd.fit(Banks[,2],thresh=0.06)
m2<-gpd.fit(Banks[,2],thresh=quantile(Banks[,2],0.83))
gpd.diag(m2)

plot(Banks[,3])
quantile(Banks[,3],0.83)
gpd.fitrange(Banks[,3],umin=0,umax=0.2,nint=20)
gpd.fit(Banks[,3],thresh=0.1)
m3<-gpd.fit(Banks[,3],thresh=quantile(Banks[,3],0.83))
gpd.diag(m3)

plot(Banks[,4])
quantile(Banks[,4],0.83)
gpd.fitrange(Banks[,4],umin=0,umax=0.2,nint=20)
gpd.fit(Banks[,4],thresh=0.06)
m4<-gpd.fit(Banks[,4],thresh=quantile(Banks[,4],0.83))
gpd.diag(m4)

#################################################################################################
# Fit margins and dependence simulataneously
#################################################################################################

fit1.4<-fit.MGPD.GumbelT(x=Banks2, u=rep(0,4), std=F,
                              dep.scale.fix=T, dep.loc.fix=T, marg.shape.ind=1:4, marg.scale.ind=1:4,
                              dep.start=fit1.3$mle[1],maxit=5000)

# Need to run the below a few times:
fit1.4<-fit.MGPD.GumbelT(x=Banks2, u=rep(0,4), std=F,
                              dep.scale.fix=T, dep.loc.fix=T, marg.shape.ind=1:4, marg.scale.ind=1:4,
                              marg.scale.start=fit1.4$mle[2:5], marg.shape.start=fit1.4$mle[6:9],
                              dep.start=fit1.4$mle[1],maxit=5000)
fit1.4

# Test for common shape

fit1.5<-fit.MGPD.GumbelT(x=Banks2, u=rep(0,4), std=F,
                              dep.scale.fix=T, dep.loc.fix=T, marg.shape.ind=rep(1,4), marg.scale.ind=1:4,
                              marg.scale.start=fit1.4$mle[2:5], marg.shape.start=mean(fit1.4$mle[6:9]),
                              dep.start=fit1.4$mle[1],maxit=5000)

fit1.5<-fit.MGPD.GumbelT(x=Banks2, u=rep(0,4), std=F,
                              dep.scale.fix=T, dep.loc.fix=T, marg.shape.ind=rep(1,4), marg.scale.ind=1:4,
                              marg.scale.start=fit1.5$mle[2:5], marg.shape.start=mean(fit1.5$mle[6]),
                              dep.start=fit1.5$mle[1],maxit=5000)
fit1.5

1-pchisq(2*(fit1.5$nll-fit1.4$nll),df=3)


# Final model: fit1.5 (Single dependence parameter, combined marginal shape, different marginal scale)

# With Hessian, to get standard errors

fit1.5.1<-fit.MGPD.GumbelT(x=Banks2, u=rep(0,4), std=F,
                                dep.scale.fix=T, dep.loc.fix=T, marg.shape.ind=rep(1,4), marg.scale.ind=1:4,
                                marg.scale.start=fit1.5$mle[2:5], marg.shape.start=mean(fit1.5$mle[6]),
                                dep.start=fit1.4$mle[1],maxit=5000,hessian=T)

solve(fit1.5.1$hess)
sqrt(diag(solve(fit1.5.1$hess)))

# rounded parameter estimates and standard errors
round(fit1.5.1$mle,4)
round(fit1.5$mle,4)

round(sqrt(diag(solve(fit1.5.1$hess))),4)

# Compare marginal and jointly estimated standard errors

sqrt(diag(solve(fit1.5.1$hess)))[2]/m1$se[1]
sqrt(diag(solve(fit1.5.1$hess)))[3]/m2$se[1]
sqrt(diag(solve(fit1.5.1$hess)))[4]/m3$se[1]
sqrt(diag(solve(fit1.5.1$hess)))[5]/m4$se[1]

###############################################################################################
# Diagnostics:
###############################################################################################

# Will print fitted 4-dim chi (caluculated using Monte Carlo)
# Rounded value is 0.4 (checked also with nsim = 100000)


MGPD.diag.GumbelT(x=Banks,a = fit1.5$mle[1],beta = rep(0,4),marg.scale = fit1.5$mle[2:5],marg.shape = fit1.5$mle[6],nsim = 10000,
                       chi=T,chiBS=T,nbs=2000,chiylabel = expression(hat(chi)[HLRB]~(q)))



###############################################################################################
# Bivariate chi plots with fitted 
###############################################################################################
pdf("ChiHL-BS.pdf")
MGPD.diag.GumbelT(x=Banks,a = fit1.5$mle[1],beta = rep(0,4),marg.scale = fit1.5$mle[2:5],marg.shape = fit1.5$mle[6],nsim = 10000,
                       chi=T, cols=c(1,2),chiBS=T,nbs=2000,chiylabel = expression(hat(chi)[HL]~(q)))
abline(v=0.83,col=4)
dev.off()

pdf("ChiHR-BS.pdf")
MGPD.diag.GumbelT(x=Banks,a = fit1.5$mle[1],beta = rep(0,4),marg.scale = fit1.5$mle[2:5],marg.shape = fit1.5$mle[6],nsim = 10000,
                       chi=T, cols=c(1,3),chiBS=T,nbs=2000,chiylabel = expression(hat(chi)[HR]~(q)))
abline(v=0.83,col=4)
dev.off()

pdf("ChiHB-BS.pdf")
MGPD.diag.GumbelT(x=Banks,a = fit1.5$mle[1],beta = rep(0,4),marg.scale = fit1.5$mle[2:5],marg.shape = fit1.5$mle[6],nsim = 10000,
                       chi=T, cols=c(1,4),chiBS=T,nbs=2000,chiylabel = expression(hat(chi)[HB]~(q)))
abline(v=0.83,col=4)
dev.off()

pdf("ChiLR-BS.pdf")
MGPD.diag.GumbelT(x=Banks,a = fit1.5$mle[1],beta = rep(0,4),marg.scale = fit1.5$mle[2:5],marg.shape = fit1.5$mle[6],nsim = 10000,
                       chi=T, cols=c(2,3),chiBS=T,nbs=2000,chiylabel = expression(hat(chi)[LR]~(q)))
abline(v=0.83,col=4)
dev.off()

pdf("ChiLB-BS.pdf")
MGPD.diag.GumbelT(x=Banks,a = fit1.5$mle[1],beta = rep(0,4),marg.scale = fit1.5$mle[2:5],marg.shape = fit1.5$mle[6],nsim = 10000,
                       chi=T, cols=c(2,4),chiBS=T,nbs=2000,chiylabel = expression(hat(chi)[LB]~(q)))
abline(v=0.83,col=4)
dev.off()

pdf("ChiRB-BS.pdf")
MGPD.diag.GumbelT(x=Banks,a = fit1.5$mle[1],beta = rep(0,4),marg.scale = fit1.5$mle[2:5],marg.shape = fit1.5$mle[6],nsim = 10000,
                       chi=T, cols=c(3,4),chiBS=T,nbs=2000,chiylabel = expression(hat(chi)[RB]~(q)))
abline(v=0.83,col=4)
dev.off()

pdf("ChiHLR-BS.pdf")
MGPD.diag.GumbelT(x=Banks,a = fit1.5$mle[1],beta = rep(0,4),marg.scale = fit1.5$mle[2:5],marg.shape = fit1.5$mle[6],nsim = 10000,
                       chi=T, cols=c(1,2,3),chiBS=T,nbs=2000,chiylabel = expression(chi[HLR]~(q)))
abline(v=0.83,col=4)
dev.off()


###############################################################################################
# Trivariate plots
###############################################################################################

pdf("ChiHRB-BS.pdf")
MGPD.diag.GumbelT(x=Banks,a = fit1.5$mle[1],beta = rep(0,4),marg.scale = fit1.5$mle[2:5],marg.shape = fit1.5$mle[6],nsim = 10000,
                       chi=T, cols=c(1,2,4),chiBS=T,nbs=2000,chiylabel = expression(chi[HRB]~(q)))
abline(v=0.83,col=4)
dev.off()

pdf("ChiHLB-BS.pdf")
MGPD.diag.GumbelT(x=Banks,a = fit1.5$mle[1],beta = rep(0,4),marg.scale = fit1.5$mle[2:5],marg.shape = fit1.5$mle[6],nsim = 10000,
                       chi=T, cols=c(1,3,4),chiBS=T,nbs=2000,chiylabel = expression(chi[HLB]~(q)))
abline(v=0.83,col=4)
dev.off()

pdf("ChiLRB-BS.pdf")
MGPD.diag.GumbelT(x=Banks,a = fit1.5$mle[1],beta = rep(0,4),marg.scale = fit1.5$mle[2:5],marg.shape = fit1.5$mle[6],nsim = 10000,
                       chi=T, cols=c(2,3,4),chiBS=T,nbs=2000,chiylabel = expression(chi[LRB]~(q)))
abline(v=0.83,col=4)
dev.off()

###############################################################################################
# Marginal GPD QQ plots
###############################################################################################
pdf("QQBank.pdf",width=14,height=3.5)
GPD.diag.th(Banks2,sig=fit1.5$mle[2:5],gam=rep(fit1.5$mle[6],4))
dev.off()

###############################################################################################
# Sum of components with common shape should also be GP above 0, with scale = sum(scale pars) and same shape
###############################################################################################

BS<-apply(Banks2,1,sum)

###############################################################################################
fitsum<-gpd.fit(BS,thresh=0.0)
gpd.diag(fitsum)

pdf("PosSum_MF.pdf")
GPD.diag.th(x=matrix(BS,ncol=1),sig=fitsum$mle[1],gam=fitsum$mle[2])
dev.off()
pdf("PosSum_JF.pdf")
GPD.diag.th(x=matrix(BS,ncol=1),sig=sum(fit1.5$mle[2:5]),gam=fit1.5$mle[6])
dev.off()

# log-likelihoods
library(evd)
sum(dgpd(BS[BS>0],scale=fitsum$mle[1],shape=fitsum$mle[2],log=T))
sum(dgpd(BS[BS>0],scale=sum(fit1.5$mle[2:5]),shape=fit1.5$mle[6],log=T))

# Compare univariate and implied multivariate MLEs
fitsum$mle
sum(fit1.5$mle[2:5])
fit1.5$mle[6]

# Standard errors

# Univariate
fitsum$se

# Multivariate (first on is via delta method)
sqrt(t(c(0,1,1,1,1,0))%*%solve(fit1.5.1$hess)%*%c(0,1,1,1,1,0))
sqrt(diag(solve(fit1.5.1$hess)))[6]


###############################################################################################
# Comparison of VaR / ES of weighted sums from model versus univariate approach 

# VaR and ES (not including threshold sum(a*u))

VaR<-function(p,sig,gam,phi)
{
  if(abs(gam)>1e-6){return(sig*((p/phi)^(-gam)-1)/gam)}
  else{return(-sig*log(p/phi))}
}

ES<-function(p,sig,gam,phi)
{
  v<-VaR(p,sig,gam,phi)
  ES<-v+(sig+gam*v)/(1-gam)
  return(ES)
}

# VaR and ES (including threshold sum(a*u))

VaR2<-function(p,sig,gam,phi,th)
{
  if(abs(gam)>1e-6){return(sig*((p/phi)^(-gam)-1)/gam+th)}
  else{return(-sig*log(p/phi)+th)}
}

ES2<-function(p,sig,gam,phi,th)
{
  v<-VaR2(p,sig,gam,phi,th=th)
  ES<-v+(sig+gam*(v-th))/(1-gam)
  return(ES)
}
##################################################
# Exceedance probability

# (BSa defined below)
phi<-(dim(Banks2)[1]/dim(Banks)[1])*mean(BSa>0)
###########
# same as:
tmp<-apply(Banks,1,ws,a=a)
u<-apply(Banks,2,quantile,0.83)
mean(tmp-sum(a*u)>0)
###########

#############################################################################
# Standard errors for VaR and SE via delta method
#############################################################################

SE.VaR<-function(cov,p,sig,gam,phi)
{
  dv.dsig<-((p/phi)^(-gam)-1)/gam
  dv.dgam<--(sig/gam^2)*((p/phi)^(-gam)-1)-(sig/gam)*(log(p/phi))*(p/phi)^(-gam)
  dv.dphi<-sig*phi^(gam-1)*p^(-gam)
  grad<-c(dv.dsig,dv.dgam,dv.dphi)
  
  return(sqrt(t(grad)%*%cov%*%grad))
}


SE.ES<-function(cov,p,sig,gam,phi)
{
  dE.dsig<-(((p/phi)^(-gam)-1)/gam)/(1-gam) + 1/(1-gam)
  dv.dgam<--(sig/gam^2)*((p/phi)^(-gam)-1)-(sig/gam)*(log(p/phi))*(p/phi)^(-gam)
  dE.dgam<- VaR(p=p,sig=sig,gam=gam,phi=phi)/(1-gam)^2 + dv.dgam/(1-gam) + (sig)/(1-gam)^2
  dE.dphi<-(sig*phi^(gam-1)*p^(-gam))/(1-gam)
  grad<-c(dE.dsig,dE.dgam,dE.dphi)
  
  return(sqrt(t(grad)%*%cov%*%grad))
}


SE.VaR.MV<-function(cov,p,sig,gam,phi,a)
{
  dv.dsig<-a*((p/phi)^(-gam)-1)/gam
  dv.dgam<--(sig/gam^2)*((p/phi)^(-gam)-1)-(sig/gam)*(log(p/phi))*(p/phi)^(-gam)
  dv.dphi<-sig*phi^(gam-1)*p^(-gam)
  grad<-c(dv.dsig,dv.dgam,dv.dphi)
  
  return(sqrt(t(grad)%*%cov%*%grad))
}


SE.ES.MV<-function(cov,p,sig,gam,phi,a)
{
  dE.dsig<-(a*((p/phi)^(-gam)-1)/gam)/(1-gam) + a/(1-gam)
  dv.dgam<--(sig/gam^2)*((p/phi)^(-gam)-1)-(sig/gam)*(log(p/phi))*(p/phi)^(-gam)
  dE.dgam<- VaR(p=p,sig=sig,gam=gam,phi=phi)/(1-gam)^2 + dv.dgam/(1-gam) + (sig)/(1-gam)^2
  # same as
  #dE.dgam<- VaR2(p=p,sig=sig,gam=gam,phi=phi)/(1-gam)^2 + dv.dgam/(1-gam) + (sig-sum(a*u))/(1-gam)^2
  dE.dphi<-sig*phi^(gam-1)*p^(-gam)/(1-gam)
  grad<-c(dE.dsig,dE.dgam,dE.dphi)
  
  return(sqrt(t(grad)%*%cov%*%grad))
}


#############################################################################################
# Analysis with portfolio weights a=(10,20,30,40)
#############################################################################################

a<-c(10,20,30,40)
ws<-function(x,a){sum(a*x)}
BSa<-apply(Banks2,1,ws,a=a)
fitsuma<-gpd.fit(BSa,thresh=0)

# add in threshold
u<-apply(Banks,2,quantile,0.83)
sum(a*u)

# Re-definition exceedance probability
phi<-(dim(Banks2)[1]/dim(Banks)[1])*mean(BSa>0)

# Covariance matrix for univariate parameters
uni.cov<-matrix(0,3,3)
uni.cov[1:2,1:2]<-fitsuma$cov
uni.cov[3,3]<-(phi*(1-phi))/dim(Banks)[1]

# Covariance matrix for multivariate parameters
multi.cov<-matrix(0,6,6)
multi.cov[1:5,1:5]<-solve(fit1.5.1$hess)[2:6,2:6]
multi.cov[6,6]<-(phi*(1-phi))/dim(Banks)[1]


#############################################################################################
# plot VaR(p) vs -logp with CIs
#############################################################################################

pseq<-seq(0.15,0.001,len=100)

VaRMV<-VaR2(p=pseq,sig=sum(a*fit1.5$mle[2:5]),gam=fit1.5$mle[6],phi=phi,th=sum(a*u))
SEVaRMV<-sapply(pseq,SE.VaR.MV,sig=sum(a*fit1.5$mle[2:5]),gam=fit1.5$mle[6],phi=phi,a=a,cov=multi.cov)

VaRUV<-VaR2(p=pseq,sig=fitsuma$mle[1],gam=fitsuma$mle[2],phi=phi,th=sum(a*u))
SEVaRUV<-sapply(pseq,SE.VaR,sig=fitsuma$mle[1],gam=fitsuma$mle[2],phi=phi,cov=uni.cov)


pdf("VaRMV1234.pdf",)
par(mar=c(5,5,4,2))
plot(-log(pseq),VaRMV,typ="l",axes=F,xlab="p",ylab="VaR(p)",ylim=c(0,70),cex.lab=2)
axis(1,at=-log(c(0.15,0.1,0.01,0.001)),labels=(c(0.15,0.1,0.01,0.001)),cex.axis=2)
axis(2,at=c(0,10,20,30,40,50,60),cex.axis=2)
lines( -log(pseq),VaRMV-1.96*SEVaRMV,col=2,lty=2)
lines( -log(pseq),VaRMV+1.96*SEVaRMV,col=2,lty=2)
points(-log(phi*(length(BSa[BSa>0])+1-rank(BSa[BSa>0]))/(length(BSa[BSa>0])+1)), (BSa[BSa>0])+sum(a*u))
title("a=(10,20,30,40)",cex.main=2)
dev.off()


pdf("VaRUV1234.pdf")
par(mar=c(5,5,4,2))
plot(-log(pseq),VaRUV,typ="l",axes=F,xlab="p",ylab="VaR(p)",ylim=c(0,70),cex.lab=2)
axis(1,at=-log(c(0.15,0.1,0.01,0.001)),labels=(c(0.15,0.1,0.01,0.001)),cex.axis=2)
axis(2,at=c(0,10,20,30,40,50,60),cex.axis=2)
lines( -log(pseq),VaRUV-1.96*SEVaRUV,col=2,lty=2)
lines( -log(pseq),VaRUV+1.96*SEVaRUV,col=2,lty=2)
points(-log(phi*(length(BSa[BSa>0])+1-rank(BSa[BSa>0]))/(length(BSa[BSa>0])+1)), (BSa[BSa>0])+sum(a*u))
title("a=(10,20,30,40)",cex.main=2)
dev.off()


#############################################################################################
# Expected Shortfall plots
#############################################################################################

ESMV<-ES2(p=pseq,sig=sum(a*fit1.5$mle[2:5]),gam=fit1.5$mle[6],phi=phi,th=sum(a*u))
SEESMV<-sapply(pseq,SE.ES.MV,sig=sum(a*fit1.5$mle[2:5]),gam=fit1.5$mle[6],phi=phi,a=a,cov=multi.cov)

ESUV<-ES2(p=pseq,sig=fitsuma$mle[1],gam=fitsuma$mle[2],phi=phi,th=sum(a*u))
SEESUV<-sapply(pseq,SE.ES,sig=fitsuma$mle[1],gam=fitsuma$mle[2],phi=phi,cov=uni.cov)


Emp.ES<-NULL
BSa0<-BSa[BSa>0]
pseq2<-(1:length(BSa0))/(length(BSa0)+1)
Emp.ES<-cumsum(rev(sort(BSa0)))/(1:length(BSa0))+sum(a*u)


pdf("ESMV1234.pdf")
par(mar=c(5,5,4,2))
plot(-log(pseq),ESMV,typ="l",axes=F,xlab="p",ylab="ES(p)",ylim=c(0,150),cex.lab=2)
axis(1,at=-log(c(0.15,0.1,0.01,0.001)),labels=(c(0.15,0.1,0.01,0.001)),cex.axis=2)
axis(2,at=2.5*c(0,10,20,30,40,50,60),cex.axis=2)
lines( -log(pseq),ESMV-1.96*SEESMV,col=2,lty=2)
lines( -log(pseq),ESMV+1.96*SEESMV,col=2,lty=2)
points(-log(pseq2*phi), Emp.ES)
title("a=(10,20,30,40)",cex.main=2)
dev.off()


pdf("ESUV1234.pdf")
par(mar=c(5,5,4,2))
plot(-log(pseq),ESUV,typ="l",axes=F,xlab="p",ylab="ES(p)",ylim=c(0,150),cex.lab=2)
axis(1,at=-log(c(0.15,0.1,0.01,0.001)),labels=(c(0.15,0.1,0.01,0.001)),cex.axis=2)
axis(2,at=2.5*c(0,10,20,30,40,50,60),cex.axis=2)
lines( -log(pseq),ESUV-1.96*SEESUV,col=2,lty=2)
lines( -log(pseq),ESUV+1.96*SEESUV,col=2,lty=2)
points(-log(pseq2*phi), Emp.ES)
title("a=(10,20,30,40)",cex.main=2)
dev.off()




#############################################################################################
# Analysis with portfolio weights a=(40,30,20,10)
#############################################################################################


a<-c(40,30,20,10)
ws<-function(x,a){sum(a*x)}
BSa<-apply(Banks2,1,ws,a=a)
fitsuma<-gpd.fit(BSa,thresh=0)

# add in threshold
u<-apply(Banks,2,quantile,0.83)
sum(a*u)

# Re-definition of exc probability
phi<-(dim(Banks2)[1]/dim(Banks)[1])*mean(BSa>0)


uni.cov<-matrix(0,3,3)
uni.cov[1:2,1:2]<-fitsuma$cov
uni.cov[3,3]<-(phi*(1-phi))/dim(Banks)[1]

multi.cov<-matrix(0,6,6)
multi.cov[1:5,1:5]<-solve(fit1.5.1$hess)[2:6,2:6]
multi.cov[6,6]<-(phi*(1-phi))/dim(Banks)[1]


#############################################################################################
# plot VaR(p) vs -logp with CIs
#############################################################################################

pseq<-seq(0.15,0.001,len=100)

VaRMV<-VaR2(p=pseq,sig=sum(a*fit1.5$mle[2:5]),gam=fit1.5$mle[6],phi=phi,th=sum(a*u))
SEVaRMV<-sapply(pseq,SE.VaR.MV,sig=sum(a*fit1.5$mle[2:5]),gam=fit1.5$mle[6],phi=phi,a=a,cov=multi.cov)

VaRUV<-VaR2(p=pseq,sig=fitsuma$mle[1],gam=fitsuma$mle[2],phi=phi,th=sum(a*u))
SEVaRUV<-sapply(pseq,SE.VaR,sig=fitsuma$mle[1],gam=fitsuma$mle[2],phi=phi,cov=uni.cov)


pdf("VaRMV4321.pdf")
par(mar=c(5,5,4,2))
plot(-log(pseq),VaRMV,typ="l",axes=F,xlab="p",ylab="VaR(p)",ylim=c(0,70),cex.lab=2)
axis(1,at=-log(c(0.15,0.1,0.01,0.001)),labels=(c(0.15,0.1,0.01,0.001)),cex.axis=2)
axis(2,at=c(0,10,20,30,40,50,60),cex.axis=2)
lines( -log(pseq),VaRMV-1.96*SEVaRMV,col=2,lty=2)
lines( -log(pseq),VaRMV+1.96*SEVaRMV,col=2,lty=2)
points(-log(phi*(length(BSa[BSa>0])+1-rank(BSa[BSa>0]))/(length(BSa[BSa>0])+1)), (BSa[BSa>0])+sum(a*u))
title("a=(40,30,20,10)",cex.main=2)
dev.off()


pdf("VaRUV4321.pdf")
par(mar=c(5,5,4,2))
plot(-log(pseq),VaRUV,typ="l",axes=F,xlab="p",ylab="VaR(p)",ylim=c(0,70),cex.lab=2)
axis(1,at=-log(c(0.15,0.1,0.01,0.001)),labels=(c(0.15,0.1,0.01,0.001)),cex.axis=2)
axis(2,at=c(0,10,20,30,40,50,60),cex.axis=2)
lines( -log(pseq),VaRUV-1.96*SEVaRUV,col=2,lty=2)
lines( -log(pseq),VaRUV+1.96*SEVaRUV,col=2,lty=2)
points(-log(phi*(length(BSa[BSa>0])+1-rank(BSa[BSa>0]))/(length(BSa[BSa>0])+1)), (BSa[BSa>0])+sum(a*u))
title("a=(40,30,20,10)",cex.main=2)
dev.off()


#############################################################################################
# Expected Shortfall Plots
#############################################################################################


ESMV<-ES2(p=pseq,sig=sum(a*fit1.5$mle[2:5]),gam=fit1.5$mle[6],phi=phi,th=sum(a*u))
SEESMV<-sapply(pseq,SE.ES.MV,sig=sum(a*fit1.5$mle[2:5]),gam=fit1.5$mle[6],phi=phi,a=a,cov=multi.cov)

ESUV<-ES2(p=pseq,sig=fitsuma$mle[1],gam=fitsuma$mle[2],phi=phi,th=sum(a*u))
SEESUV<-sapply(pseq,SE.ES,sig=fitsuma$mle[1],gam=fitsuma$mle[2],phi=phi,cov=uni.cov)


Emp.ES<-NULL
BSa0<-BSa[BSa>0]
pseq2<-(1:length(BSa0))/(length(BSa0)+1)
Emp.ES<-cumsum(rev(sort(BSa0)))/(1:length(BSa0)) + sum(a*u)

pdf("ESMV4321.pdf")
par(mar=c(5,5,4,2))
plot(-log(pseq),ESMV,typ="l",axes=F,xlab="p",ylab="ES(p)",ylim=c(0,150),cex.lab=2)
axis(1,at=-log(c(0.15,0.1,0.01,0.001)),labels=(c(0.15,0.1,0.01,0.001)),cex.axis=2)
axis(2,at=2.5*c(0,10,20,30,40,50,60),cex.axis=2)
lines( -log(pseq),ESMV-1.96*SEESMV,col=2,lty=2)
lines( -log(pseq),ESMV+1.96*SEESMV,col=2,lty=2)
points(-log(pseq2*phi), Emp.ES)
title("a=(40,30,20,10)",cex.main=2)
dev.off()

pdf("ESUV4321.pdf")
par(mar=c(5,5,4,2))
plot(-log(pseq),ESUV,typ="l",axes=F,xlab="p",ylab="ES(p)",ylim=c(0,150),cex.lab=2)
axis(1,at=-log(c(0.15,0.1,0.01,0.001)),labels=(c(0.15,0.1,0.01,0.001)),cex.axis=2)
axis(2,at=2.5*c(0,10,20,30,40,50,60),cex.axis=2)
lines( -log(pseq),ESUV-1.96*SEESUV,col=2,lty=2)
lines( -log(pseq),ESUV+1.96*SEESUV,col=2,lty=2)
points(-log(pseq2*phi), Emp.ES)
title("a=(40,30,20,10)",cex.main=2)
dev.off()



####################################################################################################
####################################################################################################
# Showing the effect of different portfolios:
# a1=10
# a2, a3 change
# a4=90-(a2+a3) as long as it's positive
####################################################################################################
####################################################################################################

# Use a large sample from fitted distribution to estimate phi

Xnew<-sim.GumbelT.MGPD(100000,d=4,a=fit1.5$mle[1],beta=rep(0,4),sig=fit1.5$mle[2:5],gamma=rep(fit1.5$mle[6],4))

# V: VaR using MV model including model-based phi (from Xnew above)
# V2: VaR using MV model but empirical phi
# V.uv: VaR using UV model

a2<-seq(1,88,by=1)
a3<-seq(1,88,by=1)

phi.tilde<-dim(Banks2)[1]/dim(Banks)[1]

V<-V2<-V.uv<-matrix(0,length(a2),length(a3))

for(i in 1:length(a2))
{
  print(i)
  for(j in 1:length(a3))
  {
    if(a2[i]+a3[j]<90)
    {
      a<-c(10,a2[i],a3[j],90-(a2[i]+a3[j]))
      BSa1<-apply(Banks2,1,ws,a=a)
      phi.a<-mean(BSa1>0) * phi.tilde
      
      fitsuma1<-gpd.fit(BSa1,thresh=0,show=F)
      
      Xnewa<-apply(Xnew,1,ws,a=a)
      phi.sim<-mean(Xnewa>0) * phi.tilde
      
      V[i,j]<-VaR2(p=0.001,sig=sum(a*fit1.5$mle[2:5]),gam=fit1.5$mle[6],phi=phi.sim, th=sum(a*u))   
      V2[i,j]<-VaR2(p=0.001,sig=sum(a*fit1.5$mle[2:5]),gam=fit1.5$mle[6],phi=phi.a, th=sum(a*u))
      V.uv[i,j]<-VaR2(p=0.001,sig=fitsuma1$mle[1],gam=fitsuma1$mle[2],phi=phi.a, th=sum(a*u))
      
    }
    else{V[i,j]<-NA; V2[i,j]<-NA;V.uv[i,j]<-NA}
  }
}


image(a2,a3,V)
image(a2,a3,V2)
image(a2,a3,V.uv)



pdf("Portfolio_MV.pdf",width=5,height=5)
filled.contour(a2,a3,V,color.palette=function(x)rev(heat.colors(x)),xlab=expression(a[L]),ylab=expression(a[R]),cex.lab=1.5,cex.axis=2,nlevels=50)
dev.off()
pdf("Portfolio_MV_empphi.pdf",width=5,height=5)
filled.contour(a2,a3,V2,color.palette=function(x)rev(heat.colors(x)),xlab=expression(a[L]),ylab=expression(a[R]),cex.lab=1.5,nlevels=50)
dev.off()
pdf("Portfolio_UV.pdf",width=5,height=5)
filled.contour(a2,a3,V.uv,color.palette=function(x)rev(heat.colors(x)),xlab=expression(a[L]),ylab=expression(a[R]),cex.lab=1.5,nlevels=50)
dev.off()

#####################################################################################################################
#####################################################################################################################
