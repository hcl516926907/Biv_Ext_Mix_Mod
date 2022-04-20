source("func_tools.R")

# data directory
dir.dat <- '/home/pgrad2/2448355h/My_PhD_Project/00_Dataset'

# name of Ruerto Rico rive dataset 
riv.dat.fit.nam <- 'Puerto Rico Fit Dataset.dat'
riv.dat.test.nam <- 'Puerto Rico Test Data.dat'

# name of Leeds air pollution dataset
pol.dat.fit.nam <- 'Leeds Fit Dataset.dat'
pol.dat.test.nam <- 'Leeds Test Dataset.dat'

# load all datasets
riv.dat.fit <- read.table(file.path(dir.dat,riv.dat.fit.nam), header=TRUE)
riv.dat.test <- read.table(file.path(dir.dat,riv.dat.test.nam), header=TRUE)
riv.dat <- rbind(riv.dat.fit,riv.dat.test)


#------------------------------------------------------------------
# sample code
source("KRSW/RevExp_U_Functions.r")
source("KRSW/RevExp_T_Functions.r")
source("KRSW/MVGaussian_T_Functions.r")
source("KRSW/CommonFunctions.r")
source("KRSW/Gumbel_U_Functions.r")
source("KRSW/Gumbel_T_Functions.r")

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
fitS<-fit.MGPD.RevExpU(x=X$Z,u=c(0,0),std=T, hessian=TRUE)
# Parameter order: (a1,a2,beta1)
#$mle
#[1] 1.9095402 3.8381261 0.5267898

# Uncensored estimation
u1<-min(X$Z)-0.01
fitS2<-fit.MGPD.RevExpU(x=X$Z, u=c(u1,u1), std=T, hessian=TRUE)
# Parameter order: (a1,a2,beta1)
#$mle
#[1] 2.1742281 3.8906472 0.4348918

# Fit GP form; censored likelihood censoring at (0,0)
# marg.scale.ind = 1:2 states two scale parameters should be fitted (for common scale set marg.scale.ind = c(1,1))
# marg.shape.ind = 1:2 states two shape parameters should be fitted ( " " )

fitGP<-fit.MGPD.RevExpU(x=X$X, u=c(0,0), std=F, marg.scale.ind = c(1,1), marg.shape.ind = 1:2, maxit=5000, hessian=TRUE)
# Parameter order: (a1,a2,beta1,sig,gamma1,gamma2)
#$mle
#[1]  1.89320730  3.64727336  0.51113618  2.09892248  0.04513193 -0.14167273

#------------------------------------------------------------------

# Transform margins to uniform scale by empirical probability integral transformation
x1.uni <- my.ecdf(riv.dat$V1)
hist(x1.uni)

x2.uni <- my.ecdf(riv.dat$V2)
hist(x2.uni)

xpb <- cbind(1/(1-x1.uni), 1/(1-x2.uni))
q <- 0.95
cond <- xpb[,1] > quantile(xpb[,1], q) | xpb[,2] > quantile(xpb[,2], q)
xpb2 <- cbind(xpb[cond, 1]/quantile(xpb[,1], q), xpb[cond, 2]/quantile(xpb[,2], q))
plot(log(xpb2))
dim(xpb2)

# exponential margins
xeb2<-log(xpb2)

# Same but on MGPD scale (i.e. scale of the observations)
X.obs <- cbind(riv.dat[cond,'V1'] - quantile(riv.dat[,'V1'], q),
               riv.dat[cond,'V2'] - quantile(riv.dat[,'V2'], q))
plot(X.obs)
dim(X.obs)

# Same but on MGPD scale (i.e. scale of the log observations)
X.obs.lg <- cbind(log(riv.dat[cond,'V1']) - quantile(log(riv.dat[,'V1']), q),
               log(riv.dat[cond,'V2']) - quantile(log(riv.dat[,'V2']), q))
plot(X.obs.lg)
dim(X.obs.lg)

# Free scale parameter, Free location parameter
fit1<-fit.MGPD.RevExpU(x=xeb2, u=rep(0,2), std=T, dep.scale.fix=F)
fit1
fit1<-fit.MGPD.RevExpU(x=xeb2, u=rep(0,2), std=T, dep.scale.fix=F,dep.start = fit1$mle)
fit1

# Constrained scale parameter, free location parameters
fit1.1<-fit.MGPD.RevExpU(x=xeb2, u=rep(0,2), std=T, dep.scale.fix=T)
fit1.1
fit1.1<-fit.MGPD.RevExpU(x=xeb2, u=rep(0,2), std=T, dep.scale.fix=T,dep.start = fit1.1$mle)
fit1.1

# Free scale parameter, constrained location parameters
fit1.2<-fit.MGPD.RevExpU(x=xeb2, u=rep(0,2), std=T, dep.scale.fix=F, dep.loc.fix=T)
fit1.2
fit1.2<-fit.MGPD.RevExpU(x=xeb2, u=rep(0,2), std=T, dep.scale.fix=F, dep.loc.fix=T, dep.start = fit1.2$mle)
fit1.2

# Single scale parameter, fixed location
fit1.3<-fit.MGPD.RevExpU(x=xeb2, u=rep(0,2), std=T, dep.scale.fix=T, dep.loc.fix=T)
fit1.3
fit1.3<-fit.MGPD.RevExpU(x=xeb2, u=rep(0,2), std=T, dep.scale.fix=T, dep.loc.fix=T, dep.start = fit1.3$mle)
fit1.3

#Wald test to simplify parameters in the dependence model
1 - pchisq(2*(fit1.1$nll - fit1$nll), df=1)
1 - pchisq(2*(fit1.2$nll - fit1$nll), df=1)
1 - pchisq(2*(fit1.3$nll - fit1.1$nll), df=1)
1 - pchisq(2*(fit1.3$nll - fit1.2$nll), df=1)

#################################################################################################
# Marginal parameter stability plots to understand if 0.83 quantile might be OK
#################################################################################################

library(ismev)

plot(log(riv.dat[,'V1']))
quantile(log(riv.dat[,'V1']),q)
gpd.fitrange(log(riv.dat[,'V1']),umin=5.5,umax=7.5,nint=20)
gpd.fit(log(riv.dat[,'V1']),thresh=6.3)
m1<-gpd.fit(log(riv.dat[,'V1']),thresh=quantile(log(riv.dat[,'V1']),q))
gpd.diag(m1)

plot(log(riv.dat[,'V2']))
quantile(log(riv.dat[,'V2']),q)
gpd.fitrange(log(riv.dat[,'V2']),umin=4,umax=5.5,nint=20)
gpd.fit(log(riv.dat[,'V2']),thresh=6.4)
m2<-gpd.fit(log(riv.dat[,'V2']),thresh=quantile(log(riv.dat[,'V2']),q))
gpd.diag(m2)


#################################################################################################
# Fit margins and dependence simulataneously
#################################################################################################

fit1.4<-fit.MGPD.RevExpU(x=X.obs.lg, u=rep(0,2), std=F,
                         dep.scale.fix=T, dep.loc.fix=T, marg.shape.ind=1:2, marg.scale.ind=1:2,
                         dep.start=fit1.3$mle[1],maxit=5000)

fit1.4<-fit.MGPD.RevExpU(x=X.obs.lg, u=rep(0,2), std=F,
                         dep.scale.fix=T, dep.loc.fix=T, marg.shape.ind=1:2, marg.scale.ind=1:2,
                         marg.scale.start=fit1.4$mle[2:3], marg.shape.start=fit1.4$mle[4:5],
                         dep.start=fit1.4$mle[1],maxit=5000)

# Test for common shape
fit1.5<-fit.MGPD.RevExpU(x=X.obs.lg, u=rep(0,2), std=F,
                            dep.scale.fix=T, dep.loc.fix=T, marg.shape.ind=rep(1,2), marg.scale.ind=1:2,
                            dep.start=fit1.4$mle[1],maxit=5000)

fit1.5<-fit.MGPD.RevExpU(x=X.obs.lg, u=rep(0,2), std=F,
                            dep.scale.fix=T, dep.loc.fix=T, marg.shape.ind=rep(1,2), marg.scale.ind=1:2,
                            marg.scale.start=fit1.5$mle[2:3], marg.shape.start=fit1.5$mle[4],
                            dep.start=fit1.5$mle[1],maxit=5000)

1-pchisq(2*(fit1.5$nll-fit1.4$nll),df=1)


# Final model: fit1.4 (Single dependence parameter, different marginal scale and shape)

# With Hessian, to get standard errors
fit1.4.1<-fit.MGPD.RevExpU(x=X.obs.lg, u=rep(0,2), std=F,
                                   dep.scale.fix=T, dep.loc.fix=T, marg.shape.ind=1:2, marg.scale.ind=1:2,
                                   marg.scale.start=fit1.4$mle[2:3], marg.shape.start=fit1.4$mle[4:5],
                                   dep.start=fit1.4$mle[1],maxit=5000,hessian=TRUE)

solve(fit1.4.1$hess)
sqrt(diag(solve(fit1.4.1$hess)))

# Compare marginal and jointly estimated standard errors

sqrt(diag(solve(fit1.4.1$hess)))[2]/m1$se[1]
sqrt(diag(solve(fit1.4.1$hess)))[3]/m2$se[1]

###############################################################################################
# Diagnostics:
###############################################################################################

# Will print fitted 4-dim chi (caluculated using Monte Carlo)
# Rounded value is 0.4 (checked also with nsim = 100000)


# MGPD.diag.RevExpU(x=log(riv.dat),a = fit1.4$mle[1],beta = rep(0,2),marg.scale = fit1.4$mle[2:3],marg.shape = fit1.4$mle[4:5],nsim = 10000,
#                   chi=T,chiBS=T,nbs=2000,chiylabel = expression(hat(chi)[HLRB]~(q)))
pdf("/home/pgrad2/2448355h/My_PhD_Project/01_Output/Biv_Ext_Mix_Mod/ChiRF.pdf")
MGPD.diag.RevExpU(x=log(riv.dat),a = fit1.4$mle[1],beta = rep(0,2),marg.scale = fit1.4$mle[2:3],marg.shape = fit1.4$mle[2:3],nsim = 10000,
                  chi=T, cols=c(1,2),chiBS=T,nbs=2000,chiylabel = expression(hat(chi)[V1V2]~(q)))
abline(v=q,col=4)
dev.off()


###############################################################################################
# Marginal GPD QQ plots
###############################################################################################
GPD.diag.th(X.obs.lg,sig=fit1.4$mle[2:3],gam=fit1.4$mle[4:5])





