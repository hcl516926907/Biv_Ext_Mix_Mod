#running command "R CMD BATCH thres_influ.R logfile &"
#running command "nohup R CMD BATCH thres_influ.R &"
#running command "nohup R CMD BATCH thres_influ.R > output.txt &"  #print output to output.txt
#running command "nohup R CMD BATCH thres_influ.R > output.txt 2>&1"  & #redirect stderr to stdout. & indicates that what follows and precedes is a file descriptor, and not a filename.     

source("KRSW/RevExp_U_Functions.r")
source("KRSW/CommonFunctions.r")
source("KRSW/ModelDiagnosticsNewNames.r")


dir.in <- '/home/pgrad2/2448355h/My_PhD_Project/01_Output/Biv_Ext_Mix_Mod/biv_ext_mix_mod_simdat'

dir.out <- '/home/pgrad2/2448355h/My_PhD_Project/01_Output/Biv_Ext_Mix_Mod/thres_influ'

load(file=file.path(dir.in,'simulation_data.RData'))

#-------------------Compare the censored and uncensored estimators----------------------------
n <- 1000
d<-2
a<-c(1.656,1.656) # NB this is 1/\alpha, where \alpha is the parameterization in Kiriliouk, Rootzen, Segers and Wadsworth. 
beta<-c(0,0)
sig<-c(0.571,0.451)
gamma<-c(0.253,0.035)


t1 <- Sys.time()
mle.cen.mat <- c()
mle.ucen.mat <- c()
sd.cen.mat <- c()
sd.ucen.mat <- c()
X.all <- list()
cnt.issue <- 1
i <- 1
while (i <= n){
  
  sim.dat <-sim.RevExpU.MGPD(n=200, d=d, a=a, beta=beta, sig=sig, gamma=gamma, MGPD = T,std=T)
  X <- sim.dat$X

  # fit without censoring
  result1.1 <- try({
    fit1.1<-fit.MGPD.RevExpU(x=X, u=apply(X,2,min)-0.01, std=F,
                             dep.scale.fix=T, dep.loc.fix=T, marg.shape.ind=1:2, marg.scale.ind=1:2,
                             maxit=5000)
    
    for (j in 1:5){
      if (j!=5){
        fit1.1<-fit.MGPD.RevExpU(x=X, u=apply(X,2,min)-0.01, std=F,
                                 dep.scale.fix=T, dep.loc.fix=T, marg.shape.ind=1:2, marg.scale.ind=1:2,
                                 marg.scale.start=fit1.1$mle[2:3], marg.shape.start=fit1.1$mle[4:5],
                                 dep.start=fit1.1$mle[1], maxit=5000)
      }else{
        fit1.1<-fit.MGPD.RevExpU(x=X, u=apply(X,2,min)-0.01, std=F,
                                 dep.scale.fix=T, dep.loc.fix=T, marg.shape.ind=1:2, marg.scale.ind=1:2,
                                 marg.scale.start=fit1.1$mle[2:3], marg.shape.start=fit1.1$mle[4:5],
                                 dep.start=fit1.1$mle[1], maxit=5000, hessian=TRUE)
      } 
    }
  }, silent=T)
  
  if (inherits(result1.1, 'try-error')){next}
  
  # fit with censoring
  result <- try({
    fit1<-fit.MGPD.RevExpU(x=X, u=rep(0,2), std=F,
                         dep.scale.fix=T, dep.loc.fix=T, marg.shape.ind=1:2, marg.scale.ind=1:2,
                         maxit=5000)
    for (j in 1:5){
      if (j!=5){
        fit1<-fit.MGPD.RevExpU(x=X, u=rep(0,2), std=F,
                               dep.scale.fix=T, dep.loc.fix=T, marg.shape.ind=1:2, marg.scale.ind=1:2,
                               marg.scale.start=fit1$mle[2:3], marg.shape.start=fit1$mle[4:5],
                               dep.start=fit1$mle[1], maxit=5000)
      }else{
        fit1<-fit.MGPD.RevExpU(x=X, u=rep(0,2), std=F,
                               dep.scale.fix=T, dep.loc.fix=T, marg.shape.ind=1:2, marg.scale.ind=1:2,
                               marg.scale.start=fit1$mle[2:3], marg.shape.start=fit1$mle[4:5],
                               dep.start=fit1$mle[1], maxit=5000, hessian=TRUE)
      }
    }
  }, silent=T)

  
  if (inherits(result,'try-error')) {next}
  
  sd1 <- sqrt(diag(solve(fit1$hess)))
  sd1.1 <- sqrt(diag(solve(fit1.1$hess)))
  
  if (!((anyNA(sd1))|(anyNA(sd1.1)))){
    mle.cen.mat <- rbind(mle.cen.mat, fit1$mle)
    sd.cen.mat <- rbind(sd.cen.mat, sd1)
    
    mle.ucen.mat <- rbind(mle.ucen.mat, fit1.1$mle)
    sd.ucen.mat <- rbind(sd.ucen.mat, sd1.1)
    
    X.all[[i]] <- X
    if (! i %% 25){
      print(paste('Done',i,'simulations.'))
      t2 <- Sys.time()
      print(t2-t1)
    }
    i <- i+1
  }

}

t2 <- Sys.time()

print(t2-t1)

save(mle.cen.mat, sd.cen.mat, mle.ucen.mat, sd.ucen.mat, file=file.path(dir.out,'mle_sd.RData'))

load(file.path(dir.out,'mle_sd.RData'))
para.true <- c(a[1], sig, gamma)

# calculate the bias
bias.cen <- colMeans(mle.cen.mat) - para.true
bias.ucen <- colMeans(mle.ucen.mat) - para.true

# ??? calculate the average of sd
# ??? var of mle
sd.cen <- colMeans(sd.cen.mat, na.rm=TRUE)
sd.ucen <- colMeans(sd.ucen.mat, na.rm=TRUE)

sd.cen.1 <- apply(mle.cen.mat, 2, sd)
sd.ucen.1 <- apply(mle.ucen.mat, 2, sd)

# calculate the rmse--root mean squared error
rmse.cen <- sqrt(bias.cen^2 + sd.cen^2)
rmse.ucen <- sqrt(bias.ucen^2 + sd.ucen^2)

rmse.cen.1 <- sqrt(bias.cen^2 + sd.cen.1^2)
rmse.ucen.1 <- sqrt(bias.ucen^2 + sd.ucen.1^2)

plot(density(mle.cen.mat[,1]))
abline(v=a[1],col='red')

plot(density(mle.ucen.mat[,1]))
abline(v=a[1],col='red')

plot(density(mle.cen.mat[,2]))
abline(v=sig[1],col='red')

plot(density(mle.ucen.mat[,2]))
abline(v=sig[1],col='red')

plot(density(mle.cen.mat[,3]))
abline(v=sig[2],col='red')

plot(density(mle.ucen.mat[,3]))
abline(v=sig[2],col='red')

plot(density(mle.cen.mat[,4]))
abline(v=gamma[1],col='red')

plot(density(mle.ucen.mat[,4]))
abline(v=gamma[1],col='red')

plot(density(mle.cen.mat[,5]))
abline(v=gamma[2],col='red')

plot(density(mle.ucen.mat[,5]))
abline(v=gamma[2],col='red')

#calculate the rate of confidence intervals covering the true parameters
#wrong
upp.cen <- mle.cen.mat + 1.96*sd.cen.mat
low.cen <- mle.cen.mat - 1.96*sd.cen.mat

#upp.cen <- mle.cen.mat + 1.96*sd.cen
#low.cen <- mle.cen.mat - 1.96*sd.cen

upp.cen.1 <- t(t(mle.cen.mat) + 1.96*sd.cen.1)
low.cen.1 <- t(t(mle.cen.mat) - 1.96*sd.cen.1)

upp.ucen <- mle.ucen.mat + 1.96*sd.ucen.mat
low.ucen <- mle.ucen.mat - 1.96*sd.ucen.mat

upp.ucen.1 <- t(t(mle.ucen.mat) + 1.96*sd.ucen.1)
low.ucen.1 <- t(t(mle.ucen.mat) - 1.96*sd.ucen.1)

cover.cen.ind <- t(t(upp.cen) >= para.true) & t(t(low.cen) <= para.true)
cover.ucen.ind <- t(t(upp.ucen) >= para.true) & t(t(low.ucen) <= para.true)

cover.cen.ind.1 <- t(t(upp.cen.1) >= para.true) & t(t(low.cen.1) <= para.true)
cover.ucen.ind.1 <- t(t(upp.ucen.1) >= para.true) & t(t(low.ucen.1) <= para.true)

colMeans(cover.cen.ind, na.rm=TRUE)
colMeans(cover.ucen.ind, na.rm=TRUE)

colMeans(cover.cen.ind.1, na.rm=TRUE)
colMeans(cover.ucen.ind.1, na.rm=TRUE)
#--------------------------------check threshold for each margin---------------------------------

u.x <- u.x.p10
u.z <- u.z.p10

X <- X.p10
Z <- Z.p10


q1.x <- sum((X[,1]> u.x[1]))/length(X[,1]) 
q2.x <- sum((X[,2]> u.x[2]))/length(X[,2]) 

q1.z <- sum((Z[,1]> u.z[1]))/length(Z[,1]) 
q2.z <- sum((Z[,2]> u.z[2]))/length(Z[,2]) 
library(ismev)
q <- 0.95
dat <- Z[,2]
u <- quantile(dat,q)
plot(dat)
gpd.fitrange(dat, umin=8,umax=12,nint=20)
m1<-gpd.fit(dat,thresh=u)
gpd.diag(m1)

