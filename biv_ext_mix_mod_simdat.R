# sample code
source("KRSW/RevExp_U_Functions.r")
source("KRSW/RevExp_T_Functions.r")
source("KRSW/MVGaussian_T_Functions.r")
source("KRSW/CommonFunctions.r")
source("KRSW/Gumbel_U_Functions.r")
source("KRSW/Gumbel_T_Functions.r")
source("KRSW/ModelDiagnosticsNewNames.r")
library(evd)
library(tmvtnorm)

dir.out <- '/home/pgrad2/2448355h/My_PhD_Project/01_Output/Biv_Ext_Mix_Mod/biv_ext_mix_mod_simdat'


############################################################################################
#------------------------sample size ratio bulk:tail = 4:1----------------------------------
############################################################################################

set.seed(1234)

d<-2
a<-c(1.656,1.656) # NB this is 1/\alpha, where \alpha is the parameterization in Kiriliouk, Rootzen, Segers and Wadsworth. 
beta<-c(0,0)
sig<-c(0.571,0.451)
gamma<-c(-0.253,-0.035)

# Final beta parameter fixed at zero in estimation, so needs to be relative to this

# Simulate with (conditionally) GP(sig,gamma) and (conditionally) exponential margins 
X.tail<-sim.RevExpU.MGPD(n=500,d=d, a=a, beta=beta, sig=sig, gamma=gamma, MGPD = T,std=T)

#  GP(sig,gamma)
plot(X.tail$X, pch=20, xlab=expression(X[1]),ylab=expression(X[2]),xlim=c(-10,2))

# Exponential
plot(X.tail$Z, pch=20, xlab=expression(Z[1]),ylab=expression(Z[2]))

# difference is that, chiplot in the evd package plots the Chi(u) plot in Cole's paper, 
# while chiPlot in ModelDiagnosticsNewNames.r directly plot the Chi across u, 
# where Chi is the limit of Chi(u) as u approaches 1
chiplot(X.tail$X, xlim=c(0.5,1), qlim=c(0.5,0.99), ask=F)
chiplot(X.tail$Z, xlim=c(0.5,1), qlim=c(0.5,0.99), ask=F)
#chiPlot(data=X.tail$X, ylabel=expression(chi[HLRB]~(q)), chimod=NULL, nsim=1000, nq = 35, qmin = 0.5, qmax = 0.99)


#--------------simulate bulk data and combine them with the tail data---------------
u.x <- c(-min(X.tail$X[,1]),-min(X.tail$X[,2]))
u.z <- c(-min(X.tail$Z[,1]),-min(X.tail$Z[,2]))

rho=0.65
sigma <- matrix(c(1,rho,rho,0.8),ncol=2)

# GP scale tail data combined with the bulk data
set.seed(1234)
X.bulk <- rtmvnorm(2000, mean=u.x-c(1,1), sigma=1*sigma, lower=c(0,0), upper=u.x)
plot(X.bulk)
X <- rbind(X.bulk, sweep(X.tail$X,2,u.x,"+"))
plot(X,xlab='X1', ylab='X2', main='Simulated tNormal + mGP Data')
plot(X,pch=20, xlim=c(16,23),ylim=c(0,10),xlab='X1', ylab='X2', main='Simulated tNormal + mGP Data')
plot(X,pch=20, xlim=c(16,23),ylim=c(0,10),
     xlab='X1', ylab='X2', main='Simulated tNormal + mGP Data',col=c(rep('grey50',2000),rep('black',500)))
segments(u.x[1], 0, u.x[1], u.x[2], col='red', lty=2)
segments(0, u.x[2], u.x[1], u.x[2], col='red', lty=2)
legend('bottomright',   legend = c("Tail Data", "Bulk Data", "Threshold"), lty=c(NA,NA,2),col = c('black', 'grey50', 'red'),         
       pch=c(20,20,NA)) 
chiplot(X, qlim=c(0.02,0.99), ask=F)

# Exponential scale tail data combined with the bulk data
set.seed(1234)
Z.bulk <- rtmvnorm(2000, mean=u.z-c(2,2), sigma=6*sigma, lower=c(0,0), upper=u.z)
Z <- rbind(Z.bulk, sweep(X.tail$Z,2,u.z,"+"))
plot(Z,pch=20,xlab='Z1', ylab='Z2', main='Simulated tNormal + std GP Data')
plot(Z,pch=20, xlab='Z1', ylab='Z2', main='Simulated tNormal + std GP Data',col=c(rep('grey50',2000),rep('black',500)))
segments(u.z[1], 0, u.z[1], u.z[2], col='red', lty=2)
segments(0, u.z[2], u.z[1], u.z[2], col='red', lty=2)
legend('bottomright',   legend = c("Tail Data", "Bulk Data", "Threshold"), lty=c(NA,NA,2),col = c('black', 'grey50', 'red'),         
       pch=c(20,20,NA)) 
chiplot(Z, qlim=c(0.02,0.99), ask=F)


X.4to1 <- X
Z.4to1 <- Z
u.x.4to1 <- u.x
u.z.4to1 <- u.z


#############################################################################################
#------------------------sample size ratio bulk:tail = 10:1----------------------------------
#############################################################################################
set.seed(1234)

d<-2
a<-c(1.656,1.656) # NB this is 1/\alpha, where \alpha is the parameterization in Kiriliouk, Rootzen, Segers and Wadsworth. 
beta<-c(0,0)
sig<-c(0.571,0.451)
gamma<-c(-0.253,-0.035)

# Final beta parameter fixed at zero in estimation, so needs to be relative to this

# Simulate with (conditionally) GP(sig,gamma) and (conditionally) exponential margins 
X.tail<-sim.RevExpU.MGPD(n=200,d=d, a=a, beta=beta, sig=sig, gamma=gamma, MGPD = T,std=T)

#  GP(sig,gamma)
plot(X.tail$X, pch=20, xlab=expression(X[1]),ylab=expression(X[2]),xlim=c(-10,2))

# Exponential
plot(X.tail$Z, pch=20, xlab=expression(Z[1]),ylab=expression(Z[2]))


#--------------simulate bulk data and combine them with the tail data---------------
u.x <- c(-min(X.tail$X[,1]),-min(X.tail$X[,2]))
u.z <- c(-min(X.tail$Z[,1]),-min(X.tail$Z[,2]))

rho=0.65
sigma <- matrix(c(1,rho,rho,0.8),ncol=2)

# GP scale tail data combined with the bulk data
set.seed(1234)
X.bulk <- rtmvnorm(2000, mean=u.x-c(1,1), sigma=1*sigma, lower=c(0,0), upper=u.x)
plot(X.bulk)
X <- rbind(X.bulk, sweep(X.tail$X,2,u.x,"+"))
plot(X,xlab='X1', ylab='X2', main='Simulated tNormal + mGP Data')
plot(X,pch=20, xlim=c(16,23),ylim=c(0,10),xlab='X1', ylab='X2', main='Simulated tNormal + mGP Data')
plot(X,pch=20, xlim=c(16,23),ylim=c(0,10),
     xlab='X1', ylab='X2', main='Simulated tNormal + mGP Data',col=c(rep('grey50',2000),rep('black',200)))
segments(u.x[1], 0, u.x[1], u.x[2], col='red', lty=2)
segments(0, u.x[2], u.x[1], u.x[2], col='red', lty=2)
legend('bottomright',   legend = c("Tail Data", "Bulk Data", "Threshold"), lty=c(NA,NA,2),col = c('black', 'grey50', 'red'),         
       pch=c(20,20,NA)) 
chiplot(X, qlim=c(0.02,0.99), ask=F)

# Exponential scale tail data combined with the bulk data
set.seed(1234)
Z.bulk <- rtmvnorm(2000, mean=u.z-c(3,3), sigma=4*sigma, lower=c(0,0), upper=u.z)
Z <- rbind(Z.bulk, sweep(X.tail$Z,2,u.z,"+"))
plot(Z,pch=20,xlab='Z1', ylab='Z2', main='Simulated tNormal + std GP Data')
plot(Z,pch=20, xlab='Z1', ylab='Z2', main='Simulated tNormal + std GP Data',col=c(rep('grey50',2000),rep('black',200)))
segments(u.z[1], 0, u.z[1], u.z[2], col='red', lty=2)
segments(0, u.z[2], u.z[1], u.z[2], col='red', lty=2)
legend('bottomright',   legend = c("Tail Data", "Bulk Data", "Threshold"), lty=c(NA,NA,2),col = c('black', 'grey50', 'red'),         
       pch=c(20,20,NA)) 
chiplot(Z, qlim=c(0.02,0.99), ask=F)


X.10to1 <- X
Z.10to1 <- Z
u.x.10to1 <- u.x
u.z.10to1 <- u.z

#############################################################################################
#------------------------sample size ratio bulk:tail = 20:1----------------------------------
#############################################################################################
set.seed(1234)

d<-2
a<-c(1.656,1.656) # NB this is 1/\alpha, where \alpha is the parameterization in Kiriliouk, Rootzen, Segers and Wadsworth. 
beta<-c(0,0)
sig<-c(0.571,0.451)
gamma<-c(-0.253,-0.035)

# Final beta parameter fixed at zero in estimation, so needs to be relative to this

# Simulate with (conditionally) GP(sig,gamma) and (conditionally) exponential margins 
X.tail<-sim.RevExpU.MGPD(n=100,d=d, a=a, beta=beta, sig=sig, gamma=gamma, MGPD = T,std=T)

#  GP(sig,gamma)
plot(X.tail$X, pch=20, xlab=expression(X[1]),ylab=expression(X[2]),xlim=c(-10,2))

# Exponential
plot(X.tail$Z, pch=20, xlab=expression(Z[1]),ylab=expression(Z[2]))


#--------------simulate bulk data and combine them with the tail data---------------
u.x <- c(-min(X.tail$X[,1]),-min(X.tail$X[,2]))
u.z <- c(-min(X.tail$Z[,1]),-min(X.tail$Z[,2]))

rho=0.65
sigma <- matrix(c(1,rho,rho,0.8),ncol=2)

# GP scale tail data combined with the bulk data
set.seed(1234)
X.bulk <- rtmvnorm(2000, mean=u.x-c(1.3,1.1), sigma=0.6*sigma, lower=c(0,0), upper=u.x)
plot(X.bulk)
X <- rbind(X.bulk, sweep(X.tail$X,2,u.x,"+"))
plot(X,xlab='X1', ylab='X2', main='Simulated tNormal + mGP Data')
plot(X,pch=20, xlim=c(0,7),ylim=c(0,5),xlab='X1', ylab='X2', main='Simulated tNormal + mGP Data')
plot(X,pch=20, xlim=c(0,7),ylim=c(0,5),
     xlab='X1', ylab='X2', main='Simulated tNormal + mGP Data',col=c(rep('grey50',2000),rep('black',100)))
segments(u.x[1], 0, u.x[1], u.x[2], col='red', lty=2)
segments(0, u.x[2], u.x[1], u.x[2], col='red', lty=2)
legend('bottomright',   legend = c("Tail Data", "Bulk Data", "Threshold"), lty=c(NA,NA,2),col = c('black', 'grey50', 'red'),         
       pch=c(20,20,NA)) 
chiplot(X, qlim=c(0.02,0.99), ask=F)

# Exponential scale tail data combined with the bulk data
set.seed(1234)
Z.bulk <- rtmvnorm(2000, mean=u.z-c(2,2), sigma=1.5*matrix(c(1,0.65,0.65,0.8),ncol=2), lower=c(0,0), upper=u.z)
Z <- rbind(Z.bulk, sweep(X.tail$Z,2,u.z,"+"))
plot(Z,pch=20,xlab='Z1', ylab='Z2', main='Simulated tNormal + std GP Data')
plot(Z,pch=20, xlab='Z1', ylab='Z2', main='Simulated tNormal + std GP Data',col=c(rep('grey50',2000),rep('black',500)))
segments(u.z[1], 0, u.z[1], u.z[2], col='red', lty=2)
segments(0, u.z[2], u.z[1], u.z[2], col='red', lty=2)
legend('bottomright',   legend = c("Tail Data", "Bulk Data", "Threshold"), lty=c(NA,NA,2),col = c('black', 'grey50', 'red'),         
       pch=c(20,20,NA)) 
chiplot(Z, qlim=c(0.02,0.99), ask=F)

X.20to1 <- X
Z.20to1 <- Z
u.x.20to1 <- u.x
u.z.20to1 <- u.z

save(X.4to1,Z.4to1,u.x.4to1,u.z.4to1, 
     X.10to1,Z.10to1,u.x.10to1,u.z.10to1,
     X.20to1,Z.20to1,u.x.20to1,u.z.20to1,
     file=file.path(dir.out,'simulation_data.RData'))
