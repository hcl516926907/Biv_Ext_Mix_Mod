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
#------------------------tail prop: 20%-------------------------------
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
X.bulk <- rtmvnorm(2500, mean=u.x-c(1,1), sigma=1*sigma, lower=c(0,0))
cond.x <- t(t(X.bulk) < u.x)
X.bulk2 <- X.bulk[rowSums(cond.x)==2,]
plot(X.bulk2)
print(dim(X.bulk2)[1])
X <- rbind(X.bulk2, sweep(X.tail$X,2,u.x,"+"))
print(dim(X)[1])
plot(X,xlab='X1', ylab='X2', main='Simulated tNormal + mGP Data')
plot(X,pch=20, xlim=c(16,23),ylim=c(0,10),xlab='X1', ylab='X2', main='Simulated tNormal + mGP Data')
plot(X,pch=20, xlim=c(16,23),ylim=c(0,10),
     xlab='X1', ylab='X2', main='Simulated tNormal + mGP Data',col=c(rep('grey50',dim(X.bulk2)[1]),rep('black',500)))
segments(u.x[1], 0, u.x[1], u.x[2], col='red', lty=2)
segments(0, u.x[2], u.x[1], u.x[2], col='red', lty=2)
legend('bottomright',   legend = c("Tail Data", "Bulk Data", "Threshold"), lty=c(NA,NA,2),col = c('black', 'grey50', 'red'),         
       pch=c(20,20,NA)) 
chiplot(X, qlim=c(0.02,0.99), ask=F)

# Exponential scale tail data combined with the bulk data
set.seed(1234)
Z.bulk <- rtmvnorm(2700, mean=u.z-c(2,2), sigma=6*sigma, lower=c(0,0))
cond.z <- t(t(Z.bulk) < u.z)
Z.bulk2 <- Z.bulk[rowSums(cond.z)==2,]
plot(Z.bulk2)
print(dim(Z.bulk2)[1])
Z <- rbind(Z.bulk2, sweep(X.tail$Z,2,u.z,"+"))
print(dim(Z)[1])
plot(Z,pch=20,xlab='Z1', ylab='Z2', main='Simulated tNormal + std GP Data')
plot(Z,pch=20, xlab='Z1', ylab='Z2', main='Simulated tNormal + std GP Data',col=c(rep('grey50',dim(Z.bulk2)[1]),rep('black',500)))
segments(u.z[1], 0, u.z[1], u.z[2], col='red', lty=2)
segments(0, u.z[2], u.z[1], u.z[2], col='red', lty=2)
legend('bottomright',   legend = c("Tail Data", "Bulk Data", "Threshold"), lty=c(NA,NA,2),col = c('black', 'grey50', 'red'),         
       pch=c(20,20,NA)) 
chiplot(Z, qlim=c(0.02,0.99), ask=F)

print(500/dim(X)[1])
print(500/dim(Z)[1])
X.p20 <- X
Z.p20 <- Z
u.x.p20 <- u.x
u.z.p20 <- u.z


#############################################################################################
#------------------------tail prop: 10%----------------------------------
#############################################################################################
set.seed(1234)

d<-2
a<-c(1.656,1.656) # NB this is 1/\alpha, where \alpha is the parameterization in Kiriliouk, Rootzen, Segers and Wadsworth. 
beta<-c(0,0)
sig<-c(0.571,0.451)
gamma<-c(-0.253,-0.035)

# Final beta parameter fixed at zero in estimation, so needs to be relative to this

# Simulate with (conditionally) GP(sig,gamma) and (conditionally) exponential margins 
X.tail<-sim.RevExpU.MGPD(n=250,d=d, a=a, beta=beta, sig=sig, gamma=gamma, MGPD = T,std=T)

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
X.bulk <- rtmvnorm(2900, mean=u.x-c(1,1), sigma=1*sigma, lower=c(0,0))
cond.x <- t(t(X.bulk) < u.x)
X.bulk2 <- X.bulk[rowSums(cond.x)==2,]
plot(X.bulk2)
print(dim(X.bulk2)[1])
X <- rbind(X.bulk2, sweep(X.tail$X,2,u.x,"+"))
print(dim(X)[1])
plot(X,xlab='X1', ylab='X2', main='Simulated tNormal + mGP Data')
plot(X,pch=20, xlim=c(16,23),ylim=c(0,10),xlab='X1', ylab='X2', main='Simulated tNormal + mGP Data')
plot(X,pch=20, xlim=c(16,23),ylim=c(0,10),
     xlab='X1', ylab='X2', main='Simulated tNormal + mGP Data',col=c(rep('grey50',dim(X.bulk2)[1]),rep('black',250)))
segments(u.x[1], 0, u.x[1], u.x[2], col='red', lty=2)
segments(0, u.x[2], u.x[1], u.x[2], col='red', lty=2)
legend('bottomright',   legend = c("Tail Data", "Bulk Data", "Threshold"), lty=c(NA,NA,2),col = c('black', 'grey50', 'red'),         
       pch=c(20,20,NA)) 
chiplot(X, qlim=c(0.02,0.99), ask=F)

# Exponential scale tail data combined with the bulk data
set.seed(1234)
Z.bulk <- rtmvnorm(2500, mean=u.z-c(3,3), sigma=4*sigma, lower=c(0,0))
cond.z <- t(t(Z.bulk) < u.z)
Z.bulk2 <- Z.bulk[rowSums(cond.z)==2,]
plot(Z.bulk2)
print(dim(Z.bulk2)[1])
Z <- rbind(Z.bulk2, sweep(X.tail$Z,2,u.z,"+"))
print(dim(Z)[1])
plot(Z,pch=20,xlab='Z1', ylab='Z2', main='Simulated tNormal + std GP Data')
plot(Z,pch=20, xlab='Z1', ylab='Z2', main='Simulated tNormal + std GP Data',col=c(rep('grey50',dim(Z.bulk2)[1]),rep('black',250)))
segments(u.z[1], 0, u.z[1], u.z[2], col='red', lty=2)
segments(0, u.z[2], u.z[1], u.z[2], col='red', lty=2)
legend('bottomright',   legend = c("Tail Data", "Bulk Data", "Threshold"), lty=c(NA,NA,2),col = c('black', 'grey50', 'red'),         
       pch=c(20,20,NA)) 
chiplot(Z, qlim=c(0.02,0.99), ask=F)


print(250/dim(X)[1])
print(250/dim(Z)[1])
X.p10 <- X
Z.p10 <- Z
u.x.p10 <- u.x
u.z.p10 <- u.z

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
X.tail<-sim.RevExpU.MGPD(n=125,d=d, a=a, beta=beta, sig=sig, gamma=gamma, MGPD = T,std=T)

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
X.bulk <- rtmvnorm(2600, mean=u.x-c(1.3,1.1), sigma=0.6*sigma, lower=c(0,0))
cond.x <- t(t(X.bulk) < u.x)
X.bulk2 <- X.bulk[rowSums(cond.x)==2,]
plot(X.bulk2)
print(dim(X.bulk2)[1])
X <- rbind(X.bulk2, sweep(X.tail$X,2,u.x,"+"))
print(dim(X)[1])
plot(X,xlab='X1', ylab='X2', main='Simulated tNormal + mGP Data')
plot(X,pch=20, xlim=c(0,7),ylim=c(0,5),xlab='X1', ylab='X2', main='Simulated tNormal + mGP Data')
plot(X,pch=20, xlim=c(0,7),ylim=c(0,5),
     xlab='X1', ylab='X2', main='Simulated tNormal + mGP Data',col=c(rep('grey50',dim(X.bulk2)[1]),rep('black',125)))
segments(u.x[1], 0, u.x[1], u.x[2], col='red', lty=2)
segments(0, u.x[2], u.x[1], u.x[2], col='red', lty=2)
legend('bottomright',   legend = c("Tail Data", "Bulk Data", "Threshold"), lty=c(NA,NA,2),col = c('black', 'grey50', 'red'),         
       pch=c(20,20,NA)) 
chiplot(X, qlim=c(0.02,0.99), ask=F)

# Exponential scale tail data combined with the bulk data
set.seed(1234)
Z.bulk <- rtmvnorm(2570, mean=u.z-c(1.7,1.7), sigma=1.2*matrix(c(1,0.65,0.65,0.8),ncol=2), lower=c(0,0))
cond.z <- t(t(Z.bulk) < u.z)
Z.bulk2 <- Z.bulk[rowSums(cond.z)==2,]
plot(Z.bulk2)
print(dim(Z.bulk2)[1])
Z <- rbind(Z.bulk2, sweep(X.tail$Z,2,u.z,"+"))
print(dim(Z)[1])
plot(Z,pch=20,xlab='Z1', ylab='Z2', main='Simulated tNormal + std GP Data')
plot(Z,pch=20, xlab='Z1', ylab='Z2', main='Simulated tNormal + std GP Data',col=c(rep('grey50',dim(Z.bulk2)[1]),rep('black',125)))
segments(u.z[1], 0, u.z[1], u.z[2], col='red', lty=2)
segments(0, u.z[2], u.z[1], u.z[2], col='red', lty=2)
legend('bottomright',   legend = c("Tail Data", "Bulk Data", "Threshold"), lty=c(NA,NA,2),col = c('black', 'grey50', 'red'),         
       pch=c(20,20,NA)) 
chiplot(Z, qlim=c(0.02,0.99), ask=F)


X.p05 <- X
Z.p05 <- Z
u.x.p05 <- u.x
u.z.p05 <- u.z

save(X.p05,Z.p05,u.x.p05,u.z.p05, 
     X.p10,Z.p10,u.x.p10,u.z.p10, 
     X.p20,Z.p20,u.x.p20,u.z.p20, 
     file=file.path(dir.out,'simulation_data.RData'))

#############################################################################################
#-----------------use positive gamma to solve the left endpoint issue
#############################################################################################

set.seed(1234)

d<-2
a<-c(1.656,1.656) # NB this is 1/\alpha, where \alpha is the parameterization in Kiriliouk, Rootzen, Segers and Wadsworth. 
beta<-c(0,0)
sig<-c(0.571,0.451)
gamma<-c(0.253,0.035)

# Final beta parameter fixed at zero in estimation, so needs to be relative to this

# Simulate with (conditionally) GP(sig,gamma) and (conditionally) exponential margins 
X.tail<-sim.RevExpU.MGPD(n=250,d=d, a=a, beta=beta, sig=sig, gamma=gamma, MGPD = T,std=T)

#  GP(sig,gamma)
plot(X.tail$X, pch=20, xlab=expression(X[1]),ylab=expression(X[2]))

plot(X.tail$Z, pch=20, xlab=expression(X[1]),ylab=expression(X[2]))

# C is to adjust the position of the tail distribution
C <- c(3,1)
u.x <- c(-min(X.tail$X[,1]),-min(X.tail$X[,2])) + C
u.x.1 <- sig/gamma
 
rho=0.65
sigma <- matrix(c(1,rho,rho,0.8),ncol=2)

# GP scale tail data combined with the bulk data
set.seed(1234)
X.bulk <- rtmvnorm(2250, mean=u.x-c(1.5,1.5), sigma=1.5*sigma, lower=c(0,0),upper=u.x)

X <- rbind(X.bulk, sweep(X.tail$X,2,u.x,"+"))
plot(X,pch=20)
plot(X,pch=20, xlab='Z1', ylab='Z2', main='Simulated tNormal + GP Data',
     col=c(rep('grey50',dim(X.bulk)[1]),rep('black',250)))
segments(u.x[1], 0, u.x[1], u.x[2], col='red', lty=2)
segments(0, u.x[2], u.x[1], u.x[2], col='red', lty=2)
legend('bottomright',   legend = c("Tail Data", "Bulk Data", "Threshold"), lty=c(NA,NA,2),col = c('black', 'grey50', 'red'),         
       pch=c(20,20,NA)) 

chiplot(X, qlim=c(0.02,0.99), ask=F)
chiPlot(data=X, ylabel=expression(chi(q)~"prop=0.10"), chimod=NULL, nsim=1000, nq = 50, qmin = 0.5, qmax = 0.99)

X.p10 <- X
u.x.p10 <- u.x
#------------------------------tail proportaiton 5%-----------------


set.seed(1234)

d<-2
a<-c(1.656,1.656) # NB this is 1/\alpha, where \alpha is the parameterization in Kiriliouk, Rootzen, Segers and Wadsworth. 
beta<-c(0,0)
sig<-c(0.571,0.451)
gamma<-c(0.253,0.035)

# Final beta parameter fixed at zero in estimation, so needs to be relative to this

# Simulate with (conditionally) GP(sig,gamma) and (conditionally) exponential margins 
X.tail<-sim.RevExpU.MGPD(n=125,d=d, a=a, beta=beta, sig=sig, gamma=gamma, MGPD = T,std=T)

#  GP(sig,gamma)
plot(X.tail$X, pch=20, xlab=expression(X[1]),ylab=expression(X[2]))

plot(X.tail$Z, pch=20, xlab=expression(X[1]),ylab=expression(X[2]))

# C is to adjust the position of the tail distribution
C <- c(4,4)
u.x <- c(-min(X.tail$X[,1]),-min(X.tail$X[,2])) + C
u.x.1 <- sig/gamma

rho=0.65
sigma <- matrix(c(1,rho,rho,0.8),ncol=2)

# GP scale tail data combined with the bulk data
set.seed(1234)
X.bulk <- rtmvnorm(2375, mean=u.x-c(2,1.8), sigma=1.5*sigma, lower=c(0,0),upper=u.x)

X <- rbind(X.bulk, sweep(X.tail$X,2,u.x,"+"))
plot(X,pch=20)
plot(X,pch=20, xlab='Z1', ylab='Z2', main='Simulated tNormal + GP Data',
     col=c(rep('grey50',dim(X.bulk)[1]),rep('black',250)))
segments(u.x[1], 0, u.x[1], u.x[2], col='red', lty=2)
segments(0, u.x[2], u.x[1], u.x[2], col='red', lty=2)
legend('bottomright',   legend = c("Tail Data", "Bulk Data", "Threshold"), lty=c(NA,NA,2),col = c('black', 'grey50', 'red'),         
       pch=c(20,20,NA)) 

chiplot(X, qlim=c(0.02,0.99), ask=F)
chiPlot(data=X, ylabel=expression(chi(q)~"prop=0.05"), chimod=NULL, nsim=2000, nq = 50, qmin = 0.5, qmax = 0.99)


#----------------regenerate the samples with correct sample size-------------------
d<-2
a<-c(1.656,1.656)
beta<-c(0,0)
sig<-c(0.571,0.451)
gamma<-c(0.253,0.035)
n <- 2500
mu <- c(3.549946,4.412749)
u.x <- c(5.549946,6.212749)
rho=0.65
sigma <- 1.5* matrix(c(1,rho,rho,0.8),ncol=2)
p <- pmvnorm(lower=rep(0,2), upper=u.x, mean=mu, sigma=sigma, keepAttr = F)

set.seed(1234)
X.tail<-sim.RevExpU.MGPD(n=n-floor(n*p),d=d, a=a, beta=beta, sig=sig, gamma=gamma, MGPD = T,std=T)

# GP scale tail data combined with the bulk data
set.seed(1234)
X.bulk <- rtmvnorm(floor(n*p), mean=mu, sigma=sigma, lower=c(0,0),upper=u.x)

X <- rbind(X.bulk, sweep(X.tail$X,2,u.x,"+"))
X.p09 <- X
u.x.p09 <- u.x
save(X.p09, u.x.p09,
     file=file.path(dir.out,'simulation_data_pos_shape_fixed.RData'))

#-------------------generate data with non-stationarity in the threshold---------
set.seed(1234)
x1 <- rnorm(2500,0,1)
x2 <- rnorm(2500,0,2)
x3 <- rnorm(2500,0,3)
x4 <- rnorm(2500,0,4)
X <- cbind(x1,x2,x3,x4)
beta1 <- c(1,2,3,4)
beta2 <- c(-2,-3,4,5)
BETA <- cbind(beta1,beta2)
eta <- X %*% BETA
upper <- 8
lower <- 5
U <- lower + (upper-lower)*sigmoid(0.1*eta)
plot(U)

d<-2
a<-c(1.656,1.656)
beta<-c(0,0)
sig<-c(0.571,0.451)
gamma<-c(0.253,0.035)
N <- 2500
p <- 0.05

Y <- matrix(NA, nrow=N,ncol=2)

mu <- c(5,5.41)
rho=0.5
sigma <- 1.5* matrix(c(1,rho,rho,0.8),ncol=2)

p.vec <- c()
for (i in 1:2500){
    p.vec[i] <- pmvnorm(lower=rep(0,2), upper=U[i,], mean=mu, sigma=sigma, keepAttr = F)
}

rv.uni <- runif(2500)
Y.bulk <- c()
Y.tail <- c()
for (i in 1:N){
    if (rv.uni[i] < p.vec[i]){
        y<- rtmvnorm(1, mean=mu, sigma=sigma, lower=c(0,0),upper=U[i,])
        Y[i,] <- y
        Y.bulk <- rbind(Y.bulk, y)
    }else{
        y <- sim.RevExpU.MGPD(n=1,d=d, a=a, beta=beta, sig=sig, gamma=gamma, MGPD = T,std=T)$X
        Y[i,] <- y + U[i,]   
        Y.tail <- rbind(Y.tail,y + U[i,]  )
    }
}
plot(Y)
plot(Y.tail)
# lines(U, col='red')

save(X, Y, U,
     file=file.path(dir.out,'simulation_data_non_stationary.RData'))
