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
library(Rcpp)
library(pracma)

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
N <- 1000
x1 <- rnorm(N,0,1)
X1 <- cbind(rep(1,N),x1)


x2 <- rnorm(N,0,1)
X2 <- cbind(rep(1,N), x2)

beta1 <- c(0.1, 0.2)
beta2 <- c(0.3, -0.4)

eta <- cbind(X1%*%beta1, X2%*%beta2)
lower <- c(6,6)
upper <- c(8,8)
U <- (lower + sweep(sigmoid(eta),2, upper-lower, "*"))
# U <- sweep(sigmoid(eta),2, upper, "*")
# U <- sweep(pexp(exp(eta)),2, upper, "*")
plot(U)

d<-2
a<-c(1.656,1.656)
beta<-c(0,0)
sig<-c(0.571,0.451)
gamma<-c(0.253,0.035)
p <- 0.05

Y <- matrix(NA, nrow=N,ncol=2)

mu <- c(5,5.41)
rho=0.5
sigma <- 1.5* matrix(c(1,rho,rho,0.8),ncol=2)

t1 <- Sys.time()
p.vec <- c()
for (i in 1:N){
    p.vec[i] <- pmvnorm(lower=rep(0,2), upper=U[i,], mean=mu, sigma=sigma, keepAttr = F)
}
t2 <- Sys.time()
print(t2-t1)

rv.uni <- runif(N)
Y.bulk <- c()
Y.tail <- c()
Y.tail.raw <- c()
for (i in 1:N){
    if (rv.uni[i] < p.vec[i]){
        y<- rtmvnorm(1, mean=mu, sigma=sigma, lower=c(0,0),upper=U[i,])
        Y[i,] <- y
        Y.bulk <- rbind(Y.bulk, y)
    }else{
        y <- sim.RevExpU.MGPD(n=1,d=d, a=a, beta=beta, sig=sig, gamma=gamma, MGPD = T,std=T)$X
        Y[i,] <- y + U[i,]   
        Y.tail.raw <- rbind(Y.tail.raw, y)
        Y.tail <- rbind(Y.tail,y + U[i,]  )
    }
}
plot(Y)
# plot(Y.tail)
# plot(Y.tail.raw)
# lines(U, col='red')

save(X1,X2, Y, U,
     file=file.path(dir.out,'simulation_data_non_stationary.RData'))

#------------------------------Rcpp-------------------------------------
src1 <-
'
double rcpp_sum(NumericVector v){
  double sum = 0;
  for(int i=0; i<v.length(); ++i){
    sum += v[i];
  }

  NumericVector d = {1,2,3};
  Rcout << d;
  return(sum);
}
'

Rcpp::cppFunction(code = src1)
rcpp_sum(1:10000)

src2 <- '
NumericVector my_fun(){
    // calling rnorm()
    Function f("rnorm");   

    // Next code is interpreted as rnorm(n=5, mean=10, sd=2)
    return f(5, Named("sd")=2, Named("mean")=10);
}
'
Rcpp::cppFunction(code = src2)
my_fun()


src3 <-
'
NumericVector row_sum(NumericMatrix M){
  NumericVector v (M.nrow());
  Function f("sum");
  for(int i=0; i<M.nrow(); ++i){
    v[i] = sum(M(i,_));
  }
  return(v);
}
'

Rcpp::cppFunction(code = src3)
row_sum(matrix(c(1,2,3,4),ncol=2))


src4 <-
'
double my_function(NumericVector ub, NumericMatrix M){
  NumericVector lb {0,0};
  Environment pkg = Environment::namespace_env("mvtnorm");
  Function pmvnorm = pkg["pmvnorm"];
  SEXP val = pmvnorm(Named("lower",lb), Named("upper",ub),Named("mean",ub),
                Named("sigma",M), Named("keepAttr", false));
  double res =  Rcpp::as<double>(val);
  return(res);
}
'

Rcpp::cppFunction(code = src4)
my_function(c(2,2),matrix(c(1,0,0,4),ncol=2))



src5 <-
'
NumericVector mix_prob(NumericMatrix UB, NumericVector mu, NumericMatrix M){
  NumericVector lb (mu.length());
  NumericVector prob_vec (UB.nrow());
  
  Environment pkg = Environment::namespace_env("mvtnorm");
  Function pmvnorm = pkg["pmvnorm"];
  
  for(int i=0; i<UB.nrow(); ++i){
    SEXP val = pmvnorm(Named("lower",lb), Named("upper",UB(i,_)),Named("mean",mu),
                Named("sigma",M), Named("keepAttr", false));
    double res =  Rcpp::as<double>(val);
    prob_vec[i] = res;
  }
  return(prob_vec);
}
'

Rcpp::cppFunction(code = src5)
t1 <- Sys.time()
res <- mix_prob(U,mu,sigma)
t2 <- Sys.time()
print(t2-t1)


#----------generate data with non-stationarity in the extremeral dependence---------
seed <- 1235
# set.seed(1234)
set.seed(seed)
N <- 2500
x1 <- rnorm(N,0,1)
X1 <- cbind(rep(1,N), x1)
X1 <-  sweep(X1, 2, c(0,mean(X1[,2])), '-')

x2 <- rnorm(N,0,1)
X2 <- cbind(rep(1,N), x2)
X2 <-  sweep(X2, 2, c(0,mean(X2[,2])), '-')

beta.a1 <- c(0.5, 0.2)
beta.a2 <- c(0.5, -0.4)
beta.b1 <- c(0.25,-0.25)

# beta.a1 <- c(0.5, 0)
# beta.a2 <- c(-0.4, 0)
# beta.b1 <- c(0.25,0)

d<-2
a<- cbind(exp(X1%*%beta.a1), exp(X1%*%beta.a2))
b <- cbind(X1%*%beta.b1, 0)
sig<-c(0.571,0.451)
gamma<-c(0.253,0.135)

u.x <- c(7,7.2)

mu <- c(5,5.41)
rho=0.5
sigma <- 1.5* matrix(c(1,rho,rho,0.8),ncol=2)
p <- pmvnorm(lower=rep(0,2), upper=u.x, mean=mu, sigma=sigma, keepAttr = F)

n.tail <- N-floor(N*p)
Y.tail <- matrix(NA, nrow=n.tail,ncol=2)
set.seed(seed)
for (i in 1:n.tail){
  Y.tail[i,] <- sim.RevExpU.MGPD(n=1,d=d, a=a[i,], beta=b[i,], sig=sig, gamma=gamma, MGPD = T,std=T)$X
}

# GP scale tail data combined with the bulk data
set.seed(seed)
Y.bulk <- rtmvnorm(floor(N*p), mean=mu, sigma=sigma, lower=c(0,0),upper=u.x)
Y <- rbind(sweep(Y.tail,2,u.x,"+"), Y.bulk)
plot(Y)
# plot(Y.tail)
# plot(Y.tail.raw)
# lines(U, col='red')

save(X1,X2, Y, u.x,
     file=file.path(dir.out,'simulation_data_non_stationary_extr_dep_1235.RData'))


#---------generate data with stationarity in the extremeral dependence by MCMC---------
set.seed(1234)
N <- 2500
x1 <- rnorm(N,0,1)
X1 <- cbind(rep(1,N), x1)
X1 <-  sweep(X1, 2, c(0,mean(X1[,2])), '-')

x2 <- rnorm(N,0,1)
X2 <- cbind(rep(1,N), x2)
X2 <-  sweep(X2, 2, c(0,mean(X2[,2])), '-')

# beta.a1 <- c(0.5, 0.2)
# beta.a2 <- c(0.5, -0.4)
# beta.b1 <- c(0.25,-0.25)

beta.a1 <- c(0.5, 0)
beta.a2 <- c(-0.4, 0)
beta.b1 <- c(0.25,0)

d<-2
a<- cbind(exp(X1%*%beta.a1), exp(X1%*%beta.a2))
b <- cbind(X1%*%beta.b1, 0)
sig<-c(0.571,0.451)
gamma<-c(0.253,0.135)

u.x <- c(7,7.2)

mu <- c(5,5.41)
rho=0.5
sigma <- 1.5* matrix(c(1,rho,rho,0.8),ncol=2)
p <- pmvnorm(lower=rep(0,2), upper=u.x, mean=mu, sigma=sigma, keepAttr = F)

n.tail <- N-floor(N*p)
Y.tail <- matrix(NA, nrow=n.tail,ncol=2)


ll.powunif.GPD<-function(theta,x,u,a.ind,lam.ind,sig.ind,gamma.ind, lamfix=FALSE, balthresh=FALSE, marg.scale.ind,marg.shape.ind)
{ 
  x.mat.ind <- 1
  if ((x[1]<=0)&(x[2]<=0)) return(-10e10)
  
  if (is.null(dim(x))){
    d <- length(x)
    x.mat.ind <- 0
  }else{
    d<-dim(x)[2]
  }
  
  a<-theta[a.ind]
  if(length(a)==1)
  {
    a<-rep(a,d)
  }
  
  if(lamfix){
    lam<-rep(1,d)
  }else{
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
  
  if(any(lam<0.01)||any(a<=0.01)||any(sig<=0.001)||any(rej)){return(-10e10)}
  
  ll.uc <- 0
  ll.pc <- 0
  if (!x.mat.ind){
    uc <- comp.gt(x, u)
    if (uc){
      L <- fX.powunif(x=x, a=a, lam=lam, sig=sig, gamma=gamma)
      ll.uc <- log(L)
    }else{
      L2 <- fX.powunif.cens(x=x, u=u, lam=lam, a=a, sig=sig, gamma=gamma)
      ll.pc <- log(L2)
    }
  }else{
    uc<-apply(x,1,comp.gt,u=u)
    
    x.uc<-x[uc,]
    x.pc<-x[!uc,]
    
    L<-apply(x.uc,1,fX.powunif,a=a,lam=lam,sig=sig,gamma=gamma)
    ll.uc<-sum(log(L))
    
    if(sum(!uc)>0)
    {
      L2<-apply(x.pc,1,fX.powunif.cens,a=a,lam=lam,u=u,sig=sig,gamma=gamma)
      ll.pc<- sum(log(L2))
    }
  }
  if (is.nan(ll.uc)|is.nan(ll.pc)|(ll.uc==-Inf)|(ll.uc==Inf)){
    return(-10e10)
  }
  ll<- ll.uc+ll.pc
  
  
  return(ll)
}

theta <- c(a[1,], exp(b[1,1]), sig, gamma)
a.ind <- c(1, 2)
b.ind <- 3
sig.ind <- c(4, 5)
gamma.ind <- c(6, 7)
marg.scale.ind <- c(1, 2)
marg.shape.ind <- c(1, 2)

y <- c(1,0)
ll.powunif.GPD(theta=theta,x=y,u=min(y)-0.01,a.ind=a.ind,
               lam.ind=lam.ind,sig.ind=sig.ind,gamma.ind=gamma.ind,
               marg.scale.ind=marg.scale.ind, marg.shape.ind=marg.shape.ind )


nim_ll.powunif.GPD <- nimbleRcall(function(theta=double(1), x=double(1), u=double(0), a.ind=double(1),
                                            lam.ind=double(0), sig.ind=double(1), gamma.ind=double(1), 
                                            lamfix=logical(0, default = 0), balthresh=logical(0, default = 0), 
                                            marg.scale.ind=double(1), marg.shape.ind=double(1)){}, 
                                   Rfun = 'll.powunif.GPD',
                                   returnType = double(0))


dpowunif.GPD <- nimbleFunction(
  run = function(x=double(1), y=double(1),
                 a.ind=double(1), lam.ind=double(0), lamfix=logical(0, default = 0), 
                 sig.ind=double(1), gamma.ind=double(1),
                 log = logical(0, default = 0)) {
    returnType(double(0))
    
    ll <- nim_ll.powunif.GPD(x=y, theta=x, u=min(y)-0.01, a.ind=a.ind,
                                                lam.ind=lam.ind, sig.ind=sig.ind, gamma.ind=gamma.ind, 
                                                lamfix=lamfix, balthresh=FALSE, 
                                                marg.scale.ind=1:2, marg.shape.ind=1:2)
    if (log) {
      totalProb <- ll
    }else{
      totalProb <- exp(ll)
    }
    return(totalProb)
  })

dpowunif.GPD(x=theta, y=y, a.ind=a.ind, lam.ind=lam.ind, lamfix=FALSE,
             sig.ind=sig.ind, gamma.ind=gamma.ind,log=1)

rpowunif.GPD <- nimbleFunction(
  run = function(n=integer(0), y=double(1), a.ind=double(1), lam.ind=double(0), lamfix=logical(0, default = 0), 
                 sig.ind=double(1), gamma.ind=double(1)) {
    returnType(double(1))
    
    totalProb <- 1:7
    return(totalProb)
  })

registerDistributions(list(
  dpowunif.GPD = list(
    BUGSdist = "dpowunif.GPD(y, a.ind, lam.ind, lamfix, sig.ind, gamma.ind)",
    types = c('value = double(1)', 'y = double(1)', 'a.ind = double(1)', 
              'lam.ind = double(0)', 'lamfix = logical(0)', 'sig.ind = double(1)',
              'gamma.ind = double(1)')
  )))



Simu.GPD.code <- nimbleCode({
  
  theta[1:7] ~ dpowunif.GPD(y=y[1:2], a.ind=1:2,
                        lam.ind=3, lamfix=0, 
                        sig.ind=4:5, gamma.ind=6:7)
  y[1] ~ dunif(-100, 100)
  y[2] ~ dunif(-100, 100)
})



Simu.GPD.model <- nimbleModel(Simu.GPD.code, constants = list(theta=theta),
                                                            check = FALSE)

Simu.GPD.model$setData(list(theta = theta))

cSimu.GPD.model <- compileNimble(Simu.GPD.model)
Simu.GPD.model.confg <- configureMCMC(Simu.GPD.model)

Simu.GPD.model.MCMC<- buildMCMC(Simu.GPD.model.confg)

cSimu.GPD.model.MCMC <- compileNimble(Simu.GPD.model.MCMC, project = Simu.GPD.model)

t1 <- Sys.time()
results <- runMCMC(cSimu.GPD.model.MCMC, niter = 60000,nburnin=20000,thin=4,
                   summary = TRUE, WAIC = FALSE,setSeed = 1234)
t2 <- Sys.time()
print(t2-t1)

plot(results$samples[,'y[2]'],type='l')
plot(results$samples[,c('y[1]','y[2]')])
# plot(density(results$samples[,c('y[1]')]))

cond1 <- results$samples[,'y[1]']>0
plot(density(results$samples[, ][cond1, 'y[1]']))

cond2 <- results$samples[,'y[2]']>0
plot(density(results$samples[, ][cond2, 'y[2]']))

x <- seq(0,10,0.1)
plot(x, dgpd(x, loc=0, scale=sig[1], shape=gamma[1]),type='l')
sp.mg1 <- rgpd(1000, loc=0, scale=sig[1], shape=gamma[1])
sp.mg2 <- rgpd(1000, loc=0, scale=sig[2], shape=gamma[2])
qqplot(sp.mg1,results$samples[, ][cond1, 'y[1]'])
abline(a=0, b=1)
qqplot(sp.mg2,results$samples[, ][cond2, 'y[2]'])
abline(a=0, b=1)

plot(x, dgpd(x, loc=0, scale=sig[2], shape=gamma[2]),type='l')

Y.tail.1 <- sim.RevExpU.MGPD(n=200,d=d, a=a[1,], beta=b[1,], sig=sig, gamma=gamma, MGPD = T,std=T)$X
plot(Y.tail.1)
# plot(density(Y.tail.1[,1]))

cond1 <- Y.tail.1[,1]>0
plot(density(Y.tail.1[cond1, 1]))
qqplot(sp.mg1,Y.tail.1[cond1, 1])
abline(a=0, b=1)

cond2 <- Y.tail.1[,2]>0
plot(density(Y.tail.1[cond2, 2]))
qqplot(sp.mg2,Y.tail.1[cond2, 2])
abline(a=0, b=1)

x <- Y.tail.1
u <- min(x)-0.01
par <- theta
opt1<-optim(nll.powunif.GPD, par=par, u=u, x=x, lamfix=F, a.ind=a.ind, lam.ind=lam.ind,
           sig.ind=sig.ind, gamma.ind=gamma.ind, marg.scale.ind=marg.scale.ind, marg.shape.ind=marg.shape.ind, control=list(maxit=10000,reltol=1e-6),balthresh=FALSE)
par <- opt1$par
opt1.1 <- optim(nll.powunif.GPD, par=par, u=u, x=x, lamfix=F, a.ind=a.ind, lam.ind=lam.ind,
                sig.ind=sig.ind, gamma.ind=gamma.ind, marg.scale.ind=marg.scale.ind, marg.shape.ind=marg.shape.ind, control=list(maxit=10000,reltol=1e-6),balthresh=FALSE)
par <- opt1.1$par
opt1.2 <- optim(nll.powunif.GPD, par=par, u=u, x=x, lamfix=F, a.ind=a.ind, lam.ind=lam.ind,
                sig.ind=sig.ind, gamma.ind=gamma.ind, marg.scale.ind=marg.scale.ind, marg.shape.ind=marg.shape.ind, control=list(maxit=10000,reltol=1e-6),balthresh=FALSE)
par <- opt1.2$par
opt1.3 <- optim(nll.powunif.GPD, par=par, u=u, x=x, lamfix=F, a.ind=a.ind, lam.ind=lam.ind,
                sig.ind=sig.ind, gamma.ind=gamma.ind, marg.scale.ind=marg.scale.ind, marg.shape.ind=marg.shape.ind, control=list(maxit=10000,reltol=1e-6),balthresh=FALSE,hessian=TRUE)

opt1.3$par
1.96*sqrt(diag(solve(opt1.3$hessian)))

0.970 0.960 0.697 0.773 0.838 0.776 0.765

x <- results$samples[,c('y[1]','y[2]')]
u <- min(x)-0.01
par <- theta
opt2<-optim(nll.powunif.GPD, par=par, u=u, x=x, lamfix=F, a.ind=a.ind, lam.ind=lam.ind,
            sig.ind=sig.ind, gamma.ind=gamma.ind, marg.scale.ind=marg.scale.ind, marg.shape.ind=marg.shape.ind, control=list(maxit=10000,reltol=1e-6),balthresh=FALSE)
par <- opt2$par
opt2.1 <- optim(nll.powunif.GPD, par=par, u=u, x=x, lamfix=F, a.ind=a.ind, lam.ind=lam.ind,
                sig.ind=sig.ind, gamma.ind=gamma.ind, marg.scale.ind=marg.scale.ind, marg.shape.ind=marg.shape.ind, control=list(maxit=10000,reltol=1e-6),balthresh=FALSE)
par <- opt2.1$par
opt2.2 <- optim(nll.powunif.GPD, par=par, u=u, x=x, lamfix=F, a.ind=a.ind, lam.ind=lam.ind,
                sig.ind=sig.ind, gamma.ind=gamma.ind, marg.scale.ind=marg.scale.ind, marg.shape.ind=marg.shape.ind, control=list(maxit=10000,reltol=1e-6),balthresh=FALSE)

par <- opt2.2$par
opt2.3 <- optim(nll.powunif.GPD, par=par, u=u, x=x, lamfix=F, a.ind=a.ind, lam.ind=lam.ind,
                sig.ind=sig.ind, gamma.ind=gamma.ind, marg.scale.ind=marg.scale.ind, marg.shape.ind=marg.shape.ind, control=list(maxit=10000,reltol=1e-6),balthresh=FALSE,hessian=TRUE)

0.056 0.030 0.032 0.013 0.015 0.007 0.015
1.96*sqrt(diag(solve(opt2.3$hessian)))

opt2.3$par
-nll.powunif.GPD(theta=theta, x=x,u=u, lamfix=F, a.ind=a.ind, lam.ind=lam.ind,
                sig.ind=sig.ind, gamma.ind=gamma.ind, marg.scale.ind=marg.scale.ind, marg.shape.ind=marg.shape.ind )

-nll.powunif.GPD(theta=opt2.1$par, x=x,u=u, lamfix=F, a.ind=a.ind, lam.ind=lam.ind,
                 sig.ind=sig.ind, gamma.ind=gamma.ind, marg.scale.ind=marg.scale.ind, marg.shape.ind=marg.shape.ind )



set.seed(1234)
for (i in 1:n.tail){
  Y.tail[i,] <- sim.RevExpU.MGPD(n=1,d=d, a=a[i,], beta=b[i,], sig=sig, gamma=gamma, MGPD = T,std=T)$X
}

# GP scale tail data combined with the bulk data
set.seed(1234)
Y.bulk <- rtmvnorm(floor(N*p), mean=mu, sigma=sigma, lower=c(0,0),upper=u.x)
Y <- rbind(sweep(Y.tail,2,u.x,"+"), Y.bulk)
plot(Y)
# plot(Y.tail)
# plot(Y.tail.raw)
# lines(U, col='red')
# 
# save(X1,X2, Y, u.x,
#      file=file.path(dir.out,'simulation_data_non_stationary_extr_dep.RData'))
