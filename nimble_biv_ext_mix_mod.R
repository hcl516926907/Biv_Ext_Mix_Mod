source("KRSW/RevExp_U_Functions.r")
source("KRSW/CommonFunctions.r")
library(nimble, warn.conflicts = F)
library(mvtnorm)

R_pmnorm_chol <- function(lower, upper, mean, cholesky){
  sigma <- t(cholesky) %*% cholesky
  return(pmvnorm(lower=lower, upper=upper, mean=mean, sigma=sigma, keepAttr = F))
}



pmnorm_chol <- nimbleRcall(function(lower = double(1), upper = double(1), 
                                        mean=double(1),cholesky=double(2)){}, 
                               Rfun = 'R_pmnorm_chol',
                               returnType = double(0))


nim_nll_powunif_GPD <- nimbleRcall(function(x=double(2), theta=double(1), u=double(0), a.ind=double(1),
                                            lam.ind=double(1), sig.ind=double(1), gamma.ind=double(1), 
                                            lamfix=logical(0, default = 0), balthresh=logical(0, default = 0), 
                                            marg.scale.ind=double(1), marg.shape.ind=double(1)){}, 
                                   Rfun = 'nll.powunif.GPD',
                                   returnType = double(0))

thres <- c(5.55,6.213)
theta <- c(1.656, 0.571, 0.451, 0.253, 0.035)
a.ind <- c(1)
lam.ind <- c(1)
sig.ind <- c(2,3)
gamma.ind <- c(4,5)
marg.scale.ind <- c(1,2)
marg.shape.ind <- c(1,2)

cond <- (X[,1]>thres[1]) | (X[,2]>thres[2])
y.tail <- cbind(X[cond,1] - thres[1],
                X[cond,2] - thres[2])
y.single <- y.tail[1,]

rbind(y.single,y.single)
nll.powunif.GPD(theta=theta, x=rbind(y.single,y.single), u=min(y.tail)-0.01, a.ind=a.ind, lam.ind=lam.ind, sig.ind=sig.ind,
                gamma.ind=gamma.ind, marg.scale.ind=marg.scale.ind, 
                marg.shape.ind = marg.shape.ind)

nim_nll_powunif_GPD(x=rbind(y.single,y.single), theta=theta, u=min(y.tail)-0.01, a.ind=a.ind, lam.ind=lam.ind, sig.ind=sig.ind,
                gamma.ind=gamma.ind, marg.scale.ind=marg.scale.ind, 
                marg.shape.ind = marg.shape.ind)


dbiextmix <- nimbleFunction(
  run = function(x=double(1), theta=double(1), thres=double(1), mu=double(1), 
                 chol=double(2), d=double(0, default=2),
                 a.ind=double(1), lam.ind=double(1), lamfix=logical(0, default = 0), 
                 sig.ind=double(1), gamma.ind=double(1),
                 log = logical(0, default = 0)) {
    returnType(double())
    
    tail.ind <- any(x>thres)
    
    pi <- pmnorm_chol(lower=rep(0,d), upper=thres, mean=mu, cholesky = chol)
    
    sig <- theta[sig.ind]
    gamma <- theta[gamma.ind]
    eta <- -sig/gamma
    
    dtail <- 0
    dbulk <- 0
    
    if(tail.ind){
      if (all((x-thres)>eta)){
        y.tail <- t(matrix(c(x-thres, x-thres),nrow=d))
        print(y.tail)
        twollt <- -nim_nll_powunif_GPD(x=y.tail, theta=theta, u=min(y.tail)-0.01, a.ind=a.ind,
                                      lam.ind=lam.ind, sig.ind=sig.ind, gamma.ind=gamma.ind, 
                                      lamfix=lamfix, balthresh=FALSE, 
                                      marg.scale.ind=1:2, marg.shape.ind=1:2)
        print(twollt)
        dtail <- exp(twollt/2)
      }
    }else{
    
    dbulk <- dmnorm_chol(x, mean=mu, cholesky = chol, prec_param = FALSE)
    print(1)
    }
    
    totalProb <- dtail + dbulk
    if (log) return(log(totalProb))
    return(totalProb)
  })

dbiextmix(y.single+thres,theta=theta, thres=thres, mu=c(3,4), 
          chol=diag(2),
          a.ind=a.ind, lam.ind=lam.ind, lamfix=TRUE, 
          sig.ind=sig.ind, gamma.ind=gamma.ind)




dSS <- nimbleFunction(
  run = function(x = double(1), mu = double(1), cov=double(2), theta=double(0), log = logical(0, default = 0)) {
    returnType(double())
    cov.cho <- chol(cov)
    comp1 <- theta*dmnorm_chol(x, mean=mu[1:2], cholesky=cov.cho, prec_param = FALSE)
    comp2 <- (1-theta)*dmnorm_chol(x, mean=mu[3:4], cholesky=cov.cho, prec_param = FALSE)
    
    totalProb <- comp1 + comp2
    if (log) return(log(totalProb))
    return(totalProb)
  })


test <- nimbleFunction(
  run = function(lower = double(1), upper = double(1), mean=double(1),cholesky=double(2)) {
    returnType(double(0))
    return(pmnorm_chol(lower=lower, upper=upper, mean=mean, cholesky = cholesky))
  })

pmvnorm(lower=rep(0,2), upper=c(5,5), mean=c(3,3), sigma=diag(2), keepAttr = F)
test(lower=rep(0,2), upper=c(5,5), mean=c(3,3), cholesky=diag(2))
## Not run: 
## Say we want an R function that adds 2 to every value in a vector
add2 <- function(x) {
  x + 2 
}
Radd2 <- nimbleRcall(function(x = double(1)){}, Rfun = 'add2',
                     returnType = double(1))
demoCode <- nimbleCode({
  for(i in 1:4) {x[i] ~ dnorm(0,1)} 
  z[1:4] <- Radd2(x[1:4])
})
demoModel <- nimbleModel(demoCode, inits = list(x = rnorm(4)),
                         check = FALSE, calculate = FALSE)
CdemoModel <- compileNimble(demoModel)