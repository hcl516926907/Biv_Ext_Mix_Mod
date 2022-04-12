######################################################################################
# Censored likelihood, estimates the dependence parameters and sigma jointly
# Note that this is for d = 3 only, for constant threshold u
nll.Z.censor.gamma0 <- function(par, data, lambda1){
  d <- ncol(data)
  n <- nrow(data)
  lambda <- c(lambda1, par[1:(d-1)],0); sigma <- par[d]
  if(any(par < 0.01) || any(abs(c(lambda[1]-lambda[2],lambda[1]-lambda[3],lambda[2]-lambda[3])) < 1e-03)){
    return(10^10)
  } else{
      D <- n*log(sum(1/lambda[1:d])) - n*sum(log(lambda[1:d])) - n*log(factorial(d))
      cvec <- -diff(lambda)
      final <- sapply(c(1:n), function(i){
        censor <- which(data[i,] <= 0)
        if(length(censor) == 0){ #u < z1 < z2 < z3
          return((d+1)*log(sum(sapply(c(1:d), function(j) cvec[j]*exp(data[i,j]/sigma)))) - sum(data[i,]/sigma) + 3*log(sigma))
        } else if(length(censor) == 1){ #z1 < u < z2 < z3
          temp1 <- sum(cvec[2:3]*exp(data[i,2:3]/sigma))
          return(log((d*cvec[1]*(sigma^2))/(temp1^(-d) - (cvec[1] + temp1)^(-d))) - sum(data[i,2:3]/sigma))
        } else if(length(censor) == 2){ #z1 < z2 < u < z3
          temp1 <- ((cvec[2] + cvec[3]*exp(data[i,3]/sigma))^(1-d))/cvec[2]
          temp2 <- ((sum(cvec[1:2]) + cvec[3]*exp(data[i,3]/sigma))^(1-d))/sum(cvec[1:2])
          temp3 <- (cvec[1]*(cvec[3]^(1-d))*exp(data[i,3]*(1-d)/sigma))/(cvec[2]*sum(cvec[1:2]))
          return(log((d*(1-d)*cvec[1]*sigma)/(temp1 - temp2 - temp3)) - data[i,3]/sigma)
        }})
      return(D + sum(final))
  }
}


nll.Z.censor <- function(par, data, lambda1){
    d <- ncol(data)
    n <- nrow(data)
    lambda <- c(lambda1, par[1:(d-1)],0); sigma <- par[d]; gamma <- par[d+1]
    maxdat <- max(data[,3])
    if(any(par[1:d] < 0.01) || any(abs(c(lambda[1]-lambda[2],lambda[1]-lambda[3],lambda[2]-lambda[3])) < 1e-03)){
        return(10^10)
    } else if(maxdat > -sigma/gamma && gamma < 0){
      return(10^10)
    } else{
        D <- n*log(sum(1/lambda[1:d])) - n*sum(log(lambda[1:d])) - n*log(factorial(d))
        cvec <- -diff(lambda)
        final <- sapply(c(1:n), function(i){
            censor <- which(data[i,] <= 0)
            if(length(censor) == 0){ #u < z1 < z2 < z3
                return((d+1)*log(sum(sapply(c(1:d), function(j) cvec[j]*((1+(data[i,j]*gamma)/sigma)^(1/gamma))))) -
                           (1/gamma-1)*sum(log(1+(data[i,]*gamma)/sigma)) + 3*log(sigma))
            } else if(length(censor) == 1){ #z1 < u < z2 < z3
                temp1 <- sum(cvec[2:3]*((1+(data[i,2:3]*gamma)/sigma)^(1/gamma)))
                return(log((d*cvec[1]*(sigma^2))/(temp1^(-d) - (cvec[1] + temp1)^(-d))) -
                           (1/gamma-1)*sum(log(1+(data[i,2:3]*gamma)/sigma)))
            } else if(length(censor) == 2){ #z1 < z2 < u < z3
                temp1 <- ((cvec[2] + cvec[3]*((1+(data[i,3]*gamma)/sigma)^(1/gamma)))^(1-d))/cvec[2]
                temp2 <- ((sum(cvec[1:2]) + cvec[3]*((1+(data[i,3]*gamma)/sigma)^(1/gamma)))^(1-d))/sum(cvec[1:2])
                temp3 <- (cvec[1]*(cvec[3]^(1-d))*((1+(data[i,3]*gamma)/sigma)^((1-d)/gamma)))/(cvec[2]*sum(cvec[1:2]))
                return(log((d*(1-d)*cvec[1]*sigma)/(temp1 - temp2 - temp3)) - (1/gamma-1)*log(1+(data[i,3]*gamma)/sigma))
            }})
        return(D + sum(final))
    }
}

# gam0 = TRUE means that gamma = 0 
EstimationOrdered <- function(data, gam0 = TRUE, start, lambda1){
  d <- ncol(data)
  if(gam0 == FALSE){
        result <- optim(start, fn = nll.Z.censor, data = data, lambda1 = lambda1, hessian = TRUE, 
                        control = list(reltol = 1e-15,maxit=10000))
  } else{
      result <- optim(start, fn = nll.Z.censor.gamma0, data = data, lambda1 = lambda1, hessian = TRUE, 
                      control = list(reltol = 1e-15,maxit=10000))
  }
  return(list(result = result$par, minimum = result$value, covar = solve(result$hessian)))
}

#### Yearly probability of a landslide: density functions
dens1 <- function(t2, t3, x1, par){
  if(t2 < t3){
    temp <- (par[2]-par[3])*t2 + par[3]*t3
    return(temp^(-3) - (temp + (par[1]-par[2])*min(exp(x1),t2))^(-3))
  } else{
    return(0)
  }
}
dens2 <- function(t3, x1, x2, par){
  dens1V <- Vectorize(dens1,"t2")
  return(integrate(dens1V,lower=0,upper=exp(x2),t3=t3,x1=x1,par=par)$value)
}
dens2V <- Vectorize(dens2,"t3")

densToInt <- function(t, x, par){
  if(t[1] < t[2]){
    temp <- (par[2]-par[3])*t[1] + par[3]*t[2]
    return(temp^(-3) - (temp + (par[1]-par[2])*min(exp(x),t[1]))^(-3))
  } else{
    return(0)
  }
}

######### Marginal QQ-plots ###########
q.gpd<-function(q,sig,xi){
  if(abs(xi)>1e-6){return(sig*((1-q)^(-xi)-1)/xi)}
  else{return(-sig*log(1-q))}
}

GPD.diag <- function(x,sig,gam){
  d <- dim(x)[2]
  n <- dim(x)[1]
  par(mfrow=c(1,d))
  for(j in 1:d){
    n[j]<-sum(x[,j]>0)
    gpdQ <- q.gpd((1:n[j])/(n[j]+1),sig=sig[j],xi=gam[j])
    par(cex.lab=1.5,cex.axis=1.5,cex.main=1.5,mar=c(5,4.4,4,2)+0.6)
    plot(gpdQ,sort(x[x[,j]>0,j]),ylab="Empirical",xlab="Model",cex.lab=2,cex.axis=2)
    abline(a=0,b=1,lwd=2)
    lowUnif <- sapply(c(1:n[j]), function(i) qbeta(0.025,i,n[j]+1-i))
    highUnif <- sapply(c(1:n[j]), function(i) qbeta(0.975,i,n[j]+1-i))
    lines(gpdQ,q.gpd(lowUnif, sig=sig[j],xi=gam[j]),lty=2,lwd=2,col="gray")
    lines(gpdQ,q.gpd(highUnif, sig=sig[j],xi=gam[j]),lty=2,lwd=2,col="gray")
  }
}

######### Diagnostic based on a set A #########################
diagA.plot<- function(t, u, dataY, I){
  dataSt <- dataY - u
  dataAll <- dataSt[dataSt[,3] > 0,]
  dataSt2 <- dataY - u - c(7.36, 9.18, 9.55)*log(t)
  dataAll2 <- dataSt2[dataSt[,3] > 0,]
  length(which(dataAll[,I] > 0))/(t*length(which(dataAll2[,I] > 0)))
}

############### Model-based chi ###############################
chi12 <- function(par){return(1 - par[1]/(2*(par[1]+par[2])))}
chi13 <- function(par){
  temp1 <- (par[3]+2*par[2])*(par[2]+2*par[3])*(par[2]*par[3] + par[1]*par[3] + par[1]*par[2])
  return(1 - (par[1]*((par[2]+par[3])^3))/temp1)
}
chi23 <- function(par){
  temp1 <- (par[1] + 2*par[2])*(par[2] + 2*par[1])*(par[2]*par[3] + par[1]*par[3] + par[1]*par[2])
  return(1 - ((par[1]*par[2]*((par[1]+par[2])^2))/temp1))
}
chi123 <- function(par){ 
  temp <- ((par[1]*par[2])*(4*par[1]*par[2] + par[1]*par[3] + 3*(par[2]^2) + par[2]*par[3]))/(3*(2*par[1]+par[2])*(2*par[2]+par[3])*(par[1]*par[2] + par[1]*par[3] + par[2]*par[3]))
  return(1 - par[1]/(2*(par[1]+par[2])) - temp)
}

################## Empirical chi #################################
chiEmp <-function(data,nq,qmin,qmax){
  n<-nrow(data)
  datatr<-(apply(data,2,rank) - 0.5)/n   
  qlim<-c(qmin,qmax)
  u <- seq(qlim[1], qlim[2], length = nq)
  cu<-sapply(c(1:nq), function(i) mean(apply(datatr,1,min) >= u[i]))
  return(cbind(u,cu/(1-u)))
}

chiPlot <- function(data, ylabel, chimod, nsim, nq = 35, qmin = 0.5, qmax = 0.95){
  tmp <- matrix(,nrow=nsim,ncol=nq)
  n <- nrow(data)
  for(j in 1:nsim){
    nsample <- sample(1:n,size=n,replace=T)
    newdata <- data[nsample,]
    tmp[j,] <- chiEmp(newdata,nq=nq,qmin=qmin,qmax=qmax)[,2]
  }
  CIlow <- apply(tmp,2,quantile,0.025)
  CIhigh <- apply(tmp,2,quantile,0.975)
  chi <-chiEmp(data,nq=nq,qmin=qmin,qmax=qmax)
  
  par(cex.lab=2,cex.axis=2,cex.main=1.5,mar=c(5,4.4,4,2)+0.9)
  plot(chi[,1],chi[,2],ylim=c(0,1),xlab="q",
       ylab=ylabel,lwd=2)
  lines(chi[,1],CIlow,lty=3,lwd=2)
  lines(chi[,1],CIhigh,lty=3,lwd=2)
  abline(h = chimod, lwd = 2)
}

### Goodness-if-fit test from Einmahl et al. (2016)

ell12 <- function(par){return(1 + par[1]/(2*(par[1]+par[2])))}
ell13 <- function(par){
  temp1 <- (par[3]+2*par[2])*(par[2]+2*par[3])*(par[2]*par[3] + par[1]*par[3] + par[1]*par[2])
  return(1 + (par[1]*((par[2]+par[3])^3))/temp1)
}
ell23 <- function(par){
  temp1 <- (par[1] + 2*par[2])*(par[2] + 2*par[1])*(par[2]*par[3] + par[1]*par[3] + par[1]*par[2])
  return(1 + ((par[1]*par[2]*((par[1]+par[2])^2))/temp1))
}
ell123 <- function(par){ #lambdas
  temp <- ((par[1]*par[2])*(4*par[1]*par[2] + par[1]*par[3] + 3*(par[2]^2) + par[2]*par[3]))/(3*(2*par[1]+par[2])*(2*par[2]+par[3])*(par[1]*par[2] + par[1]*par[3] + par[2]*par[3]))
  return(1 + par[1]/(2*(par[1]+par[2])) + temp)
}

stdf12<- function(x,l){
  a <- c(l[1]*x[1], (l[1]*l[2]*x[2])/(l[1] + l[2]))
  R1 <- ((a[1]>a[2])*a[1]*l[1])*(1/(l[1]^2) - (a[2]/(a[2]*(l[1]-l[2]) + a[1]*l[2]))^2)
  R2 <- ((a[2]*l[1]*l[2])/(l[1]-l[2]))*(1/(l[2]^2) - (a[1]/(min(a[1],a[2])*(l[1]-l[2]) + a[1]*l[2]))^2)
  return(R1 + R2)
}

stdf13<- function(x,l){
  a <- c(l[1]*x[1], (prod(l)*x[2])/(l[2]*l[3] + l[1]*l[3] + l[1]*l[2]))
  temp1 <- 1/(l[1]^2) - (a[2]/(a[2]*(l[1]-l[3]) + a[1]*l[3]))^2
  temp2 <- 1/(l[1]^2) - (a[2]/(a[2]*(l[1]-l[2]) + a[1]*l[2]))^2
  R1 <- ((a[1]*l[1]*(a[1]>a[2]))/(l[2]-l[3]))*(l[2]*temp1 - l[3]*temp2)
  temp3 <- 1/(l[3]^2) - (a[1]/(a[1]*l[3] + min(a[1],a[2])*(l[1]-l[3])))^2
  temp4 <- 1/(l[2]^2) - (a[1]/(a[1]*l[2] + min(a[1],a[2])*(l[1]-l[2])))^2
  R2 <- ((a[2]*prod(l))/(l[2]-l[3]))*(temp3/(l[1]-l[3]) - temp4/(l[1]-l[2]))
  return(R1 + R2)
}

stdf23<- function(x,l){
  a <- c((l[1]*l[2]*x[1])/(l[1] + l[2]), (prod(l)*x[2])/(l[2]*l[3] + l[1]*l[3] + l[1]*l[2]))
  temp1 <- (a[2]/(a[2]*(l[2]-l[3]) + a[1]*l[3]))^2
  temp2 <- (a[2]/(a[2]*(l[1]-l[3]) + a[1]*l[3]))^2
  R1 <- ((a[1]*l[1]*l[2]*(a[1]>a[2]))/(l[1]-l[2]))*(1/(l[2]^2) - 1/(l[1]^2) - temp1 + temp2)
  temp3 <- 1/(l[3]^2) - (a[1]/(a[1]*l[3] + min(a[1],a[2])*(l[2]-l[3])))^2
  temp4 <- 1/(l[3]^2) - (a[1]/(a[1]*l[3] + min(a[1],a[2])*(l[1]-l[3])))^2
  R2 <- ((a[2]*prod(l))/(l[1]-l[2]))*(temp3/(l[2]-l[3]) - temp4/(l[1]-l[3]))
  return(R1 + R2)
}


stdf <- function(x, l){
  a <- c(l[1]*x[1], (l[1]*l[2]*x[2])/(l[1] + l[2]), (prod(l)*x[3])/(l[2]*l[3] + l[1]*l[3] + l[1]*l[2]))
  #### part 1
  I1anpt1 <- a[1]/(l[2]*a[1] + a[2]*(l[1]-l[2]))
  I1anpt2 <- (a[1]*a[3])/(a[1]*a[3]*(l[2]-l[3]) + a[2]*a[3]*(l[1]-l[2]) + a[1]*a[2]*l[3])
  I1an <- ((1*(a[1] > a[2])*(a[2] > a[3])*a[2]*prod(l))/((l[1]-l[2])*l[3]))*(I1anpt1^2 - I1anpt2^2)
  I2anpt <- a[3]/(a[3]*(l[1]-l[3]) + a[2]*l[3])
  I2an <- ((1*(a[1] > a[2])*(a[2] > a[3])*a[1]*prod(l))/((l[1]-l[2])*l[3]))*(1/(l[1]^2) - I2anpt^2)
  I3an <- ((1*(a[1] > a[2])*(a[2] > a[3])*a[1]*prod(l))/(((l[1]-l[2])^2)*l[3]))*(I1anpt1 - I1anpt2)
  I4an <- ((1*(a[1] > a[2])*(a[2] > a[3])*a[1]*prod(l))/(((l[1]-l[2])^2)*l[3]))*(1/l[1] - I2anpt)
  I5anpt1 <- a[1]/(a[1]*l[3] + a[3]*(l[1]-l[3]))
  I5anpt2 <- (a[1]*a[2])/(a[1]*a[2]*l[3] + a[2]*a[3]*(l[1]-l[2]) + a[1]*min(a[2],a[3])*(l[2]-l[3]))
  I5anpt3 <- ((a[1]*prod(l))/(l[3]*(l[1]-l[2])))*((a[3]/(max(a[1],a[2],a[3])*l[3] + a[3]*(l[1]-l[3])))^2)
  I5an <- ((1*(a[1] > max(a[2],a[3]))*a[3]*prod(l))/((l[1]-l[2])*(l[2]-l[3])))*(I5anpt1^2 - I5anpt2^2) + I5anpt3
  I6an <- ((a[1]*prod(l))/(l[3]*(l[1]-l[2])))*((a[3]/(a[3]*(l[1]-l[3]) + max(a[2],a[3])*l[3]))^2)
  I7anpt1 <- ((a[1]*prod(l))/(((l[1]-l[2])^2)*(l[1]-l[3])))*(1/l[3] - 1/(l[3] + min(1,a[3]/a[2],a[3]/a[1])*(l[1]-l[3])))
  I7an <- ((1*(a[1] > max(a[2],a[3]))*a[1]*prod(l))/(((l[1]-l[2])^2)*(l[2]-l[3])))*(I5anpt1 - I5anpt2) + I7anpt1
  I8an <- ((a[1]*prod(l))/((l[1]-l[3])*((l[1]-l[2])^2)))*(1/l[3] - a[2]/(a[2]*l[3] + min(a[2],a[3])*(l[1]-l[3])))
  Final1 <- I1an - I2an + I3an - I4an + I5an - I6an + I7an - I8an
  ### part 2
  J1an <-(((a[2]>a[3])*a[2]*prod(l))/(l[3]*(l[1]-l[2])))*(1/(l[2]^2) - (a[3]/(a[3]*(l[2]-l[3]) + a[2]*l[3]))^2)
  J2anpt1 <- a[1]/(a[1]*l[2] + a[3]*(l[1]-l[2]))
  J2anpt2 <- (a[1]*a[3])/(a[1]*a[3]*(l[2]-l[3]) + min(a[1],a[2])*(a[1]*l[3] + a[3]*(l[1]-l[2])))
  J2an <- ((a[1]*a[2]*prod(l)*(a[3] < a[1])*(a[3] < a[2]))/((l[1]-l[2])*(a[1]*l[3] + a[3]*(l[1]-l[2]))))*(J2anpt1^2 - J2anpt2^2)
  J3anpt1 <- a[3]/(a[3]*(l[1]-l[3]) + max(a[1],a[3])*l[3])
  J3anpt2 <- a[3]/(a[3]*(l[1]-l[3]) + a[2]*l[3])
  J3an <- (((a[2] > a[3])*(a[2] > a[1])*a[2]*prod(l))/(l[3]*(l[1]-l[2])))*(J3anpt1^2 - J3anpt2^2)
  J4anpt1 <- a[1]/(a[1]*l[2] + a[3]*(l[1]-l[2]))
  J4anpt2 <- (a[1])/(a[1]*l[2] + min(a[1],a[2])*(l[1]-l[2]))
  J4an <- ((a[2]*l[1]*l[2]*(a[3]<a[1])*(a[3] < a[2]))/(l[1]-l[2]))*(J4anpt1^2 - J4anpt2^2)
  J5anpt1 <- (a[1]*a[3])/(a[1]*a[3]*(l[2]-l[3]) + min(a[1],a[2],a[3])*(a[1]*l[3] + a[3]*(l[1]-l[2])))
  J5an <- ((a[2]*l[1]*l[2]*a[3])/(a[3]*(l[1]-l[2]) + a[1]*l[3]))*((1/(l[2]-l[3]))^2 - J5anpt1^2)
  J6anpt1 <- (a[1]*a[3])/(a[1]*a[3]*(l[2]-l[3]) + min(a[1],a[2])*(a[1]*l[3] + a[3]*(l[1]-l[2])))
  J6an <- ((a[2]*l[1]*l[2]*a[3])/(a[3]*(l[1]-l[2]) + a[1]*l[3]))*((1/(l[2]-l[3]))^2 - J6anpt1^2)
  Final2 <- J1an - J2an - J3an + J4an + J5an - J6an
  ### part 3
  Ktemp <- (a[1]*(l[2]-l[3]) + a[2]*(l[1]-l[2]))
  Ktemp2 <- (a[1]*(l[2]-l[3]) + min(a[1],a[2])*(l[1]-l[2]))
  Kanpt1 <- (1/(l[1]-l[3]))*(1/(l[3]^2) - (a[1]/(a[1]*l[3] + min(a[1],a[3])*(l[1]-l[3])))^2)
  Kanpt2 <- (a[2]/(Ktemp))*(1/(l[3]^2) - ((a[1]*a[2])/(a[1]*a[2]*l[3] + Ktemp*min(a[2],a[3])))^2)
  Kanpt3 <- (1/(l[1]-l[2]))*((a[1]/(a[1]*l[2] + (l[1]-l[2])*min(a[2],a[3])))^2 - ((a[1]/(a[1]*l[2] + (l[1]-l[2])*min(a[1],a[3])))^2))
  Kanpt4 <-(1/(l[2]-l[3]))*(1/(l[3]^2) - (a[2]/(a[2]*l[3] + min(a[2],a[3])*(l[2]-l[3])))^2)
  Kanpt5 <- (a[1]/Ktemp2)*(1/(l[3]^2) - ((a[1]*a[2])/(a[1]*a[2]*l[3] + Ktemp2*min(a[2],a[3])))^2)
  Final3 <- ((a[3]*prod(l)*(a[1]>a[2]))/((l[2]-l[3])))*(Kanpt1 - Kanpt2 - Kanpt3) + 
    ((a[3]*prod(l))/((l[1]-l[2])))*(Kanpt4 - Kanpt5)
  return(Final1 + Final2 + Final3)
}

fec<-function(par, cst){
  dim <- length(which(cst > 1e-05))
  if(dim == 2){
    if(cst[1] > 0 && cst[2] > 0){
      return(stdf12(cst[1:2], par))
    } else if(cst[1] > 0 && cst[3] > 0){
      return(stdf13(cst[c(1,3)], par))
    } else{
      return(stdf23(cst[2:3], par))
    }
  } else if(dim == 1){
    return(1)
  } else if(dim == 0){
    return(0)
  } else{
    return(stdf(cst, par))
  }
}
ecder <- function(par,cst){
  dim <- length(which(cst > 1e-05))
  if(dim == 2){
    if(cst[1] > 0 && cst[2] > 0){
      return(c(grad(stdf12,x=cst[1:2],l=par),0))
    } else if(cst[1] > 0 && cst[3] > 0){
      temp <- grad(stdf13,x=cst[c(1,3)],l=par)
      return(c(temp[1],0,temp[2]))
    } else{
      return(c(0,grad(stdf23,x=cst[2:3],l=par)))
    }
  } else{
    return(grad(stdf,x=cst,l=par))
  }
}

Asym <- function(C1,C2,theta,d){
  ec1 <- fec(theta, C1)
  ec2 <- fec(theta, C2)
  der1 <- ecder(theta,C1)
  der2 <- ecder(theta,C2)
  T1 <- ec1 + ec2 - fec(theta, pmax(C1,C2))
  T2 <- sum(sapply(c(1:d), function(k) {
    Ctemp <- rep(0,d)
    Ctemp[k] <- C2[k]
    der2[k]*(ec1 + C2[k] - fec(theta,pmax(C1,Ctemp)))
  }))
  T3 <- sum(sapply(c(1:d), function(k) {
    Ctemp <- rep(0,d)
    Ctemp[k] <- C1[k]
    der1[k]*(ec2 + C1[k] - fec(theta,pmax(Ctemp,C2)))
  }))
  T4 <- sum(sapply(c(1:d), function(k) sapply(c(1:d), function(l){
    Ctemp1 <- Ctemp2 <- rep(0,d)
    Ctemp1[k] <- C1[k]
    Ctemp2[l] <- C2[l]
    der1[k]*der2[l]*(C1[k] + C2[l] - fec(theta,pmax(Ctemp1,Ctemp2)))
  })))
  return(T1 - T2 - T3 + T4)
}

AsymVar <- function(theta,indices){
  qd <- dim(indices)
  resmat <- matrix(0,nrow=qd[1],ncol=qd[1])
  resmat[lower.tri(resmat, diag = TRUE)] <- unlist(sapply(c(1:qd[1]), function(i) sapply(c(i:qd[1]),
                                                                                         function(j) Asym(indices[i,],indices[j,],theta,qd[2]))))
  if(all(dim(resmat) == c(1,1))){
    return(resmat)
  } else{
    return(t(resmat) + resmat - diag(diag(resmat)))
  }
}


psiDot<-function(par){return(rbind(grad(ell12,x=par),grad(ell13,x=par),grad(ell23,x=par),grad(ell123,x=par)))}


ECminimize<-function(par,indices, totlist,w,fix){
  if(any(par <= 1e-05)){
    return(10^6)
  } else{
    pars<- c(fix,par)
    res <- c(ell12(pars),ell13(pars),ell23(pars),ell123(pars)) - totlist
    return(t(res) %*% w %*% res)
  }
}

ECminimizecu<-function(par,indices,totlist,fix){
  if(any(par <= 1e-05)){
    return(10^6)
  } else{
    pars<- c(fix,par)
    res <- c(ell12(pars),ell13(pars),ell23(pars),ell123(pars)) - totlist
    Sigma <- AsymVar(pars, indices)
    return(t(res) %*% solve(Sigma) %*% res)
  }
}


ECestimator <- function(x, indices, k, fix = 1, covMat = FALSE, iterate = FALSE) {
  ranks <- apply(x, 2, rank)
  Omega <- diag(nrow(indices))
  totL <- apply(indices, 1, function(j) stdfEmp(ranks, k = k, cst = j))
  theta <- thetaPilot <- optim(c(5,8), fn = ECminimize, totlist = totL, indices = indices, w = Omega, fix = fix)
  if (iterate) {
    theta <- optim(thetaPilot$par, fn = ECminimizecu, totlist = totL, indices = indices, fix = fix)
  }
  covMatrix <- NULL
  if(covMat){
    phidot <- psiDot(c(fix,theta$par))
    GammaMat <- AsymVar(theta = c(fix,theta$par), indices = indices)
    if(iterate){
      covMatrix <- solve(t(phidot[,2:3]) %*% solve(GammaMat) %*% phidot[,2:3]) / k
    } else{
      temp <- solve(t(phidot[,2:3]) %*% Omega %*% phidot[,2:3]) %*% t(phidot[,2:3])
      covMatrix <- (temp %*% Omega %*% GammaMat %*% Omega %*% t(temp)) / k
    }
  }
  return(list(theta = theta$par, thetaPilot = thetaPilot$par, covMatrix = covMatrix, value = theta$value))
}

