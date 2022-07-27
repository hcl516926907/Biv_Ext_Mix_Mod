############################### Normal Distribution #############################
set.seed(1234)
sp <- rnorm(200,5,0.1)

nll <- function(theta,x){
  return(-sum(dnorm(x, mean=theta[1], sd=theta[2], log=T)))
}

fit <- optim(nll, par=runif(2), x=sp)
for (i in 1:5){
  if (i<5){
    fit <- optim(nll, par=fit$par, x=sp)
  }else{
    fit <- optim(nll, par=fit$par, x=sp, hessian=T)
  }
}

sd.hess <- sqrt(diag(solve(fit$hessian)))


#----------------------------Monte Carlo Estimation--------------------------------
mle.mat <- c()
sd.mat <- c()
hess.list <- list()
for (i in 1:10000){
  set.seed(i)
  sp <- rnorm(200,5,0.1)
  
  fit <- optim(nll, par=runif(2), x=sp)
  for (j in 1:5){
    if (j<5){
      fit <- optim(nll, par=fit$par, x=sp)
    }else{
      fit <- optim(nll, par=fit$par, x=sp, hessian=T)
    }
  }
  
  sd.hess <- sqrt(diag(solve(fit$hessian)))
  hess.list[[i]] <- fit$hessian
  mle.mat <- rbind(mle.mat, fit$par)
  sd.mat <- rbind(sd.mat, sd.hess)
}


sd.mt.hess <- colMeans(sd.mat, na.rm=TRUE)

# average the hessian matrix before doing the inverse operation.
sd.mt.hess.1 <- sqrt(diag(solve(Reduce("+", hess.list) / length(hess.list))))
sd.mt <- apply(mle.mat, 2, sd)
