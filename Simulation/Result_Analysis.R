library(evd)
library(scoringRules)
library(mvtnorm)
library(tmvtnorm)
library(HDInterval)
source("KRSW/RevExp_U_Functions.r")
source("KRSW/CommonFunctions.r")

dir.out <- "/home/pgrad2/2448355h/My_PhD_Project/01_Output/Biv_Ext_Mix_Mod/Simulation"


##################################--------
load(file=file.path(dir.out, 'Scenario_1.1_itr1_20_lamfix.RData'))
chain_out <- chain_res

para.name <- colnames(chain_out[[1]][[1]]$samples)
post.mean.mat <- matrix(NA, nrow=length(chain_out), ncol= length(para.name))
colnames(post.mean.mat) <- para.name
cover.ind.mat <- matrix(NA, nrow=length(chain_out), ncol= length(para.name))
colnames(cover.ind.mat) <- para.name

for (itr in 1:20){
  rhat.seq <- c()
  ess.seq <- c()
  for (name in para.name){
    post.sp <- cbind(chain_out[[itr]][[1]]$samples[,name],
                     chain_out[[itr]][[2]]$samples[,name],
                     chain_out[[itr]][[3]]$samples[,name]
                     )
    rhat.seq <- c(rhat.seq, rhat(post.sp))
    ess.seq <- c(ess.seq, ess_basic(post.sp))
  }
  convg.stat <- data.frame(para.name,rhat.seq,ess.seq )
  convg.stat
  
  samples.all <- rbind(chain_out[[itr]][[1]]$samples[,],
                       chain_out[[itr]][[2]]$samples[,],
                       chain_out[[itr]][[3]]$samples[,]
                       )

  
  # True parameters
  a <- c(0.5, 1.2)
  beta <- c(0, 0)
  sig <- c(0.5, 1.2)
  gamma <- c(0.3, 0.1)
  mu <- c(3.5, 4.0)
  sd1 <- 1
  sd2 <- 1.5
  rho <- 0.7
  sigma <- matrix(c(sd1^2, rho*sd1*sd2, rho*sd1*sd2, sd2^2),ncol=2)
  chol_corr <- chol(matrix(c(1,rho,rho,1),ncol=2))
  u.x <- c(5.5, 6.7)
  # u.x <- c(4.7,6)
  true.para <- c(1,0,chol_corr[1,2],chol_corr[2,2], mu, sd1, sd2, a, exp(3),sig, gamma, u.x)
  
  post.mean <- colMeans(samples.all[,])
  post.mean.mat[itr,] <- post.mean
  post.lb <- apply(samples.all, 2, quantile, 0.025)
  post.ub <- apply(samples.all, 2, quantile, 0.975)
  cover.ind <- (true.para >= post.lb)&(true.para <= post.ub)
  # print(cover.ind)
  cover.ind.mat[itr, ] <- cover.ind
}

avg.post.mean <- colMeans(post.mean.mat)

plot(samples.all[,'thres[2]'],type='l')
plot(density(samples.all[,'thres[2]']))

plot((chain_out[[2]][[2]]$samples[,'thres[2]']),type='l')








post.pred <- function(n, samples){
  Y.pred <- matrix(NA,nrow=n, ncol=2)
  idx <- sample(nrow(samples),size=n, replace=TRUE)
  d <- 2
  for (i in 1:length(idx)){
    a <-  samples[idx[i], c('theta[1]','theta[2]')]
    beta <-  c(log(samples[idx[i], 'theta[3]']), 0)
    sig <- samples[idx[i], c('theta[4]','theta[5]')]
    gamma <- samples[idx[i], c('theta[6]','theta[7]')]
    mu <- samples[idx[i], c('mu[1]','mu[2]')]
    sd1 <- samples[idx[i], 'sds[1]']
    sd2 <- samples[idx[i], 'sds[2]']
    corr.chol <- matrix(samples[idx[i],c('Ustar[1, 1]','Ustar[2, 1]',
                                    'Ustar[1, 2]','Ustar[2, 2]')],ncol=2)
    sigma <- diag(c(sd1,sd2))%*%t(corr.chol)%*%corr.chol%*%diag(c(sd1,sd2))
    thres <- samples[idx[i], c('thres[1]','thres[2]')]
    p <- pmvnorm(lower=rep(0,2), upper=thres, mean=mu, sigma=sigma, keepAttr = F)
    u <- runif(1)
    if (u<p){
      Y.pred[i,] <- rtmvnorm(1, mean=mu, sigma=sigma, lower=c(0,0),upper=thres)
    }else{
      Y.tail <- sim.RevExpU.MGPD(n=1,d=d, a=a, beta=beta, sig=sig, gamma=gamma, MGPD = T,std=T)$X
      Y.pred[i,] <- thres + Y.tail
    }
  }
  return(Y.pred)
}

Y.pred <- post.pred(2000, samples)
mu <- c(3.5, 4.5)
sd1 <- 1.22
sd2 <- 1.10
rho <- 0.72
sigma <- matrix(c(sd1^2, rho*sd1*sd2, rho*sd1*sd2, sd2^2),ncol=2)
u.x <- c(5.6, 7)

es <- 0
es_w1 <- 0
es_w2 <- 0
w1 <- rep(NA,nrow(Y.pred))
for (i in 1:nrow(Y.pred)){
  w1[i] <- pmvnorm(upper=Y.pred[i,],mean=mu,sigma=sigma , keepAttr=FALSE)
}
w1 <- w1/sum(w1)
w2 <- (Y.pred[,1]>u.x[1])|(Y.pred[,2]>u.x[2])
w2 <- w2/sum(w2)
for (i in nrow(Y)){
  es <- es + es_sample(y = Y[i,], dat = t(Y.pred))
  es_w1 <- es_w1 + es_sample(y = Y[i,], dat = t(Y.pred), w = w1)
  es_w2 <- es_w2 + es_sample(y = Y[i,], dat = t(Y.pred), w = w2)
}


print(chiplot(Y))
chiplot(Y.pred)
post.pred(10, samples)

chi.value <- function(data, nq=100,qlim=NULL,conf = 0.95){
  data <- na.omit(data)
  n <- nrow(data)
  data <- cbind(rank(data[, 1])/(n + 1), rank(data[, 2])/(n + 
                                                            1))
  rowmax <- apply(data, 1, max)
  rowmin <- apply(data, 1, min)
  eps <- .Machine$double.eps^0.5
  qlim2 <- c(min(rowmax) + eps, max(rowmin) - eps)
  if (!is.null(qlim)) {
    if (qlim[1] < qlim2[1]) 
      stop("lower quantile limit is too low")
    if (qlim[2] > qlim2[2]) 
      stop("upper quantile limit is too high")
    if (qlim[1] > qlim[2]) 
      stop("lower quantile limit is less than upper quantile limit")
  }else qlim <- qlim2
  u <- seq(qlim[1], qlim[2], length = nq)
  cu <- cbaru <- numeric(nq)
  for (i in 1:nq) cu[i] <- mean(rowmax < u[i])
  for (i in 1:nq) cbaru[i] <- mean(rowmin > u[i])
  chiu <- 2 - log(cu)/log(u)
  chibaru <- (2 * log(1 - u))/log(cbaru) - 1
  cnst <- qnorm((1 + conf)/2)
  varchi <- ((1/log(u)^2 * 1)/cu^2 * cu * (1 - cu))/n
  varchi <- cnst * sqrt(varchi)
  varchibar <- (((4 * log(1 - u)^2)/(log(cbaru)^4 * cbaru^2)) * 
                  cbaru * (1 - cbaru))/n
  varchibar <- cnst * sqrt(varchibar)
  chiu <- cbind(chilow = chiu - varchi, chi = chiu, chiupp = chiu + 
                  varchi)
  chibaru <- cbind(chiblow = chibaru - varchibar, chib = chibaru, 
                   chibupp = chibaru + varchibar)
  chiulb <- 2 - log(pmax(2 * u - 1, 0))/log(u)
  chibarulb <- 2 * log(1 - u)/log(1 - 2 * u + pmax(2 * u -  1, 0)) - 1
  return(list(chiu,chibaru,u))
}



n.rep <- 200
chi.Y.pred <- matrix(NA,nrow=n.rep, ncol=100)
chib.Y.pred <- matrix(NA,nrow=n.rep, ncol=100)
for (i in 1:n.rep){
  dat <- post.pred(2000, samples)
  chi.Y.pred[i,] <-  chi.value(dat)[[1]][,'chi']
  chib.Y.pred[i,] <- chi.value(dat)[[2]][,'chib']
}

apply(chi.Y.pred,2,mean)

chi.Y <- chi.value(Y)[[1]]
chib.Y <- chi.value(Y)[[2]]
u <- chi.value(Y)[[3]]

chi.Y <- cbind(chi.Y,
               apply(chi.Y.pred,2,quantile,0.025),
               apply(chi.Y.pred,2,mean),
               apply(chi.Y.pred,2,quantile,0.975))
colnames(chi.Y) <- c('chilow','chi','chiupp','pred_chilow','pred_chi','pred_chiupp')
matplot(u, chi.Y, type = "l",lty=c(2,1,2,2,1,2), col=c(1,1,1,2,2,2), main='chi plot')
legend("bottomleft", colnames(chi.Y),col=c(1,1,1,2,2,2),cex=0.8, lty = c(2,1,2,2,1,2))

chib.Y <- cbind(chib.Y,
               apply(chib.Y.pred,2,quantile,0.025),
               apply(chib.Y.pred,2,mean),
               apply(chib.Y.pred,2,quantile,0.975))
colnames(chib.Y) <- c('chiblow','chib','chibupp','pred_chiblow','pred_chib','pred_chibupp')
matplot(u, chib.Y, type = "l",lty=c(2,1,2,2,1,2), col=c(1,1,1,2,2,2),main='chi bar plot')
legend("bottom", colnames(chib.Y),col=c(1,1,1,2,2,2),cex=0.8, lty = c(2,1,2,2,1,2))

