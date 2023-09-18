library(evd)
library(scoringRules)
library(mvtnorm)
library(tmvtnorm)
library(HDInterval)
library(posterior)
library(ggplot2)
library(RColorBrewer)
source("KRSW/RevExp_U_Functions.r")
source("KRSW/CommonFunctions.r")
source("KRSW/ModelDiagnosticsNewNames.r")

dir.out <- "/home/pgrad2/2448355h/My_PhD_Project/01_Output/Biv_Ext_Mix_Mod/UK_Temp"
dir.data <- "/home/pgrad2/2448355h/My_PhD_Project/00_Dataset/UK_Temp"
#################################Used functions#######################
post.pred <- function(n, samples, seed=1234){
  set.seed(seed)
  Y.pred <- matrix(NA,nrow=n, ncol=2)
  idx <- sample(nrow(samples),size=n, replace=TRUE)
  d <- 2
  for (i in 1:length(idx)){
    a <-  samples[idx[i], c('theta[1]','theta[2]')]
    sig <- samples[idx[i], c('theta[4]','theta[5]')]
    gamma <- samples[idx[i], c('theta[6]','theta[7]')]
    mu <- samples[idx[i], c('mu[1]','mu[2]')]
    sd1 <- samples[idx[i], 'sds[1]']
    sd2 <- samples[idx[i], 'sds[2]']
    corr.chol <- matrix(samples[idx[i],c('Ustar[1, 1]','Ustar[2, 1]',
                                         'Ustar[1, 2]','Ustar[2, 2]')],ncol=2)
    sigma <- diag(c(sd1,sd2))%*%t(corr.chol)%*%corr.chol%*%diag(c(sd1,sd2))
    thres <- samples[idx[i], c('thres[1]','thres[2]')]
    p <- pmvnorm( upper=thres, mean=mu, sigma=sigma, keepAttr = F)
    u <- runif(1)
    if (u<p){
      Y.pred[i,] <- rtmvnorm(1, mean=mu, sigma=sigma,upper=thres)
    }else{
      Y.tail <- sim.RevExpU.MGPD(n=1,d=d, a=a, beta=c(0,0), sig=sig, gamma=gamma, MGPD = T,std=T)$X
      Y.pred[i,] <- thres + Y.tail
    }
  }
  return(Y.pred)
}



# load(file=file.path(dir.out, filename='durham_london_0.8_0.99.RData'))
# load(file=file.path(dir.data, "durham_london.RData"))

# load(file=file.path(dir.out, filename='western-isles_argyll_0.8_0.99.RData'))
# load(file=file.path(dir.data, "western-isles_argyll.RData"))

load(file=file.path(dir.out, filename='east-sussex_oxfordshire_0.8_0.99_v1.RData'))
load(file=file.path(dir.data, "east-sussex_oxfordshire.RData"))
Y.fit <- -Y.all

# load(file=file.path(dir.out, filename='manchester_edinburgh_0.8_0.99.RData'))
# load(file=file.path(dir.data, "manchester_edinburgh.RData"))
# Y.fit <- Y.all


plot(Y.fit, main='Negative Residuals',xlab='Ringmer',ylab='Shirburn')

chain_out <- chain_res
para.name <- colnames(chain_out[[1]]$samples)


rhat.seq <- c()
ess.seq <- c()
for (name in para.name){
  post.sp <- cbind(chain_out[[1]]$samples[,name],
                   chain_out[[2]]$samples[,name],
                   chain_out[[3]]$samples[,name]
  )
  rhat.seq <- c(rhat.seq, rhat(post.sp))
  ess.seq <- c(ess.seq, ess_basic(post.sp))
}
convg.stat <- data.frame(para.name,rhat.seq,ess.seq )

samples.all <- rbind(chain_out[[1]]$samples[,],
                     chain_out[[2]]$samples[,],
                     chain_out[[3]]$samples[,])



apply(Y.all,2,quantile,0.83)

plot(samples.all[,'sds[2]'],type='l')
plot(samples.all[,'thres[1]'],type='l')
abline(h=quantile(Y.fit[,1], 0.8), lty=2,col='red')
plot(samples.all[,'thres[2]'],type='l')
abline(h=quantile(Y.fit[,2], 0.8), lty=2,col='red')
plot(samples.all[,'theta[6]'],type='l')

hist(samples.all[,'mu[2]'],breaks=100)
hist(samples.all[,'thres[2]'],breaks=100)
plot(density(samples.all[,'thres[2]']))

df.samples.all <- data.frame(samples.all)
colnames(df.samples.all) <- c("Ustar_1_1","Ustar_2_1","Ustar_1_2","Ustar_2_2","mu_1",'mu_2','s_1','s_2','a_1','a_2','lam','sigma_1',
                              'sigma_2','gamma_1','gamma_2','u_1','u_2')
hist.mypalette <- brewer.pal(9,"Greys")


res <- 500

for (col in colnames(df.samples.all)){
  filename <- paste("Hist_",col,".png",sep='')
  p <- ggplot(  df.samples.all, aes(x=.data[[col]])) +
    geom_histogram( bins=50, fill=hist.mypalette[5], color=hist.mypalette[7], alpha=0.9) +
    ggtitle(col) +
    theme(
      plot.title = element_text(size=15,hjust = 0.5),
      axis.title.x=element_blank(),
      axis.title.y=element_blank()
    )
  png(filename = file.path("/home/pgrad2/2448355h/My_PhD_Project/01_Output/Biv_Ext_Mix_Mod/Plots",filename), width = 6*res, height = 5*res, res=res)

  print(p)
  dev.off()
}


#seed durham london 123
#seed moray 1234
#seed western_isles 1
#east-sussex_oxfordshire 1, 123
#manchester_edinburgh 12 3000 prediction
#17

Y.fit.pred <- post.pred(50*nrow(Y.fit), samples.all,seed=17)
# Y.fit.pred <- post.pred(3000, samples.all,seed=12)
plot(Y.fit.pred, main='Predicted Residuals')

qqplot(Y.fit.pred[,1], Y.fit[,1],xlab='Predicted Y1', ylab="Empirical Y1")
abline(a=0,b=1)

emp_quant1 <- quantile(Y.fit[,1], probs = seq(0, 1, length.out = nrow(Y.fit)))
thy_quant1 <- qnorm(seq(0, 1, length.out = nrow(Y.fit)), mean = mean(Y.fit[,1]), sd = sd(Y.fit[,1]))
qqplot(thy_quant1, emp_quant1,xlab="Theoretical Y1", ylab="Empirical Y1")
abline(a=0,b=1)

qqplot(Y.fit.pred[,2],Y.fit[,2], xlab='Predicted Y2', ylab="Empirical Y2")
abline(a=0,b=1)

emp_quant2 <- quantile(Y.fit[,2], probs = seq(0, 1, length.out = nrow(Y.fit)))
thy_quant2 <- qnorm(seq(0, 1, length.out = nrow(Y.fit)), mean = mean(Y.fit[,2]), sd = sd(Y.fit[,2]))
qqplot(thy_quant2, emp_quant2,xlab="Theoretical Y2", ylab="Empirical Y2")
abline(a=0,b=1)


chi.thy <- function(a){
  a1 <- max(1/a)
  a2 <- min(1/a)
  chi <- 1 - ((1+1/a1)/(1+1/a2))^(1+a2)*a1/a2/(1+a1+a2)
  return(chi)
}
chiEmp1<-function(data,nq=25,qmin, qmax){
  n<-nrow(data)
  datatr<-(apply(data,2,rank) - 0.5)/n   
  qlim<-c(qmin,qmax)
  u <- seq(qlim[1], qlim[2], length = nq)
  cu<-sapply(c(1:nq), function(i) mean(apply(datatr,1,min) >= u[i]))
  res <- cbind(u,cu/(1-u))
  res[which(res[,2]>1),2] <- 1
  return(res)
}

# set.seed(1)
# Y <- rmvnorm(1000000,mean=c(0,0), sigma=sigma)
# nq <- 25
# qmin=0.5
# qmax <- 0.9999

etaEmp1<-function(data,nq=25,qmin, qmax){
  n<-nrow(data)
  datatr<-(apply(data,2,rank) - 0.5)/n   
  qlim<-c(qmin,qmax)
  u <- seq(qlim[1], qlim[2], length = nq)
  cu<-sapply(c(1:nq), function(i) mean(apply(datatr,1,min) >= u[i]))
  res <- cbind(u,2*log(1-u)/log(cu)-1)
  return(res)
}

post.dsp <- function(n.sp, samples, n.pred=2000, seed=1234){
  set.seed(seed)
  idx <- sample(nrow(samples),size=n.sp, replace=TRUE)
  d <- 2
  chi.seq <- rep(NA, length(idx))
  tau.seq <- rep(NA, length(idx))
  sCor.seq <- rep(NA<length(idx))
  chi.emp.seq <- rep(NA, length(idx))
  chiu.emp <- list()
  eta.emp.seq <- rep(NA, length(idx))
  etau.emp <- list()
  for (i in 1:length(idx)){
    a <-  samples[idx[i], c('theta[1]','theta[2]')]
    sig <- samples[idx[i], c('theta[4]','theta[5]')]
    gamma <- samples[idx[i], c('theta[6]','theta[7]')]
    mu <- samples[idx[i], c('mu[1]','mu[2]')]
    sd1 <- samples[idx[i], 'sds[1]']
    sd2 <- samples[idx[i], 'sds[2]']
    corr.chol <- matrix(samples[idx[i],c('Ustar[1, 1]','Ustar[2, 1]',
                                         'Ustar[1, 2]','Ustar[2, 2]')],ncol=2)
    sigma <- diag(c(sd1,sd2))%*%t(corr.chol)%*%corr.chol%*%diag(c(sd1,sd2))
    thres <- samples[idx[i], c('thres[1]','thres[2]')]
    p <- pmvnorm( upper=thres, mean=mu, sigma=sigma, keepAttr = F)
    u <- runif(1)
    
    Y.tail<-sim.RevExpU.MGPD(n=n.pred-floor(n.pred*p),d=d, a=a, beta=c(0,0), sig=sig, gamma=gamma, MGPD = T,std=T)
    
    # GP scale tail data combined with the bulk data
    Y.bulk <- rtmvnorm(floor(n.pred*p), mean=mu, sigma=sigma, upper=thres)
    
    # The name of the dataset should be Y for further WAIC calculation.
    Y <- rbind(Y.bulk, sweep(Y.tail$X,2,thres,"+"))
    chi.seq[i] <- chi.thy(a)
    chiu.emp[[i]] <- chiEmp1(data=Y, nq=50, qmin=0.5, qmax=0.99)
    chi.emp.seq[i] <- chiu.emp[[i]][50,2]
    etau.emp[[i]] <- etaEmp1(data=Y, nq=50, qmin=0.5, qmax=0.99)
    eta.emp.seq[i] <- etau.emp[[i]][50,2]
    tau.seq[i] <- cor(Y[,1], Y[,2], method='kendall')
    sCor.seq[i] <- cor(Y[,1], Y[,2], method = "spearman")
  }
  dps.mat <- cbind(tau.seq, sCor.seq, chi.seq, chi.emp.seq, eta.emp.seq)
  return(list("dps.mat"=dps.mat, 'chiu'=chiu.emp, 'etau'=etau.emp))
}


t1 <- Sys.time()
# dept.est <- post.dsp(3000, samples.all, n.pred=nrow(Y.fit))
dept.est <- post.dsp(3000, samples.all, n.pred=25*nrow(Y.fit))
t2 <- Sys.time()
print(t2-t1)
plot(density(dept.est[[1]][,3]))


emp.dsp <- function(n.sp, samples,seed=1234){
  set.seed(seed)
  d <- 2
  chi.seq <- rep(NA, n.sp)
  tau.seq <- rep(NA, n.sp)
  sCor.seq <- rep(NA, n.sp)
  chiu.emp <- list()
  eta.emp.seq <- rep(NA, length(n.sp))
  etau.emp <- list()
  i <- 1
  while (i <= n.sp){
    idx <- sample(nrow(samples),size=nrow(samples), replace=TRUE)
    Y.bs <- samples[idx,]
    chiu.emp[[i]] <- chiEmp1(data=Y.bs, nq=50, qmin=0.5, qmax=0.99)
    chi.seq[i] <- chiu.emp[[i]][50,2]
    if ( chi.seq[i]>1) break
    etau.emp[[i]] <- etaEmp1(data=Y.bs, nq=50, qmin=0.5, qmax=0.99)
    eta.emp.seq[i] <- etau.emp[[i]][50,2]
    tau.seq[i] <- cor(Y.bs[,1], Y.bs[,2], method='kendall')
    sCor.seq[i] <- cor(Y.bs[,1], Y.bs[,2], method = "spearman")
    i <- i + 1
  }
  dps.mat <- cbind(tau.seq, sCor.seq, chi.seq, eta.emp.seq)
  return(list("dps.mat"=dps.mat, 'chiu'=chiu.emp,  'etau'=etau.emp))
}

chiEmp(data=Y.fit,nq=50,qmin=0.50,qmax=0.99)
t3 <- Sys.time()
dept.emp.300 <- emp.dsp(300, Y.fit)
dept.pred.emp.300 <- emp.dsp(300, Y.fit.pred)
dept.emp.3000 <- emp.dsp(3000, Y.fit)
t4 <- Sys.time()
print(t4-t3)


# load(file=file.path(dir.out, 'Depd_Est_east-sussex_oxfordshire.RData'))

plot(density(dept.emp.3000[[1]][,1]), lwd=2, lty=2, main='kendall tau density',col=2)
lines(density(dept.est[[1]][,1]), lwd=2)
abline(v= cor(Y.fit[,1], Y.fit[,2], method='kendall'), col=2, lwd=2, lty=2)
legend("topleft", c('empirical','posterior'),col=c(2,1),cex=1, lty = c(2,1))

plot(density(dept.emp.3000[[1]][,2]), lwd=2, lty=2, main='spearman density',col=2)
lines(density(dept.est[[1]][,2]), lwd=2)
abline(v= cor(Y.fit[,1], Y.fit[,2], method='spearman'), col=2, lwd=2, lty=2)
legend("topleft", c('empirical','posterior'),col=c(2,1),cex=1, lty = c(2,1))

plot(density(dept.emp.3000[[1]][,3]), lwd=2, lty=2, main='chi density',col=2,ylim=c(0,14))
lines(density(dept.est[[1]][,4]), lwd=2)
abline(v=chiEmp(data=Y.fit,nq=50,qmin=0.50,qmax=0.99)[50,2], col=2, lwd=2, lty=2)
legend("topleft", c('empirical','posterior'),col=c(2,1),cex=1, lty = c(2,1))

plot(density(dept.emp.3000[[1]][,3]), lwd=2, lty=2, main='chi density',col=2,ylim=c(0,14))
lines(density(dept.est[[1]][,3]), lwd=2)
abline(v=chiEmp(data=Y.fit,nq=50,qmin=0.50,qmax=0.99)[50,2], col=2, lwd=2, lty=2)
legend("topleft", c('empirical','posterior'),col=c(2,1),cex=1, lty = c(2,1))


length(dept.emp.3000[[2]])

chi_cb <- function(dept.list){
  n <- length(dept.list)
  u <- dept.list[[1]][,1]
  chiu <- c()
  for (i in 1:n){
    chiu <- cbind(chiu,dept.list[[i]][,2])    
  }
  lb <- apply(chiu,1, quantile, 0.025)
  ub <- apply(chiu,1, quantile, 0.975)
  mean <- apply(chiu,1, mean)
  return(cbind(u,lb,ub,mean))
}

chi.cb.emp <- chi_cb(dept.emp.3000[[2]])
chi.cb.est <- chi_cb(dept.est[[2]])
chi.cb.thy <- quantile(dept.est[[1]][,3],c(0.025,0.975))
chi.cb.thy <- c(chi.cb.thy, mean(dept.est[[1]][,3]))

plot(chi.cb.emp[,1],chi.cb.emp[,4],type='l',ylim=c(0.2,0.8),lwd=1)
lines(chi.cb.emp[,1],chi.cb.emp[,2],type='l',lty=2,lwd=1)
lines(chi.cb.emp[,1],chi.cb.emp[,3],type='l',lty=2,lwd=1)

lines(chi.cb.est[,1],chi.cb.est[,4],type='l',lty=3, col='red', lwd=4)
# segments(x0=chi.cb.est[nrow(chi.cb.est),1], y0=chi.cb.thy[1], 
#          x1=chi.cb.est[nrow(chi.cb.est),1],y1=chi.cb.thy[2], lty=3,col='blue')
arrows(chi.cb.est[nrow(chi.cb.est),1], chi.cb.thy[1], 
       chi.cb.est[nrow(chi.cb.est),1], chi.cb.thy[2], angle=90, code=3, length=0.05,col='blue',lwd=4)
points(chi.cb.est[nrow(chi.cb.est),1], chi.cb.thy[3],col='blue',lwd=4)
lines(chi.cb.est[,1],chi.cb.est[,2],type='l',lty=4, col='red',lwd=2)
lines(chi.cb.est[,1],chi.cb.est[,3],type='l',lty=4, col='red',lwd=2)

chi.cb.emp <- data.frame(chi.cb.emp)
chi.cb.est <- data.frame(chi.cb.est)


df <- chi.cb.emp
df$lb2 <- chi.cb.est$lb
df$ub2 <- chi.cb.est$ub
df$mean2 <- chi.cb.est$mean

df.thy <- data.frame(u=0.99,lb=chi.cb.thy[1],ub=chi.cb.thy[2],mean=chi.cb.thy[3])

display.brewer.all(colorblindFriendly = TRUE) 

emp.mypalette <- brewer.pal(9,"Greys")
pred.mypalette <- brewer.pal(9,"Blues")

ggplot(df, aes(x=u, y=mean)) +
  # geom_point(aes(y=mean), color="blue") + 
  geom_line(aes(y=mean,linetype="dashed"), color=emp.mypalette[9]) +
  geom_ribbon(aes(ymin=lb, ymax=ub), fill=emp.mypalette[4], alpha=0.7) +
  scale_linetype_manual(values=c(solid="solid", dashed="dashed",twodash='twodash'))+

  geom_line(aes(y=mean2,linetype="twodash"), color=pred.mypalette[9]) +
  geom_ribbon(aes(ymin=lb2, ymax=ub2), fill=pred.mypalette[5], alpha=0.4) +
  
  geom_errorbar(data=df.thy, aes(x=u, ymin=lb, ymax=ub), width=0.01, color=pred.mypalette[7], size=1) +
  geom_point(data=df.thy, aes(x=u, y=mean), size=2, color=pred.mypalette[7])+
  
  labs(x='v',y='chi(v)') + 
  theme(legend.position = "none")


eta.cb.emp <- data.frame(chi_cb(dept.emp.3000[[3]]))
eta.cb.est <- data.frame(chi_cb(dept.est[[3]]))

df.eta <- eta.cb.emp
df.eta$lb2 <- eta.cb.est$lb
df.eta$ub2 <- eta.cb.est$ub
df.eta$mean2 <- eta.cb.est$mean

ggplot(df.eta, aes(x=u, y=mean)) +
  geom_line(aes(y=mean,linetype="dashed"), color=emp.mypalette[9]) +
  geom_ribbon(aes(ymin=lb, ymax=ub), fill=emp.mypalette[4], alpha=0.7) +
  scale_linetype_manual(values=c(solid="solid", dashed="dashed",twodash='twodash'))+
  
  geom_line(aes(y=mean2,linetype="twodash"), color=pred.mypalette[9]) +
  geom_ribbon(aes(ymin=lb2, ymax=ub2), fill=pred.mypalette[5], alpha=0.4) +
  
  labs(x='v',y='chi_bar(v)') + 
  theme(legend.position = "none")
# save(dept.emp.300,dept.emp.3000, dept.est, file=file.path(dir.out, 'Depd_Est_east-sussex_oxfordshire.RData'))
# save(dept.emp.300,dept.emp.3000, dept.est, file=file.path(dir.out, 'Depd_Manchester_Edinburgh.RData'))
# save(dept.emp.300,dept.emp.3000, dept.est, file=file.path(dir.out, 'Depd_Est_east-sussex_oxfordshire_seed123.RData'))
