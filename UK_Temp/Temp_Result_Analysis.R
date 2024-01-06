library(evd)
library(scoringRules)
library(mvtnorm)
library(tmvtnorm)
library(HDInterval)
library(posterior)
library(ggplot2)
library(RColorBrewer)
library(latex2exp)
source("KRSW/RevExp_U_Functions.r")
source("KRSW/CommonFunctions.r")
source("KRSW/ModelDiagnosticsNewNames.r")

dir.out <- "/home/pgrad2/2448355h/My_PhD_Project/01_Output/Biv_Ext_Mix_Mod/UK_Temp"
dir.data <- "/home/pgrad2/2448355h/My_PhD_Project/00_Dataset/UK_Temp"
#################################Used functions#######################
emp.mypalette <- brewer.pal(9,"Greys")
pred.mypalette <- brewer.pal(9,"Blues")
qq.palette <- brewer.pal(8,"Accent")
res <- 500

# Function to draw samples from the posterior predictive distribution
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

# Function to generate posterior predictive replicates
post.rep <-  function(n.rep, rep.size, samples, seed=1234){
  set.seed(seed)
  Y.rep <- list()
  idx <- sample(nrow(samples),size=n.rep, replace=TRUE)
  d <- 2
  for (i in 1:n.rep){
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
    

    Y.tail<-sim.RevExpU.MGPD(n=rep.size-floor(rep.size*p),d=d, a=a, beta=c(0,0), sig=sig, gamma=gamma, MGPD = T,std=T)
    
    # GP scale tail data combined with the bulk data
    Y.bulk <- rtmvnorm(floor(rep.size*p), mean=mu, sigma=sigma, upper=thres)
    
    # The name of the dataset should be Y for further WAIC calculation.
    Y <- rbind(Y.bulk, sweep(Y.tail$X,2,thres,"+"))
    Y.rep[[i]] <- Y
  }
  return(Y.rep)
}


load(file=file.path(dir.out, filename='east-sussex_oxfordshire_0.8_0.99_v2.RData'))
load(file=file.path(dir.data, "east-sussex_oxfordshire.RData"))
Y.fit <- -Y.all


# Scatter plot of the data
plot(Y.fit, main='Negative Residuals',xlab='Ringmer',ylab='Shirburn')

p <- ggplot(data.frame(Y.fit), aes(x=X1, y=X2)) + 
  geom_point(size=2,shape=19,col=pred.mypalette[7], alpha=0.9) +
  labs(x='Ringmer', y='Shirburn')+
  theme(axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        plot.title = element_text(size=18,hjust = 0.5))
png(filename = file.path(dir.out,"Temp_Residuals.png"), width = 6*res, height = 5*res, res=res)
print(p)
dev.off()

chain_out <- chain_res
para.name <- colnames(chain_out[[1]]$samples)

# Check the convergence of the MCMC
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

################################  histogram of the posteriors ################
df.samples.all <- data.frame(samples.all)
colnames(df.samples.all) <- c("Ustar_1_1","Ustar_2_1","Ustar_1_2","Ustar_2_2","mu_1",'mu_2','s_1','s_2','a_1','a_2','lam','sigma_1',
                              'sigma_2','gamma_1','gamma_2','u_1','u_2')

org.name <- colnames(df.samples.all)
latex.name <- c(r'($U[1,1]$)', r'($U[2,1]$)', r'($U[1,2]$)', r'($U[2,2]$)', r'($\mu_1$)', r'($\mu_2$)', r'($\s_1$)',
                r'($\s_2$)', r'($a_1$)', r'($a_2$)', r'($lam$)', r'($\sigma_1$)', r'($\sigma_2$)',
                r'($\gamma_1$)', r'($\gamma_2$)', r'($u_1$)', r'($u_2$)')
for (i in 1:length(org.name)){
  filename <- paste("Hist_",org.name[i],".png",sep='')
  p <- ggplot(  df.samples.all, aes(x=.data[[org.name[i]]])) +
    geom_histogram( bins=50, fill=pred.mypalette[5], color=pred.mypalette[7], alpha=0.9) +
    labs(x=TeX(latex.name[i])) + 
    theme(axis.text.x=element_text(size=15),
          axis.text.y=element_text(size=15),
          axis.title.x=element_text(size=15),
          axis.title.y=element_blank(),
          legend.position ='None')
  png(filename = file.path(dir.out,filename), width = 6*res, height = 5*res, res=res)
  print(p)
  dev.off()
}



##########################posterior predictive check#########################
n.rep <- 3000
Y.fit.rep <- post.rep(n.rep, nrow(Y.fit), samples.all, seed=1234)

#Quantile-Quantile plot of the emprical values and posterior predictive values for margin 1
probs <- c(seq(0.01, 0.9, length.out = 100),seq(0.91,0.999,length.out=50))
emp_quant1 <- quantile(Y.fit[,1], probs = probs)
thy_quant1 <- qnorm(probs, mean = mean(Y.fit[,1]), sd = sd(Y.fit[,1]))

pred_quant1.mat <- c() 
for (i in 1:n.rep){
  pred_quant1 <-  quantile( Y.fit.rep[[i]][,1], probs = probs)
  pred_quant1.mat <- cbind(pred_quant1, pred_quant1.mat)
}

lower1 <- apply(pred_quant1.mat,1,quantile,0.025)
upper1 <- apply(pred_quant1.mat,1,quantile,0.975)
mean1 <- apply(pred_quant1.mat,1,mean)


df.qq1 <- data.frame('Predicted_Quantile'= mean1,
                     'Empirical_Quantile'=emp_quant1,
                     'Gaussian_Quantile'=thy_quant1,
                     'lower_bound'=lower1,
                     'upper_bound'=upper1)


p <- ggplot(df.qq1, aes(x=emp_quant1, y=Predicted_Quantile,colour='BEMM')) +
  geom_abline(intercept = 0, slope = 1, color = 'black') +
  geom_ribbon(aes(ymin=lower_bound, ymax=upper_bound), colour=NA,  fill=pred.mypalette[5], alpha=0.25) +
  geom_point(aes(y=Gaussian_Quantile,colour='Gaussian')) +
  geom_point() +
  scale_color_manual(values = c("Gaussian" = qq.palette[3], "BEMM" = qq.palette[5]))+
  ggtitle("Ringmer") +
  xlab("Predicted Quantiles") +
  ylab("Empirical Quantiles") +
  theme(axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        plot.title = element_text(size=18,hjust = 0.5),
        legend.position=c(0, 1),
        legend.justification = c(0, 1),
        legend.title = element_blank(),
        legend.text = element_text(size=15),
        legend.key.size = unit(0.05, 'npc'),
        legend.key.width = unit(0.1, 'npc')) 
png(filename = file.path(dir.out, "QQ_plot_X1_with_CB.png"), width = 6*res, height = 5*res, res=res)
print(p)
dev.off()

#Quantile-Quantile plot of the emprical values and posterior predictive values for margin 2
emp_quant2 <- quantile(Y.fit[,2], probs = probs)
thy_quant2 <- qnorm(probs, mean = mean(Y.fit[,2]), sd = sd(Y.fit[,2]))

pred_quant2.mat <- c() 
for (i in 1:n.rep){
  pred_quant2 <-  quantile( Y.fit.rep[[i]][,2], probs = probs)
  pred_quant2.mat <- cbind(pred_quant2, pred_quant2.mat)
}

lower2 <- apply(pred_quant2.mat,1,quantile,0.025)
upper2 <- apply(pred_quant2.mat,1,quantile,0.975)
mean2 <- apply(pred_quant2.mat,1,mean)
df.qq2 <- data.frame('Predicted_Quantile'= mean2,
                     'Empirical_Quantile'=emp_quant2,
                     'Gaussian_Quantile'=thy_quant2,
                     'lower_bound'=lower2,
                     'upper_bound'=upper2)

p <- ggplot(df.qq2, aes(x=emp_quant2, y=Predicted_Quantile,colour='BEMM')) +
  geom_abline(intercept = 0, slope = 1, color = 'black') +
  geom_ribbon(aes(ymin=lower_bound, ymax=upper_bound), colour=NA,  fill=pred.mypalette[5], alpha=0.25) +
  geom_point(aes(y=Gaussian_Quantile,colour='Gaussian')) +
  geom_point() +
  scale_color_manual(values = c("Gaussian" = qq.palette[3], "BEMM" = qq.palette[5]))+
  ggtitle("Shirburn") +
  xlab("Empirical Quantiles") +
  ylab("Predicted Quantiles") +
  theme(axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        plot.title = element_text(size=18,hjust = 0.5),
        legend.position=c(0, 1),
        legend.justification = c(0, 1),
        legend.title = element_blank(),
        legend.text = element_text(size=15),
        legend.key.size = unit(0.05, 'npc'),
        legend.key.width = unit(0.1, 'npc')) 
png(filename = file.path(dir.out, "QQ_plot_X2_with_CB.png"), width = 6*res, height = 5*res, res=res)
print(p)
dev.off()


# Comparison of chi, chi bar plot and histogram of the Kendall's tau between empirical values and 
# posterior predictive values.
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

etaEmp1<-function(data,nq=25,qmin, qmax){
  n<-nrow(data)
  datatr<-(apply(data,2,rank) - 0.5)/n   
  qlim<-c(qmin,qmax)
  u <- seq(qlim[1], qlim[2], length = nq)
  cu<-sapply(c(1:nq), function(i) mean(apply(datatr,1,min) >= u[i]))
  res <- cbind(u,2*log(1-u)/log(cu)-1)
  return(res)
}



chi.seq <- rep(NA, length(n.rep))
tau.seq <- rep(NA, length(n.rep))
chiu.emp <- list()
etau.emp <- list()
set.seed(1234)
idx <- sample(nrow(samples.all),size=n.rep, replace=TRUE)
for (i in 1:n.rep){
  print(i)
  Y <- Y.fit.rep[[i]]
  a <-  samples.all[idx[i], c('theta[1]','theta[2]')]
  chi.seq[i] <- chi.thy(a)
  chiu.emp[[i]] <- chiEmp1(data=Y, nq=50, qmin=0.5, qmax=0.99)
  etau.emp[[i]] <- etaEmp1(data=Y, nq=50, qmin=0.5, qmax=0.99)
  tau.seq[i] <- cor(Y[,1], Y[,2], method='kendall')
}

calc_cb <- function(depd.list){
  n <- length(depd.list)
  u <- depd.list[[1]][,1]
  values <- c()
  for (i in 1:n){
    values <- cbind(values,depd.list[[i]][,2])    
  }
  lb <- apply(values,1, quantile, 0.025)
  ub <- apply(values,1, quantile, 0.975)
  mean <- apply(values,1, mean)
  return(data.frame('u'=u,'lb'=lb,'ub'=ub,'mean'=mean))
}

df.chi <- calc_cb(chiu.emp)
df.chi$emp <- chiEmp1(data=Y.fit, nq=50, qmin=0.5, qmax=0.99)[,2]
df.thy <- data.frame(u=0.99,lb=quantile(chi.seq,0.025),
                            ub=quantile(chi.seq,0.975),
                            mean=mean(chi.seq))

cols <- c("empirical"=qq.palette[6],
          'posterior mean'=pred.mypalette[6], 'cb'=pred.mypalette[4],
          'theorical'=pred.mypalette[9])
p <- ggplot(df.chi, aes(x=u, y=mean)) +
  geom_ribbon(aes(ymin=lb, ymax=ub,y=mean,fill='cb'),alpha=0.9) +
  geom_line(aes(y=mean,linetype="predicted",colour='posterior mean'),size=0.8) +
  geom_line(aes(y=emp,linetype="empirical",colour='empirical'),size=0.8) +

  geom_errorbar(data=df.thy, aes(x=u, ymin=lb, ymax=ub,color='theorical'), width=0.01,  size=1) +
  geom_point(data=df.thy, aes(x=u, y=mean,color='theorical'), size=2)+
  
  labs(x='r',y=TeX(r'($\chi(r)$)'), color  = "Mean", linetype = "Mean") + 
  scale_colour_manual("",values=cols)+
  # guide = 'none' removes the legend
  scale_fill_manual(name="Uncertainty",values=cols,guide="none")+
  scale_linetype_manual("",values=c('empirical'="dashed",'predicted'='twodash'),guide='none')+
  theme(axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        legend.position=c(0, 0),
        legend.justification = c(0, 0),
        legend.title = element_text(size=0),
        legend.text = element_text(size=15),
        legend.key.size = unit(0.05, 'npc'),
        legend.key.width = unit(0.1, 'npc'))

png(filename = file.path(dir.out, "chi_ppc.png"), width = 6*res, height = 5*res, res=res)
print(p)
dev.off()





df.eta <- calc_cb(etau.emp)
df.eta$emp <- etaEmp1(data=Y.fit, nq=50, qmin=0.5, qmax=0.99)[,2]


cols <- c("empirical"=qq.palette[6],
          'posterior mean'=pred.mypalette[6], 'cb'=pred.mypalette[4],
          'theorical'=pred.mypalette[9])
p <- ggplot(df.eta, aes(x=u, y=mean)) +
  geom_ribbon(aes(ymin=lb, ymax=ub,y=mean,fill='cb'),alpha=0.9) +
  geom_line(aes(y=mean,linetype="predicted",colour='posterior mean'),size=0.8) +
  geom_line(aes(y=emp,linetype="empirical",colour='empirical'),size=0.8) +

  labs(x='r',y=TeX(r'($\bar{\chi}(r)$)')) + 
  scale_colour_manual("",values=cols)+
  # guide = 'none' removes the legend
  scale_fill_manual(name="Uncertainty",values=cols,guide="none")+
  scale_linetype_manual("",values=c('empirical'="dashed",'predicted'='twodash'),guide='none')+
  theme(axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        legend.position = "none",
        legend.justification = c(0, 0),
        legend.title = element_text(size=0),
        legend.text = element_text(size=15),
        legend.key.size = unit(0.05, 'npc'),
        legend.key.width = unit(0.1, 'npc'))

png(filename = file.path(dir.out, "chibar_ppc.png"), width = 6*res, height = 5*res, res=res)
print(p)
dev.off()


df.tau <- data.frame(tau=tau.seq)
p <- ggplot(data=df.tau, aes(x=tau))+
  geom_histogram(bins=50, fill=pred.mypalette[5], color=pred.mypalette[7], alpha=0.9)+
  labs(x=TeX(r'($\tau$)'),y='')+
  geom_vline(xintercept = cor(Y.fit[,1], Y.fit[,2], method='kendall') , linetype="dashed", 
             color = qq.palette[6], size=0.8)+
  theme(axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        legend.position = "none",
        legend.justification = c(0, 0),
        legend.title = element_text(size=0),
        legend.text = element_text(size=15),
        legend.key.size = unit(0.05, 'npc'),
        legend.key.width = unit(0.1, 'npc'),
        plot.margin = unit(c(0.02, 0.05, 0, 0), 
                           "npc"))
png(filename = file.path(dir.out, "tau_ppc.png"), width = 6*res, height = 5*res, res=res)
print(p)
dev.off()


#################################Energy Score#############################
library(scoringRules)

t1 <- Sys.time()
energy_score<- function(Y.data, Y.pred){
  
  thres <- apply(Y.data, 2, quantile, 0.9) 
  
  mu <- mean(Y.data)
  sigma=cov(Y.data)
  weight_func_gauss <- function(x) prod(pnorm(x, mu, diag(sigma)))
  chain_func_gauss <- function(x){
    (x - mu)*pnorm(x, mu, diag(sigma)) + (diag(sigma)^2)*dnorm(x, mu, diag(sigma))
  }
  

  es.mat <- matrix(0, nrow=2,ncol=3)

  for (k in 1:nrow(Y.data)){
    es.mat[1,1] = es.mat[1,1]  + es_sample(y = Y.data[k,], dat = t(Y.pred))
  }
  es.mat[2,1] = es.mat[1,1]
  
  for (k in 1:nrow(Y.data)){
    es.mat[1,2] = es.mat[1,2]  + owes_sample(y = Y.data[k,], dat = t(Y.pred), a=thres)
    es.mat[2,2] = es.mat[2,2]  + twes_sample(y = Y.data[k,], dat = t(Y.pred), a=thres)
  }

  for (k in 1:nrow(Y.data)){
    es.mat[1,3] = es.mat[1,3]  + owes_sample(y = Y.data[k,], dat = t(Y.pred), weight_func = weight_func_gauss)
    es.mat[2,3] = es.mat[2,3]  + twes_sample(y = Y.data[k,], dat = t(Y.pred), chain_func = chain_func_gauss)
  }

  return(es.mat)
}
energy_score(Y.fit,Y.fit.rep[[1]])
t2 <- Sys.time()
print(t2-t1)


Y.pred.BEMM <- post.pred(10000, samples.all,seed=1234)
set.seed(1234)
Y.pred.gauss <- rmvnorm(10000, mean = colMeans(Y.fit), sigma = cov(Y.fit))

es_BEMM <- energy_score(Y.fit, Y.pred.BEMM)
es_Gauss <- energy_score(Y.fit, Y.pred.gauss)


