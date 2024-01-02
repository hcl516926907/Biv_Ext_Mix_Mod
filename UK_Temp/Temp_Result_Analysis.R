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

# load(file=file.path(dir.out, filename='durham_london_0.8_0.99.RData'))
# load(file=file.path(dir.data, "durham_london.RData"))

# load(file=file.path(dir.out, filename='western-isles_argyll_0.8_0.99.RData'))
# load(file=file.path(dir.data, "western-isles_argyll.RData"))

load(file=file.path(dir.out, filename='east-sussex_oxfordshire_0.8_0.99_v2.RData'))
load(file=file.path(dir.data, "east-sussex_oxfordshire.RData"))
Y.fit <- -Y.all

# load(file=file.path(dir.out, filename='manchester_edinburgh_0.8_0.99.RData'))
# load(file=file.path(dir.data, "manchester_edinburgh.RData"))
# Y.fit <- Y.all


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

a <- r'()'

################################# diagnostic based on 1 replicate data #############
#seed durham london 123
#seed moray 1234
#seed western_isles 1
#east-sussex_oxfordshire 1, 123
#manchester_edinburgh 12 3000 prediction
#17

Y.fit.pred <- post.pred(nrow(Y.fit), samples.all,seed=17)
# Y.fit.pred <- post.pred(3000, samples.all,seed=12)
plot(Y.fit.pred, main='Predicted Residuals')

qqplot(Y.fit.pred[,1], Y.fit[,1],xlab='Predicted Y1', ylab="Empirical Y1")
abline(a=0,b=1)

probs <- seq(10^-6, 1-10^-6, length.out = nrow(Y.fit))

emp_quant1 <- quantile(Y.fit[,1], probs = seq(10^-6, 1-10^-6, length.out = nrow(Y.fit)) )
thy_quant1 <- qnorm(probs, mean = mean(Y.fit[,1]), sd = sd(Y.fit[,1]))

qqplot(thy_quant1, emp_quant1,xlab="Theoretical Y1", ylab="Empirical Y1")
abline(a=0,b=1)
pred_quant1 <- quantile(Y.fit.pred[,1], probs = probs)

qqplot(Y.fit.pred[,2],Y.fit[,2], xlab='Predicted Y2', ylab="Empirical Y2")
abline(a=0,b=1)

emp_quant2 <- quantile(Y.fit[,2], probs)
thy_quant2 <- qnorm(probs , mean = mean(Y.fit[,2]), sd = sd(Y.fit[,2]))
pred_quant2 <- quantile(Y.fit.pred[,2], probs = probs)
qqplot(thy_quant2, emp_quant2,xlab="Theoretical Y2", ylab="Empirical Y2")
abline(a=0,b=1)


df1.1 <- data.frame(
  theoretical = thy_quant1,
  sample = emp_quant1,
  dataset = "Gaussian"
)
df1.2 <- data.frame(
  theoretical = pred_quant1,
  sample = emp_quant1,
  dataset = "BEMM"
)


df2.1 <- data.frame(
  theoretical = thy_quant2,
  sample = emp_quant2,
  dataset = "Gaussian"
)
df2.2 <- data.frame(
  theoretical = pred_quant2,
  sample = emp_quant2,
  dataset = "BEMM"
)
combined_df.1 <- rbind(df1.1, df1.2)
combined_df.2 <- rbind(df2.1, df2.2)



ggplot(combined_df.1, aes(x = theoretical, y = sample, color = dataset)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = 'black') +
  ggtitle("Ringmer") +
  xlab("Predicted Quantiles") +
  ylab("Empirical Quantiles") +
  scale_color_manual(values = c("Gaussian" = qq.palette[3], "BEMM" = qq.palette[5]))+
  # xlim(-7, 12)+
  # ylim(-7, 12)+
  theme(axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        plot.title = element_text(size=18,hjust = 0.5),
        legend.position ='None')

# png(filename = file.path(dir.out, "QQ_plot_X1.png"), width = 6*res, height = 5*res, res=res)
# print(p)
# dev.off()

ggplot(combined_df.2, aes(x = theoretical, y = sample, color = dataset)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = 'black') +
  ggtitle("Shirburn") +
  xlab("Predicted Quantiles") +
  ylab("Empirical Quantiles") +
  scale_color_manual(values = c("Gaussian" = qq.palette[3], "BEMM" = qq.palette[5]))+
  # xlim(-7, 12)+
  # ylim(-7, 12)+
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

# png(filename = file.path(dir.out, "QQ_plot_X2.png"), width = 6*res, height = 5*res, res=res)
# print(p)
# dev.off()

##########################posterior predictive check#########################
n.rep <- 3000
Y.fit.rep <- post.rep(n.rep, nrow(Y.fit), samples.all, seed=1234)


# probs <- c(seq(0.01, 0.99, length.out = 100),0.999)
# probs <- c(seq(0.01, 0.99, length.out = 100),0.995)
# probs <- c(seq(0.01, 0.99, length.out = 100))
probs <- c(seq(0.01, 0.9, length.out = 100),seq(0.91,0.999,length.out=50))
# probs <- seq(10^-6, 1-10^-6, length.out = nrow(Y.fit))
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

qqplot(mean1,emp_quant1)


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
##########################################################
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
dept.est <- post.dsp(3000, samples.all, n.pred=nrow(Y.fit))
# dept.est <- post.dsp(3000, samples.all, n.pred=25*nrow(Y.fit))
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

# save(dept.emp.300,dept.emp.3000, dept.est, file=file.path(dir.out, 'Depd_Est_east-sussex_oxfordshire.RData'))
# save(dept.emp.300,dept.emp.3000, dept.est, file=file.path(dir.out, 'Depd_Manchester_Edinburgh.RData'))
# save(dept.emp.300,dept.emp.3000, dept.est, file=file.path(dir.out, 'Depd_Est_east-sussex_oxfordshire_seed123.RData'))
# save(dept.emp.300,dept.emp.3000, dept.est, file=file.path(dir.out, 'Depd_Est_east-sussex_oxfordshire_seed17.RData'))
# load(file=file.path(dir.out, 'Depd_Est_east-sussex_oxfordshire_seed17.RData'))

# plot(density(dept.emp.3000[[1]][,1]), lwd=2, lty=2, main='kendall tau density',col=2)
# lines(density(dept.est[[1]][,1]), lwd=2)
# abline(v= cor(Y.fit[,1], Y.fit[,2], method='kendall'), col=2, lwd=2, lty=2)
# legend("topleft", c('empirical','posterior'),col=c(2,1),cex=1, lty = c(2,1))
# 
# plot(density(dept.emp.3000[[1]][,2]), lwd=2, lty=2, main='spearman density',col=2)
# lines(density(dept.est[[1]][,2]), lwd=2)
# abline(v= cor(Y.fit[,1], Y.fit[,2], method='spearman'), col=2, lwd=2, lty=2)
# legend("topleft", c('empirical','posterior'),col=c(2,1),cex=1, lty = c(2,1))
# 
# plot(density(dept.emp.3000[[1]][,3]), lwd=2, lty=2, main='chi density',col=2,ylim=c(0,14))
# lines(density(dept.est[[1]][,4]), lwd=2)
# abline(v=chiEmp(data=Y.fit,nq=50,qmin=0.50,qmax=0.99)[50,2], col=2, lwd=2, lty=2)
# legend("topleft", c('empirical','posterior'),col=c(2,1),cex=1, lty = c(2,1))
# 
# plot(density(dept.emp.3000[[1]][,3]), lwd=2, lty=2, main='chi density',col=2,ylim=c(0,14))
# lines(density(dept.est[[1]][,3]), lwd=2)
# abline(v=chiEmp(data=Y.fit,nq=50,qmin=0.50,qmax=0.99)[50,2], col=2, lwd=2, lty=2)
# legend("topleft", c('empirical','posterior'),col=c(2,1),cex=1, lty = c(2,1))


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



cols <- c("emp_chi"=emp.mypalette[9],'empirical'=emp.mypalette[4],
          'pred_chi'=pred.mypalette[9], 'predicted'=pred.mypalette[5],
          'theorical'=pred.mypalette[7])
p <- ggplot(df, aes(x=u, y=mean)) +
  geom_line(aes(y=mean,linetype="empirical",colour='emp_chi')) +
  
  geom_ribbon(aes(ymin=lb, ymax=ub,y=mean,fill='empirical'),alpha=0.7) +

  geom_line(aes(y=mean2,linetype="predicted",color='pred_chi') ) +
  geom_ribbon(aes(ymin=lb2, ymax=ub2, fill='predicted'),alpha=0.4 ) +
  
  geom_errorbar(data=df.thy, aes(x=u, ymin=lb, ymax=ub,color='theorical'), width=0.01,  size=1) +
  geom_point(data=df.thy, aes(x=u, y=mean,color='theorical'), size=2)+
  
  labs(x='r',y=TeX(r'($\chi(r)$)')) + 
  scale_colour_manual(name="Mean",values=cols,guide = "none")+
  scale_fill_manual(name="Uncertainty",values=cols)+
  scale_linetype_manual(name='Mean',values=c('empirical'="dashed",'predicted'='twodash'))+
  theme(axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        legend.position=c(0, 0),
        legend.justification = c(0, 0),
        legend.title = element_text(size=15),
        legend.text = element_text(size=15),
        legend.key.size = unit(0.05, 'npc'),
        legend.key.width = unit(0.1, 'npc'))

png(filename = file.path(dir.out, "chi.png"), width = 6*res, height = 5*res, res=res)
print(p)
dev.off()

eta.cb.emp <- data.frame(chi_cb(dept.emp.3000[[3]]))
eta.cb.est <- data.frame(chi_cb(dept.est[[3]]))

df.eta <- eta.cb.emp
df.eta$lb2 <- eta.cb.est$lb
df.eta$ub2 <- eta.cb.est$ub
df.eta$mean2 <- eta.cb.est$mean

p <- ggplot(df.eta, aes(x=u, y=mean)) +
  geom_line(aes(y=mean,linetype="dashed"), color=emp.mypalette[9]) +
  geom_ribbon(aes(ymin=lb, ymax=ub), fill=emp.mypalette[4], alpha=0.7) +
  scale_linetype_manual(values=c(solid="solid", dashed="dashed",twodash='twodash'))+
  
  geom_line(aes(y=mean2,linetype="twodash"), color=pred.mypalette[9]) +
  geom_ribbon(aes(ymin=lb2, ymax=ub2), fill=pred.mypalette[5], alpha=0.4) +
  
  labs(x='r',y=TeX(r'($\bar{\chi}(r)$)')) + 
  theme(axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        legend.position = "none",
        legend.title = element_text(size=15),
        legend.text = element_text(size=15),
        legend.key.size = unit(0.05, 'npc'),
        legend.key.width = unit(0.1, 'npc'))

png(filename = file.path(dir.out, "chibar.png"), width = 6*res, height = 5*res, res=res)
print(p)
dev.off()



df.tau1 <- data.frame(tau=dept.emp.3000[[1]][,1],dataset='Empirical')
df.tau2 <- data.frame(tau=dept.est[[1]][,1],dataset='Predicted')
df.tau <- rbind(df.tau1,df.tau2)

p <- ggplot(df.tau, aes(x=dataset,y=tau,colour=dataset))+
  geom_boxplot() +
  labs(x='model',y=TeX(r'($\tau$)'))+
  scale_colour_manual(values=c('Empirical'=emp.mypalette[4],'Predicted'=pred.mypalette[5]))+
  theme(axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=15),
        legend.title = element_blank(),
        legend.text = element_text(size=15))
png(filename = file.path(dir.out, "tau.png"), width = 6*res, height = 5*res, res=res)
print(p)
dev.off()

###########################################CRPS###############################

Y.tail <- sim.RevExpU.MGPD(n=1000000,d=2, a=c(0.5,1), beta=c(0,0), sig=c(1,1), gamma=c(-1,0.9), MGPD = T,std=T)$X
colMeans(Y.tail)
library(cubature)
intCrps <- function(d,igd,l,u){
  
  converged <- F
  reltol <- 0.01
  
  if(sum(l<u)==d){
    while(converged == F){
      hide <- capture.output(int <- vegas(igd,lowerLimit=l,upperLimit=u,relTol=reltol))
      if(converged == F){
        reltol <- reltol + 0.01
      }
      else{
        converged <- T
      }
      
      if(reltol > 0.1){
        converged <- T
      }
    }
  }else{
    int<-NULL
    int$integral <- 0
    int$returnCode <- 1
  }
  
  return(int)
  
}

emcdf <- function(x,u){
  x <- as.matrix(x)
  n <- nrow(x)
  cnt <- apply(t(x)<=u,2,all )
  return(sum(cnt)/(n+1))
}

X <-cbind( rnorm(2000),rnorm(2000))
y <- c(1,2)
emcdf(x,c(4,4))

crps <- function(X,y){
  
  # Bounds and dimensions
  X <- as.matrix(X)
  y <- c(y)
  l <- apply(X,2,min)
  u <- apply(X,2,max)
  d <- dim(X)[2]
  
  # Lower integrand
  igdL <- function(u){
    emcdf(X,u)^2
  }
  
  # Upper integrand
  igdU <- function(u){
    (emcdf(X,u)-1)^2
  }
  
  # Lower integral
  intL <- intCrps(d,igdL,l,y)
  
  # Upper integral
  intU <- intCrps(d,igdU,y,u)
  
  # CRPS result
  if((intL$returnCode== 0) && (intU$returnCode==0)){
    return(intL$integral + intU$integral)
  }else{
    return(NA)
  }
}

t1 <- Sys.time()
crps(X,c(2,2))
print(Sys.time()-t1)

crps1 <- function(X,y){
  
  # Bounds and dimensions
  X <- as.matrix(X)
  y <- c(y)
  l <- apply(X,2,min)
  u <- apply(X,2,max)
  d <- dim(X)[2]
  
  # Lower integrand
  igd <- function(u){
    emcdf(X,u)^2
  }
  
  # Upper integrand
  igdU <- function(u){
    (emcdf(X,u)-1)^2
  }
  
  # Lower integral
  intL <- intCrps(d,igdL,l,y)
  
  # Upper integral
  intU <- intCrps(d,igdU,y,u)
  
  # CRPS result
  if((intL$returnCode== 0) && (intU$returnCode==0)){
    return(intL$integral + intU$integral)
  }else{
    return(NA)
  }
}

crps(x,c(2,2))



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



Y.pred.BEMM_2k <- post.pred(2000, samples.all,seed=1234)
set.seed(1234)
Y.pred.gauss_2k <- rmvnorm(2000, mean = colMeans(Y.fit), sigma = cov(Y.fit))


es_BEMM_2k_0.9 <- energy_score(Y.fit, Y.pred.BEMM_2k)
es_Gauss_2k_0.9 <- energy_score(Y.fit, Y.pred.gauss_2k)  

