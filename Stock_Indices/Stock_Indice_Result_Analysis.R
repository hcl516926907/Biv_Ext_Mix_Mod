library(scoringRules)
library(mvtnorm)
library(tmvtnorm)
library(posterior)
library(ggplot2)
library(RColorBrewer)
library(latex2exp)
library(extraDistr)
library(ggplot2)
library(latex2exp)
library(MASS)
library(copula)
source("Simulation/RevExp_U_Functions.r")
source("Simulation/CommonFunctions.r")
source("Simulation/Gumbel_U_Functions.r")

dir.out <- "/home/pgrad2/2448355h/My_PhD_Project/01_Output/Biv_Ext_Mix_Mod/Stock_Market_Index"
dir.data <- "/home/pgrad2/2448355h/My_PhD_Project/00_Dataset/UK_Temp"
#################################Used functions#######################
emp.mypalette <- brewer.pal(9,"Greys")
pred.mypalette <- brewer.pal(9,"Blues")
qq.palette <- brewer.pal(8,"Accent")
res <- 300

# load(file=file.path("/home/pgrad2/2448355h/My_PhD_Project/01_Output/Biv_Ext_Mix_Mod/Simulation", 'daily_index_res_Gumbel_thres_0.85.RData') )
load(file=file.path("/home/pgrad2/2448355h/My_PhD_Project/01_Output/Biv_Ext_Mix_Mod/Simulation", 'daily_index_res_Gumbel_thres_0.01_sd_0.1.RData') )



p <- ggplot(data.frame(Y.fit), aes(x=log.cac.return, y=log.dax.return)) + 
  geom_point(size=2,shape=19,col=pred.mypalette[7], alpha=0.9) +
  labs(x='CAC 40', y='DAX')+
  theme(axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        plot.title = element_text(size=18,hjust = 0.5))
print(p)
png(filename = file.path(dir.out,"Log_Return_Residuals.png"), width = 6*res, height = 5*res, res=res)
print(p)
dev.off()


# Function to generate posterior predictive replicates
post.rep <-  function(n.rep, rep.size, samples, seed=1234, lbound.tail=rep(-Inf,2), ubound.tail=rep(Inf,2)){
  set.seed(seed)
  Y.rep <- list()
  idx <- sample(nrow(samples),size=n.rep, replace=TRUE)
  d <- 2
  for (i in 1:n.rep){
    a <-  rep(samples[idx[i], c('theta[1]')])
    sig <- samples[idx[i], c('theta[2]','theta[3]')]
    gamma <- samples[idx[i], c('theta[4]','theta[5]')]
    
    params.bulk <- samples[idx[i], c('params.bulk[1]','params.bulk[2]','params.bulk[3]',
                                     'params.bulk[4]','params.bulk[5]','params.bulk[6]',
                                     'params.bulk[7]')]
    
    cop <- gumbelCopula(param = params.bulk[1], dim = 2)
    
    bulk.dist <-  mvdc(cop, margins = c('lst','lst'), 
                                  
                                  paramMargins=list(list(df=params.bulk[2], mu=params.bulk[3], sigma=params.bulk[4]),
                                                    
                                                    list(df=params.bulk[5], mu=params.bulk[6], sigma=params.bulk[7])))
    

    # thres <- samples[idx[i], c('thres[1]','thres[2]')]
    thres <- c(0.85,0.85)
    
    p <- pMvdc(thres, bulk.dist)
    
    Y.bulk <- rMvdc(floor(rep.size*p), bulk.dist)
    Y.bulk <- Y.bulk[Y.bulk[,1]<thres[1] & Y.bulk[,2]<thres[2],]
    
    
    Y.tail<-sim.GumbelU.MGPD(n=rep.size-floor(rep.size*p),d=2, a=a, beta=c(0,0), sig=sig, gamma=gamma, MGPD = T,std=F)
    
    bound.cond.1 <- Y.tail[,1] > (lbound.tail-thres)[1] & Y.tail[,1] < (ubound.tail-thres)[1]
    bound.cond.2 <- Y.tail[,2] > (lbound.tail-thres)[2] & Y.tail[,2] < (ubound.tail-thres)[2]
    Y.tail <- Y.tail[bound.cond.1&bound.cond.2,]
    
    # The name of the dataset should be Y for further WAIC calculation.
    Y <- rbind(Y.bulk, sweep(Y.tail,2,thres,"+"))
    Y.rep[[i]] <- Y
  }
  return(Y.rep)
}


n.rep <- 1000
Y.fit
post.burnin <- 2000:4999

plot(density(results$samples[post.burnin, 'thres[1]' ]))
plot(results$samples[post.burnin, 'thres[2]' ],type='l')
plot(results$samples[1:4999, 'thres[1]' ],type='l')

theta=colMeans(results$samples[post.burnin,c('theta[1]','theta[2]','theta[3]','theta[4]'
                                             ,'theta[5]')])
thres=colMeans(results$samples[post.burnin,c('thres[1]','thres[2]')])
params.bulk=colMeans(results$samples[post.burnin,c('params.bulk[1]','params.bulk[2]','params.bulk[3]',
                                                   'params.bulk[4]','params.bulk[5]',
                                                   'params.bulk[6]','params.bulk[7]')])

Y.fit.rep <- post.rep(n.rep=n.rep, rep.size=nrow(Y.fit), samples=results$samples[post.burnin, ], seed=1234,
                      lbound.tail=rep(-Inf,2))


probs <- c(seq(0.01, 0.9, length.out = 100),seq(0.91,0.999,length.out=50))
emp_quant1 <- quantile(Y.fit[,1], probs = probs)

est1 <- fitdistr(Y.fit[,1], "t", start = list(m=0,s=1, df=6), lower=c(-1, 0.001,1))

thy_quant1 <- qlst(probs, df=est1$estimate[3], mu = est1$estimate[1], sigma = est1$estimate[2])

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
                     'Stuent_t_Quantile'=thy_quant1,
                     'lower_bound'=lower1,
                     'upper_bound'=upper1)


p <- ggplot(df.qq1, aes(x=emp_quant1, y=Predicted_Quantile,colour='BEMM')) +
  geom_abline(intercept = 0, slope = 1, color = 'black') +
  geom_ribbon(aes(ymin=lower_bound, ymax=upper_bound), colour=NA,  fill=pred.mypalette[5], alpha=0.25) +
  geom_point(aes(y=Stuent_t_Quantile,colour='Student_t')) +
  geom_point() +
  scale_color_manual(values = c("Student_t" = qq.palette[3], "BEMM" = qq.palette[5]))+
  ggtitle("Residuals of log return of CAC 40 ") +
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
print(p)
png(filename = file.path(dir.out, "QQ_plot_X1_with_CB.png"), width = 6*res, height = 5*res, res=res)
print(p)
dev.off()


emp_quant2 <- quantile(Y.fit[,2], probs = probs)

est2 <- fitdistr(Y.fit[,2], "t", start = list(m=0,s=1, df=3), lower=c(-1, 0.001,1))

thy_quant2 <- qlst(probs, df=est2$estimate[3], mu = est2$estimate[1], sigma = est2$estimate[2])


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
                     'Student_t_Quantile'=thy_quant2,
                     'lower_bound'=lower2,
                     'upper_bound'=upper2)

p <- ggplot(df.qq2, aes(x=emp_quant2, y=Predicted_Quantile,colour='BEMM')) +
  geom_abline(intercept = 0, slope = 1, color = 'black') +
  geom_ribbon(aes(ymin=lower_bound, ymax=upper_bound), colour=NA,  fill=pred.mypalette[5], alpha=0.25) +
  geom_point(aes(y=Student_t_Quantile,colour='Student_t')) +
  geom_point() +
  scale_color_manual(values = c("Student_t" = qq.palette[3], "BEMM" = qq.palette[5]))+
  ggtitle("Residuals of log return of DAX") +
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
print(p)

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


samples.all <-  results$samples[post.burnin, ]
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
          'posterior mean'=pred.mypalette[6], 'cb'=pred.mypalette[4])
p <- ggplot(df.chi, aes(x=u, y=mean)) +
  geom_ribbon(aes(ymin=lb, ymax=ub,y=mean,fill='cb'),alpha=0.9) +
  geom_line(aes(y=mean,linetype="predicted",colour='posterior mean'),size=0.8) +
  geom_line(aes(y=emp,linetype="empirical",colour='empirical'),size=0.8) +
  
  # geom_errorbar(data=df.thy, aes(x=u, ymin=lb, ymax=ub,color='theorical'), width=0.01,  size=1) +
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
print(p)

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

print(p)
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

print(p)
png(filename = file.path(dir.out, "tau_ppc.png"), width = 6*res, height = 5*res, res=res)
print(p)
dev.off()