library(scoringRules)
library(mvtnorm)
library(tmvtnorm)
library(posterior)
library(ggplot2)
library(RColorBrewer)
library(latex2exp)
library(cowplot)
source("Simulation/RevExp_U_Functions.r")
source("Simulation/CommonFunctions.r")

dir.out <- "/home/pgrad2/2448355h/My_PhD_Project/01_Output/Biv_Ext_Mix_Mod/UK_Temp"
dir.data <- "/home/pgrad2/2448355h/My_PhD_Project/00_Dataset/UK_Temp"
#################################Used functions#######################
emp.mypalette <- brewer.pal(9,"Greys")
pred.mypalette <- brewer.pal(9,"Blues")
qq.palette <- brewer.pal(8,"Accent")
res <- 300

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


load(file=file.path(dir.out, filename='east-sussex_oxfordshire_0.8_0.99_v3.RData'))
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
colnames(df.samples.all) <- c("Ustar_1_1","Ustar_2_1","Ustar_1_2","Ustar_2_2",
                              "delta_thres_1",
                              "delta_thres_2",
                              "mu_1",'mu_2','s_1','s_2','a_1','a_2','lam','sigma_1',
                              'sigma_2','gamma_1','gamma_2','u_1','u_2')

org.name <- colnames(df.samples.all)
latex.name <- c(r'($L[1,1]$)', r'($L[2,1]$)', r'($L[1,2]$)', r'($L[2,2]$)',
                r'($\delta_u1$)',
                r'($\delta_u2$)',
                r'($\mu_1$)', r'($\mu_2$)', r'($\s_1$)',
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
  # png(filename = file.path(dir.out,filename), width = 6*res, height = 5*res, res=res)
  print(p)
  # dev.off()
}



df.chain <- lapply(1:3, function(k) {
  df <- data.frame(chain_out[[k]]$samples[,])
  colnames(df) <- c("Ustar_1_1","Ustar_2_1","Ustar_1_2","Ustar_2_2","mu_1",'mu_2','s_1','s_2','a_1','a_2','lam','sigma_1',
                    'sigma_2','gamma_1','gamma_2','u_1','u_2')
  df
})

# --- compute common x-limits per parameter across chains ---
xlim_u1 <- range(unlist(lapply(df.chain, \(d) d$u_1)), na.rm = TRUE)
xlim_u2 <- range(unlist(lapply(df.chain, \(d) d$u_2)), na.rm = TRUE)

# --- Helper: same style + centered title + shared x-lims ---
make_hist <- function(df, var, latex_lab, title_text, xlim_common) {
  ggplot(df, aes(x = .data[[var]])) +
    geom_histogram(
      bins = 50,
      fill = pred.mypalette[5],
      color = pred.mypalette[7],
      alpha = 0.9
    ) +
    labs(x = TeX(latex_lab), title = title_text) +
    coord_cartesian(xlim = xlim_common) +
    theme(
      axis.text.x  = element_text(size = 15),
      axis.text.y  = element_text(size = 15),
      axis.title.x = element_text(size = 15),
      axis.title.y = element_blank(),
      legend.position = "none",
      plot.title = element_text(size = 15, hjust = 0.5)  # center title
    )
}

# --- Make 6 plots (3 for u_1, 3 for u_2) ---
p_u1 <- lapply(1:3, function(k)
  make_hist(df.chain[[k]], "u_1", r'($u_1$)', paste0("Chain ", k), xlim_u1)
)

p_u2 <- lapply(1:3, function(k)
  make_hist(df.chain[[k]], "u_2", r'($u_2$)', paste0("Chain ", k), xlim_u2)
)

# --- Arrange 2 × 3 ---
final_plot <- plot_grid(
  p_u1[[1]], p_u1[[2]], p_u1[[3]],
  p_u2[[1]], p_u2[[2]], p_u2[[3]],
  nrow = 2
)

ggsave(
  filename = file.path(dir.out, "hist_u1_u2_grid.png"),
  plot     = final_plot,
  width    = 12, height = 7, units = "in",
  dpi      = 300
)


################################  traceplot tcolours()################################  traceplot the posteriors ################
samples.all <- samples.all[, !colnames(samples.all) %in% c("delta_thres[1]", "delta_thres[2]")]
colnames(samples.all) <- c("Ustar_1_1","Ustar_2_1","Ustar_1_2","Ustar_2_2", "mu_1",'mu_2','s_1','s_2','a_1','a_2','lam','sigma_1',
                              'sigma_2','gamma_1','gamma_2','u_1','u_2')


org.name <- colnames(samples.all)
latex.name <- c(r'($L[1,1]$)', r'($L[2,1]$)', r'($L[1,2]$)', r'($L[2,2]$)', r'($\mu_1$)', r'($\mu_2$)', r'($\s_1$)',
                r'($\s_2$)', r'($a_1$)', r'($a_2$)', r'($lam$)', r'($\sigma_1$)', r'($\sigma_2$)',
                r'($\gamma_1$)', r'($\gamma_2$)', r'($u_1$)', r'($u_2$)')
for (i in 1:(length(org.name))){
  chain1 <- data.frame(iteration = 1:1000, value = samples.all[1:1000, org.name[i]], chain = 'Chain 1')
  chain2 <- data.frame(iteration = 1:1000, value = samples.all[1001:2000, org.name[i]], chain = 'Chain 2')
  chain3 <- data.frame(iteration = 1:1000, value = samples.all[2001:3000, org.name[i]], chain = 'Chain 3')
  mcmc_data <- rbind(chain1, chain2, chain3)
  # Generating the traceplot
  p <- ggplot(mcmc_data, aes(x = iteration, y = value, colour=chain)) +
    geom_line() +
    theme_minimal() +
    scale_color_manual(values = c("Chain 1" = qq.palette[1], "Chain 2" = qq.palette[3], "Chain 3" = qq.palette[5]))+
    ggtitle(TeX(latex.name[i]))+
    labs(x = "Iteration", y = "Parameter Value") +
    theme(axis.text.x=element_text(size=15),
          axis.text.y=element_text(size=15),
          axis.title.x=element_text(size=15),
          axis.title.y=element_blank(),
          plot.title = element_text(size=18,hjust = 0.5),
          legend.position ='None')
  filename <- paste("Traceplot_",org.name[i],".png",sep='')
  # png(filename = file.path(dir.out,filename), width = 6*res, height = 5*res, res=res)
  print(p)
  # dev.off()
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



#######################################
library(hexbin)

## --- Helper: safely fetch a scalar column from samples by its nimble name
get_col <- function(samples, name) {
  if (!name %in% colnames(samples)) stop("Missing column: ", name)
  samples[, name]
}

## --- Helper: reconstruct 2x2 upper-triangular Cholesky U
## Option A: if U[1,1], U[1,2], U[2,2] are monitored, use them
get_U_from_samples <- function(samples, s) {
  nms <- colnames(samples)
  if (all(c("U[1, 1]", "U[1, 2]", "U[2, 2]") %in% nms)) {
    U <- matrix(0, 2, 2)
    U[1,1] <- samples[s, "U[1, 1]"]
    U[1,2] <- samples[s, "U[1, 2]"]
    U[2,2] <- samples[s, "U[2, 2]"]
    return(U)
  }
  
  ## Option B: reconstruct from Ustar and sds if those are monitored
  if (all(c("sds[1]", "sds[2]", "Ustar[1, 1]", "Ustar[1, 2]", "Ustar[2, 2]") %in% nms)) {
    sds <- c(samples[s, "sds[1]"], samples[s, "sds[2]"])
    Ustar <- matrix(0, 2, 2)
    Ustar[1,1] <- samples[s, "Ustar[1, 1]"]
    Ustar[1,2] <- samples[s, "Ustar[1, 2]"]
    Ustar[2,2] <- samples[s, "Ustar[2, 2]"]
    ## U = uppertri_mult_diag(Ustar, sds) means: multiply each row by sds (diag on left)
    ## i.e. U[i,j] = Ustar[i,j] * sds[i]
    U <- Ustar
    U[1,] <- U[1,] * sds[1]
    U[2,] <- U[2,] * sds[2]
    return(U)
  }
  
  stop("Cannot build cholesky U: monitor either U[1,1],U[1,2],U[2,2] OR (sds + Ustar).")
}

## --- Compute total log-likelihood for each draw (or a subset) using dbiextmix
## NOTE: Call your nimbleFunction in BEMM_Function_Temp.R
compute_loglik_total <- function(samples, idx,
                                 y,
                                 a.ind, lam.ind, lamfix,
                                 sig.ind, gamma.ind) {
  out <- numeric(length(idx))
  for (k in seq_along(idx)) {
    s <- idx[k]
    
    thres <- c(samples[s, "thres[1]"], samples[s, "thres[2]"])
    mu    <- c(samples[s, "mu[1]"],    samples[s, "mu[2]"])
    theta <- as.numeric(samples[s, paste0("theta[", 1:7, "]")])
    U     <- get_U_from_samples(samples, s)
    
    ## your nimbleFunction signature uses x=double(2) but then x[,1], so pass matrix y
    out[k] <- dbiextmix(
      x = y,
      theta = theta,
      thres = thres,
      mu = mu,
      cholesky = U,
      D = 2L,
      a.ind = a.ind,
      lam.ind = lam.ind,
      lamfix = lamfix,
      sig.ind = sig.ind,
      gamma.ind = gamma.ind,
      log = TRUE
    )
  }
  out
}

## --- 1) Build a plotting data frame of draws (subsample for speed)
S <- nrow(samples.all)

## choose up to ~5000 draws for plotting; adjust as you like
set.seed(1)
S_plot <- min(S, 5000)
idx <- sort(sample.int(S, S_plot))

draws <- data.frame(
  u1 = samples.all[idx, "thres[1]"],
  u2 = samples.all[idx, "thres[2]"]
)

## --- 2) Compute log-likelihood per draw (total)
## You must have these objects from your model setup:
## a.ind, lam.ind, lamfix, sig.ind, gamma.ind, and y (N x 2)
draws$ll_total <- compute_loglik_total(
  samples = samples.all,
  idx = idx,
  y = Y.fit,
  a.ind = 1:2,
  lam.ind = 3,
  lamfix = TRUE,
  sig.ind = 4:5,
  gamma.ind = 6:7
)

## --- 3) Hexbin aggregation: mean log-lik per bin
hb <- hexbin(draws$u1, draws$u2, xbins = 35,IDs = TRUE)  ## tune xbins for resolution

hex_df <- data.frame(
  u1 = hb@xcm,
  u2 = hb@ycm,
  n  = hb@count,
  mean_ll = tapply(draws$ll_total, hb@cID, mean)
)


p <- ggplot(draws, aes(u1, u2)) +
  ## Hex plot: mean(ll_total) per hex bin
  stat_summary_hex(
    aes(z = ll_total, fill = after_stat(value)),
    bins = 55,
    fun = mean
  ) +
  ## Density contours with color mapped to contour level (adds legend)
  stat_density_2d(
    aes(colour = after_stat(level)),
    bins = 6,
    linewidth = 0.5
  ) +
  ## --- Hex fill colors (pick one)
  # scale_fill_viridis_c(option = "D", name = "Mean log-lik\n(per hex)") +
  ## Alternatively:
  ## scale_fill_gradient(low = "white", high = "steelblue", name = "Mean log-lik\n(per hex)")
  ##
  ## --- Contour line colors + legend
  scale_colour_viridis_c(option = "C", name = "Posterior\ndensity level") +
  ## Alternatively:
  ## scale_colour_gradient(low = "grey40", high = "black", name = "Posterior\ndensity level")
  ##
  labs(
    x = expression(u[1]),
    y = expression(u[2]),
    fill = "Mean log-lik\n(per hex)",
  ) +
  theme_bw() +
  ## Make sure legends are both shown and readable
  guides(
    fill = guide_colorbar(order = 1),
    colour = guide_colorbar(order = 2)
  )+
  theme(axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        plot.title = element_text(size=18,hjust = 0.5)) 

png(filename = file.path(dir.out, "posterior_density.png"), width = 6*res, height = 5*res, res=res)
print(p)
dev.off()