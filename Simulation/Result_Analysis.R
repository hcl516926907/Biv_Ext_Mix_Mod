library(evd)
library(scoringRules)
library(mvtnorm)
library(tmvtnorm)
library(HDInterval)
library(posterior)
library(latex2exp)
library(RColorBrewer)
library(ggplot2)
source("KRSW/RevExp_U_Functions.r")
source("KRSW/CommonFunctions.r")
source("KRSW/ModelDiagnosticsNewNames.r")

dir.out <- "/home/pgrad2/2448355h/My_PhD_Project/01_Output/Biv_Ext_Mix_Mod/Simulation"


##################################Scenario 1#########################
# The following code load the MCMC samples and calculate the posterior mean, credible interval and coverage rate.
ver <- '1.3'
all.files <- list.files(dir.out, pattern=paste('Scenario_',ver, "_itr*",sep=''))

chain_out <- c()
for (file in all.files){
  load(file=file.path(dir.out, file))
  chain_out <- c(chain_out,chain_res)
}

para.name <- colnames(chain_out[[1]][[1]]$samples)
post.mean.mat <- matrix(NA, nrow=length(chain_out), ncol= length(para.name))
colnames(post.mean.mat) <- para.name
cover.ind.mat <- matrix(NA, nrow=length(chain_out), ncol= length(para.name))
colnames(cover.ind.mat) <- para.name
ci.mat <- matrix(NA, nrow=length(chain_out), ncol= length(para.name))
max_rhat <- rep(NA, length(chain_out))
Y.pred.list <- list()
for (itr in 1:length(chain_out)){
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
  max_rhat[itr] <- max(convg.stat$rhat.seq,na.rm = T)
  
  samples.all <- rbind(chain_out[[itr]][[1]]$samples[,],
                       chain_out[[itr]][[2]]$samples[,],
                       chain_out[[itr]][[3]]$samples[,]
                       )
  # True parameters
  a <- c(0.5, 1.2)
  beta <- c(0, 0)
  sig <- c(0.5, 1.2)
  if (ver=='1.1'){
    gamma <- c(0.3, 0.1)
  }else if(ver=='1.2'){
    gamma <- c(0.2, -0.2)
  }else{
    gamma <- c(-0.1, -0.3) 
  }
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
  ci.mat[itr,] <- post.ub-post.lb
}


avg.post.mean <- colMeans(post.mean.mat)
avg.ci.length <- colMeans(ci.mat)

coverage.rate <- colSums(cover.ind.mat)/length(chain_out)
max_rhat

plot(samples.all[,'theta[7]'],type='l')
plot(density(samples.all[,'thres[2]']))

plot((chain_out[[2]][[2]]$samples[,'theta[3]']),type='l')

# boxplot of the MCMC samples
set.seed(1234)
seq <- rep(NA,1000)
for (n in 1:1000){
  seq[n] <- sum(runif(n)<0.95)/n
}
plot(seq, type='l')
abline(v=150,col='red')

para.name <- c('Ustar_1_1', 'Ustar_2_1', 'Ustar_1_2','Ustar_2_2','mu_1','mu_2','s_1','s_2','a_1','a_2','lam','sigma_1',
               'sigma_2','gamma_1','gamma_2','u_1','u_2')
group.name <- c(rep('bulk',8),rep('tail',7),rep('thres',2))
df.para.plot <- data.frame()
for (i in 1:length(para.name)){
  df.tmp <- data.frame(est=post.mean.mat[,i])
  df.tmp$name <- para.name[i]
  df.tmp$tv <- true.para[i]
  df.tmp$group <- group.name[i]
  df.para.plot <- rbind(df.para.plot,df.tmp)
}
df.para.plot$name <- factor(df.para.plot$name,     # Reorder factor levels
                            para.name)
df.para.plot <- df.para.plot[!(df.para.plot$name %in% c('lam','Ustar_1_1','Ustar_2_1')),]
facet.name <- paste(para.name," (",round(coverage.rate,2),")", sep='')
facet.name <- facet.name[which(!(para.name  %in% c('lam','Ustar_1_1','Ustar_2_1')) )]
names(facet.name) <- para.name[which(!(para.name  %in% c('lam','Ustar_1_1','Ustar_2_1')) )]

ggplot(df.para.plot, aes(x=name, y=est,col=group)) +
  geom_boxplot()+
  geom_point(aes(x=name, y=tv),size=1,col='black')+
  # facet_wrap(~ name, scales="free", nrow=3, labeller = labeller(name = custom_labeller)) +
  facet_wrap(~ name, scales="free", nrow=3, labeller = labeller(name = facet.name)) +
  labs(x='parameter',y='posterior mean')+
  theme(axis.text.x=element_blank())

##################################Scenario 2#########################
# The following code compares chi, chi bar and Kendall's tau from Lidia's model with that from ours.
source("Simulation/Functions.R")
# load(file=file.path("/home/pgrad2/2448355h/My_PhD_Project/01_Output/Biv_Ext_Mix_Mod/Simulation", filename='Lidia_model_1_80.RData'))
load(file=file.path("/home/pgrad2/2448355h/My_PhD_Project/01_Output/Biv_Ext_Mix_Mod/Simulation", filename='Lidia_model_1_100.RData'))

Kconst<-NULL
est<-matrix(NA,ncol=4,nrow=length(outputM1))

thresh<-seq(0.01,0.99,len=100)
# surv is wrong in itr 1-80. 
surv<-list()
for(i in 1:length(outputM1)){
  est[i,]<-outputM1[[i]][[1]]
  Kconst[i]<-outputM1[[i]][[5]]
  surv[[i]]<-outputM1[[i]][[6]]
}

# for(i in 1:length(outputM1)){
#   est[i,]<-outputM1[[i]][[1]]
#   Kconst[i]<-outputM1[[i]][[5]]
#   sur <- rep(NA, length(thresh))
#   parct <- est[i,1]
#   parcb <- est[i,2:3]
#   parweight <- est[i,4]
#   KM1 <- Kconst[i]
#   for (j in 1:length(thresh)){
#     sur[j] <-survivalf(x=thresh[j], parct=parct,parcb=parcb,
#                            ct="InverGumbel",cb="t",weightfun=function(u,v,theta){pifun(u,v,theta)},parweight=parweight,K=KM1) 
#   }
#   surv[[i]]<-sur
#   print(i)
# }

etau_model.new<-function(thresh,survP){
  eta_model<-c()
  for(i in 1:length(thresh)){
    eta_model[i]<-2*(log(1-thresh[i]))/(log(survP[i]))-1
  }
  return(eta_model)
}

chi0.65<-chi0.7<-chi0.75<-chi0.8<-chi0.85<-chi0.9<-chi0.95<-chi0.99<-eta0.65<-eta0.7<-eta0.75<-eta0.8<-eta0.85<-eta0.9<-eta0.95<-eta0.99<-c()
for(i in 1:length(outputM1)){
  chi0.65[i]<-chiu_model(thresh=thresh[66],survP=surv[[i]][66])
  chi0.7[i]<-chiu_model(thresh=thresh[71],survP=surv[[i]][71])
  chi0.75[i]<-chiu_model(thresh=thresh[76],survP=surv[[i]][76])
  chi0.8[i]<-chiu_model(thresh=thresh[81],survP=surv[[i]][81])
  chi0.85[i]<-chiu_model(thresh=thresh[86],survP=surv[[i]][86])
  chi0.9[i]<-chiu_model(thresh=thresh[91],survP=surv[[i]][91])
  chi0.95[i]<-chiu_model(thresh=thresh[96],survP=surv[[i]][96])
  chi0.99[i]<-chiu_model(thresh=thresh[100],survP=surv[[i]][100])
  eta0.65[i]<-etau_model.new(thresh=thresh[66],survP=surv[[i]][66])
  eta0.7[i]<-etau_model.new(thresh=thresh[71],survP=surv[[i]][71])
  eta0.75[i]<-etau_model.new(thresh=thresh[76],survP=surv[[i]][76])
  eta0.8[i]<-etau_model.new(thresh=thresh[81],survP=surv[[i]][81])
  eta0.85[i]<-etau_model.new(thresh=thresh[86],survP=surv[[i]][86])
  eta0.9[i]<-etau_model.new(thresh=thresh[91],survP=surv[[i]][91])
  eta0.95[i]<-etau_model.new(thresh=thresh[96],survP=surv[[i]][96])
  eta0.99[i]<-etau_model.new(thresh=thresh[100],survP=surv[[i]][100])
}

rho <- 0.7
chi0.65geral<-(1-2*thresh[66]+pCopula(cbind(thresh[66],thresh[66]),normalCopula(rho)))/(1-thresh[66])
chi0.7geral<-(1-2*thresh[71]+pCopula(cbind(thresh[71],thresh[71]),normalCopula(rho)))/(1-thresh[71])
chi0.75geral<-(1-2*thresh[76]+pCopula(cbind(thresh[76],thresh[76]),normalCopula(rho)))/(1-thresh[76])
chi0.8geral<-(1-2*thresh[81]+pCopula(cbind(thresh[81],thresh[81]),normalCopula(rho)))/(1-thresh[81])
chi0.85geral<-(1-2*thresh[86]+pCopula(cbind(thresh[86],thresh[86]),normalCopula(rho)))/(1-thresh[86])
chi0.9geral<-(1-2*thresh[91]+pCopula(cbind(thresh[91],thresh[91]),normalCopula(rho)))/(1-thresh[91])
chi0.95geral<-(1-2*thresh[96]+pCopula(cbind(thresh[96],thresh[96]),normalCopula(rho)))/(1-thresh[96])
chi0.99geral<-(1-2*thresh[100]+pCopula(cbind(thresh[100],thresh[100]),normalCopula(rho)))/(1-thresh[100])

eta0.65geral<-2*(log(1-thresh[66]))/(log(1-2*thresh[66]+pCopula(cbind(thresh[66],thresh[66]),normalCopula(rho))))-1
eta0.7geral<-2*(log(1-thresh[71]))/(log(1-2*thresh[71]+pCopula(cbind(thresh[71],thresh[71]),normalCopula(rho))))-1
eta0.75geral<-2*(log(1-thresh[76]))/(log(1-2*thresh[76]+pCopula(cbind(thresh[76],thresh[76]),normalCopula(rho))))-1
eta0.8geral<-2*(log(1-thresh[81]))/(log(1-2*thresh[81]+pCopula(cbind(thresh[81],thresh[81]),normalCopula(rho))))-1
eta0.85geral<-2*(log(1-thresh[86]))/(log(1-2*thresh[86]+pCopula(cbind(thresh[86],thresh[86]),normalCopula(rho))))-1
eta0.9geral<-2*(log(1-thresh[91]))/(log(1-2*thresh[91]+pCopula(cbind(thresh[91],thresh[91]),normalCopula(rho))))-1
eta0.95geral<-2*(log(1-thresh[96]))/(log(1-2*thresh[96]+pCopula(cbind(thresh[96],thresh[96]),normalCopula(rho))))-1
eta0.99geral<-2*(log(1-thresh[100]))/(log(1-2*thresh[100]+pCopula(cbind(thresh[100],thresh[100]),normalCopula(rho))))-1
taumodel<-c()
nsim<-50000
n<-100000
for(i in 1:length(outputM1)){
  cstar<-sim(n=n,nsim=nsim,parct=est[i,1],parcb=c(est[i,2],est[i,3]),ct="InverGumbel",cb="t",
             weightfun=function(u,v,theta){pifun(u,v,theta)},parweight=est[i,4])
  taumodel[i]<-cor(cstar[,1],cstar[,2],method="kendall")
  print(i)
}
# kendall's \tau for the gaussian with \rho=0.65
kendalgeral<-2/pi*asin(rho)


# load(file=file.path('/home/pgrad2/2448355h/My_PhD_Project/01_Output/Biv_Ext_Mix_Mod/Simulation', filename='Scenario_2_Depd_Param.RData'))
load(file=file.path('/home/pgrad2/2448355h/My_PhD_Project/01_Output/Biv_Ext_Mix_Mod/Simulation', filename='Scenario_2_Depd_Param_1000simu.RData'))
# The 1000simu versions wrongly included the MCMC output indexed in 21-40.
depd_res <- c(depd_res[1:20],depd_res[41:1020])
#------------------------------------------tau---------------------
get_tau <- function(depd_res){
  post.mean.seq <- c()
  for (i in 1:length(depd_res)){
    post.mean <- mean(depd_res[[i]][[1]][,1])
    post.mean.seq <- c(post.mean.seq,post.mean)
  }
  return(post.mean.seq)
}

tau.palette <- brewer.pal(8,"Accent")
df.tau1 <- data.frame(tau=taumodel,tv=kendalgeral,dataset='BMCM')
df.tau2 <- data.frame(tau=get_tau(depd_res),tv=kendalgeral,dataset='BEMM')
df.tau <- rbind(df.tau1,df.tau2)

res <- 500
p <- ggplot(df.tau, aes(x=dataset,y=tau,colour=dataset))+
  geom_boxplot() +
  geom_point(aes(x=dataset, y=tv),size=2,col='black')+
  labs(x='model',y=TeX(r'($\tau$)'))+
  scale_colour_manual(values=c('BMCM'=tau.palette[3],'BEMM'=tau.palette[5]))+
  theme(axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size=15))
print(p)
png(filename = file.path(dir.out, "Simulation_2_tau.png"), width = 6*res, height = 5*res, res=res)
print(p)
dev.off()

#---------------------------------------------chi----------------------------
lidia.chi <- data.frame('u' = c(0.65, 0.7,0.75,0.8,0.85,0.9,0.95,0.99),
                        'chi_data'=c(chi0.65geral,chi0.7geral,chi0.75geral,chi0.8geral,
                                     chi0.85geral,chi0.9geral,chi0.95geral,chi0.99geral),
                        'chi_lower'=c(quantile(chi0.65,0.025),quantile(chi0.7,0.025),
                                      quantile(chi0.75,0.025),quantile(chi0.8,0.025),
                                      quantile(chi0.85,0.025),quantile(chi0.9,0.025),
                                      quantile(chi0.95,0.025),quantile(chi0.99,0.025)),
                        'chi_upper'=c(quantile(chi0.65,0.975),quantile(chi0.7,0.975),
                                      quantile(chi0.75,0.975),quantile(chi0.8,0.975),
                                      quantile(chi0.85,0.975),quantile(chi0.9,0.975),
                                      quantile(chi0.95,0.975),quantile(chi0.99,0.975))
                        )

depd_res[[2]]

cb <- function(depd_res){
  n <- length(depd_res)
  u <- depd_res[[1]][[2]][[1]][,1]
  chiu <- c()
  for (i in 1:n){
    post.mean.mat <- c()
    for (j in 1:length(depd_res[[i]][[2]])){
      post.mean.mat <- cbind(post.mean.mat,depd_res[[i]][[2]][[j]][,2])
    }
    post.mean <- rowMeans(post.mean.mat)
    chiu <- cbind(chiu,post.mean)    
  }
  lb <- apply(chiu,1, quantile, 0.025)
  ub <- apply(chiu,1, quantile, 0.975)
  mean <- apply(chiu,1, mean)
  return(cbind(u,lb,ub,mean))
}
my.chi.cb<- data.frame(cb(depd_res))

chi_plot <- lidia.chi
chi_plot$my_ub <- my.chi.cb[as.character(my.chi.cb$u) %in% as.character(chi_plot$u),]$ub
chi_plot$my_lb <- my.chi.cb[as.character(my.chi.cb$u) %in% as.character(chi_plot$u),]$lb

simu2.palette <- brewer.pal(8,"Accent")
dodge_amount <- 0.01
p <- ggplot(chi_plot, aes(x=u, y=chi_data)) +
  geom_point(aes(x = u - dodge_amount/2),size=2,col=simu2.palette[3]) +
  geom_point(aes(x = u + dodge_amount/2),size=2,col=simu2.palette[5]) +

  geom_errorbar(
    aes(ymin=chi_lower, ymax=chi_upper, x = u - dodge_amount/2),
    width=0.005, color=simu2.palette[3],
    position=position_dodge(dodge_amount)
  ) +
  geom_errorbar(
    aes(ymin=my_lb, ymax=my_ub, x = u + dodge_amount/2),
    width=0.005, color=simu2.palette[5],
    position=position_dodge(dodge_amount)
  ) +
  scale_x_continuous(breaks = c(0.65,0.7,0.75,0.8,0.85,0.9,0.95,0.99))+
  labs(x='r',y=expression(chi(r))) +
  theme(axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15))


print(p)

png(filename = file.path(dir.out, "Simulation_2_chi.png"), width = 6*res, height = 5*res, res=res)
print(p)
dev.off()

#--------------------------------------boxplot version of chi------------------------
lidia.chi <- data.frame('u' = rep(c(0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.99),each=100),
                        'chi_est'=c(chi0.65,chi0.7, chi0.75, chi0.8, chi0.85,chi0.9,
                                    chi0.95, chi0.99))
lidia.chi$chi_data=rep(c(chi0.65geral, chi0.7geral, chi0.75geral, chi0.8geral, 
                       chi0.85geral, chi0.9geral, chi0.95geral,chi0.99geral),each=100)
lidia.chi$model <- 'BMCM'

chiboxplot <- function(depd_res){
  n <- length(depd_res)
  u <- c(0.65,0.7,0.75,0.8,0.85,0.9,0.95,0.99)
  chiu <- data.frame()
  u.flt <- which(as.character(depd_res[[1]][[2]][[1]][,'u']) %in% as.character(u))
  for (i in 1:n){
    post.mean.mat <- c()
    for (j in 1:length(depd_res[[i]][[2]])){
      post.mean.mat <- cbind(post.mean.mat,depd_res[[i]][[2]][[j]][u.flt,2])
    }
    post.mean <- rowMeans(post.mean.mat)
    df.post.mean <- data.frame('u'=u, 'chi_est'=post.mean)
    chiu <- rbind(chiu,df.post.mean)    
  }
  return(chiu)
}
my.chi<- data.frame(chiboxplot(depd_res))
my.chi$chi_data=rep(c(chi0.65geral, chi0.7geral, chi0.75geral, chi0.8geral, 
                         chi0.85geral, chi0.9geral, chi0.95geral,chi0.99geral),times=1000)
my.chi$model <- 'BEMM'

df.chi.boxplot <- rbind(lidia.chi,my.chi)
df.chi.boxplot$u <- as.character(df.chi.boxplot$u)
p <- ggplot(df.chi.boxplot, aes(x=u,y=chi_est,color=model))+
  geom_boxplot() +
  geom_point(aes(x=u, y=chi_data),size=2,col='black')+
  labs(x='model',y=TeX(r'($\tau$)'))+
  scale_colour_manual(values=c('BMCM'=tau.palette[3],'BEMM'=tau.palette[5]))+
  scale_x_discrete(breaks = c(0.65,0.7,0.75,0.8,0.85,0.9,0.95,0.99))+
  labs(x='r',y=expression(chi(r))) +
  theme(axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        legend.position = "none")
print(p)

png(filename = file.path(dir.out, "Simulation_2_chi_boxplot.png"), width = 6*res, height = 5*res, res=res)
print(p)
dev.off()


#------------------------------------------ chi bar --------------------------------------------
lidia.chib <- data.frame('u' = c(0.65, 0.7,0.75,0.8,0.85,0.9,0.95,0.99),
                        'chib_data'=c(eta0.65geral,eta0.7geral,eta0.75geral,eta0.8geral,
                                      eta0.85geral,eta0.9geral,eta0.95geral,eta0.99geral),
                        'chib_lower'=c(quantile(eta0.65,0.025),quantile(eta0.7,0.025),
                                      quantile(eta0.75,0.025),quantile(eta0.8,0.025),
                                      quantile(eta0.85,0.025),quantile(eta0.9,0.025),
                                      quantile(eta0.95,0.025),quantile(eta0.99,0.025)),
                        'chib_upper'=c(quantile(eta0.65,0.975),quantile(eta0.7,0.975),
                                      quantile(eta0.75,0.975),quantile(eta0.8,0.975),
                                      quantile(eta0.85,0.975),quantile(eta0.9,0.975),
                                      quantile(eta0.95,0.975),quantile(eta0.99,0.975))
)


chib.cb <- function(depd_res){
  n <- length(depd_res)
  u <- depd_res[[1]][[2]][[1]][,1]
  chib.u <- c()
  for (i in 1:n){
    post.mean.mat <- c()
    for (j in 1:length(depd_res[[i]][[3]])){
      post.mean.mat <- cbind(post.mean.mat,depd_res[[i]][[3]][[j]][,2])
    }
    post.mean <- rowMeans(post.mean.mat)
    chib.u <- cbind(chib.u,post.mean)    
  }
  lb <- apply(chib.u,1, quantile, 0.025)
  ub <- apply(chib.u,1, quantile, 0.975)
  mean <- apply(chib.u,1, mean)
  return(cbind(u,lb,ub,mean))
}
my.chib.cb<- data.frame(chib.cb(depd_res))


chib_plot <- lidia.chib
chib_plot$my_ub <- my.chib.cb[as.character(my.chib.cb$u) %in% as.character(chib_plot$u),]$ub
chib_plot$my_lb <- my.chib.cb[as.character(my.chib.cb$u) %in% as.character(chib_plot$u),]$lb

simu2.palette <- brewer.pal(8,"Accent")
dodge_amount <- 0.01
cols <- c("BMCM"=simu2.palette[3],'BEMM'=simu2.palette[5])
p <- ggplot(chib_plot, aes(x=u, y=chib_data)) +
  geom_point(aes(x = u - dodge_amount/2),size=2,color=simu2.palette[3]) +
  geom_point(aes(x = u + dodge_amount/2),size=2,color=simu2.palette[5]) +
  
  geom_errorbar(
    aes(ymin=chib_lower, ymax=chib_upper, x = u - dodge_amount/2,colour="BMCM"),
    width=0.005,
    position=position_dodge(dodge_amount)
  ) +
  geom_errorbar(
    aes(ymin=my_lb, ymax=my_ub, x = u + dodge_amount/2, colour='BEMM'),
    width=0.005,
    position=position_dodge(dodge_amount)
  ) +
  scale_x_continuous(breaks = c(0.65,0.7,0.75,0.8,0.85,0.9,0.95,0.99))+
  labs(x='r',y=TeX(r'($\bar{\chi}(r)$)')) +
  theme(axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        legend.position=c(0, 1),
        legend.justification = c(0, 1),
        legend.title = element_blank(),
        legend.text = element_text(size=15),
        legend.key.size = unit(0.05, 'npc'),
        legend.key.width = unit(0.1, 'npc')) + 
  scale_colour_manual(values=cols)
print(p)

png(filename = file.path(dir.out, "Simulation_2_chibar.png"), width = 6*res, height = 5*res, res=res)
print(p)
dev.off()


#----------------------------------chi bar boxplot--------------------------
lidia.chibar <- data.frame('u' = rep(c(0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.99),each=100),
                        'chibar_est'=c(eta0.65,eta0.7, eta0.75, eta0.8, eta0.85,eta0.9,
                                    eta0.95, eta0.99))
lidia.chibar$chibar_data=rep(c(eta0.65geral, eta0.7geral, eta0.75geral, eta0.8geral, 
                         eta0.85geral, eta0.9geral, eta0.95geral,eta0.99geral),each=100)
lidia.chibar$model <- 'BMCM'

chibar.boxplot <- function(depd_res){
  n <- length(depd_res)
  u <- c(0.65,0.7,0.75,0.8,0.85,0.9,0.95,0.99)
  chibaru <- data.frame()
  u.flt <- which(as.character(depd_res[[1]][[3]][[1]][,'u']) %in% as.character(u))
  for (i in 1:n){
    post.mean.mat <- c()
    for (j in 1:length(depd_res[[i]][[3]])){
      post.mean.mat <- cbind(post.mean.mat,depd_res[[i]][[3]][[j]][u.flt,2])
    }
    post.mean <- rowMeans(post.mean.mat)
    df.post.mean <- data.frame('u'=u, 'chibar_est'=post.mean)
    chibaru <- rbind(chibaru,df.post.mean)    
  }
  return(chibaru)
}
my.chibar<- data.frame(chibar.boxplot(depd_res))
my.chibar$chibar_data=rep(c(eta0.65geral, eta0.7geral, eta0.75geral, eta0.8geral, 
                      eta0.85geral, eta0.9geral, eta0.95geral,eta0.99geral),times=1000)
my.chibar$model <- 'BEMM'

df.chibar.boxplot <- rbind(lidia.chibar,my.chibar)
df.chibar.boxplot$u <- as.character(df.chibar.boxplot$u)
p <- ggplot(df.chibar.boxplot, aes(x=u,y=chibar_est,color=model))+
  geom_boxplot() +
  geom_point(aes(x=u, y=chibar_data),size=2,col='black')+
  labs(x='model',y=TeX(r'($\tau$)'))+
  scale_colour_manual(values=c('BMCM'=tau.palette[3],'BEMM'=tau.palette[5]))+
  scale_x_discrete(breaks = c(0.65,0.7,0.75,0.8,0.85,0.9,0.95,0.99))+
  labs(x='r',y=TeX(r'($\bar{\chi}(r)$)')) +
  theme(axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        legend.position=c(0, 1),
        legend.justification = c(0, 1),
        legend.title = element_blank(),
        legend.text = element_text(size=15),
        legend.key.size = unit(0.05, 'npc'),
        legend.key.width = unit(0.1, 'npc'))
print(p)

png(filename = file.path(dir.out, "Simulation_2_chibar_boxplot.png"), width = 6*res, height = 5*res, res=res)
print(p)
dev.off()

