
source("Simulation/RevExp_U_Functions.r")
source("Simulation/CommonFunctions.r")
source("Simulation/Gumbel_U_Functions.r")


load_install_packages <- function(packages) {
  for(package in packages){
    # If the package is not installed, install it
    if(!require(package, character.only = TRUE)) {
      install.packages(package, dependencies = TRUE,INSTALL_opts = '--no-lock')
      # Load the package after installation
      library(package, character.only = TRUE)
    } else {
      # If the package is already installed, just load it
      library(package, character.only = TRUE)
    }
  }
}

# List the packages you want to load
packages <- c("nimble", "foreach","doSNOW","parallel",'gsl',
              'copula','extraDistr','rugarch','tidyr','RColorBrewer',
              'MASS','ggplot2','latex2exp')  


load_install_packages(packages)

dir.out <- "/home/pgrad2/2448355h/My_PhD_Project/01_Output/Biv_Ext_Mix_Mod/Simulation"
#################################Used functions#######################
emp.mypalette <- brewer.pal(9,"Greys")
pred.mypalette <- brewer.pal(9,"Blues")
qq.palette <- brewer.pal(8,"Accent")
res <- 300

# load(file=file.path(dir.out, 'Simulation_4_3332_thres_0.01_0.90_0.99_sd_20.RData') )
load(file=file.path(dir.out, 'Simulation_4_3332_thres_0.01_0.90_0.99_sd_0.5.RData') )
# load(file=file.path(dir.out, 'Simulation_4_3332_thres_0.01_0.90_0.99_sd_1.RData') )




# Function to generate posterior predictive replicates
post.rep <-  function(n.rep, rep.size, samples, seed=1234, lbound.tail=rep(-Inf,2), ubound.tail=rep(Inf,2),
                      total.dist.name){
  set.seed(seed)
  Y.rep <- list()
  idx <- sample(nrow(samples),size=n.rep, replace=TRUE)
  d <- 2
  for (k in 1:n.rep){
    print(k)
    a <-  samples[idx[k], c('theta[1]')]
    sig <- samples[idx[k], c('theta[2]','theta[3]')]
    gamma <- samples[idx[k], c('theta[4]','theta[5]')]
    
    params <- samples[idx[k], c('params.bulk[1]','params.bulk[2]','params.bulk[3]',
                                'params.bulk[4]','params.bulk[5]')]
    
    bulk.dist.name <- total.dist.name[1:3]
    cop.name.map <- c("normal","clayton","gumbel","frank","joe","plackett")
    
    cop <- switch(cop.name.map[bulk.dist.name[1]],
                  "normal" = normalCopula(param = params[1], dim = 2),
                  "clayton" = claytonCopula(param = params[1], dim = 2),
                  "gumbel" = gumbelCopula(param = params[1], dim = 2),
                  "frank" = frankCopula(param = params[1], dim = 2),
                  "joe" = joeCopula(param = params[1], dim = 2),
                  "plackett" = plackettCopula(param = params[1]),
                  stop("Unsupported copula name"))
    
    margin.name.map <- c("norm","exp","gamma","lnorm","weibull",'lst')
    margins.name <- margin.name.map[bulk.dist.name[2:3]]
    params.margin <- params[-1]
    transformed_params <- list()
    start.idx <- 1
    for (i in 1:length(margins.name)) {
      m <- margins.name[i]
      
      if (m == "norm") {
        n.param <- 2
        p <- params.margin[ start.idx : (start.idx+1)]
        transformed_params[[i]] <- list(mean = p[1], sd = p[2]) # ensure sd > 0
        start.idx <- start.idx + n.param
        
      } else if (m == "exp") {
        n.param <- 1
        p <- params.margin[ start.idx]
        transformed_params[[i]] <- list(rate = p[1]) # ensure rate > 0
        start.idx <- start.idx + n.param
        
      } else if (m == "gamma") {
        n.param <- 2
        p <- params.margin[ start.idx : (start.idx+1)]
        transformed_params[[i]] <- list(shape = p[1], rate = p[2]) # ensure shape, rate > 0
        start.idx <- start.idx + n.param
        
      } else if (m == "lnorm") {
        n.param <- 2
        p <- params.margin[ start.idx : (start.idx+1)]
        transformed_params[[i]] <- list(meanlog = p[1], sdlog = p[2]) # ensure sdlog > 0
        start.idx <- start.idx + n.param
      } else if (m == "weibull") {
        n.param <- 2
        p <- params.margin[ start.idx : (start.idx+1)]
        transformed_params[[i]] <- list(shape = p[1], scale = p[2]) # ensure shape, scale > 0
        start.idx <- start.idx + n.param
      } else if (m == 'lst'){
        n.param <- 3
        p <- params.margin[ start.idx : (start.idx+2)]
        transformed_params[[i]] <- list(df = p[1], mu = p[2], sigma = p[3]) 
        start.idx <- start.idx + n.param
      } else {
        stop("Unsupported marginal distribution")
      }
    }
    
    
    bulk.dist <- mvdc(cop, margins = margins.name, 
                      
                      paramMargins=list(transformed_params[[1]],
                                        
                                        transformed_params[[2]])
    )
    
    thres <- samples[idx[i], c('thres[1]','thres[2]')]
    
    p <- pMvdc(thres, bulk.dist)
    print(p)
    
    n.Y.bulk <- 0
    Y.bulk.ut <- matrix(, nrow = 0, ncol = 2)
    while (n.Y.bulk< (floor(rep.size*p) )){
      Y.bulk <- rMvdc(floor(rep.size*p), bulk.dist)
      Y.bulk <- Y.bulk[Y.bulk[,1]<thres[1] & Y.bulk[,2]<thres[2],]
      n.Y.bulk <- n.Y.bulk + nrow(Y.bulk)
      Y.bulk.ut <- rbind(Y.bulk.ut, Y.bulk)
    }
    Y.bulk <- Y.bulk.ut[sample(1:n.Y.bulk, floor(rep.size*p), replace=TRUE),]
    
    n.Y.tail <- 0
    Y.tail.ut <- matrix(, nrow = 0, ncol = 2)
    while (n.Y.tail< (rep.size-floor(rep.size*p) )){
      if (total.dist.name[4]==1){
        Y.tail <- sim.RevExpU.MGPD(n=rep.size, d=2, a=rep(a,2), beta=c(0,0), sig=sig, gamma=gamma, MGPD = T,std=F)
      }else{
        Y.tail<-sim.GumbelU.MGPD(n=rep.size, d=2, a=a, beta=c(0,0), sig=sig, gamma=gamma, MGPD = T,std=F)
      }
      
      
      bound.cond.1 <- Y.tail[,1] > (lbound.tail-thres)[1] & Y.tail[,1] < (ubound.tail-thres)[1]
      bound.cond.2 <- Y.tail[,2] > (lbound.tail-thres)[2] & Y.tail[,2] < (ubound.tail-thres)[2]
      Y.tail <- Y.tail[bound.cond.1&bound.cond.2,]
      n.Y.tail <- n.Y.tail+ nrow(Y.tail)
      Y.tail.ut <- rbind(Y.tail.ut, Y.tail)
    }
    
    Y.tail <- Y.tail.ut[sample(1:n.Y.tail, rep.size-floor(rep.size*p), replace=TRUE),]
    
    
    # The name of the dataset should be Y for further WAIC calculation.
    Y <- rbind(Y.bulk, sweep(Y.tail,2,thres,"+"))
    Y.rep[[k]] <- Y
  }
  return(Y.rep)
}


n.rep <- 500
Y.fit
post.burnin <- 3000:4999

plot(density(results$samples[post.burnin, 'thres[1]' ]))

plot(density(results$samples[post.burnin, 'thres[2]' ]))
plot(results$samples[1:99, 'thres[1]' ],type='l')
plot(results$samples[post.burnin, 'thres[1]' ],type='l')

theta=colMeans(results$samples[post.burnin,c('theta[1]','theta[2]','theta[3]','theta[4]'
                                             ,'theta[5]')])
thres=colMeans(results$samples[post.burnin,c('thres[1]','thres[2]')])
params.bulk=colMeans(results$samples[post.burnin,c('params.bulk[1]','params.bulk[2]','params.bulk[3]',
                                                   'params.bulk[4]','params.bulk[5]')])

apply(Y.fit,2,quantile,0.98)
Y.fit.rep <- post.rep(n.rep=n.rep, rep.size=nrow(Y.fit), samples=results$samples[post.burnin, ], seed=1234,
                      lbound.tail=rep(0,2),total.dist.name=total.dist.name)


probs <- c(seq(0.01, 0.9, length.out = 100),seq(0.91,0.999,length.out=50))
# probs <- seq(0.01, 0.99, length.out = 100)
emp_quant1 <- quantile(Y.fit[,1], probs = probs)

est1 <- fitdistr(Y.fit[,1], "gamma")
thy_quant1 <- qgamma(probs, shape= est1$estimate[1], rate = est1$estimate[2])

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
  ggtitle("log.ftse.return") +
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


ggplot(data.frame(Y), aes(x=X1)) + 
  geom_histogram(color="black", fill="white")



ggplot(data.frame(Y), aes(x=X2)) + 
  geom_histogram(color="black", fill="white")




emp_quant2 <- quantile(Y.fit[,2], probs = probs)

est2 <- fitdistr(Y.fit[,2], "gamma")

thy_quant2 <- qgamma(probs, shape= est2$estimate[1], rate = est2$estimate[2])


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
  ggtitle("log.dax.return") +
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


df.joint.exceed <- expand.grid(seq(0,10,0.1),seq(0,10,0.1))
df.joint.exceed$cnt <- 0
for (i in 1:nrow(df.joint.exceed)){
  thres <- as.numeric(df.joint.exceed[i,c("Var1",'Var2')])
  for (j in 1:n.rep){
    df.joint.exceed[i,'cnt'] <- df.joint.exceed[i,'cnt']+ sum(Y.fit.rep[[j]][,1]>thres[1] & Y.fit.rep[[j]][,2]>thres[2])
  }
}

df.joint.exceed$prob <- df.joint.exceed$cnt/n.rep/nrow(Y.fit.rep[[1]])
library(viridis)
ggplot(df.joint.exceed, aes(x = Var1, y = Var2, fill = prob))+  geom_tile()+scale_fill_viridis(discrete=FALSE)
