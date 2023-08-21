dir.work <- '/home/pgrad2/2448355h/My_PhD_Project/Biv_Ext_Mix_Mod'
dir.data <- '/home/pgrad2/2448355h/My_PhD_Project/00_Dataset/Air_Pollution_India'
dir.out <- "/home/pgrad2/2448355h/My_PhD_Project/01_Output/Biv_Ext_Mix_Mod/Air_Pollution_India"
source("KRSW/RevExp_U_Functions.r")
source("KRSW/CommonFunctions.r")
source("KRSW/ModelDiagnosticsNewNames.r")

load_install_packages <- function(packages) {
  for(package in packages){
    # If the package is not installed, install it
    if(!require(package, character.only = TRUE)) {
      install.packages(package, dependencies = TRUE,repos='http://cran.us.r-project.org')
      # Load the package after installation
      library(package, character.only = TRUE)
    } else {
      # If the package is already installed, just load it
      library(package, character.only = TRUE)
    }
  }
}

# List the packages you want to load
packages <- c("nimble", "mvtnorm", "tmvtnorm","foreach","doSNOW","parallel", "dplyr",'posterior')  


load_install_packages(packages)


pollute_2023 <- read.csv(file.path(dir.data,"DL035_2023.csv"))
pollute_2022 <- read.csv(file.path(dir.data,"DL035_2022.csv"))
pollute_2021 <- read.csv(file.path(dir.data,"DL035_2021.csv"))
pollute_2020 <- read.csv(file.path(dir.data,"DL035_2020.csv"))
pollute_2019 <- read.csv(file.path(dir.data,"DL035_2019.csv"))
pollute <- rbind(pollute_2019,pollute_2020,pollute_2021,pollute_2022,pollute_2023)

pollute$To.Date <- substring(pollute$To.Date,1,10)
for (col in c('PM2.5','Ozone','NO','NO2','NOx','NH3',"PM10",'WS',"WD")){
  pollute[,col] <- as.numeric(pollute[,col])
}

pollute[is.na(pollute[,'PM2.5']),]


pollute.max <- pollute %>% group_by(To.Date)%>% 
  summarise(PM2.5=max(PM2.5, na.rm=T),
            Ozone=max(Ozone, na.rm=T),
            NO=max(NO, na.rm=T),
            NO2=max(NO2, na.rm=T),
            NOx=max(NOx, na.rm=T),
            NH3=max(NH3, na.rm=T),
            PM10=max(PM10, na.rm=T))

pm2.5_no2 <- pollute.max %>% 
           filter_at(vars(PM2.5,NO2), all_vars(!is.infinite(.))) %>% 
           select(PM2.5,NO2)
print(dim(pm2.5_no2))
plot(pm2.5_no2)

pollute.mean <- pollute %>% group_by(To.Date)%>% 
  summarise(PM2.5=mean(PM2.5, na.rm=T),
            Ozone=mean(Ozone, na.rm=T),
            NO=mean(NO, na.rm=T),
            NO2=mean(NO2, na.rm=T),
            NOx=mean(NOx, na.rm=T),
            NH3=mean(NH3, na.rm=T),
            PM10=mean(PM10, na.rm=T))

pm2.5_no2_mean <- pollute.mean %>% 
  filter_at(vars(PM2.5,NO2), all_vars(!is.na(.))) %>% 
  select(PM2.5,NO2)
print(dim(pm2.5_no2_mean))
plot(pm2.5_no2_mean)


# load(file=file.path(dir.out, filename='Air_pollution_mvn_0.8_0.99.RData'))
# load(file=file.path(dir.out, filename='Air_pollution_mvtn_0.5_0.99.RData'))
# load(file=file.path(dir.out, filename='Air_pollution_mvtn_0.6_0.99_AFSlice.RData'))
load(file=file.path(dir.out, filename='Air_pollution_mvtn_0.6_0.99_AFSlice_daily_mean.RData'))

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



plot(samples.all[,'sds[2]'],type='l')
plot(samples.all[,'thres[1]'],type='l')
plot(samples.all[,'thres[2]'],type='l')
plot(samples.all[,'theta[6]'],type='l')



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
      Y.pred[i,] <- rtmvnorm(1, mean=mu, sigma=sigma, upper=thres)
    }else{
      Y.tail <- sim.RevExpU.MGPD(n=1,d=d, a=a, beta=c(0,0), sig=sig, gamma=gamma, MGPD = T,std=T)$X
      Y.pred[i,] <- thres + Y.tail
    }
  }
  return(Y.pred)
}

Y.pred <- post.pred(1600, samples.all,seed=2)
plot(Y.pred)
plot(Y.pred[(Y.pred[,1]<600)&(Y.pred[,1]>0),])

plot(pm2.5_no2_mean)




chi.thy <- function(a){
  a1 <- max(1/a)
  a2 <- min(1/a)
  chi <- 1 - ((1+1/a1)/(1+1/a2))^(1+a2)*a1/a2/(1+a1+a2)
  return(chi)
}


post.dsp <- function(n.sp, samples, n.pred=2000, seed=1234){
  set.seed(seed)
  idx <- sample(nrow(samples),size=n.sp, replace=TRUE)
  d <- 2
  chi.seq <- rep(NA, length(idx))
  tau.seq <- rep(NA, length(idx))
  sCor.seq <- rep(NA<length(idx))
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
    tau.seq[i] <- cor(Y[,1], Y[,2], method='kendall')
    sCor.seq[i] <- cor(Y[,1], Y[,2], method = "spearman")
  }
  dps.mat <- cbind(tau.seq, sCor.seq, chi.seq)
  return(dps.mat)
}



t1 <- Sys.time()
dept.est <- post.dsp(3000, samples.all, n.pred=1639)
t2 <- Sys.time()


tau.emp <- cor(pm2.5_no2_mean$PM2.5,pm2.5_no2_mean$NO2, method='kendall')
sCor.emp <- cor(pm2.5_no2_mean$PM2.5,pm2.5_no2_mean$NO2, method='spearman')
chi.emp <- chiEmp(data=pm2.5_no2_mean,nq=100,qmin=0.50,qmax=0.99)

plot(density(dept.est[,1]), main='predicted Kendall tau')
abline(v=tau.emp, col='red')
legend('topright', legend=c("Predicted", "Empirical"),
       col=c("black", "red"), lty=1:1, cex=2)


plot(density(dept.est[,2]),main='predicted Spearman correlation')
abline(v=sCor.emp,col='red')
legend('topright', legend=c("Predicted", "Empirical"),
       col=c("black", "red"), lty=1:1, cex=2)

plot(density(dept.est[,3]), main='predicted chi')
plot(chi.emp,col='red', main='predicted chi')
abline(h=mean(dept.est[3]),col='black')
legend('topright', legend=c("Posterior mean", "Empirical"),
       col=c("black", "red"), lty=1:1, cex=2)

chiEmp(Y.pred[(Y.pred[,1]<600)&(Y.pred[,1]>0)&(Y.pred[,2]<100),], nq=100,qmin=0.50,qmax=0.99)

chiEmp(pm2.5_no2_mean, nq=100,qmin=0.50,qmax=0.99)


plot(density(Y.pred[,1]), main='margin PM2.5')
lines(density(pm2.5_no2_mean$PM2.5), col='red')
abline(v=mean(samples.all[,'thres[1]']), lty=2,col='blue')
legend('topright', legend=c("Predicted", "Empirical","predicted threshold"),
       col=c("black", "red",'blue'), lty=1:1:2, cex=2)



plot(density(Y.pred[,2]), main='margin NO2')
lines(density(pm2.5_no2_mean$NO2), col='red')
abline(v=mean(samples.all[,'thres[2]']), lty=2,col='blue')
legend('topright', legend=c("Predicted", "Empirical","predicted threshold"),
       col=c("black", "red",'blue'), lty=1:1:2, cex=2)


