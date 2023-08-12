library(scoringRules)
library(HDInterval)

source("KRSW/RevExp_U_Functions.r")
source("KRSW/CommonFunctions.r")
source("KRSW/ModelDiagnosticsNewNames.r")

dir.out <- "/home/pgrad2/2448355h/My_PhD_Project/01_Output/Biv_Ext_Mix_Mod/Simulation"

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
packages <- c("nimble", "mvtnorm", "tmvtnorm","foreach","doSNOW","parallel",'scoringRules')  


load_install_packages(packages)
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

chi.thy <- function(a){
  a1 <- max(1/a)
  a2 <- min(1/a)
  chi <- 1 - ((1+1/a1)/(1+1/a2))^(1+a2)*a1/a2/(1+a1+a2)
  return(chi)
}

load(file=file.path(dir.out, 'Scenario_1.1_itr1_20_lamfix.RData'))
chain_out1 <- chain_res

load(file=file.path(dir.out, 'Scenario_1.1_itr21_100_lamfix.RData'))
chain_out2 <- chain_res

load(file=file.path(dir.out, 'Scenario_1.1_itr101_150_lamfix.RData'))
chain_out3 <- chain_res
chain_out <- c(chain_out1,chain_out2,chain_out3)

dim(chain_out[[1]]$samples)
itr <- 1
samples.all <- rbind(chain_out[[itr]][[1]]$samples[,],
                     chain_out[[itr]][[2]]$samples[,],
                     chain_out[[itr]][[3]]$samples[,])

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


NumberOfCluster <- 40
cl <- makeCluster(NumberOfCluster)
registerDoSNOW(cl)

t1 <- Sys.time()
chain_res <-
  foreach(itr = 1:3, .packages = c('nimble','mvtnorm','tmvtnorm')) %dopar%{
    samples.all <- rbind(chain_out[[itr]][[1]]$samples[,],
                         chain_out[[itr]][[2]]$samples[,],
                         chain_out[[itr]][[3]]$samples[,])
    
    post.dsp(3000, samples.all)
  }
t2 <- Sys.time()
print(t2-t1)
stopCluster(cl)

save(chain_res,  file=file.path(dir.out, filename='Scenario_1.1_Depd_Param.RData'))