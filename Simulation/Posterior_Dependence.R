library(scoringRules)

dir.work <- '/home/pgrad2/2448355h/My_PhD_Project/Biv_Ext_Mix_Mod'
dir.out <- "/home/pgrad2/2448355h/My_PhD_Project/01_Output/Biv_Ext_Mix_Mod/Simulation"

source(file.path(dir.work, "Simulation/RevExp_U_Functions.r"))
source(file.path(dir.work, "Simulation/CommonFunctions.r"))

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

# Calculate theoretical chi of the bivariate extreme mixture model.
chi.thy <- function(a){
  a1 <- max(1/a)
  a2 <- min(1/a)
  chi <- 1 - ((1+1/a1)/(1+1/a2))^(1+a2)*a1/a2/(1+a1+a2)
  return(chi)
}

ver <- '2'

all.files <- list.files(dir.out, pattern=paste('Scenario_',ver, "_itr*",sep=''))

chain_out <- c()
for (file in all.files){
  load(file=file.path(dir.out, file))
  chain_out <- c(chain_out,chain_res)
}




# Calculate the empirical chi
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

# Calculate the empirical chi bar
etaEmp1<-function(data,nq=25,qmin, qmax){
  n<-nrow(data)
  datatr<-(apply(data,2,rank) - 0.5)/n   
  qlim<-c(qmin,qmax)
  u <- seq(qlim[1], qlim[2], length = nq)
  cu<-sapply(c(1:nq), function(i) mean(apply(datatr,1,min) >= u[i]))
  res <- cbind(u,2*log(1-u)/log(cu)-1)
  return(res)
}

# Calculate the empirical chi, chi bar and Kendall's tau on the replicate data sampled from the 
# posterior predictive distribution.
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
    

    Y.bulk <- rtmvnorm(floor(n.pred*p), mean=mu, sigma=sigma, upper=thres)
    
    # GP scale tail data combined with the bulk data
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



chain_out.mat <- list()
for (itr in 1:length(chain_out)){
  tmp <- rbind(chain_out[[itr]][[1]]$samples[,],
               chain_out[[itr]][[2]]$samples[,],
               chain_out[[itr]][[3]]$samples[,])
  chain_out.mat[[itr]] <- tmp
}
rm(chain_out1,chain_out2,chain_res)

#################################Parallelly run the code#######################
NumberOfCluster <- 40
cl <- makeCluster(NumberOfCluster)
registerDoSNOW(cl)


t1 <- Sys.time()
depd_res <-
  foreach(itr = 1:length(chain_out.mat), .packages = c('mvtnorm','tmvtnorm')) %dopar%{
    samples.all <- chain_out.mat[[itr]]
    post.dsp(3000, samples.all)
  }
t2 <- Sys.time()
print(t2-t1)
stopCluster(cl)

# save(depd_res,  file=file.path(dir.out, filename='Scenario_2_Depd_Param.RData'))
save(depd_res,  file=file.path(dir.out, filename='Scenario_2_Depd_Param_1000simu.RData'))