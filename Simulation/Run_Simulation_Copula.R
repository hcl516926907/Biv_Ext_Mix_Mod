dir.work <- '/home/pgrad2/2448355h/My_PhD_Project/Biv_Ext_Mix_Mod'
dir.out <- "/home/pgrad2/2448355h/My_PhD_Project/01_Output/Biv_Ext_Mix_Mod/Simulation"
source(file.path(dir.work, "KRSW/RevExp_U_Functions.r"))
source(file.path(dir.work, "KRSW/CommonFunctions.r"))


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
packages <- c("nimble","foreach","doSNOW","parallel",'copula','posterior')  


load_install_packages(packages)

print(detectCores())


NumberOfCluster <- 3
cl <- makeCluster(NumberOfCluster)
registerDoSNOW(cl)

source(file.path(dir.work, 'Simulation/BEMM_Functions_Copula.R'))

t1 <- Sys.time()
chain_res <-
  foreach(i = 2:2) %:%
  foreach(j = 1:3, .packages = c('nimble', 'copula')) %dopar%{
    seed <- i
    d <- 2
    a <- c(0.5, 1.2)
    beta <- c(0.5, 0)
    sig <- c(0.5, 1.2)
    gamma <- c(0.3, 0.1)
    n <- 2000
    mu <- c(3.5, 4.0)
    sd1 <- 1
    sd2 <- 1.5
    rho <- 0.7
    sigma <- matrix(c(sd1^2, rho*sd1*sd2, rho*sd1*sd2, sd2^2),ncol=2)
    
    u.x <- c(5.5, 6.8)
    norm.cop <- ellipCopula("normal", param = rho, dim = 2)
    gaussian <- mvdc(norm.cop, margins = c('norm','norm'), 
                     paramMargins=list(list(mean=mu[1],sd=sd1),
                                  list(mean=mu[2],sd=sd2)))
    p <- pMvdc(u.x, gaussian)
    
    set.seed(seed)
    Y.tail<-sim.RevExpU.MGPD(n=n-floor(n*p),d=d, a=a, beta=beta, sig=sig, gamma=gamma, MGPD = T,std=T)
    
    # GP scale tail data combined with the bulk data
    Y.bulk <- matrix(NA, nrow=floor(n*p),ncol=2)
    k <- 1
    while (k<= floor(n*p)){
      rsp <- rMvdc(1, gaussian)
      if (all(rsp<=u.x)){
        Y.bulk[k,] <- rsp
        k <- k + 1
      }
    }
    plot(Y.bulk)
    
    # The name of the dataset should be Y for further WAIC calculation.
    Y <- rbind(Y.bulk, sweep(Y.tail$X,2,u.x,"+"))
    
    run_MCMC_parallel(seed=j, dat=Y, niter=20000, nburnin = 10000, thin=10)
  }

t2 <- Sys.time()
print(t2-t1)

para.name <- colnames(chain_res[[1]][[1]]$samples)
rhat.seq <- c()
ess.seq <- c()
for (name in para.name){
  post.sp <- cbind(chain_res[[1]][[1]]$samples[,name],
                   chain_res[[1]][[2]]$samples[,name],
                   chain_res[[1]][[3]]$samples[,name])
  rhat.seq <- c(rhat.seq, rhat(post.sp))
  ess.seq <- c(ess.seq, ess_basic(post.sp))
}
convg.stat <- data.frame(para.name,rhat.seq,ess.seq )


samples.all <- rbind(chain_res[[1]][[1]]$samples,
                     chain_res[[1]][[2]]$samples,
                     chain_res[[1]][[3]]$samples)

plot(chain_res[[1]][[2]]$samples[,'thres[2]'],type='l')
