dir.work <- '/home/pgrad2/2448355h/My_PhD_Project/Biv_Ext_Mix_Mod'
dir.out <- "/home/pgrad2/2448355h/My_PhD_Project/01_Output/Biv_Ext_Mix_Mod/Simulation"
source(file.path(dir.work, "KRSW/RevExp_U_Functions.r"))
source(file.path(dir.work, "KRSW/CommonFunctions.r"))

# install_load_packages <- function(packages) {
#   nodename <- Sys.info()['nodename']
#   lib_path <- file.path("/home/pgrad2/2448355h/R/library", nodename)
#   .libPaths(lib_path)
#   # Create a separate library folder for each node, if it doesn't exist
#   if (!dir.exists(lib_path)) {
#     dir.create(lib_path, recursive = TRUE)
#   }
#   
#   for (package in packages) {
#     # Check if the package is installed in the specific folder
#     if(!require(package, character.only = TRUE, lib.loc = lib_path)) {
#       install.packages(package, lib = lib_path, dependencies = TRUE)
#       library(package, character.only = TRUE, lib.loc = lib_path)
#     }
#     
#     # Load the package
#     library(package, character.only = TRUE, lib.loc = lib_path)
#   }
# }

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
packages <- c("nimble", "mvtnorm", "tmvtnorm","foreach","doSNOW","parallel")  


load_install_packages(packages)

print(detectCores())
# t3 <- Sys.time()
# for (i in 1:1) {
#   print(i)
#   seed <- i
#   d <- 2
#   a <- c(1.5, 2)
#   beta <- c(0.5, 0)
#   sig <- c(0.5, 0.4)
#   gamma <- c(0.2, 0.1)
#   n <- 2000
#   mu <- c(3.5, 4.5)
#   sd1 <- 1.22
#   sd2 <- 1.10
#   rho <- 0.72
#   sigma <- matrix(c(sd1^2, rho*sd1*sd2, rho*sd1*sd2, sd2^2),ncol=2)
# 
#   u.x <- c(5.6, 7)
#   p <- pmvnorm(lower=rep(0,2), upper=u.x, mean=mu, sigma=sigma, keepAttr = F)
# 
#   set.seed(seed)
#   Y.tail<-sim.RevExpU.MGPD(n=n-floor(n*p),d=d, a=a, beta=beta, sig=sig, gamma=gamma, MGPD = T,std=T)
# 
#   # GP scale tail data combined with the bulk data
#   set.seed(seed)
#   Y.bulk <- rtmvnorm(floor(n*p), mean=mu, sigma=sigma, lower=c(0,0),upper=u.x)
# 
#   # The name of the dataset should be Y for further WAIC calculation.
#   Y <- rbind(Y.bulk, sweep(Y.tail$X,2,u.x,"+"))
#   # plot(Y)
# 
#   source(file.path(dir.work, 'Simulation/BEMM_Functions.R'))
# 
#   detectCores()
#   this_cluster <- makeCluster(4)
# 
#   t1 <- Sys.time()
#   chain_output <- parLapply(cl = this_cluster, X = 1:4,
#                             fun = run_MCMC_parallel,
#                              dat = Y, niter = 200, nburnin=0, thin=10)
#   t2 <- Sys.time()
#   print(t2-t1)
# 
# 
#   stopCluster(this_cluster)
# 
#   para.name <- colnames(chain_output[[1]]$samples)
#   rhat.seq <- c()
#   ess.seq <- c()
#   for (name in para.name){
#     post.sp <- cbind(chain_output[[1]]$samples[,name],
#                      chain_output[[2]]$samples[,name],
#                      chain_output[[3]]$samples[,name])
#     rhat.seq <- c(rhat.seq, rhat(post.sp))
#     ess.seq <- c(ess.seq, ess_basic(post.sp))
#   }
#   convg.stat <- data.frame(para.name,rhat.seq,ess.seq )
# 
# 
# 
#   samples.all <- rbind(chain_output[[1]]$samples,
#                        chain_output[[2]]$samples,
#                        chain_output[[3]]$samples)
#   # waic <- calculateWAIC(samples.all, BivExtMixmodel)
#   filename <- paste('Scenario1.1_seed',i,'.RData',sep='')
#   # save(Y, chain_output, convg.stat, file=file.path(dir.out, filename))
# }
# t4 <- Sys.time()
# print(t4-t3)
#36mins for 200 iter on euclid 01, 2.15 mins for 200 iters on RStudio server
#22mins for 100 iter on euclid 01, 1.45 mins for 100 iters on RStudio server


NumberOfCluster <- 15
cl <- makeCluster(NumberOfCluster)
registerDoSNOW(cl)

source(file.path(dir.work, 'Simulation/BEMM_Functions.R'))

t1 <- Sys.time()
chain_res <-
foreach(i = 1:5) %:%
  foreach(j = 1:3, .packages = c('nimble','mvtnorm','tmvtnorm')) %dopar%{
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
    p <- pmvnorm(lower=rep(0,2), upper=u.x, mean=mu, sigma=sigma, keepAttr = F)

    set.seed(seed)
    Y.tail<-sim.RevExpU.MGPD(n=n-floor(n*p),d=d, a=a, beta=beta, sig=sig, gamma=gamma, MGPD = T,std=T)

    # GP scale tail data combined with the bulk data
    set.seed(seed)
    Y.bulk <- rtmvnorm(floor(n*p), mean=mu, sigma=sigma, lower=c(0,0),upper=u.x)

    # The name of the dataset should be Y for further WAIC calculation.
    Y <- rbind(Y.bulk, sweep(Y.tail$X,2,u.x,"+"))
    plot(Y)
    run_MCMC_parallel(seed=j, dat=Y, niter=20000, nburnin = 10000, thin=10)
  }

save(chain_res,  file=file.path(dir.out, filename='Scenario_1.1_itr1_15.RData'))
t2 <- Sys.time()
print(t2-t1)


# 4 cores, 1:1,1:3, 100 iterations, 1.45 min  on RStudio server
# 4 cores, 1:2,1:2, 100 iterations, 1.48 min  on RStudio server
#30 cores, 1:10,1:3, 100 iterations, 10 mins on Euclid 01

# 4 cores, 1:1,1:1, 200 iterations, 1.86 min  on RStudio server
# 4 cores, 1:1,1:3, 200 iterations, 2.14 min  on RStudio server
# 4 cores, 1:2,1:2, 200 iterations, 2.19 min  on RStudio server


# 4 cores, 1:1,1:3, 200 iterations, 5 min  on Euclid 01
# 4 cores, 1:1,1:1, 200 iterations, 2.12 min  on Euclid 01
# 4 cores, 1:2,1:2, 200 iterations, 12 min on Euclid 01
# 4 cores, 1:2,1:2, 200 iterations, 4.5 min on Euclid 27 #CPU utilization is full
# 4 cores, 1:2,1:2, 200 iterations, 4.5 min on Euclid 28 #CPU utilization is full
# 4 cores, 1:2,1:2, 200 iterations, 25.5 min on Euclid 29 #CPU utilization is full

# 4 cores, 1:2,1:2, 200 iterations, 6.7 min on Euclid 30
# 7 cores, 1:2,1:3, 200 iterations, 7.7 min on Euclid 30
# 20 cores, 1;6,1:3, 200 iterations, 14.9 min on Euclid 30 #potentially not utilize all 20 cores because of full utilization.
# 20 cores, 1;6,1:3, 200 iterations, 15.5 min on Euclid 29 #potentially not utilize all 20 cores because of full utilization.

# 15 cores, 1;5,1:3, 200 iterations, 13.9 min on Euclid 29 
# 15 cores, 1;5,1:3, 200 iterations, 5.6 min on Euclid 28 
# 15 cores, 1;5,1:3, 200 iterations, 5.3 min on Euclid 27


#15 cores, 1:5, 1:3, 200 iterations. 13 mins on Euclid 01
#12 cores, 1:4, 1:3, 200 iterations. 11 mins on Euclid 01
#15 cores, 1:4, 1:3, 200 iterations. 11 mins on Euclid 01

#30 cores, 1:10,1:3, 200 iterations, 21 mins on Euclid 01

# full CPU utilization on Euclid 01
#https://stackoverflow.com/questions/24145693/r-foreach-not-using-multiple-cores