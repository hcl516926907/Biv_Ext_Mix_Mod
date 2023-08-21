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
packages <- c("nimble", "mvtnorm", "tmvtnorm","foreach","doSNOW","parallel")  


load_install_packages(packages)

print(detectCores())

NumberOfCluster <- 40
cl <- makeCluster(NumberOfCluster)
registerDoSNOW(cl)

# functions used in the parallelization is in BEMM_Functions.R
source(file.path(dir.work, 'Simulation/BEMM_Functions.R'))

#Use foreach to carry out a nested loop parallelization
t1 <- Sys.time()
chain_res <-
  foreach(i = 1:10) %:%
  foreach(j = 1:3, .packages = c('nimble','mvtnorm','tmvtnorm')) %dopar%{

    # Generate some data
    seed <- i
    d <- 2
    a <- c(1, 10)
    beta <- c(0, 0)
    sig <- c(0.5, 1.2)
    gamma <- c(0.3, 0.1)
    n <- 2000
    mu <- c(3.5, 4.0)
    sd1 <- 1
    sd2 <- 3
    rho <- 0.7
    sigma <- matrix(c(sd1^2, rho*sd1*sd2, rho*sd1*sd2, sd2^2),ncol=2)
    u.x <- c(3.5, 6.7)
    p <- pmvnorm(upper=u.x, mean=mu, sigma=sigma, keepAttr = F)
    set.seed(seed)
    Y.tail<-sim.RevExpU.MGPD(n=n-floor(n*p),d=d, a=a, beta=beta, sig=sig, gamma=gamma, MGPD = T,std=T)
    set.seed(seed)
    Y.bulk <- rtmvnorm(floor(n*p), mean=mu, sigma=sigma, upper=u.x)
    Y <- rbind(Y.bulk, sweep(Y.tail$X,2,u.x,"+"))

    # function run_MCMC_parallel is called in each parallelled loop
    run_MCMC_parallel(seed=j, dat=Y, niter=200, nburnin = 0, thin=10)
  }
stopCluster(cl)

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
# 15 cores, 1;5,1:3, 200 iterations, 12.3 min on Euclid 30

#15 cores, 1:5, 1:3, 200 iterations. 13 mins on Euclid 01
#12 cores, 1:4, 1:3, 200 iterations. 11 mins on Euclid 01
#15 cores, 1:4, 1:3, 200 iterations. 11 mins on Euclid 01

#30 cores, 1:10,1:3, 200 iterations, 21 mins on Euclid 01

# full CPU utilization on Euclid 01
#https://stackoverflow.com/questions/24145693/r-foreach-not-using-multiple-cores