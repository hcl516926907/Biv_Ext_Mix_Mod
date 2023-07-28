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

source(file.path(dir.work, 'Simulation/BEMM_Functions.R'))

t1 <- Sys.time()
chain_res <-
  foreach(i = c(1:20)) %:%
  foreach(j = 1:3, .packages = c('nimble','mvtnorm','tmvtnorm')) %dopar%{
    seed <- i
    d <- 2
    mu <- c(3.5, 4.0)
    sd1 <- 1
    sd2 <- 1.5
    rho <- 0.7
    sigma <- matrix(c(sd1^2, rho*sd1*sd2, rho*sd1*sd2, sd2^2),ncol=2)
    n <- 2000
    Y <- rmvnorm(n, mean=mu, sigma=sigma)

    run_MCMC_parallel(seed=j, dat=Y, niter=30000, nburnin = 20000, thin=10)
  }

save(chain_res,  file=file.path(dir.out, filename='Scenario_2_itr1_20_lamfix.RData'))
t2 <- Sys.time()
print(t2-t1)

