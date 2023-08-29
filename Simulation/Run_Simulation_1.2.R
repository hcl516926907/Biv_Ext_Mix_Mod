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
  foreach(i = c(251:300)) %:%
  foreach(j = 1:3, .packages = c('nimble','mvtnorm','tmvtnorm')) %dopar%{
    seed <- i
    d <- 2
    a <- c(0.5, 1.2)
    beta <- c(0, 0)
    sig <- c(0.5, 1.2)
    gamma <- c(0.2, -0.2)
    n <- 2000
    mu <- c(3.5, 4.0)
    sd1 <- 1
    sd2 <- 1.5
    rho <- 0.7
    sigma <- matrix(c(sd1^2, rho*sd1*sd2, rho*sd1*sd2, sd2^2),ncol=2)
    
    u.x <- c(5.5, 6.7)
    # u.x <- c(4.7,6)
    p <- pmvnorm(upper=u.x, mean=mu, sigma=sigma, keepAttr = F)
    
    set.seed(seed)
    Y.tail<-sim.RevExpU.MGPD(n=n-floor(n*p),d=d, a=a, beta=beta, sig=sig, gamma=gamma, MGPD = T,std=T)
    
    # GP scale tail data combined with the bulk data
    set.seed(seed)
    Y.bulk <- rtmvnorm(floor(n*p), mean=mu, sigma=sigma, upper=u.x)
    
    # The name of the dataset should be Y for further WAIC calculation.
    Y <- rbind(Y.bulk, sweep(Y.tail$X,2,u.x,"+"))
    run_MCMC_parallel(seed=j, dat=Y, niter=30000, nburnin = 20000, thin=10)
  }
stopCluster(cl)

# save(chain_res,  file=file.path(dir.out, filename='Scenario_1.2_itr101_150_lamfix.RData'))
# save(chain_res,  file=file.path(dir.out, filename='Scenario_1.2_itr151_200_lamfix.RData'))
# save(chain_res,  file=file.path(dir.out, filename='Scenario_1.2_itr201_250_lamfix.RData'))
save(chain_res,  file=file.path(dir.out, filename='Scenario_1.2_itr251_300_lamfix.RData'))
t2 <- Sys.time()
print(t2-t1)

