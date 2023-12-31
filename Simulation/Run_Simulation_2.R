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

itr <- 900

cl <- makeCluster(NumberOfCluster)
registerDoSNOW(cl)

source(file.path(dir.work, 'Simulation/BEMM_Functions_AFSlice.R'))

t1 <- Sys.time()
chain_res <-
  foreach(i = (itr+1):(itr+25)) %:%
  foreach(j = 1:3, .packages = c('nimble','mvtnorm','tmvtnorm')) %dopar%{
    seed <- i
    d <- 2
    mu <- c(3.5, 4.0)
    sd1 <- 1
    sd2 <- 1.5
    rho <- 0.7
    sigma <- matrix(c(sd1^2, rho*sd1*sd2, rho*sd1*sd2, sd2^2),ncol=2)
    n <- 2000
    set.seed(seed)
    Y <- rmvnorm(n, mean=mu, sigma=sigma)

    run_MCMC_parallel(seed=j, dat=Y, nburnin=20000, niter=30000, thin=10)
  }
stopCluster(cl)

# save(chain_res,  file=file.path(dir.out, filename='Scenario_2_itr1_20_lamfix.RData'))
# save(chain_res,  file=file.path(dir.out, filename='Scenario_2_lamfix_APT.RData'))
# save(chain_res,  file=file.path(dir.out, filename='Scenario_2_itr1_20_lamfix_AFSlice.RData'))
# save(chain_res,  file=file.path(dir.out, filename='Scenario_2_itr21_50_lamfix_AFSlice.RData'))
# save(chain_res,  file=file.path(dir.out, filename='Scenario_2_itr51_75_lamfix_AFSlice.RData'))
# save(chain_res,  file=file.path(dir.out, filename='Scenario_2_itr76_100_lamfix_AFSlice.RData'))
# save(chain_res,  file=file.path(dir.out, filename='Scenario_2_itr101_125_lamfix_AFSlice.RData'))
# save(chain_res,  file=file.path(dir.out, filename='Scenario_2_itr126_150_lamfix_AFSlice.RData'))
# save(chain_res,  file=file.path(dir.out, filename='Scenario_2_itr151_175_lamfix_AFSlice.RData'))
# save(chain_res,  file=file.path(dir.out, filename='Scenario_2_itr176_200_lamfix_AFSlice.RData'))
# save(chain_res,  file=file.path(dir.out, filename='Scenario_2_itr201_225_lamfix_AFSlice.RData'))
# save(chain_res,  file=file.path(dir.out, filename='Scenario_2_itr226_250_lamfix_AFSlice.RData'))
# save(chain_res,  file=file.path(dir.out, filename='Scenario_2_itr251_275_lamfix_AFSlice.RData'))
# save(chain_res,  file=file.path(dir.out, filename='Scenario_2_itr276_300_lamfix_AFSlice.RData'))
# save(chain_res,  file=file.path(dir.out, filename='Scenario_2_itr301_325_lamfix_AFSlice.RData'))
# save(chain_res,  file=file.path(dir.out, filename='Scenario_2_itr326_350_lamfix_AFSlice.RData'))
# save(chain_res,  file=file.path(dir.out, filename='Scenario_2_itr351_375_lamfix_AFSlice.RData'))
# save(chain_res,  file=file.path(dir.out, filename='Scenario_2_itr376_400_lamfix_AFSlice.RData'))
# save(chain_res,  file=file.path(dir.out, filename='Scenario_2_itr401_425_lamfix_AFSlice.RData'))
# save(chain_res,  file=file.path(dir.out, filename='Scenario_2_itr426_450_lamfix_AFSlice.RData'))
# save(chain_res,  file=file.path(dir.out, filename='Scenario_2_itr451_475_lamfix_AFSlice.RData'))
# save(chain_res,  file=file.path(dir.out, filename='Scenario_2_itr476_500_lamfix_AFSlice.RData'))
save(chain_res,  file=file.path(dir.out, filename=paste('Scenario_2_itr',itr+1,'_',itr+25,'_lamfix.RData',sep='')))
t2 <- Sys.time()
print(t2-t1)

#itr=8 could be an example to show the power of AFSlice sampler.


cl <- makeCluster(NumberOfCluster)
registerDoSNOW(cl)

t3 <- Sys.time()
chain_res <-
  foreach(i = (itr+26):(itr+50)) %:%
  foreach(j = 1:3, .packages = c('nimble','mvtnorm','tmvtnorm')) %dopar%{
    seed <- i
    d <- 2
    mu <- c(3.5, 4.0)
    sd1 <- 1
    sd2 <- 1.5
    rho <- 0.7
    sigma <- matrix(c(sd1^2, rho*sd1*sd2, rho*sd1*sd2, sd2^2),ncol=2)
    n <- 2000
    set.seed(seed)
    Y <- rmvnorm(n, mean=mu, sigma=sigma)
    
    run_MCMC_parallel(seed=j, dat=Y, nburnin=20000, niter=30000, thin=10)
  }
stopCluster(cl)
save(chain_res,  file=file.path(dir.out, filename=paste('Scenario_2_itr',itr+26,'_',itr+50,'_lamfix.RData',sep='')))
t4 <- Sys.time()
print(t4-t3)


cl <- makeCluster(NumberOfCluster)
registerDoSNOW(cl)


t5 <- Sys.time()
chain_res <-
  foreach(i = (itr+51):(itr+75)) %:%
  foreach(j = 1:3, .packages = c('nimble','mvtnorm','tmvtnorm')) %dopar%{
    seed <- i
    d <- 2
    mu <- c(3.5, 4.0)
    sd1 <- 1
    sd2 <- 1.5
    rho <- 0.7
    sigma <- matrix(c(sd1^2, rho*sd1*sd2, rho*sd1*sd2, sd2^2),ncol=2)
    n <- 2000
    set.seed(seed)
    Y <- rmvnorm(n, mean=mu, sigma=sigma)
    
    run_MCMC_parallel(seed=j, dat=Y, nburnin=20000, niter=30000, thin=10)
  }
stopCluster(cl)


save(chain_res,  file=file.path(dir.out, filename=paste('Scenario_2_itr',itr+51,'_',itr+75,'_lamfix.RData',sep='')))
t6 <- Sys.time()
print(t6-t5)



cl <- makeCluster(NumberOfCluster)
registerDoSNOW(cl)


t7 <- Sys.time()
chain_res <-
  foreach(i = (itr+76):(itr+100)) %:%
  foreach(j = 1:3, .packages = c('nimble','mvtnorm','tmvtnorm')) %dopar%{
    seed <- i
    d <- 2
    mu <- c(3.5, 4.0)
    sd1 <- 1
    sd2 <- 1.5
    rho <- 0.7
    sigma <- matrix(c(sd1^2, rho*sd1*sd2, rho*sd1*sd2, sd2^2),ncol=2)
    n <- 2000
    set.seed(seed)
    Y <- rmvnorm(n, mean=mu, sigma=sigma)
    
    run_MCMC_parallel(seed=j, dat=Y, nburnin=20000, niter=30000, thin=10)
  }
stopCluster(cl)


save(chain_res,  file=file.path(dir.out, filename=paste('Scenario_2_itr',itr+76,'_',itr+100,'_lamfix.RData',sep='')))
t8 <- Sys.time()
print(t8-t7)
print(t8-t1)