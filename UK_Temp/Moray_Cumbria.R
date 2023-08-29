dir.work <- '/home/pgrad2/2448355h/My_PhD_Project/Biv_Ext_Mix_Mod'
dir.data <- '/home/pgrad2/2448355h/My_PhD_Project/00_Dataset/UK_Temp'
dir.out <- "/home/pgrad2/2448355h/My_PhD_Project/01_Output/Biv_Ext_Mix_Mod/UK_Temp"

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

load(file=file.path(dir.data, "moray_cumbria.RData"))

Y1 <- data.frame(time=1:length(temp.monthly.max.all$air_temp_city1),Month=temp.monthly.max.all$Month, Temp=temp.monthly.max.all$air_temp_city1)
m1 <- lm(Temp ~ time+ Month, data=Y1)

Y2 <- data.frame(time=1:length(temp.monthly.max.all$air_temp_city2),Month=temp.monthly.max.all$Month, Temp=temp.monthly.max.all$air_temp_city2)
m2 <- lm(Temp ~ time+ Month, data=Y2)

Y <- cbind(Y1$Temp-predict(m1), Y2$Temp-predict(m2))

Y.fit <- Y

NumberOfCluster <- 3
cl <- makeCluster(NumberOfCluster)
registerDoSNOW(cl)

source(file.path(dir.work, 'UK_Temp/BEMM_Function_Temp.R'))


t1 <- Sys.time()
chain_res <-
  foreach(j = 1:3, .packages = c('nimble','mvtnorm','tmvtnorm')) %dopar%{
    
    run_MCMC_parallel(seed=j, dat=Y.fit, niter=20000, nburnin = 10000, thin=10)
  }
stopCluster(cl)
t2 <- Sys.time()
print(t2-t1)

save(chain_res, file=file.path(dir.out, filename='moray_cumbria_0.8_0.99.RData'))