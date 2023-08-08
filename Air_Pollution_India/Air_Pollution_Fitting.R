dir.work <- '/home/pgrad2/2448355h/My_PhD_Project/Biv_Ext_Mix_Mod'
dir.data <- '/home/pgrad2/2448355h/My_PhD_Project/00_Dataset/Air_Pollution_India'
dir.out <- "/home/pgrad2/2448355h/My_PhD_Project/01_Output/Biv_Ext_Mix_Mod/Air_Pollution_India"

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


NumberOfCluster <- 3
cl <- makeCluster(NumberOfCluster)
registerDoSNOW(cl)

source(file.path(dir.work, 'Air_Pollution_India/BEMM_Functions_Air_Pollution.R'))

t1 <- Sys.time()
chain_res <-
  foreach(j = 1:3, .packages = c('nimble','mvtnorm','tmvtnorm')) %dopar%{
    
    run_MCMC_parallel(seed=j, dat=pm2.5_no2, niter=30000, nburnin = 20000, thin=10)
  }
stopCluster(cl)
t2 <- Sys.time()
print(t2-t1)

save(chain_res, file=file.path(dir.out, filename='Air_pollution_mvtn_0.6_0.99_AFSlice.RData'))
