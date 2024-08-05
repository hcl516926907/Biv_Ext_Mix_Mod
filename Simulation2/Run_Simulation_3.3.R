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
packages <- c("nimble", "mvtnorm", "tmvtnorm","foreach","doSNOW","parallel")  


load_install_packages(packages)

print(detectCores())


#########################################
rep.size <- 1000
set.seed(1234)

a <-  1.5
sig <- c(1.2,1.5)
gamma <- c(0.1,0.2)

params <- c(-0.4, 0, 2, 0.5, 2)
thres <- c(2.5,3.5)

lbound.tail <- c(-Inf,-Inf)
ubound.tail <- c(Inf,Inf)

total.dist.name <- c(1,1,1,2)

bulk.dist.name <- total.dist.name[1:3]
cop.name.map <- c("normal","clayton","gumbel","frank","joe","plackett")

cop <- switch(cop.name.map[bulk.dist.name[1]],
              "normal" = normalCopula(param = params[1], dim = 2),
              "clayton" = claytonCopula(param = params[1], dim = 2),
              "gumbel" = gumbelCopula(param = params[1], dim = 2),
              "frank" = frankCopula(param = params[1], dim = 2),
              "joe" = joeCopula(param = params[1], dim = 2),
              "plackett" = plackettCopula(param = params[1]),
              stop("Unsupported copula name"))

margin.name.map <- c("norm","exp","gamma","lnorm","weibull",'lst')
margins.name <- margin.name.map[bulk.dist.name[2:3]]
params.margin <- params[-1]
transformed_params <- list()
start.idx <- 1
for (i in 1:length(margins.name)) {
  m <- margins.name[i]
  
  if (m == "norm") {
    n.param <- 2
    p <- params.margin[ start.idx : (start.idx+1)]
    transformed_params[[i]] <- list(mean = p[1], sd = p[2]) # ensure sd > 0
    start.idx <- start.idx + n.param
    
  } else if (m == "exp") {
    n.param <- 1
    p <- params.margin[ start.idx]
    transformed_params[[i]] <- list(rate = p[1]) # ensure rate > 0
    start.idx <- start.idx + n.param
    
  } else if (m == "gamma") {
    n.param <- 2
    p <- params.margin[ start.idx : (start.idx+1)]
    transformed_params[[i]] <- list(shape = p[1], rate = p[2]) # ensure shape, rate > 0
    start.idx <- start.idx + n.param
    
  } else if (m == "lnorm") {
    n.param <- 2
    p <- params.margin[ start.idx : (start.idx+1)]
    transformed_params[[i]] <- list(meanlog = p[1], sdlog = p[2]) # ensure sdlog > 0
    start.idx <- start.idx + n.param
  } else if (m == "weibull") {
    n.param <- 2
    p <- params.margin[ start.idx : (start.idx+1)]
    transformed_params[[i]] <- list(shape = p[1], scale = p[2]) # ensure shape, scale > 0
    start.idx <- start.idx + n.param
  } else if (m == 'lst'){
    n.param <- 3
    p <- params.margin[ start.idx : (start.idx+2)]
    transformed_params[[i]] <- list(df = p[1], mu = p[2], sigma = p[3]) 
    start.idx <- start.idx + n.param
  } else {
    stop("Unsupported marginal distribution")
  }
}


bulk.dist <- mvdc(cop, margins = margins.name, 
                  
                  paramMargins=list(transformed_params[[1]],
                                    
                                    transformed_params[[2]])
)


p <- pMvdc(thres, bulk.dist)
print(p)

n.Y.bulk <- 0
Y.bulk.ut <- matrix(, nrow = 0, ncol = 2)
while (n.Y.bulk< (floor(rep.size*p) )){
  Y.bulk <- rMvdc(floor(rep.size*p), bulk.dist)
  Y.bulk <- Y.bulk[Y.bulk[,1]<thres[1] & Y.bulk[,2]<thres[2],]
  n.Y.bulk <- n.Y.bulk + nrow(Y.bulk)
  Y.bulk.ut <- rbind(Y.bulk.ut, Y.bulk)
}
Y.bulk <- Y.bulk.ut[sample(1:n.Y.bulk, floor(rep.size*p), replace=TRUE),]

n.Y.tail <- 0
Y.tail.ut <- matrix(, nrow = 0, ncol = 2)
while (n.Y.tail< (rep.size-floor(rep.size*p) )){
  if (total.dist.name[4]==1){
    Y.tail <- sim.RevExpU.MGPD(n=rep.size, d=2, a=rep(a,2), beta=c(0,0), sig=sig, gamma=gamma, MGPD = T,std=F)
  }else{
    Y.tail<-sim.GumbelU.MGPD(n=rep.size, d=2, a=a, beta=c(0,0), sig=sig, gamma=gamma, MGPD = T,std=F)
  }
  
  
  bound.cond.1 <- Y.tail[,1] > (lbound.tail-thres)[1] & Y.tail[,1] < (ubound.tail-thres)[1]
  bound.cond.2 <- Y.tail[,2] > (lbound.tail-thres)[2] & Y.tail[,2] < (ubound.tail-thres)[2]
  Y.tail <- Y.tail[bound.cond.1&bound.cond.2,]
  n.Y.tail <- n.Y.tail+ nrow(Y.tail)
  Y.tail.ut <- rbind(Y.tail.ut, Y.tail)
}

Y.tail <- Y.tail.ut[sample(1:n.Y.tail, rep.size-floor(rep.size*p), replace=TRUE),]


# The name of the dataset should be Y for further WAIC calculation.
Y <- rbind(Y.bulk, sweep(Y.tail,2,thres,"+"))
dim(Y)
plot(Y)
#########################################



######################################################################################
# Parallelled code for running the simulation multiple times
# The reason of dividing the code into two parts instead of one is because large number 
# iterations in the foreach cound cause some unexpected error.
######################################################################################
NumberOfCluster <- 40
cl <- makeCluster(NumberOfCluster)
registerDoSNOW(cl)

# source(file.path(dir.work, 'Simulation/BEMM_Functions.R'))

# itr <- 500
# itr <- 600
# itr <- 700
# itr <- 800
itr <- 900

t1 <- Sys.time()
chain_res <-
  foreach(i = (itr+1):(itr+50)) %:%
  foreach(j = 1:3, .packages = c('nimble','mvtnorm','tmvtnorm')) %dopar%{
    seed <- i
    d <- 2
    a <- c(0.5, 1.2)
    beta <- c(0, 0)
    sig <- c(0.5, 1.2)
    gamma <- c(0.3, 0.1)
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
# save(chain_res,file=file.path(dir.out, filename='Scenario_1.1_seed_1_9_p84.RData'))
# save(chain_res,  file=file.path(dir.out, filename='Scenario_1.1_itr1_20_lamfix.RData'))
# save(chain_res,  file=file.path(dir.out, filename='Scenario_1.1_itr21_100_lamfix.RData'))
# save(chain_res,  file=file.path(dir.out, filename='Scenario_1.1_itr101_150_lamfix.RData'))
# save(chain_res,  file=file.path(dir.out, filename='Scenario_1.1_itr151_200_lamfix.RData'))
# save(chain_res,  file=file.path(dir.out, filename='Scenario_1.1_itr201_250_lamfix.RData'))
# save(chain_res,  file=file.path(dir.out, filename='Scenario_1.1_itr251_300_lamfix.RData'))
# save(chain_res,  file=file.path(dir.out, filename='Scenario_1.1_itr301_350_lamfix.RData'))
# save(chain_res,  file=file.path(dir.out, filename='Scenario_1.1_itr351_400_lamfix.RData'))
# save(chain_res,  file=file.path(dir.out, filename='Scenario_1.1_itr401_450_lamfix.RData'))
# save(chain_res,  file=file.path(dir.out, filename='Scenario_1.1_itr451_500_lamfix.RData'))
save(chain_res,  file=file.path(dir.out, filename=paste('Scenario_1.1_itr',itr+1,'_',itr+50,'_lamfix.RData',sep='')))
t2 <- Sys.time()
print(t2-t1)


cl <- makeCluster(NumberOfCluster)
registerDoSNOW(cl)
t1 <- Sys.time()
chain_res <-
  foreach(i = (itr+51):(itr+100)) %:%
  foreach(j = 1:3, .packages = c('nimble','mvtnorm','tmvtnorm')) %dopar%{
    seed <- i
    d <- 2
    a <- c(0.5, 1.2)
    beta <- c(0, 0)
    sig <- c(0.5, 1.2)
    gamma <- c(0.3, 0.1)
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
# save(chain_res,file=file.path(dir.out, filename='Scenario_1.1_seed_1_9_p84.RData'))
# save(chain_res,  file=file.path(dir.out, filename='Scenario_1.1_itr1_20_lamfix.RData'))
# save(chain_res,  file=file.path(dir.out, filename='Scenario_1.1_itr21_100_lamfix.RData'))
# save(chain_res,  file=file.path(dir.out, filename='Scenario_1.1_itr101_150_lamfix.RData'))
# save(chain_res,  file=file.path(dir.out, filename='Scenario_1.1_itr151_200_lamfix.RData'))
# save(chain_res,  file=file.path(dir.out, filename='Scenario_1.1_itr201_250_lamfix.RData'))
# save(chain_res,  file=file.path(dir.out, filename='Scenario_1.1_itr251_300_lamfix.RData'))
# save(chain_res,  file=file.path(dir.out, filename='Scenario_1.1_itr301_350_lamfix.RData'))
# save(chain_res,  file=file.path(dir.out, filename='Scenario_1.1_itr351_400_lamfix.RData'))
# save(chain_res,  file=file.path(dir.out, filename='Scenario_1.1_itr401_450_lamfix.RData'))
# save(chain_res,  file=file.path(dir.out, filename='Scenario_1.1_itr451_500_lamfix.RData'))
save(chain_res,  file=file.path(dir.out, filename=paste('Scenario_1.1_itr',itr+51,'_',itr+100,'_lamfix.RData',sep='')))
t2 <- Sys.time()

