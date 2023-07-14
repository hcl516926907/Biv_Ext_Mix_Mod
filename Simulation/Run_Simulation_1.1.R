dir.work <- '/home/pgrad2/2448355h/My_PhD_Project/Biv_Ext_Mix_Mod'
dir.out <- "/home/pgrad2/2448355h/My_PhD_Project/01_Output/Biv_Ext_Mix_Mod/Simulation"
source(file.path(dir.work, "KRSW/RevExp_U_Functions.r"))
source(file.path(dir.work, "KRSW/CommonFunctions.r"))



install.packages("pacman",repos = "http://cran.us.r-project.org")

pacman::p_load(nimble, mvtnorm, parallel,tmvtnorm,posterior,foreach,doSNOW)

# t3 <- Sys.time()
# for (i in 1:10) {
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
#                              dat = Y, niter = 15000, nburnin=10000, thin=5)
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
#   save(Y, chain_output, convg.stat, file=file.path(dir.out, filename))
# }
# t4 <- Sys.time()
# print(t4-t3)



NumberOfCluster <- 10
cl <- makeCluster(NumberOfCluster) 
registerDoSNOW(cl) 

source(file.path(dir.work, 'Simulation/BEMM_Functions.R'))

t1 <- Sys.time()
chain_res <- 
foreach(i = 1:10) %:% 
  foreach(j = 1:3, .packages = c('nimble','mvtnorm','tmvtnorm')) %dopar%{
    seed <- i
    d <- 2
    a <- c(1.5, 2)
    beta <- c(0.5, 0)
    sig <- c(0.5, 0.4)
    gamma <- c(0.2, 0.1)
    n <- 2000
    mu <- c(3.5, 4.5)
    sd1 <- 1.22
    sd2 <- 1.10
    rho <- 0.72
    sigma <- matrix(c(sd1^2, rho*sd1*sd2, rho*sd1*sd2, sd2^2),ncol=2)
    
    u.x <- c(5.6, 7)
    p <- pmvnorm(lower=rep(0,2), upper=u.x, mean=mu, sigma=sigma, keepAttr = F)
    
    set.seed(seed)
    Y.tail<-sim.RevExpU.MGPD(n=n-floor(n*p),d=d, a=a, beta=beta, sig=sig, gamma=gamma, MGPD = T,std=T)
    
    # GP scale tail data combined with the bulk data
    set.seed(seed)
    Y.bulk <- rtmvnorm(floor(n*p), mean=mu, sigma=sigma, lower=c(0,0),upper=u.x)
    
    # The name of the dataset should be Y for further WAIC calculation.
    Y <- rbind(Y.bulk, sweep(Y.tail$X,2,u.x,"+"))
    run_MCMC_parallel(seed=j, dat=Y, niter=20000, nburnin = 10000, thin=10)
  }

save(chain_res,  file=file.path(dir.out, filename='Scenario_1.1_itr1_30.RData'))
t2 <- Sys.time()
print(t2-t1)
