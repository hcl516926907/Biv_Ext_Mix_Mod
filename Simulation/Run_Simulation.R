source("KRSW/RevExp_U_Functions.r")
source("KRSW/CommonFunctions.r")
library(nimble, warn.conflicts = F)
library(mvtnorm)
library(parallel)
library(tmvtnorm)
library(posterior)

dir.out <- "/home/pgrad2/2448355h/My_PhD_Project/01_Output/Biv_Ext_Mix_Mod/Simulation"

seed <- 1234
d <- 2
a <- c(1.5, 2)
beta <- c(0.5, 0)
sig <- c(0.5, 0.4)
gamma <- c(-0.2, -0.1)
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

plot(Y.tail$X)
# GP scale tail data combined with the bulk data
set.seed(seed)
Y.bulk <- rtmvnorm(floor(n*p), mean=mu, sigma=sigma, lower=c(0,0),upper=u.x)

# The name of the dataset should be Y for further WAIC calculation.
Y <- rbind(Y.bulk, sweep(Y.tail$X,2,u.x,"+"))
plot(Y)


source('Simulation/BEMM_Functions.R')
detectCores()
this_cluster <- makeCluster(3)


t1 <- Sys.time()
chain_output <- parLapply(cl = this_cluster, X = 1:3, 
                          fun = run_MCMC_parallel, 
                          dat = Y, niter = 20000, nburnin=10000, thin=5)
t2 <- Sys.time()

print(t2-t1)


stopCluster(this_cluster)

para.name <- colnames(chain_output[[1]]$samples)
rhat.seq <- c()
ess.seq <- c()
for (name in para.name){
  post.sp <- cbind(chain_output[[1]]$samples[,name],
                   chain_output[[2]]$samples[,name],
                   chain_output[[3]]$samples[,name])
  rhat.seq <- c(rhat.seq, rhat(post.sp))
  ess.seq <- c(ess.seq, ess_basic(post.sp))
}
convg.stat <- data.frame(para.name,rhat.seq,ess.seq )



samples.all <- rbind(chain_output[[1]]$samples,
                     chain_output[[2]]$samples,
                     chain_output[[3]]$samples)
waic <- calculateWAIC(samples.all, BivExtMixmodel)

save(chain_output, convg.stat, waic, file=file.path(dir.out, 'Scenario1.3.RData'))

plot(chain_output[[1]]$samples[,'theta[6]'],type='l')
plot(density(chain_output[[1]]$samples[,'theta[6]']))
plot(density(samples.all[,'theta[6]']))
