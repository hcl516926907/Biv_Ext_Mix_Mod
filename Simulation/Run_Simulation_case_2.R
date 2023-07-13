source("KRSW/RevExp_U_Functions.r")
source("KRSW/CommonFunctions.r")
library(nimble, warn.conflicts = F)
library(mvtnorm)
library(parallel)
library(tmvtnorm)
library(posterior)
library(copula)

dir.out <- "/home/pgrad2/2448355h/My_PhD_Project/01_Output/Biv_Ext_Mix_Mod/Simulation"

#----------------------------------Gaussian Bulk-----------------------------------------
# #1234 has some issues
# #1235 works
seed <- 1234

n <- 2000
mu <- c(3.5, 4.5)
sd1 <- 1.22
sd2 <- 1.10
rho <- 0.72
sigma <- matrix(c(sd1^2, rho*sd1*sd2, rho*sd1*sd2, sd2^2),ncol=2)

set.seed(seed)
Y <- rtmvnorm(n, mean=mu, sigma=sigma, lower=c(0,0))

# The name of the dataset should be Y for further WAIC calculation.
plot(Y)

#----------------------------------Gumbel Bulk-----------------------------------------
# seed <- 1234
# 
# n <- 2000
# mu <- c(3.5, 4.5)
# sd1 <- 1.22
# sd2 <- 1.10
# rho <- 0.72
# sigma <- matrix(c(sd1^2, rho*sd1*sd2, rho*sd1*sd2, sd2^2),ncol=2)
# 
# set.seed(seed)
# Y <- rtmvnorm(n, mean=mu, sigma=sigma, lower=c(0,0))

#---------------------------------------------------------------------------------------
source('Simulation/BEMM_Functions.R')
detectCores()
this_cluster <- makeCluster(3)


t1 <- Sys.time()
chain_output <- parLapply(cl = this_cluster, X = 1:3, 
                          fun = run_MCMC_parallel, 
                          dat = Y, niter = 110000, nburnin=10000, thin=10)
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
# waic <- calculateWAIC(samples.all, BivExtMixmodel)

# save(chain_output, convg.stat, file=file.path(dir.out, 'Scenario2.1.RData')) #seed 1234 org
# save(chain_output, convg.stat, file=file.path(dir.out, 'Scenario2.1_1234_AF_slice.RData')) #seed 1234
# save(chain_output, convg.stat, file=file.path(dir.out, 'Scenario2.1_1235.RData')) #seed 1235
# load(file=file.path(dir.out, 'Scenario1.3.RData'))
# save(chain_output, convg.stat, file=file.path(dir.out, 'Scenario2.1_1234_0.99upper.RData'))
# save(chain_output, convg.stat, file=file.path(dir.out, 'Scenario2.1_1234_0.999upper.RData'))

plot(chain_output[[3]]$samples[1:10000,'thres[2]'],type='l')
plot(density(chain_output[[1]]$samples[,'thres[2]']))
plot(density(samples.all[,'theta[6]']))

pairs(chain_output[[1]]$samples[4000:6000,c('theta[1]','thres[2]','theta[4]','theta[6]')])


# load(file=file.path(dir.out, 'Scenario2.1_1234_0.99upper.RData'))
# samples.all <- rbind(chain_output[[1]]$samples,
#                      chain_output[[2]]$samples,
#                      chain_output[[3]]$samples)
# plot(samples.all[,'thres[1]'],type='l',main='Traceplot of thres1')
# abline(h=4.5 ,col='red')