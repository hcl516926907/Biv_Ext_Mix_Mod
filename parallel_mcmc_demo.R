library(adaptMCMC)
library(posterior)
library(RColorBrewer)

seed(1234)
y <- rnorm(1000,4,1)

p.log <- function(theta){
  ll <- sum(dnorm(y, mean=theta[1], sd=theta[2], log=T)) + 
  dunif(theta[1], -10, 10, log=T) + 
  dunif(theta[1], 0, 10, log=T)
  return(ll)
}

t1 <- Sys.time()
itr <-  500000
burnin <- 500
res <- MCMC.parallel(p.log, n=itr,
                     init=c(runif(1,0,5),runif(1,0,5)),n.chain=3,n.cpu=6,
                     scale=c(0.001,0.001),adapt=FALSE)
t2 <- Sys.time()

print(t2-t1)

sims.a1 <- cbind(res[[1]]$samples[burnin:itr,1],res[[2]]$samples[burnin:itr,1],
                 res[[3]]$samples[burnin:itr,1])
sims.a2 <- cbind(res[[1]]$samples[burnin:itr,2],res[[2]]$samples[burnin:itr,2],
                 res[[3]]$samples[burnin:itr,2])

my_traceplot <- function(samples, para.name){
  n.chain <- ncol(samples)
  n.sample <- nrow(samples)
  # only 8 colors in Dark2
  palette <- brewer.pal(n = n.chain, name = "Dark2")
  for (i in 1:n.chain){
    if (i==1){
      plot(1:n.sample, samples[,1], type='l', xlab='iterations', ylab='sample',
           main=paste('Traceplot of',para.name), col=palette[i])
    }else{
      lines(1:n.sample,samples[,i], type='l', col=palette[i])
    }
    
  }
}

my_traceplot(sims.a1,'mu')
my_traceplot(sims.a2,'sigma')

rhat(sims.a1)
ess_basic(sims.a1)
rhat(sims.a2)
ess_basic(sims.a2)


library(doParallel) 

registerDoParallel(cl <- makeCluster(3))
t3 <- Sys.time() 
result <- foreach(x = 1:3, .packages="adaptMCMC") %dopar% {
  set.seed(x)
  MCMC.parallel(p.log, n=itr,
                init=c(runif(1,0,5),runif(1,0,5)),n.chain=1,n.cpu=3,
                scale=c(0.001,0.001),adapt=FALSE)
}
t4 <- Sys.time() 
stopCluster(cl)

print(t4-t3)

sims.a1 <- cbind(result[[1]][[1]]$samples[burnin:itr,1],result[[2]][[1]]$samples[burnin:itr,1],
                 result[[3]][[1]]$samples[burnin:itr,1])
sims.a2 <- cbind(result[[1]][[1]]$samples[burnin:itr,2],result[[2]][[1]]$samples[burnin:itr,2],
                 result[[3]][[1]]$samples[burnin:itr,2])

my_traceplot(sims.a1,'mu')
my_traceplot(sims.a2,'sigma')

foreach(param1 = tuning_grid$param1, param2 = tuning_grid$param2) %dopar% {
  sink("my_mcmc.log", append=TRUE) # Write to log file
  
  # Your individual MCMC wrapped inside a function
  f_mcmc(data = my_data, max_iter = 20000,
         param1=param1, param2=param2)
}

