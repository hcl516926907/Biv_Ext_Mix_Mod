
# data directory
dir.dat <- '/home/pgrad2/2448355h/My_PhD_Project/00_Dataset'

# name of Ruerto Rico rive dataset 
riv.dat.fit.nam <- 'Puerto Rico Fit Dataset.dat'
riv.dat.test.nam <- 'Puerto Rico Test Data.dat'

# name of Leeds air pollution dataset
pol.dat.fit.nam <- 'Leeds Fit Dataset.dat'
pol.dat.test.nam <- 'Leeds Test Dataset.dat'

# load all datasets
riv.dat.fit <- read.table(file.path(dir.dat,riv.dat.fit.nam), header=TRUE)
riv.dat.test <- read.table(file.path(dir.dat,riv.dat.test.nam), header=TRUE)
riv.dat <- rbind(riv.dat.fit,riv.dat.test)

pol.dat.fit <- read.table(file.path(dir.dat,pol.dat.fit.nam), header=TRUE)
pol.dat.test <- read.table(file.path(dir.dat,pol.dat.test.nam), header=TRUE)
pol.dat <- rbind(pol.dat.fit,pol.dat.test)

# exploratory plots
riv1 <- riv.dat$V1
riv2 <- riv.dat$V2

# river1 is Fajardo and river2 is Espiritu Santu, by comparing the 
# plot with that in Leonelli's paper.

# plot in original scale
plot(riv1, riv2, main = "Weekly River Flow, Original Scale",
     xlab = "Fajardo", ylab = "Espiritu Santu",
     pch = 19, frame = FALSE)

# plot in log scale
plot(log(riv1), log(riv2), main = "Weekly River Flow, log Scale",
     xlab = "Fajardo", ylab = "Espiritu Santu",
     pch = 19, frame = FALSE)

# find the quantile of the ecdf, i.e. applying F(X) to X
my.ecdf <- function(x){
  #adjust to avoid having F(x)=1, according to p36 in cole's book.
  ecdf <- rank(x,ties.method='max')/(length(x)+1)
  return(ecdf)
}

seq.test <- c(5,6,5,6,7)
my.ecdf(seq.test)
#0.4 0.8 0.4 0.8 1.0

quantile(seq.test,type=1)
#0%  25%  50%  75% 100% 
#5    5    6    6    7 

riv1.unif <- my.ecdf(riv1)
riv2.unif <- my.ecdf(riv2)

# qqplot to compare two margins
qqplot(riv1.unif, riv2.unif)
# compare margins with a standard uniform distribution
qqplot(runif(length(riv1.unif)), riv1.unif)
qqplot(runif(length(riv2.unif)), riv2.unif)

# plot with margins transfromed into standard uniform distribution
plot(riv1.unif, riv2.unif, main = "Weekly River Flow, mariginally uniform",
     xlab = "Fajardo", ylab = "Espiritu Santu",
     pch = 19, frame = FALSE)

step.size <- 0.01
u.seq <- seq(step.size, 1-step.size, step.size)
u <- 0.99

tail.coef <- function(x,y,u){
  n <- length(x)
  
  x.0 <- which(x<=u)
  x.1 <- which(x>u)
  y.0 <- which(y<=u)
  y.1 <- which(y>u)
  
  x.p.0 <- length(x.0)/n
  x.p.1 <- length(x.1)/n
  y.p.0 <- length(y.0)/n
  y.p.1 <- length(y.1)/n
  
  #conditional distribution
  p.0c0 <- length(intersect(x.0,y.0))/length(y.0)
  p.0c1 <- length(intersect(x.0,y.1))/length(y.1)
  p.1c0 <- length(intersect(x.1,y.0))/length(y.0)
  p.1c1 <- length(intersect(x.1,y.1))/length(y.1)
  
  res <- c(x.p.0, x.p.1, y.p.0, y.p.1,
           p.0c0, p.0c1, p.1c0, p.1c1)
  
  names(res) <- c('x.p.0', 'x.p.1', 'y.p.0', 'y.p.1',
                  'p.0c0', 'p.0c1', 'p.1c0', 'p.1c1')
  return(res)
}

tail.coef(riv1.unif,riv2.unif,u)
coef.tab <- sapply(u.seq, tail.coef, x=riv1.unif, y=riv2.unif)

conf <- 0.95
p.0j0 <- coef.tab['p.0c0',]*coef.tab['y.p.0',]

chi <- 2 - log(p.0j0)/log(u.seq)
# detla method
chi.var <- p.0j0*(1-p.0j0)*(1/(p.0j0*log(u.seq)))^2/length(riv1.unif)
chi.sd <- sqrt(chi.var)
upp.bond <- chi + qnorm((1+conf)/2)*chi.sd
low.bond <- chi + qnorm((1-conf)/2)*chi.sd
plot(u.seq, chi, type='l', main="Chi Plot", xlab='u', ylab='chi', ylim=c(0,1),frame = FALSE)
# confidence interval
lines(u.seq, upp.bond, type = "l", lty = 2, pch = 18)
lines(u.seq, low.bond, type = "l", lty = 2, pch = 18)

p.1j1 <- coef.tab['p.1c1',]*coef.tab['y.p.1',]
chi.bar <- 2*log(1-u.seq)/(log(p.1j1))-1
plot(u.seq, chi.bar, type='l', main="Chi Bar Plot", xlab='u', ylab='chi bar',ylim=c(-1,1))


library(evd)
chiplot(cbind(riv1,riv2))

