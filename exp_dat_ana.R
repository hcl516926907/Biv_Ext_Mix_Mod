
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
qecdf <- function(x){
  x.sort <- sort(x)
  # ecdf without consideration of repeated observations
  pse.ecdf <- 1:length(x)/length(x)
  
  # Apply ecdf(X) to a random variable X
  # For repeated observations, ecdf remains constant is equal to the sequence's
  # last pse.ecdf
  quan <- function(y){pse.ecdf[tail(which(x.sort==y),1)]} 
  
  return(sapply(x,quan))
}

seq.test <- c(5,6,5,6,7)
qecdf(seq.test)
#0.4 0.8 0.4 0.8 1.0

quantile(seq.test,type=1)
#0%  25%  50%  75% 100% 
#5    5    6    6    7 


riv1.unif <- qecdf(riv1)
riv2.unif <- qecdf(riv2)

# plot transfromed into standard uniform distribution
plot(riv1.unif, riv2.unif, main = "Weekly River Flow, mariginally uniform",
     xlab = "Fajardo", ylab = "Espiritu Santu",
     pch = 19, frame = FALSE)



u.seq <- seq(0,1,0.01)
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
  
  p.0.0 <- length(intersect(x.0,y.0))/length(y.0)
  p.0.1 <- length(intersect(x.0,y.1))/length(y.1)
  p.1.0 <- length(intersect(x.1,y.0))/length(y.0)
  p.1.1 <- length(intersect(x.1,y.1))/length(y.1)
  
  res <- c(x.p.0, x.p.1, y.p.0, y.p.1,
           p.0.0, p.0.1, p.1.0, p.1.1)
  names(res) <- c('x.p.0', 'x.p.1', 'y.p.0', 'y.p.1',
                  'p.0.0', 'p.0.1', 'p.1.0', 'p.1.1')
  return(res)
}

tail.coef(riv1.unif,riv2.unif,u)


coef.tab <- sapply(u.seq, tail.coef, x=riv1.unif, y=riv2.unif)

chi <- coef.tab['p.1.1',]
plot(u.seq, chi, type='l', main="Chi Plot", xlab='u', ylab='chi')

chi.bar <- 2*log(coef.tab['x.p.1',])/(log(coef.tab['p.1.1',]*coef.tab['y.p.1',]))-1
plot(u.seq, chi.bar, type='l', main="Chi Bar Plot", xlab='u', ylab='chi bar')

ddd
