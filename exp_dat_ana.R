
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





