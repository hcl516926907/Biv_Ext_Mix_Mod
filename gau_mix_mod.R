library(mixtools)

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




set.seed(100)
x.1 <- rmvnorm(40, c(0, 0))
x.2 <- rmvnorm(60, c(3, 4))
plot(X.1[,1],X.1[,2])
X.1 <- rbind(x.1, x.2)
mu <- list(c(0, 0), c(3, 4))
out.1 <- mvnormalmixEM(X.1, arbvar = FALSE, mu = mu,
                       epsilon = 1e-02)
out.1[2:5]
##Fitting randomly generated data with a 2-component scale mixture of bivariate normals.
x.3 <- rmvnorm(40, c(0, 0), sigma =
                 matrix(c(200, 1, 1, 150), 2, 2))
x.4 <- rmvnorm(60, c(0, 0))
X.2 <- rbind(x.3, x.4)
plot(X.2[,1],X.2[,2])
lambda <- c(0.40, 0.60)
sigma <- list(diag(1, 2), matrix(c(200, 1, 1, 150), 2, 2))
out.2 <- mvnormalmixEM(X.2, arbmean = FALSE,
                       sigma = sigma, lambda = lambda,
                       epsilon = 1e-02)
out.2[2:5]

summary(out.2)
plot(out.2)


