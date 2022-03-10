library(mixtools)
library(plotmm)
library(ggplot2)
library(mvtnorm)

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


#--------------------fit on original scale------------------------
mu.2com <- list(c(84,80),c(162,183))
sigma.2com <- list(diag(c(200,250)),diag(c(500,550)))

mod.org.2com <- mvnormalmixEM(riv.dat, mu = mu.2com,
              sigma = sigma.2com, epsilon = 1e-02)

mod.org.2com.dm <- mvnormalmixEM(riv.dat, sigma = sigma.2com, epsilon = 1e-02)
mod.org.2com.ds <- mvnormalmixEM(riv.dat, mu = mu.2com, epsilon = 1e-02)

mod.org.2com[2:5]

mod.org.3com <- mvnormalmixEM(riv.dat, k=3, epsilon = 1e-02)

vis.mm <- function(mod){
component_colors <- c("red", "blue", "green", "yellow", "orange", 
                      "purple", "darksalmon", "goldenrod2", "dodgerblue", "darkorange3", 
                      "burlywood4", "darkmagenta", "firebrick", "deeppink2", 
                      "darkseagreen1")
m <- mod
k <- length(m$lambda)
x <- data.frame(m$x)
mean <- m$mu
sigma <- m$sigma
X_1 <- x[, 1]
X_2 <- x[, 2]
post <- m$posterior
out_plot <- ggplot2::qplot(x = X_1, y = X_2) + 
            ggplot2::geom_point(colour = "darkgray", fill = "lightgray", size = 0.7) + 
            ggplot2::theme_minimal()

for (i in 1:k) {
  if(is.list(mean)){
    # mean is different
    p <- data.frame(t(data.frame(mean[[i]])))    
  }else{
    # mean is equal
    p <- data.frame(t(data.frame(mean)))  
  }
  
  if(is.list(mean)&is.list(sigma)){
    # mean and sigma are both different
    e <- data.frame(mixtools::ellipse(mean[[i]], sigma[[i]], newplot = FALSE, npoints = 500))
  } else if(is.list(mean)&(!is.list(sigma))){
    # mean is different while sigma is equal
    e <- data.frame(mixtools::ellipse(mean[[i]], sigma, newplot = FALSE, npoints = 500))  
  }else if((!is.list(mean))&is.list(sigma)){
    # mean is equal while sigma is different
    e <- data.frame(mixtools::ellipse(mean, sigma[[i]], newplot = FALSE, npoints = 500))    
  }else{print('error')}
  out_plot <- out_plot + ggplot2::geom_point(data = p, 
                                             ggplot2::aes(x = X1, y = X2), colour = "black", 
                                             size = 0.7) +
    ggplot2::geom_point(data = e, 
                        ggplot2::aes(x = X1, y = X2), 
                        colour = component_colors[i], size = 0.3) +
    ggplot2::theme_minimal()
}
return(out_plot)
}


gmm.sp <- function(mod, n=1000){
  mean <- mod$mu
  sigma <- mod$sigma
  lambda <- mod$lambda
  k <- length(mod$lambda)
  # realization of lambda by mutinomial distribution
  lbd.sp <- t(rmultinom(n, 1, lambda))
  mvn.sp <- c()
  for (i in 1:k){
    mvn.sp.k <- mvtnorm::rmvnorm(n, mean=mean[[i]],sigma=sigma[[i]])
    lbd.sp.k <- lbd.sp[,i]
    mvn.sp <- rbind(mvn.sp, mvn.sp.k[lbd.sp.k==1,]) 
  }
  # labels for the group from which the sample are generated 
  lab.sp <- apply(lbd.sp, 1, function(x) which(x %in% 1))
return(list('mvn.sp'=mvn.sp, 'lab.sp'=lab.sp))
}

apply(lbd.sp, 1, function(x) which(x %in% 1))

vis.mm(mod.org.2com)
mod.org.2com.sp <- gmm.sp(mod.org.2com,n=5000)
plot(mod.org.2com.sp$mvn.sp, col=mod.org.2com.sp$lab.sp,cex=2,pch=".")
ellipse(mod.org.2com$mu[[1]], mod.org.2com$sigma[[1]])
ellipse(mod.org.2com$mu[[2]], mod.org.2com$sigma[[2]])


vis.mm(mod.org.2com.dm)

vis.mm(mod.org.2com.ds)
plot(gmm.sp(mod.org.2com.ds))

vis.mm(mod.org.3com)
mod.org.3com.sp <- gmm.sp(mod.org.3com)
plot(mod.org.3com.sp$mvn.sp, col=mod.org.3com.sp$lab.sp)

#--------------------fit on log scale------------------------

mu.log.2com <- lapply(mu.2com,log)
mod.log.2com <- mvnormalmixEM(log(riv.dat), k=2,
                              epsilon = 1e-02)

mod.log.2com[2:5]
vis.mm(mod.log.2com)
plot(gmm.sp(mod.log.2com))

#mod.log.3com <- mvnormalmixEM(log(riv.dat), k=3, epsilon = 1e-02)

#mod.log.3com[2:5]
#vis.mm(mod.log.3com)
