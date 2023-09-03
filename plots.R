dir.work <- '/home/pgrad2/2448355h/My_PhD_Project/Biv_Ext_Mix_Mod'
dir.out <- "/home/pgrad2/2448355h/My_PhD_Project/01_Output/Biv_Ext_Mix_Mod/Plots"
source(file.path(dir.work, "KRSW/RevExp_U_Functions.r"))
source(file.path(dir.work, "KRSW/CommonFunctions.r"))

# install_load_packages <- function(packages) {
#   nodename <- Sys.info()['nodename']
#   lib_path <- file.path("/home/pgrad2/2448355h/R/library", nodename)
#   .libPaths(lib_path)
#   # Create a separate library folder for each node, if it doesn't exist
#   if (!dir.exists(lib_path)) {
#     dir.create(lib_path, recursive = TRUE)
#   }
#   
#   for (package in packages) {
#     # Check if the package is installed in the specific folder
#     if(!require(package, character.only = TRUE, lib.loc = lib_path)) {
#       install.packages(package, lib = lib_path, dependencies = TRUE)
#       library(package, character.only = TRUE, lib.loc = lib_path)
#     }
#     
#     # Load the package
#     library(package, character.only = TRUE, lib.loc = lib_path)
#   }
# }

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


packages <- c( "mvtnorm", "tmvtnorm",'RColorBrewer')  


load_install_packages(packages)

seed <- 1
d <- 2
a <- c(1.4, 1.2)
beta <- c(0, 0)
sig <- c(0.5, 0.4)
gamma <- c(0.1, 0.2)
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

eta <- -sig/gamma + u.x
mypalette<-brewer.pal(7,"Greys")
color <- c(rep(mypalette[7], nrow(Y.bulk)), rep(mypalette[5], nrow(Y.tail$X)))

res <- 500

png(filename = file.path(dir.out, "BEMM_Sample_postive_gamma.png"), width = 6*res, height = 5*res, res=res)
par(mar = c(4.2, 4.2, 2, 1.5))
plot(Y, col=color,xlab='X1', ylab='X2', pch=19, cex=0.7)
segments(x0=u.x[1],y0=min(Y)-10, x1=u.x[1],y1=u.x[2])
segments(x0=min(Y)-10, y0= u.x[2], x1=u.x[1], y1=u.x[2])

segments(x0=eta[1],y0=u.x[2], x1=eta[1],y1=max(Y)+10, lty=2)
segments(x0=u.x[1], y0= eta[2], x1=max(Y)+10, y1=eta[2],lty=2)

dev.off()



seed <- 1
d <- 2
a <- c(1.4, 1.2)
beta <- c(0, 0)
sig <- c(0.5, 0.4)
gamma <- c(-0.1, -0.2)
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

eta <- -sig/gamma + u.x
mypalette<-brewer.pal(7,"Greys")
color <- c(rep(mypalette[7], nrow(Y.bulk)), rep(mypalette[5], nrow(Y.tail$X)))

png(filename = file.path(dir.out, "BEMM_Sample_negative_gamma.png"), width = 6*res, height = 5*res, res=res)
par(mar = c(4.2, 4.2, 2, 1.5))
plot(Y, col=color,xlab='X1', ylab='X2', pch=19, cex=0.7, xlim=c(min(Y),max(Y)+3),ylim=c(min(Y),max(Y)+2) )
segments(x0=u.x[1],y0=min(Y)-10, x1=u.x[1],y1=u.x[2])
segments(x0=min(Y)-10, y0= u.x[2], x1=u.x[1], y1=u.x[2])

segments(x0=eta[1],y0=min(Y)-10, x1=eta[1],y1=eta[2], lty=2)
segments(x0=min(Y)-10, y0= eta[2], x1=eta[1], y1=eta[2],lty=2)
dev.off()
