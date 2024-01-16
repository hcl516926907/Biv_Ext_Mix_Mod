dir.work <- '/home/pgrad2/2448355h/My_PhD_Project/Biv_Ext_Mix_Mod'
dir.out <- "/home/pgrad2/2448355h/My_PhD_Project/01_Output/Biv_Ext_Mix_Mod/Plots"
source(file.path(dir.work, "Simulation/RevExp_U_Functions.r"))
source(file.path(dir.work, "Simulation/CommonFunctions.r"))


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


#######################################################################################################
# Scatter plot of the bivaratie extreme mixture model

# Parameters are notated same in the paper, except beta being the location parameter in the reverse 
# exponential distribution. Setting beta to be 0 aligns with the representation in the paper.
#######################################################################################################
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



#######################################################################################################
# Contour plot of the bivaratie extreme mixture model
#######################################################################################################

library(ggplot2)
library(nimble)

#Define the density function via Nimble
R_pmnorm_chol <- function(lower, upper, mean, cholesky){
  sigma <- t(cholesky) %*% cholesky
  return(pmvnorm(lower=lower, upper=as.vector(upper), mean=as.vector(mean), sigma=sigma, keepAttr = F))
}


R_dmvnorm_chol <- function(x, mean, cholesky, log){
  sigma <- t(cholesky) %*% cholesky
  dvect <- dmvnorm(x, mean=mean, sigma=sigma, log=log)
  if (log) {
    return(sum(dvect))
  }else{
    return(prod(dvect))
  }
  
}


pmnorm_chol <- nimbleRcall(function(lower = double(1), upper = double(1), 
                                    mean=double(1),cholesky=double(2)){}, 
                           Rfun = 'R_pmnorm_chol',
                           returnType = double(0))



dmvnorm_chol <- nimbleRcall(function(x = double(2), mean = double(1), 
                                     cholesky=double(2),log=logical(0, default = 0)){}, 
                            Rfun = 'R_dmvnorm_chol',
                            returnType = double(0))



nll.powunif.GPD.1<-function(theta,x,u,a.ind,lam.ind,sig.ind,gamma.ind, lamfix=FALSE, balthresh=FALSE, marg.scale.ind,marg.shape.ind)
{ 
  x.mat.ind <- 1
  if (is.null(dim(x))){
    d <- length(x)
    x.mat.ind <- 0
  }else{
    d<-dim(x)[2]
  }
  
  a<-theta[a.ind]
  if(length(a)==1)
  {
    a<-rep(a,d)
  }
  
  if(lamfix){
    lam<-rep(1,d)
  }else{
    lam<-c(theta[lam.ind],1)
  }
  
  if(balthresh){
    lam<-1/(1+a)
  }
  
  sig<-theta[sig.ind]
  gamma<-theta[gamma.ind]
  
  sig<-sig[marg.scale.ind]
  gamma<-gamma[marg.shape.ind]
  
  rej<-NULL
  # upper bound when xi is greater than 0
  if(x.mat.ind){
    for(j in 1:d)
    {
      rej[j]<- gamma[j]<0 && any(x[,j]>-sig[j]/gamma[j])
    }
  }else{
    for(j in 1:d)
    {
      rej[j]<- gamma[j]<0 && any(x[j]>-sig[j]/gamma[j])
    }
  }
  
  if(any(lam<0.01)||any(a<=0.01)||any(sig<=0.001)||any(rej)){return(10e7)}
  
  nll.uc <- 0
  nll.pc <- 0
  if (!x.mat.ind){
    uc <- comp.gt(x, u)
    if (uc){
      L <- fX.powunif(x=x, a=a, lam=lam, sig=sig, gamma=gamma)
      nll.uc <- -log(L)
    }else{
      L2 <- fX.powunif.cens(x=x, u=u, lam=lam, a=a, sig=sig, gamma=gamma)
      nll.pc <- -log(L2)
    }
  }else{
    uc<-apply(x,1,comp.gt,u=u)
    
    x.uc<-x[uc,]
    x.pc<-x[!uc,]
    
    L<-apply(x.uc,1,fX.powunif,a=a,lam=lam,sig=sig,gamma=gamma)
    nll.uc<--sum(log(L))
    
    if(sum(!uc)>0)
    {
      L2<-apply(x.pc,1,fX.powunif.cens,a=a,lam=lam,u=u,sig=sig,gamma=gamma)
      nll.pc<--sum(log(L2))
    }
  }
  if (is.nan(nll.uc)|is.nan(nll.pc)|(nll.uc==-Inf)|(nll.uc==Inf)){
    return(10e7)
  }
  nll<-nll.uc+nll.pc
  
  return(nll)
}





nim_nll_powunif_GPD_SG <- nimbleRcall(function(theta=double(1), x=double(1), u=double(0), a.ind=double(1),
                                               lam.ind=double(0), sig.ind=double(1), gamma.ind=double(1), 
                                               lamfix=logical(0, default = 0), balthresh=logical(0, default = 0), 
                                               marg.scale.ind=double(1), marg.shape.ind=double(1)){}, 
                                      Rfun = 'nll.powunif.GPD.1',
                                      returnType = double(0))

nim_nll_powunif_GPD_MAT <- nimbleRcall(function(theta=double(1), x=double(2), u=double(0), a.ind=double(1),
                                                lam.ind=double(0), sig.ind=double(1), gamma.ind=double(1), 
                                                lamfix=logical(0, default = 0), balthresh=logical(0, default = 0), 
                                                marg.scale.ind=double(1), marg.shape.ind=double(1)){}, 
                                       Rfun = 'nll.powunif.GPD.1',
                                       returnType = double(0))
dbiextmix <- nimbleFunction(
  run = function(x=double(2), theta=double(1), thres=double(1), mu=double(1), 
                 cholesky=double(2), D=integer(0, default=2),
                 a.ind=double(1), lam.ind=double(0), lamfix=logical(0, default = 0), 
                 sig.ind=double(1), gamma.ind=double(1),
                 log = logical(0, default = 0)) {
    returnType(double(0))
    
    cond <- (x[,1]>thres[1]) | (x[,2]>thres[2])
    n.tail <- sum(cond)
    n.bulk <- sum(!cond)
    y.tail <- matrix(c(x[cond,1] - thres[1], x[cond,2] - thres[2]), ncol=D)
    y.bulk <- x[!cond,]
    
    pi <- pmnorm_chol(lower=rep(-Inf,D), upper=thres, mean=mu, cholesky = cholesky)
    
    sig <- theta[sig.ind]
    gamma <- theta[gamma.ind]
    eta <- -sig/gamma
    eta[which(gamma<=0)] <- -Inf
    
    dtail <- 0
    dbulk <- 0
    
    if (n.tail>0){
      y.min <- eta
      for (i in 1:D){
        y.min[i] <- min(y.tail[,i])
      }
      if (all(y.min>eta)){
        llt <- -nim_nll_powunif_GPD_MAT(x=y.tail, theta=theta, u=min(y.tail)-0.01, a.ind=a.ind,
                                        lam.ind=lam.ind, sig.ind=sig.ind, gamma.ind=gamma.ind, 
                                        lamfix=lamfix, balthresh=FALSE, 
                                        marg.scale.ind=1:2, marg.shape.ind=1:2)
        if (log){
          dtail <- llt
        }else{
          dtail <- exp(llt)
        }
      }else{
        if (log) dtail <- -10^10
      }
    }
    
    if (n.bulk>0){
      dbulk <- dmvnorm_chol(y.bulk, mean=mu, cholesky = cholesky, log = log )
    }
    
    if (log) {
      totalProb <- n.tail*log(1-pi) + dtail + dbulk
    }else{
      totalProb <- (1-pi)^n.tail *dtail*dbulk
    }
    
    return(totalProb)
  })


mu <- c(3.5, 4.0)
sd1 <- 1
sd2 <- 1.5
rho <- 0.7
sigma <- matrix(c(sd1^2, rho*sd1*sd2, rho*sd1*sd2, sd2^2),ncol=2)

dbiextmix(matrix(c(1,2),nrow=1), theta=c(2,2,1,0.5,0.5,0.1,0.1),thres=c(5.5,6.7),mu=c(3.5,4),
          cholesky = chol(sigma),D=2,a.ind=c(1,2),lam.ind=3,lamfix=T,
          sig.ind=4:5,gamma.ind=6:7, log=T)


#log-likelihood with standard GPD tail, i.e. sigma=c(1,1), gamma=c(0,0)
log_likelihood <- function(x){
  Y <- rbind(x,x)
  a <- c(1, 2)
  sig <- c(1, 1)
  gamma <- c(0, 0)
  theta <- c(a,1,sig,gamma)
  thres <- c(5.7, 7.5)
  mu <- c(3.5, 3.7)
  sd1 <- 1
  sd2 <- 1.5
  rho <- 0.7
  sigma <- matrix(c(sd1^2, rho*sd1*sd2, rho*sd1*sd2, sd2^2),ncol=2)
  ll <- dbiextmix(Y, theta=theta,thres=thres,mu=mu,
            cholesky = chol(sigma),D=2,a.ind=1:2,lam.ind=3,lamfix=T,
            sig.ind=4:5,gamma.ind=6:7, log=T)
  return(ll/2)
  
}


# Step 2: Create a grid of x and y values
x_seq <- seq(-2, 12, length.out = 100)
y_seq <- seq(-2, 12, length.out = 100)
grid <- expand.grid(x = x_seq, y = y_seq)

# Step 3: Evaluate the log-likelihood over the grid
grid$ll <- apply(grid,1,log_likelihood)

# Step 4: Exponentiate to get the density
grid$density <- exp(grid$ll)


bin1 <- cut(grid$density, breaks = 20)
interval1 <- unique(bin1)
all_edges1 <- unlist(lapply(interval1, function(interval) {
  string <- gsub("\\[|\\]|\\(|\\)", "", interval)
  as.numeric(unlist(strsplit(string, ",")))
}))
unique_edges1 <- unique(all_edges1)
print(unique_edges1)

bin2 <- cut(grid$density[which(grid$density<unique_edges1[2])], breaks=20)
interval2 <- unique(bin2)
all_edges2 <- unlist(lapply(interval2, function(interval) {
  string <- gsub("\\[|\\]|\\(|\\)", "", interval)
  as.numeric(unlist(strsplit(string, ",")))
}))
unique_edges2 <- unique(all_edges2)

contour_break <- c(unique_edges1[2:length(unique_edges1)],unique_edges2[2:length(unique_edges2)])
# Step 5: Create the contour plot

contour.palette <- brewer.pal(8,"Accent")
p <- ggplot(grid, aes(x = x, y = y, z = density)) +
  geom_contour(aes(color = ..level..),breaks=contour_break) +
  scale_color_continuous(low = contour.palette[5], high = contour.palette[7]) +
  geom_segment(aes(x = 5.7, y = -2, xend = 5.7, yend = 7.5), color = "black", linetype = "dashed") +
  geom_segment(aes(x = -0.5, y = 7.5, xend = 5.7, yend = 7.5), color = "black", linetype = "dashed") + 
  annotate("rect", xmin = 5.7, xmax = 10, ymin = -2, ymax = 12,
           alpha = .05,fill = contour.palette[5])+
  annotate("rect", xmin = -0.5, xmax = 5.7, ymin = 7.5, ymax = 12,
           alpha = .05,fill = contour.palette[5])+
  annotate("rect", xmin = -0.5, xmax = 5.7, ymin = -2, ymax = 7.5,
           alpha = .05,fill = contour.palette[7])+
  labs(x = "X1",
       y = "X2",
       color = "Density") +
  xlim(-0.5,10)+
  ylim(-2,12)+
  theme(axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        legend.position = "none")

res <- 300
png(filename = file.path("/home/pgrad2/2448355h/My_PhD_Project/01_Output/Biv_Ext_Mix_Mod/Plots", "BEMM_Contour_plot_1.png"), width = 6*res, height = 5*res, res=res)
print(p)
dev.off()

#---------------------------------two positive gammas--------------------------------

#log-likelihood with sigma=c(1,1.2), gamma=c(0.2,0.3)
log_likelihood <- function(x){
  Y <- rbind(x,x)
  a <- c(1, 2)
  sig <- c(1, 1.2)
  gamma <- c(.2, 0.3)
  theta <- c(a,1,sig,gamma)
  thres <- c(5.7, 7.5)
  mu <- c(3.5, 3.7)
  sd1 <- 1
  sd2 <- 1.5
  rho <- 0.7
  sigma <- matrix(c(sd1^2, rho*sd1*sd2, rho*sd1*sd2, sd2^2),ncol=2)
  ll <- dbiextmix(Y, theta=theta,thres=thres,mu=mu,
                  cholesky = chol(sigma),D=2,a.ind=1:2,lam.ind=3,lamfix=T,
                  sig.ind=4:5,gamma.ind=6:7, log=T)
  return(ll/2)
  
}


# Step 2: Create a grid of x and y values
x_seq <- seq(-2, 12, length.out = 100)
y_seq <- seq(-2, 12, length.out = 100)
grid <- expand.grid(x = x_seq, y = y_seq)

# Step 3: Evaluate the log-likelihood over the grid
grid$ll <- apply(grid,1,log_likelihood)

# Step 4: Exponentiate to get the density
grid$density <- exp(grid$ll)


bin1 <- cut(grid$density, breaks = 20)
interval1 <- unique(bin1)
all_edges1 <- unlist(lapply(interval1, function(interval) {
  string <- gsub("\\[|\\]|\\(|\\)", "", interval)
  as.numeric(unlist(strsplit(string, ",")))
}))
unique_edges1 <- unique(all_edges1)
print(unique_edges1)

bin2 <- cut(grid$density[which(grid$density<unique_edges1[2])], breaks=20)
interval2 <- unique(bin2)
all_edges2 <- unlist(lapply(interval2, function(interval) {
  string <- gsub("\\[|\\]|\\(|\\)", "", interval)
  as.numeric(unlist(strsplit(string, ",")))
}))
unique_edges2 <- unique(all_edges2)

contour_break <- c(unique_edges1[2:length(unique_edges1)],unique_edges2[2:length(unique_edges2)])
# Step 5: Create the contour plot

contour.palette <- brewer.pal(8,"Accent")
p <- ggplot(grid, aes(x = x, y = y, z = density)) +
  geom_contour(aes(color = ..level..),breaks=contour_break) +
  scale_color_continuous(low = contour.palette[5], high = contour.palette[7]) +
  geom_segment(aes(x = 5.7, y = -2, xend = 5.7, yend = 7.5), color = "black", linetype = "dashed") +
  geom_segment(aes(x = -0.5, y = 7.5, xend = 5.7, yend = 7.5), color = "black", linetype = "dashed") + 
  
  geom_segment(aes(x =  0.7, y = 7.5, xend = 0.7, yend = 12), color = "black", linetype = "twodash") + 
  geom_segment(aes(x =  5.7, y = 3.5, xend = 10, yend = 3.5), color = "black", linetype = "twodash") + 
  annotate("rect", xmin = 5.7, xmax = 10, ymin = 3.5, ymax = 12,
           alpha = 0.1,fill = contour.palette[5])+
  annotate("rect", xmin = 0.7, xmax = 5.7, ymin = 7.5, ymax = 12,
           alpha = 0.1,fill = contour.palette[5])+
  annotate("rect", xmin = -0.5, xmax = 5.7, ymin = -2, ymax = 7.5,
           alpha = 0.1,fill = contour.palette[7])+
  labs(x = "X1",
       y = "X2",
       color = "Density") +
  xlim(-0.5,10)+
  ylim(-2,12)+
  theme(axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        legend.position=c(1, 0),
        legend.justification = c(1, 0),
        legend.title = element_text(size=15),
        legend.text = element_text(size=15),
        legend.key.size = unit(0.05, 'npc'),
        legend.key.width = unit(0.05, 'npc'))
print(p)


png(filename = file.path("/home/pgrad2/2448355h/My_PhD_Project/01_Output/Biv_Ext_Mix_Mod/Plots", "BEMM_Contour_plot_2.png"), width = 6*res, height = 5*res, res=res)
print(p)
dev.off()