dir.work <- '/home/pgrad2/2448355h/My_PhD_Project/Biv_Ext_Mix_Mod'
dir.out <- "/home/pgrad2/2448355h/My_PhD_Project/01_Output/Biv_Ext_Mix_Mod/Simulation"
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

# List the packages you want to load
packages <- c("nimble", "foreach","doSNOW","parallel",'copula')  


load_install_packages(packages)

print(detectCores())

######################################################################################
# Parallelled code for running the simulation multiple times
# The reason of dividing the code into two parts instead of one is because large number 
# iterations in the foreach cound cause some unexpected error.
######################################################################################
i <- 1
seed <- i
d <- 2
a <- c(2, 2.5)
beta <- c(0, 0)
sig <- c(1, 0.8)
gamma <- c(-0.1, -0.2)

n <- 2000

par <- c(1.5,  7, 2, 2, 4)

u.x <- c(6, 6.5)

lbound.tail=c(0,0)
ubound.tail=c(Inf,Inf)

gum.cop <- gumbelCopula(param = 1.5, dim = 2)

gum.dist <- mvdc(gum.cop, margins = c('gamma','weibull'), 
                 
                 paramMargins=list(list(shape=par[2], rate=par[3]),
                                   
                                   list(shape=par[4],scale=par[5])))

p <- pMvdc(u.x, gum.dist)


set.seed(1111)
Y.tail<-sim.RevExpU.MGPD(n=n-floor(n*p),d=d, a=a, beta=beta, sig=sig, gamma=gamma, MGPD = T,std=T)$X

# GP scale tail data combined with the bulk data

Y.bulk <- rMvdc(floor(n*p), gum.dist)
Y.bulk <- Y.bulk[Y.bulk[,1]<u.x[1] & Y.bulk[,2]<u.x[2],]

Y <- rbind(Y.bulk, sweep(Y.tail,2,u.x,"+"))
apply(Y.bulk,2,min) + u.x
plot(Y)
bound.cond.1 <- Y.tail[,1] > (lbound.tail-u.x)[1] & Y.tail[,1] < (ubound.tail-u.x)[1]
bound.cond.2 <- Y.tail[,2] > (lbound.tail-u.x)[2] & Y.tail[,2] < (ubound.tail-u.x)[2]
Y.tail <- Y.tail[bound.cond.1&bound.cond.2,]
Y <- rbind(Y.bulk, sweep(Y.tail,2,u.x,"+"))
plot(Y)







jointY1.revexp <- function(x,y1,a, sig, gamma){
  if (length(x)==1){
    y <- c(y1,x)
    res <- fX.powunif(y, lam=c(1,1), a=a, sig=sig, gamma=gamma)
    if (is.na(res)) res <- 0
  }else{
    y <- cbind(rep(y1,length(x)),x)
    res <- apply(y, MARGIN = 1, fX.powunif, lam=c(1,1),a=a, sig=sig, gamma=gamma)
    res[is.na(res)] <- 0
  }
  return(res)
}

marginY1.revexp <- function(y1, y2.low=-Inf, y2.upp=Inf, a, sig, gamma){
  upp.para <- if (gamma[2] <0) -sig[2]/gamma[2] else Inf  
  upper <- min(y2.upp, upp.para)
  
  low.para <- if( gamma[2]>0) -sig[2]/gamma[2] else -Inf  
  lower <- max(y2.low, low.para )
  
  if(y1 < 0){
    res <- integrate(jointY1.revexp, lower = 0, upper = upper, y1=y1, a=a,sig=sig, gamma=gamma)$value
  }else{
    res <- integrate(jointY1.revexp, lower = lower, upper = upper, y1=y1, a=a,sig=sig, gamma=gamma)$value
  }
  
  return(res)
}


marginY1.revexp(-10, y2.upp=Inf, a=a, sig=sig, gamma=gamma)

jointY2.revexp <- function(x,y2,a, sig, gamma){
  if (length(x)==1){
    y <- c(x,y2)
    res <- fX.powunif(y, lam=c(1,1), a=a, sig=sig, gamma=gamma)
    if (is.na(res)) res <- 0
  }else{
    y <- cbind(x, rep(y2,length(x)))
    res <- apply(y, MARGIN = 1, fX.powunif, lam=c(1,1),a=a, sig=sig, gamma=gamma)
    res[is.na(res)] <- 0
  }
  return(res)
}

marginY2.revexp <- function(y2, y1.low=-Inf, y1.upp=Inf, a, sig, gamma){
  upp.para <- if (gamma[1] <0) -sig[1]/gamma[1] else Inf  
  upper <- min(y1.upp, upp.para)
  
  low.para <- if( gamma[1]>0) -sig[1]/gamma[1] else -Inf  
  lower <- max(y1.low, low.para )
  
  if(y2 < 0){
    res <- integrate(jointY2.revexp, lower = 0, upper = upper, y2=y2, a=a,sig=sig, gamma=gamma)$value
  }else{
    res <- integrate(jointY2.revexp, lower = lower, upper = upper, y2=y2, a=a,sig=sig, gamma=gamma)$value
  }
  
  return(res)
}

marginY2.revexp(-0.1,y1.upp=Inf, a=a, sig=sig, gamma=gamma)



# standarized Y1,Y1<0 with sig=1 and gamma=0
fY1_lt0.margin <- function(y,a,ub){
  EM <- EM.pu(a=a,lam=c(1,1))
  b <- 1/a
  den <- EM*(1+sum(b))
  num <- prod(b)*exp(b[1]*y)*(1-exp(-(1+b[1])*ub))/(1+b[1])
  return(num/den)
}

fY1_lt0.margin(0,a=a,ub=Inf)

marginY1.revexp(-500, y2.upp=Inf, a=a, sig=c(1,1), gamma=c(0,0))

fY1_gt0.margin <- function(y,a,lb, ub){
  EM <- EM.pu(a=a,lam=c(1,1))
  b <- 1/a
  num1 <- b[1]*(exp(-y)-exp(b[2]*lb-(1+b[2])*y))
  den1 <- EM*(1+sum(b))
  
  num2 <- prod(b)*exp(-y)-prod(b)*exp(b[1]*y-(1+b[1])*ub)
  den2 <- EM*(1+sum(b))*(1+b[1])
  return(num1/den1+num2/den2)
}

fY1_gt0.margin(5,a=a, lb=-3, ub=Inf)
marginY1.revexp(5, y2.low= -3 ,y2.upp=Inf, a=a, sig=c(1,1), gamma=c(0,0))

tail.margin1 <- function(y, a,sig, gamma, tail.low, tail.upp){
  upp.para1 <- if (gamma[1] <0) -sig[1]/gamma[1] else Inf  
  upp.para2 <- if (gamma[2] <0) -sig[2]/gamma[2] else Inf 
  upper <- pmin(tail.upp, c(upp.para1, upp.para2))
  low.para1 <- if( gamma[1]>0) -sig[1]/gamma[1] else -Inf  
  low.para2 <- if( gamma[2]>0) -sig[2]/gamma[2] else -Inf  
  lower <- pmax(tail.low, c(low.para1,low.para2) )
  
  dtail <- 0
  if ( y> lower[1] & y< upper[1]){
    
  for (i in 1:length(sig)){
    if (abs(gamma[i])<10^-6){
      lower[i] <- lower[i]/sig[i]
      upper[i] <- upper[i]/sig[i]
    }else{
      if(is.na(log(gamma[i]/sig[i]*upper[i]+1+ 10^-6))) print(list(gamma=gamma,sig=sig,upper=upper))
      lower[i] <- log(gamma[i]/sig[i]*lower[i]+1 + 10^-6)/gamma[i]
      upper[i] <- log(gamma[i]/sig[i]*upper[i]+1 + 10^-6)/gamma[i]
      
    }
  }

    if (gamma[1]==0){
      z <- y/sig[1]
      jacob <- 1/sig[1]
    }else{
      z <- log(gamma[1]*y/sig[1]+1)/gamma[1]
      jacob <- 1/(y*gamma[1]+sig[1])
    }
    
    if (z<0){
      dtail <- fY1_lt0.margin(z,a=a,ub=upper[2])*abs(jacob)
    }else{
      dtail <- fY1_gt0.margin(z,a=a,lb=lower[2], ub=upper[2])*abs(jacob)
    }
  }
  return(dtail)
}

marginY1.revexp(1, y2.low= -3 ,y2.upp=Inf, a=a, sig=c(3,1), gamma=c(0,0))

tail.margin1(1, a=a, sig=c(3,1), gamma=c(0,0), tail.low=c(-3,-3), tail.upp=c(Inf,Inf))

marginY1.revexp(10, y2.low= -3 ,y2.upp=Inf, a=a, sig=c(3,1), gamma=c(0.5,0))
tail.margin1(10, a=a, sig=c(3,1), gamma=c(0.5,0), tail.low=c(-3,-3), tail.upp=c(Inf,Inf))


x.seq <- seq(-0.5,10,0.1)
y <- sapply(x.seq, tail.margin1, a=a, sig=c(2,1), gamma=c(2,0), tail.low=c(-Inf,-3), tail.upp=c(Inf,Inf))
plot(x.seq,y,type='l')

y1 <- sapply(x.seq, marginY1.revexp, a=a, sig=c(2,1), gamma=c(2,0), y2.low= -3 ,y2.upp=Inf)
plot(x.seq,y1,type='l')



# standarized Y1,Y1<0 with sig=1 and gamma=0
fY2_lt0.margin <- function(y,a,ub){
  EM <- EM.pu(a=a,lam=c(1,1))
  b <- 1/a
  den <- EM*(1+sum(b))
  num <- prod(b)*exp(b[2]*y)*(1-exp(-(1+b[2])*ub))/(1+b[2])
  return(num/den)
}

fY2_lt0.margin(-5,a=a,ub=5)
marginY2.revexp(-5, y1.upp=5, a=a, sig=c(1,1), gamma=c(0,0))

fY2_gt0.margin <- function(y, a, lb, ub){
  EM <- EM.pu(a=a,lam=c(1,1))
  b <- 1/a
  num1 <- b[2]*(exp(-y)-exp(b[1]*lb-(1+b[1])*y))
  den1 <- EM*(1+sum(b))
  
  num2 <- prod(b)*exp(-y)-prod(b)*exp(b[2]*y-(1+b[2])*ub)
  den2 <- EM*(1+sum(b))*(1+b[2])
  
  return(num1/den1+num2/den2)
}
fY2_gt0.margin(0.1,a=a, lb=-1,ub=5)
marginY2.revexp(0.1, y1.low=-1, y1.upp=5, a=a, sig=c(1,1), gamma=c(0,0))


trunc.prob <- function(lower=rep(-Inf,2),upper=rep(Inf,2),a,sig,gamma){
  upp.para1 <- if (gamma[1] <0) -sig[1]/gamma[1] else Inf  
  upp.para2 <- if (gamma[2] <0) -sig[2]/gamma[2] else Inf 
  upper <- pmin(upper, c(upp.para1, upp.para2))
  
  low.para1 <- if( gamma[1]>0) -sig[1]/gamma[1] else -Inf  
  low.para2 <- if( gamma[2]>0) -sig[2]/gamma[2] else -Inf  
  lower <- pmax(lower, c(low.para1,low.para2) )
  
  double.int.Y1 <- function(x, y2.low=lower[2], y2.upp=upper[2], a=a, sig=sig, gamma=gamma){
    sapply(x, marginY1.revexp, y2.low = y2.low, y2.upp=y2.upp, a=a, sig=sig, gamma=gamma)
  }
  
  p <- integrate(double.int.Y1, lower =lower[1], upper = upper[1] , a=a,sig=sig, gamma=gamma)$value
  return(p)
}



trunc.norm.const <- function(a,lower,upper){
  lb1 <- lower[1]
  lb2 <- lower[2]
  ub1 <- upper[1]
  ub2 <- upper[2]
  EM <- EM.pu(a=a,lam=c(1,1))
  b <- 1/a
  
  den1 <- EM*(1+sum(b))*(1+b[2])
  num1 <- b[1]*(1-exp(-(1+b[2])*ub1))*(1-exp(b[2]*lb2))
  
  den2 <- EM*(1+sum(b))*(1+b[1])
  num2 <- (1+b[1])*b[2]*(1-exp(-ub2))-b[2]*exp(b[1]*lb1)*(1-exp(-(1+b[1])*ub2))
  
  den3 <- EM*(1+sum(b))*(1+b[2])
  num3 <- prod(b)*(1-exp(-ub2))-b[1]*(exp(-ub1)-exp(-(1+b[2])*ub1))
  
  return(num1/den1+num2/den2 + num3/den3)
  
}


trunc.norm.const.GPD <- function(a, sig, gamma, lower, upper ){
  upp.para1 <- if (gamma[1] <0) -sig[1]/gamma[1] else Inf  
  upp.para2 <- if (gamma[2] <0) -sig[2]/gamma[2] else Inf 
  upper <- pmin(upper, c(upp.para1, upp.para2))
  
  low.para1 <- if( gamma[1]>0) -sig[1]/gamma[1] else -Inf  
  low.para2 <- if( gamma[2]>0) -sig[2]/gamma[2] else -Inf  
  lower <- pmax(lower, c(low.para1,low.para2) )
  # can't do it vectorized because the formula for dim 1 and dim 2 may differ depending on the value of gamma
  for (i in 1:length(sig)){
    if (abs(gamma[i])<10^-6){
      lower[i] <- lower[i]/sig[i]
      upper[i] <- upper[i]/sig[i]
    }else{
      if(is.na(log(gamma[i]/sig[i]*upper[i]+1+ 10^-6))) print(list(gamma=gamma,sig=sig,upper=upper))
      lower[i] <- log(gamma[i]/sig[i]*lower[i]+1 + 10^-6)/gamma[i]
      upper[i] <- log(gamma[i]/sig[i]*upper[i]+1 + 10^-6)/gamma[i]
      
    }
  }
  return(trunc.norm.const(a,lower,upper))
}

trunc.norm.const.GPD(a=c(4,1), sig=c(0.1,2), gamma=c(0.5,0.5),lower=c(-0.2,-0.2), upper=c(4,5))
trunc.prob(lower=c(-0.2,-0.2),upper=c(4,5),a=c(4,1),sig=c(0.1,2),gamma=c(0.5,0.5))




x.seq <- seq(-2.499,10,0.1) 
y.tail <- sapply(x.seq, marginY1.revexp,y2.low=-5, y2.upp=Inf, a=a, sig=c(5,5), gamma=c(2,2) )
plot(x.seq,y.tail,type='l')





dMvdc(c(1,1), gum.dist)

bulk.margin.1 <- function(x, x2.low=0, x2.upp=u.x[2], dist=gum.dist){

  marg.den <- function(x2){
    point <- cbind(rep(x,length(x2)), x2)
    return(dMvdc(point, gum.dist))
  }
  p <- integrate(marg.den, lower =x2.low, upper = x2.upp)$value
  return(p)
}

bulk.margin.1(1)

bulk.margin.2 <- function(x, x1.low=0, x1.upp=u.x[1], dist=gum.dist){
  
  marg.den <- function(x1){
    point <- cbind(x1,rep(x,length(x1)))
    return(dMvdc(point, gum.dist))
  }
  p <- integrate(marg.den, lower =x1.low, upper = x1.upp)$value
  return(p)
}

bulk.margin.2(1)


BEMM.margin1 <- function(x, dist=gum.dist, thres=u.x, a.=a, sig.=sig, gamma.=gamma,
                         bulk.low=rep(0,2), tail.low=rep(0,2),tail.upp=rep(Inf,2)){
  # transfrom to mGPD scale 
  tail.low <- tail.low - u.x
  upp.para1 <- if (gamma.[1] <0) -sig.[1]/gamma.[1] else Inf  
  upp.para2 <- if (gamma.[2] <0) -sig.[2]/gamma.[2] else Inf 
  tail.upp <- pmin(tail.upp, c(upp.para1, upp.para2))
  
  low.para1 <- if( gamma.[1]>0) -sig.[1]/gamma.[1] else -Inf  
  low.para2 <- if( gamma.[2]>0) -sig.[2]/gamma.[2] else -Inf  
  tail.low <- pmax(tail.low, c(low.para1,low.para2) )
  
  
  pi <- pMvdc(thres, dist)
  M <- trunc.norm.const.GPD(a., sig., gamma., tail.low,tail.upp)
  x.tail <- x-thres[1]
  dtail <- 0
  if (x.tail > tail.low[1] & x.tail < tail.upp[1]){
    # dtail <- marginY1.revexp(x-thres[1], y2.low=tail.low[2], y2.upp=tail.upp[2], a=a., sig=sig., gamma=gamma.)
    dtail <- tail.margin1(x.tail, a=a., sig=sig., gamma=gamma., tail.low=tail.low, tail.upp=tail.upp)
  }
  tdtail <- dtail/M
  
  dbulk <- bulk.margin.1(x, x2.low=bulk.low[2], x2.upp=thres[2], dist=dist)
  
  if (x>=bulk.low[1] & x<thres[1]){
    dBEMM <-  dbulk + (1-pi)*tdtail
  }else if( x > thres[1] & x < (tail.upp[1]+thres[1]) ){
    dBEMM <- (1-pi)*tdtail
  }else dBEMM <- 0
  
  return(dBEMM)
}


x.seq <- seq(0.001,10,0.01) 
y <- sapply(x.seq, BEMM.margin1, dist=gum.dist, thres=u.x, a.=a, sig.=c(2,1), gamma.=c(0,0),
                 bulk.low=rep(0,2), tail.low=rep(0,2),tail.upp=rep(Inf,2) )

y.raw <- dgamma(x.seq, shape=par[2],rate=par[3])
plot(x.seq,y,type='l')
lines(x.seq, y.raw, col='red')

y2 <- sapply(x.seq, BEMM.margin1, dist=gum.dist, thres=u.x, a.=a, sig.=sig, gamma.=c(0.5,0.3),
            bulk.low=rep(0,2), tail.low=rep(0,2),tail.upp=rep(Inf,2) )
plot(x.seq,y2,type='l')
lines(x.seq, y.raw, col='red')

test1 <- function(x){
  sapply(x, BEMM.margin1, dist=gum.dist, thres=u.x, a.=a, sig.=c(2,1), gamma.=c(0.5,-0.8),
         bulk.low=rep(0,2), tail.low=rep(0,2),tail.upp=rep(Inf,2) )
}

integrate(test1, lower =0, upper = Inf)$value

#######################################plot for margin1#################
library(ggplot2)
library(RColorBrewer)
res <- 500
df.plot <- data.frame(x = x.seq,
                      y.BEMM = y,
                      y.BEMM2 = y2,
                      y.bulk = y.raw)
discon <- as.numeric(df.plot[df.plot$x==u.x[1]+0.001, c('x','y.BEMM')])
df.plot[df.plot$x==u.x[1]+0.001,'y.BEMM'] <- NA

discon2 <- as.numeric(df.plot[df.plot$x==4+0.001, c('x','y.BEMM2')])
df.plot[df.plot$x==u.x[1]+0.001,'y.BEMM2'] <- NA
df.plot[df.plot$x==4+0.001,'y.BEMM2'] <- NA

pal <- brewer.pal(8,"Accent")

colors <- c('BEMM'=pal[7], 'Gamma'=pal[5])
p <- ggplot(df.plot) +
  geom_line( aes(x = x, y = y.BEMM,col='BEMM', linetype='BEMM') ) +
  geom_point(aes(x=discon[1], y= discon[2] , col='BEMM'),size=1)+
  geom_line( aes(x = x, y = y.bulk,col='Gamma',linetype="Gamma"))+
  scale_color_manual(values = colors)+
  scale_linetype_manual(values=c('BEMM'=1,'Gamma'=2))+
  labs(color  = "Density", linetype = "Density")+
  theme(axis.title.y=element_blank())


png(filename = file.path("/home/pgrad2/2448355h/My_PhD_Project/01_Output/Biv_Ext_Mix_Mod/Plots", "BEMM_margin1.png"), width = 6*res, height = 5*res, res=res)
print(p)
dev.off()

p <- ggplot(df.plot) +
  geom_line( aes(x = x, y = y.BEMM2,col='BEMM', linetype='BEMM') ) +
  geom_point(aes(x=discon[1], y= discon[2] , col='BEMM'),size=1)+
  geom_point(aes(x=discon2[1], y= discon2[2] , col='BEMM'),size=1)+
  geom_line( aes(x = x, y = y.bulk,col='Gamma',linetype="Gamma"))+
  scale_color_manual(values = colors)+
  scale_linetype_manual(values=c('BEMM'=1,'Gamma'=2))+
  labs(color  = "Density", linetype = "Density")+
  theme(axis.title.y=element_blank())

png(filename = file.path("/home/pgrad2/2448355h/My_PhD_Project/01_Output/Biv_Ext_Mix_Mod/Plots", "BEMM_margin2.png"), width = 6*res, height = 5*res, res=res)
print(p)
dev.off()

ggplot(df_discontinuous, aes(x = x, y = y.BEMM)) +
  geom_line() +
  geom_point(data = df_discontinuous[!is.na(df_discontinuous$y.BEMM), ], aes(x = x, y = y.BEMM))


###################################contour density plot


R_dbulk <- function(x,params, bulk.dist.name, lbound.bulk, ubound.bulk ){
  cop.name.map <- c("normal","clayton","gumbel","frank","joe","plackett")
  cop <- switch(cop.name.map[bulk.dist.name[1]],
                "normal" = normalCopula(param = params[1], dim = 2),
                "clayton" = claytonCopula(param = params[1], dim = 2),
                "gumbel" = gumbelCopula(param = params[1], dim = 2),
                "frank" = frankCopula(param = params[1], dim = 2),
                "joe" = joeCopula(param = params[1], dim = 2),
                "plackett" = plackettCopula(param = params[1]),
                stop("Unsupported copula name"))
  
  margin.name.map <- c("norm","exp","gamma","lnorm","weibull")
  margins.name <- margin.name.map[bulk.dist.name[2:3]]
  params.margin <- params[2:5]
  transformed_params <- list()
  for (i in 1:length(margins.name)) {
    m <- margins.name[i]
    p <- params.margin[(2*i-1) : (2*i)]
    
    if (m == "norm") {
      transformed_params[[i]] <- list(mean = p[1], sd = p[2]) # ensure sd > 0
    } else if (m == "exp") {
      transformed_params[[i]] <- list(rate = p[1]) # ensure rate > 0
    } else if (m == "gamma") {
      transformed_params[[i]] <- list(shape = p[1], rate = p[2]) # ensure shape, rate > 0
    } else if (m == "lnorm") {
      transformed_params[[i]] <- list(meanlog = p[1], sdlog = p[2]) # ensure sdlog > 0
    } else if (m == "weibull") {
      transformed_params[[i]] <- list(shape = p[1], scale = p[2]) # ensure shape, scale > 0
    } else {
      stop("Unsupported marginal distribution")
    }
  }
  
  
  joint.dist <- mvdc(cop, margins = margins.name, 
                     
                     paramMargins=list(transformed_params[[1]],
                                       
                                       transformed_params[[2]])
  )
  
  p1 <- pMvdc(ubound.bulk, joint.dist)
  p2 <- pMvdc(lbound.bulk, joint.dist)
  p <- p1 - p2 
  ll.bulk <- sum(dMvdc(x, joint.dist,log=T))
  return(c(p, ll.bulk))
  
}

# gumbel, gamma, weibull
R_dbulk(Y.bulk,c(1.5,  7, 2, 2, 4),c(3,3,5),c(0,0),u.x)



dbulk <- nimbleRcall(function(x = double(2), params = double(1),bulk.dist.name = double(1),
                              lbound.bulk = double(1), ubound.bulk = double(1)){}, 
                     Rfun = 'R_dbulk',
                     returnType = double(1))

dbulk(Y.bulk,c(1.5,  7, 2, 2, 4),c(3,3,5),c(0,0),u.x)


nll.powunif.GPD.1<-function(theta,x)
{ 
  
  u <- min(x)-0.01
  x.mat.ind <- 1
  if (is.null(dim(x))){
    d <- length(x)
    x.mat.ind <- 0
  }else{
    d<-dim(x)[2]
  }
  
  a<-theta[1:2]
  
  lam<-rep(1,d)
  
  
  sig<-theta[3:4]
  gamma<-theta[5:6]
  
  
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





nim_nll_powunif_GPD_SG <- nimbleRcall(function(theta=double(1), x=double(1)){}, 
                                      Rfun = 'nll.powunif.GPD.1',
                                      returnType = double(0))

nim_nll_powunif_GPD_MAT <- nimbleRcall(function(theta=double(1), x=double(2)){}, 
                                       Rfun = 'nll.powunif.GPD.1',
                                       returnType = double(0))

trunc.norm.const <- function(a,lower,upper){
  lb1 <- lower[1]
  lb2 <- lower[2]
  ub1 <- upper[1]
  ub2 <- upper[2]
  EM <- EM.pu(a=a,lam=c(1,1))
  b <- 1/a
  
  den1 <- EM*(1+sum(b))*(1+b[2])
  num1 <- b[1]*(1-exp(-(1+b[2])*ub1))*(1-exp(b[2]*lb2))
  
  den2 <- EM*(1+sum(b))*(1+b[1])
  num2 <- (1+b[1])*b[2]*(1-exp(-ub2))-b[2]*exp(b[1]*lb1)*(1-exp(-(1+b[1])*ub2))
  
  den3 <- EM*(1+sum(b))*(1+b[2])
  num3 <- prod(b)*(1-exp(-ub2))-b[1]*(exp(-ub1)-exp(-(1+b[2])*ub1))
  
  return(num1/den1+num2/den2 + num3/den3)
  
}

trunc.norm.const.GPD <- function(a, sig, gamma, lower, upper ){
  upp.para1 <- if (gamma[1] <0) -sig[1]/gamma[1] else Inf  
  upp.para2 <- if (gamma[2] <0) -sig[2]/gamma[2] else Inf 
  upper <- pmin(upper, c(upp.para1, upp.para2))
  
  low.para1 <- if( gamma[1]>0) -sig[1]/gamma[1] else -Inf  
  low.para2 <- if( gamma[2]>0) -sig[2]/gamma[2] else -Inf  
  lower <- pmax(lower, c(low.para1,low.para2) )
  # can't do it vectorized because the formula for dim 1 and dim 2 may differ depending on the value of gamma
  for (i in 1:length(sig)){
    if (abs(gamma[i])<10^-6){
      lower[i] <- lower[i]/sig[i]
      upper[i] <- upper[i]/sig[i]
    }else{
      if(is.na(log(gamma[i]/sig[i]*upper[i]+1+ 10^-6))) print(list(gamma=gamma,sig=sig,upper=upper))
      lower[i] <- log(gamma[i]/sig[i]*lower[i]+1 + 10^-6)/gamma[i]
      upper[i] <- log(gamma[i]/sig[i]*upper[i]+1 + 10^-6)/gamma[i]
      
    }
  }
  return(trunc.norm.const(a,lower,upper))
}

nim_trunc.norm.const.GPD <- nimbleRcall(function(a=double(1), sig=double(1), gamma=double(1),
                                                 lower=double(1), upper=double(1)){}, 
                                        Rfun = 'trunc.norm.const.GPD',
                                        returnType = double(0))

dbiextmix.sg <- nimbleFunction(
  run = function(x=double(2), thres=double(1), params.bulk=double(1), bulk.dist.name=double(1),
                 theta=double(1),  
                 lower=double(1),upper=double(1),
                 log = logical(0, default = 0)) {
    returnType(double(0))
    
    cond <- (x[,1]>thres[1]) | (x[,2]>thres[2])
    n.tail <- sum(cond)
    n.bulk <- sum(!cond)
    y.tail <- matrix(c(x[cond,1] - thres[1], x[cond,2] - thres[2]), ncol=2)
    y.bulk <- x[!cond,]
    
    lbound.bulk <- lower[1:2]
    lbound.tail <- lower[3:4]
    
    ubound.bulk <- upper[1:2]
    ubound.tail <- upper[3:4]
    
    dbulk.res <- dbulk(y.bulk, params.bulk, bulk.dist.name ,lbound.bulk,thres)
    pi <- dbulk.res[1]
    
    sig <- theta[3:4]
    gamma <- theta[5:6]
    eta <- -sig/gamma
    eta[which(gamma<=0)] <- -Inf
    
    ll.tail <- -10^10
    ll.bulk <- -10^10
    
    if (n.tail>0){
      y.min <- eta
      for (i in 1:2){
        y.min[i] <- min(y.tail[,i])
      }
      if (all(y.min>eta)){
        ll.tail <- -nim_nll_powunif_GPD_MAT(x=y.tail, theta=theta)
        
        den <- nim_trunc.norm.const.GPD(a=theta[1:2], sig=theta[3:4], gamma=theta[5:6],
                                        lower=lbound.tail-thres,upper=ubound.tail-thres)
        
        ll.tail <- ll.tail - n.tail*log(den)
      }else{
        ll.tail <- -10^10
      }
      
      ll.all <- n.tail*log(1-pi) + ll.tail
    }
    
    if (n.bulk>0){
      ll.bulk <- max(dbulk.res[2],ll.bulk)
      ll.all <-  ll.bulk 
    }
    if (log) {
      totalProb <- ll.all
    }else{
      totalProb <- exp(ll.all)
    }
    
    return(totalProb)
  })

x <- Y[1:2,]
theta=c(a,sig,gamma)
thres=u.x
params.bulk=par
bulk.dist.name=c(3,3,5)
lower=rep(0,4)
upper=rep(Inf,4)
log = TRUE


dbiextmix.sg(x=Y, thres=thres, params.bulk=par, bulk.dist.name=c(3,3,5),
          theta=theta,  
          lower=rep(0,4),upper=rep(Inf,4),
          log = TRUE)

#log-likelihood with standard GPD tail, i.e. sigma=c(1,1), gamma=c(0,0)
log_likelihood <- function(x){
  Y <- rbind(x,x)
  theta=c(a,sig,gamma)
  thres=u.x
  params.bulk=par
  bulk.dist.name=c(3,3,5)
  lower=rep(0,4)
  upper=rep(Inf,4)
  log = TRUE
  ll <- dbiextmix.sg(x=Y, thres=thres, params.bulk=par, bulk.dist.name=c(3,3,5),
                  theta=theta,  
                  lower=rep(0,4),upper=rep(Inf,4),
                  log = TRUE)
  return(ll/2)
}

log_likelihood(Y[1,])



# Step 2: Create a grid of x and y values
x_seq <- seq(0, 12, 0.1)
y_seq <- seq(0, 12, 0.1)
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
  geom_segment(aes(x = u.x[1], y = 0, xend = u.x[1], yend = u.x[2]), color = "black", linetype = "dashed") +
  geom_segment(aes(x = 0, y = u.x[2], xend = u.x[1], yend = u.x[2]), color = "black", linetype = "dashed") +
  # annotate("rect", xmin = 5.7, xmax = 10, ymin = -2, ymax = 12,
  #          alpha = .05,fill = contour.palette[5])+
  # annotate("rect", xmin = -0.5, xmax = 5.7, ymin = 7.5, ymax = 12,
  #          alpha = .05,fill = contour.palette[5])+
  # annotate("rect", xmin = -0.5, xmax = 5.7, ymin = -2, ymax = 7.5,
  #          alpha = .05,fill = contour.palette[7])+
  labs(x = "X1",
       y = "X2",
       color = "Density") +
  xlim(0,11)+
  ylim(0,11)
  # theme(axis.text.x=element_text(size=15),
  #       axis.text.y=element_text(size=15),
  #       axis.title.x=element_text(size=15),
  #       axis.title.y=element_text(size=15))
print(p)

png(filename = file.path("/home/pgrad2/2448355h/My_PhD_Project/01_Output/Biv_Ext_Mix_Mod/Plots", "BEMM_contour.png"), width = 6*res, height = 5*res, res=res)
print(p)
dev.off()