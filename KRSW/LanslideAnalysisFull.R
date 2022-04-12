source("Functions/Landslide-analysis-functions.R")
library(ismev)
library(numDeriv)
library(MASS)
library(cubature)
library(tailDepFun)

### Read data
dataAll <- read.csv("Data/precipitation.csv")
data <- as.vector(dataAll[,2])
data[is.na(data)]<-0
dataTS <- ts(data,start=1913,frequency=365)
n <- length(data)
mypar <- function() par(cex.lab=2,cex.axis=2,cex.main=1.5,mar=c(5,4.4,4,2)+0.1)

### Plot three-day rainfall amounts
mypar()
data3 <- sapply(c(1:(length(dataTS)-2)), function(i) sum(dataTS[i:(i+2)]))
dataTS3 <- ts(data3,start=1913,frequency=365)
plot(dataTS3, ylab = "rainfall in mm", main = "Three-day rainfall amounts in Abisko")
abline(h = 12, col = "red", lwd = 2)

### Construct the three-dimensional data set
u <- 12
Y <- c()
indx <- indx1 <- indx2 <- indx3 <- 0
r <- 5
i <- 2
while(i < n){
  i <- i + 1
  if(data[i] > u || sum(data[(i-1):i]) > u || sum(data[(i-2):i]) > u){
     if(data[i] > u){imax <- i}
     if(sum(data[(i-1):i]) > u){imax <- i - 3 + which(data[(i-1):i] == max(data[(i-1):i]))[1]}
     if(sum(data[(i-2):i]) > u){imax <- i - 3 + which(data[(i-2):i] == max(data[(i-2):i]))[1]}
     if(max(indx) > (imax-r)){
       cluster <- data[(max(indx)+3):(imax+r)]
     } else{
       cluster <- data[(imax-r):(imax+r)]
     }
    cluster2 <- sapply(c(1:(length(cluster)-1)), function(j) sum(cluster[j:(j+1)]))
    cluster3 <- sapply(c(1:(length(cluster)-2)), function(j) sum(cluster[j:(j+2)]))
    indx1 <- append(indx1,imax-r-1+which(cluster==max(cluster))[1])
    indx2 <- append(indx2,imax-r-1+which(cluster2==max(cluster2)))
    indx3 <- append(indx3,imax-r-1+which(cluster3==max(cluster3)))
    Y <- rbind(Y, c(max(cluster),max(cluster2),max(cluster3)))
    indx <- append(indx,imax)
    i <- i + r
  }
}
indx1 <- indx1[-1]
indxfinal <- indx1[intersect(which(Y[,1] < Y[,2]),which(Y[,2] < Y[,3]))]
dataY <- Y[intersect(which(Y[,1] < Y[,2]),which(Y[,2] < Y[,3])),] # remove observations that are not strictly increasing

##############################################################################################################
################# Test for linear trend #############################
u1 <- 12; u2 <- 13.5; u3 <- 14
X1 <- dataY[,1] - u1 ; X2 <- dataY[,2] - u2 ; X3 <- dataY[,3] - u3
data1 <- X1[X1 > 0] ; data2 <- X2[X2 > 0] ; data3 <- X3[X3 > 0]

timet <- matrix(indxfinal,ncol=1)/length(data)
timex1 <- matrix(timet[which(X1 > 0),],ncol=1)
timex2 <- matrix(timet[which(X2 > 0),],ncol=1)
timex3 <- matrix(timet[which(X3 > 0),],ncol=1)
fit1trend <- gpd.fit(data1, threshold = 0, show = FALSE, ydat = timex1, sigl = 1, siglink = exp)
fit2trend <- gpd.fit(data2, threshold = 0, show = FALSE, ydat = timex2, sigl = 1, siglink = exp)
fit3trend <- gpd.fit(data3, threshold = 0, show = FALSE, ydat = timex3, sigl = 1, siglink = exp)

fit1x <- gpd.fit(data1, threshold = 0, show = FALSE)
fit2x <- gpd.fit(data2, threshold = 0, show = FALSE)
fit3x <- gpd.fit(data3, threshold = 0, show = FALSE)
### the estimates found in Table 3
gamma <- c(fit1x$mle[2],fit2x$mle[2],fit3x$mle[2],fit1x$se[2],fit2x$se[2],fit3x$se[2])
eta <- c(fit1x$mle[1],fit2x$mle[1],fit3x$mle[1],fit1x$se[1],fit2x$se[1],fit3x$se[1])

### Deviances in Table 1 of the supplement
2*(fit1x$nllh-fit1trend$nllh)
2*(fit2x$nllh-fit2trend$nllh)
2*(fit3x$nllh-fit3trend$nllh)

### Test for gamma = 0
exp1 <- fitdistr(data1, "exponential")
exp2 <- fitdistr(data2, "exponential")
exp3 <- fitdistr(data3, "exponential")

2*(fit1x$nllh+exp1$loglik)
2*(fit2x$nllh+exp2$loglik)
2*(fit3x$nllh+exp3$loglik)

##################################################################################################
u <- 24
dataSt <- dataY - u
dataAll <- dataSt[dataSt[,3] > 0,]

### the results from Table 4
result1 <- EstimationOrdered(dataAll, gam0 = TRUE, start = c(1,2,8),lambda1=1)
result2 <- EstimationOrdered(dataAll, gam0 = FALSE, start = c(1,2,8,0.1),lambda1=1)
error1 <- c(0,sqrt(diag(result1$covar)))
error2 <- c(0,sqrt(diag(result2$covar)))

#### Yearly probability of a landslide
y <- c(39.5,56.6,69.9)
sigma <- 10.17
u <- 24
lambda <- c(1,0.84,1.08)
x <- (y - u)/sigma
### takes a minute:
res <- adaptIntegrate(densToInt, lower = c(0,1), upper = exp(x[2:3]), x = x[1], par = lambda, tol = 1e-04)$integral
cnst <- (2*prod(lambda))/((lambda[1]-lambda[2])*sum(1/lambda))

mu <- (nrow(dataAll)/102)*(1-cnst*res)
mu*exp(-mu) # prob of extreme rainfall in any given year
1 - exp(-mu) # prob of at least one extreme rainfall

### marginal qq-plots (supplementary material)
GPD.diag(dataAll,rep(result1$result[3],3),rep(0,3))
GPD.diag(dataAll,rep(result2$result[3],3),rep(result2$result[4],3))

### Figure 6
t <- seq(1,2,length=50)
nsim <- 1000
tmp1 <- tmp2 <- tmp3 <- matrix(0,nrow=nsim,ncol=length(t))
for(j in 1:nsim){ 
  nsample <- sample(1:nrow(dataY),size=nrow(dataY),replace=T)
  newdata <- dataY[nsample,]
  tmp1[j,] <- sapply(t, function(i) diagA.plot(i, u, newdata, 1))
  tmp2[j,] <- sapply(t, function(i) diagA.plot(i, u, newdata, 2))
  tmp3[j,] <- sapply(t, function(i) diagA.plot(i, u, newdata, 3))
}
CIlow1 <- apply(tmp1,2,quantile,0.025)
CIhigh1 <- apply(tmp1,2,quantile,0.975)
CIlow2 <- apply(tmp2,2,quantile,0.025)
CIhigh2 <- apply(tmp2,2,quantile,0.975)
CIlow3 <- apply(tmp3,2,quantile,0.025)
CIhigh3 <- apply(tmp3,2,quantile,0.975)

dev.off()
plot(t, sapply(t, function(i) diagA.plot(i, u, dataY, 1)), ylab = "", ylim = c(0.75,1.7),
     main = expression(paste("Goodness-of-fit diagnostic for ", A[1])))
lines(t, CIlow1, lty = 3, lwd = 2)
lines(t, CIhigh1, lty = 3, lwd = 2)
abline(h = 1, lwd = 2)

plot(t, sapply(t, function(i) diagA.plot(i, u, dataY, 2)), ylab = "", ylim = c(0.75,1.25),
     main = expression(paste("Goodness-of-fit diagnostic for ", A[2])))
lines(t, CIlow2, lty = 3, lwd = 2)
lines(t, CIhigh2, lty = 3, lwd = 2)
abline(h = 1, lwd = 2)

plot(t, sapply(t, function(i) diagA.plot(i, u, dataY, 3)), ylab = "", ylim = c(0.75,1.25), 
     main = expression(paste("Goodness-of-fit diagnostic for ", A[3])))
lines(t, CIlow3, lty = 3, lwd = 2)
lines(t, CIhigh3, lty = 3, lwd = 2)
abline(h = 1, lwd = 2)

# Model-based versus empirical exceedances probabilities
u <- 24
param1 <- 1/(1 + sum(result1$result[1:2]))
param2 <- (1+result1$result[1])/(1 + sum(result1$result[1:2]))
c(param1, param2)
emp1 <- length(which(dataY[dataY[,3] > u, 1] > u))/length(which(dataY[,3] > u))
emp2 <- length(which(dataY[dataY[,3] > u, 2] > u))/length(which(dataY[,3] > u))
c(emp1, emp2)

sigma <- result1$covar[1:2,1:2]
par <- result1$result[1:2]
gradf1 <- function(l2,l3){return(c(1/((l2 + 1 + l2/l3)^2), 1/((l3 + 1 + l3/l2)^2)))}
gradf2 <- function(l2,l3){return(c(-1/(l3*(l2 + 1 + l2/l3)^2), (1 + 1/l2)/((l3 + l3/l2 + 1)^2)))}
lerr1 <- sqrt(t(gradf1(par[1],par[2])) %*% sigma %*% gradf1(par[1],par[2])) #std error of emp1
lerr2 <- sqrt(t(gradf2(par[1],par[2])) %*% sigma %*% gradf2(par[1],par[2])) #std error of emp2


# Model-based versus empirical chi (takes a couple of minutes to plot)
chimod <- c(chi12(c(1,result1$result[1:2])), chi13(c(1,result1$result[1:2])), 
            chi23(c(1,result1$result[1:2])), chi123(c(1,result1$result[1:2])))

chiPlot(dataY[,1:2], ylabel = expression(paste(hat(chi)[12] (q))),
        chimod = chimod[1], nsim = 1000, qmax = 0.99)
chiPlot(dataY[,c(1,3)], ylabel = expression(paste(hat(chi)[13] (q))),
        chimod = chimod[2], nsim = 1000, qmax = 0.99)
chiPlot(dataY[,2:3], ylabel = expression(paste(hat(chi)[23] (q))),
        chimod = chimod[3], nsim = 1000, qmax = 0.99)
chiPlot(dataY, ylabel = expression(paste(hat(chi)[123] (q))),
        chimod = chimod[4], nsim = 1000, qmax = 0.99)

### Goodness-of-fit test from Einmahl et al. (2016)
k <- c(50,75,100,125,150)
n <- nrow(dataY)
indices <- rbind(c(1,1,0),c(1,0,1),c(0,1,1),c(1,1,1))
data <- apply(dataY, 2, function(col) n/(n + 0.5 - rank(col)))
result <- vector('list', length = length(k))
for(i in 1:length(k)){
  result[[i]] <- ECestimator(data, indices, k[i], iterate = TRUE, covMat = TRUE)
  print(k[i]*result[[i]]$value)
}

