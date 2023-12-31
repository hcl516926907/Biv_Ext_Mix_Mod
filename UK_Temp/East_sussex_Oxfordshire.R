dir.work <- '/home/pgrad2/2448355h/My_PhD_Project/Biv_Ext_Mix_Mod'
dir.data <- '/home/pgrad2/2448355h/My_PhD_Project/00_Dataset/UK_Temp'
dir.out <- "/home/pgrad2/2448355h/My_PhD_Project/01_Output/Biv_Ext_Mix_Mod/UK_Temp"

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
packages <- c("nimble", "mvtnorm", "tmvtnorm","foreach","doSNOW","parallel", "dplyr",'posterior','tseries')  


load_install_packages(packages)

load(file=file.path(dir.data, "east-sussex_oxfordshire.RData"))

temp.daily.max.all.flt <- temp.daily.max.all[(temp.daily.max.all$Year>=2016),]
Y1 <- data.frame(temp=temp.daily.max.all.flt$air_temp_city1, 
                 time= 1:nrow(temp.daily.max.all.flt),
                 date = temp.daily.max.all.flt$Date,
                 month = temp.daily.max.all.flt$Month,
                 lag1=temp.daily.max.all.flt$city1_lag1,
                 lag2=temp.daily.max.all.flt$city1_lag2,
                 lag3=temp.daily.max.all.flt$city1_lag3,
                 lag4=temp.daily.max.all.flt$city1_lag4,
                 lag5=temp.daily.max.all.flt$city1_lag5)
acf(Y1$temp)

model_1 <- lm(temp ~  sin(2*pi*time/365) + cos(2*pi*time/365)  + lag1 , data=Y1)
model_1.1 <- lm(temp ~ poly(time,2) +  sin(2*pi*time/365) + cos(2*pi*time/365)  + lag1 + lag2 + lag3 + lag4 + lag5, data=Y1)

y1 <- model_1$residuals
n1 <- length(y1)
emp_quant1 <- quantile(y1, probs = seq(0, 1, length.out = n1))
thy_quant1 <- qnorm(seq(0, 1, length.out = n1), mean = mean(y1), sd = sd(y1))
qqplot(thy_quant1, emp_quant1,xlab="Theoretical Quantiles", ylab="Sample Quantiles", main="QQ Plot of Y1")
abline(a=0,b=1)


Y2 <- data.frame(temp = temp.daily.max.all.flt$air_temp_city2, 
                 time = 1:nrow(temp.daily.max.all.flt),
                 date = temp.daily.max.all.flt$Date,
                 month = temp.daily.max.all.flt$Month,
                 lag1=temp.daily.max.all.flt$city2_lag1)


model_2 <- lm(temp ~ sin(2*pi*time/365) + cos(2*pi*time/365)  + lag1 , data=Y2)
summary(model_2)
model_2.1 <- lm(temp ~ poly(time,2) + sin(2*pi*time/365) + cos(2*pi*time/365)  + lag1 , data=Y2)
summary(model_2.1)
plot(density(model_2$residuals))

y2 <- model_2$residuals
n2 <- length(y2)
emp_quant2 <- quantile(y2, probs = seq(0, 1, length.out = n2))
thy_quant2 <- qnorm(seq(0, 1, length.out = n2), mean = mean(y2), sd = sd(y2))
qqplot(thy_quant2, emp_quant2,xlab="Theoretical Quantiles", ylab="Sample Quantiles", main="QQ Plot of Y2")
abline(a=0,b=1)

acf(model_2$residuals)
adf.test(model_2$residuals)
kpss.test(model_2$residuals)
Box.test(model_2$residuals, lag=log(length(model_2$residuals)))

Y <- cbind(model_1$residuals,model_2$residuals)


thres1 <- quantile(Y[,1],0.025)
thres2 <- quantile(Y[,2],0.025)
extr.cond.joint <- (Y[,1] < thres1)&(Y[,2] < thres1)
extr.cond.par1 <- Y[,1]<thres1
extr.cond.par2 <- Y[,2]<thres2
Y1[extr.cond.joint,'month']
month.order <- c("Jan",'Feb','Mar','Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec')

Tab.joint <- table(factor(Y1[extr.cond.joint,'month'], levels=month.order))
barplot(Tab.joint,las=2)

Tab.par1 <- table(factor(Y1[extr.cond.par1,'month'], levels=month.order))
barplot(Tab.par1,las=2)

Tab.par2 <- table(factor(Y1[extr.cond.par2,'month'], levels=month.order))
barplot(Tab.par2,,las=2)

Y.fit <- -Y

NumberOfCluster <- 3
cl <- makeCluster(NumberOfCluster)
registerDoSNOW(cl)

# source(file.path(dir.work, 'UK_Temp/BEMM_Function_Temp.R'))
# run_MCMC_parallel(seed=3, dat=Y.fit, niter=20000, nburnin = 10000, thin=10)
# 

t1 <- Sys.time()
chain_res <-
  foreach(j = 1:3, .packages = c('nimble','mvtnorm','tmvtnorm')) %dopar%{
    
    run_MCMC_parallel(seed=j, dat=Y.fit, niter=20000, nburnin = 10000, thin=10)
  }
stopCluster(cl)
t2 <- Sys.time()
print(t2-t1)

# save(chain_res, file=file.path(dir.out, filename='east-sussex_oxfordshire_0.8_0.99.RData'))
## v1 corrects the prior of the mu.
# save(chain_res, file=file.path(dir.out, filename='east-sussex_oxfordshire_0.8_0.99_v1.RData')) 
## v2 add the constraints to obtain a finite marginal expectation
save(chain_res, file=file.path(dir.out, filename='east-sussex_oxfordshire_0.8_0.99_v2.RData')) 