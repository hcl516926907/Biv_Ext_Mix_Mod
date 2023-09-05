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
packages <- c( "mvtnorm", "tmvtnorm","foreach","doSNOW","parallel")  


load_install_packages(packages)

print(detectCores())


NumberOfCluster <- 40
cl <- makeCluster(NumberOfCluster)
registerDoSNOW(cl)


t1 <- Sys.time()
outputM1 <-
  foreach(i = 1:80, .packages = c('mvtnorm','tmvtnorm',"copula","pracma","evd")) %dopar%{
    dir.work <- '/home/pgrad2/2448355h/My_PhD_Project/Biv_Ext_Mix_Mod'
    source(file.path(dir.work, 'Simulation/Functions.R'))
    seed <- i
    d <- 2
    mu <- c(3.5, 4.0)
    sd1 <- 1
    sd2 <- 1.5
    rho <- 0.7
    sigma <- matrix(c(sd1^2, rho*sd1*sd2, rho*sd1*sd2, sd2^2),ncol=2)
    n <- 2000
    set.seed(seed)
    Y <- rmvnorm(n, mean=mu, sigma=sigma)
    Fn1 <- ecdf(Y[,1])
    Fn2 <- ecdf(Y[,2])
    u1 <- Fn1(Y[,1])*n/(n+1)
    u2 <- Fn2(Y[,2])*n/(n+1)
    U <- cbind(u1,u2)
    
    thresh<-seq(0.01,0.99,len=100)
    start<-Sys.time()
    optM1<-optim(loglik,par=c(3,0,3,1),x=U,ct="InverGumbel",cb="t",
                 weightfun=function(u,v,theta){pifun(u,v,theta)})
    end<-Sys.time()-start
    KM1<-K(parct=optM1$par[1],parcb=optM1$par[2:3],
           ct="InverGumbel",cb="t",weightfun=function(u,v,theta){pifun(u,v,theta)},parweight=optM1$par[4])
    survM1<-survivalf(x=U, parct=optM1$par[1],parcb=optM1$par[2:3],
                      ct="InverGumbel",cb="t",weightfun=function(u,v,theta){pifun(u,v,theta)},parweight=optM1$par[4],K=KM1)
    aicM1<-AICm(optM1$par,-optM1$value)
    return(list("Estimates M1"=optM1$par,"Estimated Loglikelihood M1"=optM1$value,"AIC M1"=aicM1,"Time"=end,"KM1"=KM1,"Survival M1"=survM1))
  }
stopCluster(cl)
t2 <- Sys.time()
print(t2-t1)

dir.out <- "/home/pgrad2/2448355h/My_PhD_Project/01_Output/Biv_Ext_Mix_Mod/Simulation"
save(outputM1,  file=file.path(dir.out, filename='Lidia_model_1_80.RData'))

# Kconst<-NULL
# est<-matrix(NA,ncol=4,nrow=length(outputM1))
# surv<-list()
# for(i in 1:length(outputM1)){
#   est[i,]<-outputM1[[i]][[1]]
#   Kconst[i]<-outputM1[[i]][[5]]
#   surv[[i]]<-outputM1[[i]][[6]]
# }
# 
# 
# chi0.65<-chi0.7<-chi0.75<-chi0.8<-chi0.85<-chi0.9<-chi0.95<-chi0.99<-eta0.65<-eta0.7<-eta0.75<-eta0.8<-eta0.85<-eta0.9<-eta0.95<-eta0.99<-c()
# for(i in 1:length(outputM1)){
#   chi0.65[i]<-chiu_model(thresh=thresh[66],survP=surv[[i]][66])
#   chi0.7[i]<-chiu_model(thresh=thresh[71],survP=surv[[i]][71])
#   chi0.75[i]<-chiu_model(thresh=thresh[76],survP=surv[[i]][76])
#   chi0.8[i]<-chiu_model(thresh=thresh[81],survP=surv[[i]][81])
#   chi0.85[i]<-chiu_model(thresh=thresh[86],survP=surv[[i]][86])
#   chi0.9[i]<-chiu_model(thresh=thresh[91],survP=surv[[i]][91])
#   chi0.95[i]<-chiu_model(thresh=thresh[96],survP=surv[[i]][96])
#   chi0.99[i]<-chiu_model(thresh=thresh[100],survP=surv[[i]][100])
#   eta0.65[i]<-etau_model(thresh=thresh[66],survP=surv[[i]][66])
#   eta0.7[i]<-etau_model(thresh=thresh[71],survP=surv[[i]][71])
#   eta0.75[i]<-etau_model(thresh=thresh[76],survP=surv[[i]][76])
#   eta0.8[i]<-etau_model(thresh=thresh[81],survP=surv[[i]][81])
#   eta0.85[i]<-etau_model(thresh=thresh[86],survP=surv[[i]][86])
#   eta0.9[i]<-etau_model(thresh=thresh[91],survP=surv[[i]][91])
#   eta0.95[i]<-etau_model(thresh=thresh[96],survP=surv[[i]][96])
#   eta0.99[i]<-etau_model(thresh=thresh[100],survP=surv[[i]][100])
# }
# 
# 
# chi0.65geral<-(1-2*thresh[66]+pCopula(cbind(thresh[66],thresh[66]),normalCopula(rho)))/(1-thresh[66])
# chi0.7geral<-(1-2*thresh[71]+pCopula(cbind(thresh[71],thresh[71]),normalCopula(rho)))/(1-thresh[71])
# chi0.75geral<-(1-2*thresh[76]+pCopula(cbind(thresh[76],thresh[76]),normalCopula(rho)))/(1-thresh[76])
# chi0.8geral<-(1-2*thresh[81]+pCopula(cbind(thresh[81],thresh[81]),normalCopula(rho)))/(1-thresh[81])
# chi0.85geral<-(1-2*thresh[86]+pCopula(cbind(thresh[86],thresh[86]),normalCopula(rho)))/(1-thresh[86])
# chi0.9geral<-(1-2*thresh[91]+pCopula(cbind(thresh[91],thresh[91]),normalCopula(rho)))/(1-thresh[91])
# chi0.95geral<-(1-2*thresh[96]+pCopula(cbind(thresh[96],thresh[96]),normalCopula(rho)))/(1-thresh[96])
# chi0.99geral<-(1-2*thresh[100]+pCopula(cbind(thresh[100],thresh[100]),normalCopula(rho)))/(1-thresh[100])
# 
# taumodel<-c()
# nsim<-50000
# n<-100000
# for(i in 1:length(outputM1)){
#   cstar<-sim(n=n,nsim=nsim,parct=est[i,1],parcb=c(est[i,2],est[i,3]),ct="InverGumbel",cb="t",
#              weightfun=function(u,v,theta){pifun(u,v,theta)},parweight=est[i,4])
#   taumodel[i]<-cor(cstar[,1],cstar[,2],method="kendall")
# }
# kendall's \tau for the gaussian with \rho=0.65
# kendalgeral<-2/pi*asin(rho)