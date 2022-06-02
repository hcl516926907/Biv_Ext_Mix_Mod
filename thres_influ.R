source("KRSW/RevExp_U_Functions.r")
source("KRSW/CommonFunctions.r")
source("KRSW/ModelDiagnosticsNewNames.r")


dir.in <- '/home/pgrad2/2448355h/My_PhD_Project/01_Output/Biv_Ext_Mix_Mod/biv_ext_mix_mod_simdat'

load(file=file.path(dir.in,'simulation_data.RData'))


#--------------------------------check threshold for each margin---------------------------------

u.x <- u.x.p10
u.z <- u.z.p10

X <- X.p10
Z <- Z.p10


q1.x <- sum((X[,1]> u.x[1]))/length(X[,1]) 
q2.x <- sum((X[,2]> u.x[2]))/length(X[,2]) 

q1.z <- sum((Z[,1]> u.z[1]))/length(Z[,1]) 
q2.z <- sum((Z[,2]> u.z[2]))/length(Z[,2]) 
library(ismev)
q <- 0.95
dat <- Z[,2]
u <- quantile(dat,q)
plot(dat)
gpd.fitrange(dat, umin=8,umax=12,nint=20)
m1<-gpd.fit(dat,thresh=u)
gpd.diag(m1)
