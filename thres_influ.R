source("KRSW/RevExp_U_Functions.r")
source("KRSW/CommonFunctions.r")
source("KRSW/ModelDiagnosticsNewNames.r")


dir.in <- '/home/pgrad2/2448355h/My_PhD_Project/01_Output/Biv_Ext_Mix_Mod/biv_ext_mix_mod_simdat'

load(file=file.path(dir.in,'simulation_data.RData'))


#--------------------------------check threshold for each margin---------------------------------

library(ismev)
q <- 0.98
dat <- Z.10to1[,2]
u <- quantile(dat,q)
plot(dat)
gpd.fitrange(dat, umin=6,umax=8,nint=20)
m1<-gpd.fit(dat,thresh=u)
gpd.diag(m1)
