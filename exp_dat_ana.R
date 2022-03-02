
# data directory
dir.dat <- '/home/pgrad2/2448355h/My_PhD_Project/00_Dataset'

# name of Ruerto Rico rive dataset 
riv.dat.fit.nam <- 'Puerto Rico Fit Dataset.dat'
riv.dat.tes.nam <- 'Puerto Rico Test Data.dat'

# name of Leeds air pollution dataset
pol.dat.fit.nam <- 'Leeds Fit Dataset.dat'
pol.dat.tes.nam <- 'Leeds Test Dataset.dat'


# load all datasets
riv.dat.fit <- read.table(file.path(dir.dat,riv.dat.fit.nam), header=TRUE)
riv.dat.tes <- read.table(file.path(dir.dat,riv.dat.tes.nam), header=TRUE)
riv.dat <- rbind(riv.dat.fit,riv.dat.tes)

pol.dat.fit <- read.table(file.path(dir.dat,pol.dat.fit.nam), header=TRUE)
pol.dat.tes <- read.table(file.path(dir.dat,pol.dat.tes.nam), header=TRUE)
pol.dat <- rbind(pol.dat.fit,pol.dat.tes)

