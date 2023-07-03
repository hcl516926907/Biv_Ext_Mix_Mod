library(qgam)

dir.data <- "/home/pgrad2/2448355h/My_PhD_Project/00_Dataset"


df <- read.csv(file.path(dir.data, 'Amaurot.csv'))
df$u<-df$WindSpeed*cos(df$WindDirection)
df$v<-df$WindSpeed*sin(df$WindDirection)
df.dropna <- na.omit(df)
# 
# fit <- qgam(Y~s(V1, k=20, bs="cr") + 
#               s(V2, k=20, bs='cr') +
#               s(V3, k=20, bs='cr') + 
#               Season +
#               s(Atmosphere, k=20, bs='cr') + 
#               s(u, k = 20, bs = "cr") + 
#               s(v, k = 20, bs = "cr"),
#             data = df.dropna, 
#             qu = 0.999)
# 
# pred <- predict(fit, newdata=df.dropna)
# mean(df.dropna$Y<pred)
# 

Quan.Est <- function(Y,p){
  Y.sorted <- sort(Y,decreasing =TRUE)
  n <- length(Y.sorted)
  n.prime <- n + 1
  eps <- n.prime*p - floor(n.prime*p) 
  if (p >= n/n.prime){
    print("Extrapolation")
    x.q <- Y.sorted[1] - (Y.sorted[1] -  Y.sorted[2])*log(n.prime*(1-p))
  }else{
    print("Interpolation")
    k <- floor(n.prime*p)
    x.q <- (1-eps)*Y.sorted[n-k+2] + eps*(Y.sorted[n-k+1])
  }
  return(x.q)
}
# Quan.Est(df.dropna$Y, 0.84)

Interpolation <- function(Y.large,Y.small,n,p){
  n.prime <- n+1
  eps <- n.prime*p - floor(n.prime*p) 
  x.q <- (1-eps)*Y.small + eps*Y.large
  return(x.q)
}

Extrapolation <- function(Y.large,Y.small,n,p){
  n.prime <- n+1
  x.q <- Y.large - (Y.large -  Y.small)*log(n.prime*(1-p))
  return(x.q)
}

Quan.reg.int <- function(n,p, df.train, df.pred, dir.out, seed){
  n.prime <- n + 1
  k <- ceiling(n.prime*p)
  q1 <- (k-1)/n.prime
  q2 <- k/n.prime
  
  fit1 <- qgam(Y~s(V1, k=20, bs="cr") + 
                 s(V2, k=20, bs='cr') +
                 s(V3, k=20, bs='cr') + 
                 Season +
                 s(Atmosphere, k=20, bs='cr') + 
                 s(u, k = 20, bs = "cr") + 
                 s(v, k = 20, bs = "cr"),qu = q1)
  fit2 <- qgam(Y~Y~s(V1, k=20, bs="cr") + 
                 s(V2, k=20, bs='cr') +
                 s(V3, k=20, bs='cr') + 
                 Season +
                 s(Atmosphere, k=20, bs='cr') + 
                 s(u, k = 20, bs = "cr") + 
                 s(v, k = 20, bs = "cr"),
               data = df.train, 
               qu = q2)
  pred1 <- predict(fit1, newdata=df.pred)
  pred2 <- predict(fit2, newdata=df.pred)
  
  save(fit1,fit2, file=paste(dir.out, '/int_model_',seed,'.RData',sep=''))
  
  if(p>0.5){
    cond <- which(pred2<pred1)
    pred2[cond] <- pred1[cond]
  }else{
    cond <- which[pred1>pred2]
    pred1[cond] <- pred2[cond]
  }
  
  pred.int <- Interpolation(pred2, pred1, n, p)
  return(pred.int)
}

Quan.reg.ext <- function(n, p, df.train, df.pred, dir.out,seed){
  n.prime <- n + 1
  k <- ceiling(n.prime*p)
  q1 <- (n-1)/n.prime
  q2 <- n/n.prime
  
  fit1 <- qgam(Y~s(V1, k=20, bs="cr") + 
                 s(V2, k=20, bs='cr') +
                 s(V3, k=20, bs='cr') + 
                 Season +
                 s(Atmosphere, k=20, bs='cr') + 
                 s(u, k = 20, bs = "cr") + 
                 s(v, k = 20, bs = "cr"),
               data = df.train, 
               qu = q1)
  fit2 <- qgam(Y~s(V1, k=20, bs="cr") + 
                 s(V2, k=20, bs='cr') +
                 s(V3, k=20, bs='cr') + 
                 Season +
                 s(Atmosphere, k=20, bs='cr') + 
                 s(u, k = 20, bs = "cr") + 
                 s(v, k = 20, bs = "cr"),
               data = df.train, 
               qu = q2)
  save(fit1,fit2, file=paste(dir.out, '/ext_model_',seed,'.RData',sep=''))
  
  pred1 <- predict(fit1, newdata=df.pred)
  pred2 <- predict(fit2, newdata=df.pred)
  
  if(p>0.5){
    cond <- which(pred2<pred1)
    pred2[cond] <- pred1[cond]
  }else{
    cond <- which[pred1>pred2]
    pred1[cond] <- pred2[cond]
  }
  
  pred.ext <- Extrapolation(pred2, pred1, n, p)
}


Quan.reg.est <- function(n, p, df.train, df.pred, dir.out, seed){
  n.prime <- n + 1
  if (p >= n/n.prime){
    print("Extrapolation")
    x.q <- Quan.reg.ext(n, p, df.train, df.pred, dir.out, seed)
  }else{
    print("Interpolation")
    x.q <- Quan.reg.int(n,p, df.train, df.pred, dir.out, seed)
  }
  return(x.q)
}




################################################

# 
# q.seq <- c(0.75,0.9,0.95,0.99,0.999)
# n.seq <- c(10,30,50,100)
# res.mat <- matrix(NA, nrow=length(n.seq), ncol=length(q.seq))
# for (i in 1:length(n.seq)){
#   for (j in 1:length(q.seq)){
#     pred <- Quan.reg.est(n.seq[i], q.seq[j])
#     res.mat[i, j] <- mean(df.dropna$Y<pred)
#     print(c('Done',n.seq[i],q.seq[j]))
#   }
# }
# 
# res.bm <- rep(NA, length(q.seq))
# for (i in 1:length(q.seq)){
#   fit <-  qgam(Y~s(V1, k=20, bs="cr") + 
#                  s(V2, k=20, bs='cr') +
#                  s(V3, k=20, bs='cr') + 
#                  Season +
#                  s(Atmosphere, k=20, bs='cr'),
#                data = df.dropna, 
#                qu = q.seq[i])
#   pred <- predict(fit, newdata=df.dropna)
#   res.bm[i] <- mean(df.dropna$Y<pred)
# }
# 
# plot(NULL, xlab='Theoretical probability',
#            ylab='Model-based probablity')
# for (i in 1:length(n.seq)){
#   lines(q.seq, res.mat[i,], color=i)
# }
# 
# save(res.mat, res.bm, file = file.path(dir.data, 'Output','Challenge1','res_int_ext_quantile.RData'))
###########################################################


#bootstrap average
data <- df
B = 5 # number of bootstrap samples
#nblock = 100 # number of blocks
bsize = 150 # number of elements in each block


nblock = ceiling(nrow(data)/bsize)
blocks = split(1:nrow(data), ceiling(seq_along(1:nrow(data))/bsize))
n.seq <- c(10,30,90)
weight <- rep(0, length(n.seq))
p <- 0.9999
score.mat <- matrix(NA, nrow=B, ncol=length(n.seq))
t1 <- Sys.time()
pred.list <- list()
for(b in 1:B){
  set.seed(b)
  idtest = sample(1:nblock, nblock, replace = T)
  data_test <- data[as.numeric(unlist(blocks[idtest])), ]
  pred1 <- Quan.reg.est(10, p, data, data_test)
  pred2 <- Quan.reg.est(30, p, data, data_test)
  pred3 <- Quan.reg.est(90, p, data, data_test)
  pred.list[[b]] <- cbind(pred1, pred2, pred3)
  
  score1 <- abs(mean(data_test$Y<pred1,na.rm = TRUE)-p)
  score2 <- abs(mean(data_test$Y<pred2,na.rm = TRUE)-p)
  score3 <- abs(mean(data_test$Y<pred3,na.rm = TRUE)-p)
  score.all <- c(score1, score2, score3)
  cond <- which(score.all==min(score.all))
  weight[cond] <- weight[cond] + 1
  score.mat[b,] <- score.all
  print(c('Done ',b,' bootstrap'))
}
t2 <- Sys.time()
print(t2-t1)



#################### Bootstrap to get the prediction and confidence interval ########
library(parallel)
library(qgam)
library(foreach)
library(doSNOW)
library(matrixStats)
library(plotrix)


data_train <- read.csv(file.path(dir.data, 'Amaurot.csv'))
data_train$u<-data_train$WindSpeed*cos(data_train$WindDirection)
data_train$v<-data_train$WindSpeed*sin(data_train$WindDirection)

data_score <- read.csv(file.path(dir.data, 'AmaurotTestSet.csv'))
data_score$u<-data_score$WindSpeed*cos(data_score$WindDirection)
data_score$v<-data_score$WindSpeed*sin(data_score$WindDirection)

dir.out <- '/home/pgrad2/2448355h/My_PhD_Project/01_Output/Data_Challenge/Q1'

bootstrap.fit <- function(x, data_train, data_score, dir.out){
  bsize = 150 # number of elements in each block
  nblock = ceiling(nrow(data_train)/bsize)
  blocks = split(1:nrow(data_train), ceiling(seq_along(1:nrow(data_train))/bsize))
  n.seq <- c(30,90)
  p <- 0.9999
  
  set.seed(x)
  idtest = sample(1:nblock, nblock, replace = T)
  data_bss <- data_train[as.numeric(unlist(blocks[idtest])), ]
  
  pred1 <- Quan.reg.est(30, p, data_bss, data_score, dir.out, seed=x)
  pred2 <- Quan.reg.est(90, p, data_bss, data_score, dir.out, seed=x)
  pred.all <- cbind(pred1, pred2)
  return(pred.all)
}
# 
# t1 <- Sys.time()
# for (x in 1:2){
#   bootstrap.fit(x, data_train, data_score, dir.out)
# }
# t2 <- Sys.time()


t3 <- Sys.time()

cl <- makeCluster(detectCores())
registerDoSNOW(cl)


iterations <- 300
pred.list <- list()
pb <- txtProgressBar(max = iterations, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

result <- foreach(i = 1:iterations,.packages='qgam',.options.snow = opts) %dopar% {
  pred <- bootstrap.fit(i, data_train, data_score, dir.out)
}

# Stop the cluster
stopCluster(cl)

t4 <- Sys.time()
print(t4-t3)

# Print the results
print(result)

save(result, file=file.path(dir.out, 'bootstrap_result.RData'))

load(file=file.path(dir.out, 'bootstrap_result.RData'))
bs.sp.1 <- c()
bs.sp.2 <- c()
for (i in 1:iterations){
  bs.sp.1 <- cbind(result[[i]][,1], bs.sp.1)
  bs.sp.2 <- cbind(result[[i]][,2], bs.sp.2)
}
bs.sp.3 <- 0.5*bs.sp.1 + 0.5*bs.sp.2

ci.low.1 <- rowQuantiles(bs.sp.1, probs=0.25)
ci.mid.1 <- rowQuantiles(bs.sp.1, probs=0.5)
ci.upp.1 <- rowQuantiles(bs.sp.1, probs=0.75)
avg.1 <- mean(ci.upp.1-ci.low.1)

ci.low.2 <- rowQuantiles(bs.sp.2, probs=0.25)
ci.mid.2 <- rowQuantiles(bs.sp.2, probs=0.5)
ci.upp.2 <- rowQuantiles(bs.sp.2, probs=0.75)
avg.2 <- mean(ci.upp.2 - ci.low.2)

ci.low.3 <- rowQuantiles(bs.sp.3, probs=0.25)
ci.mid.3 <- rowQuantiles(bs.sp.3, probs=0.5)
ci.upp.3 <- rowQuantiles(bs.sp.3, probs=0.75)
avg.3 <- mean(ci.upp.3 - ci.low.3)

par(mfrow=c(3,3))
for (i in 1:9){
  plotCI(x=1:3, 
         y=c(ci.mid.1[i], ci.mid.2[i], ci.mid.3[i]), 
         ui=c(ci.upp.1[i], ci.upp.2[i], ci.upp.3[i]), 
         li=c(ci.low.1[i], ci.low.2[i], ci.low.3[i]),
         xaxt = "n", xlab='',ylab='CI', main=paste('CI plot for point',i))
  axis(1, at=1:3, labels=c('pred1', 'pred2', 'pred_avg'))
}



df1 <- data.frame(
  idx = 1:100,
  estimate = apply(bs.sp.1, 1, mean),
  lower = apply(bs.sp.1, 1, quantile, 0.25),
  upper = apply(bs.sp.1, 1, quantile, 0.75),
  model = '30_sample_extrapolation'
)

df2 <- data.frame(
  idx = 1:100,
  estimate = apply(bs.sp.2, 1, mean),
  lower = apply(bs.sp.2, 1, quantile, 0.25),
  upper = apply(bs.sp.2, 1, quantile, 0.75),
  model = '90_sample_extrapolation'
)

df3 <- data.frame(
  idx = 1:100,
  estimate = apply(bs.sp.3, 1, mean),
  lower = apply(bs.sp.3, 1, quantile, 0.25),
  upper = apply(bs.sp.3, 1, quantile, 0.75),
  model = 'weighted_avg_extrapolation'
)

df.all <- rbind(df1,df2,df3)

ggplot(df.all, aes(x = idx, y = estimate, color = model)) +
  geom_point(position = position_dodge(0.5)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, position = position_dodge(0.5)) +
  ylab("Estimate") +
  xlab("Index") +
  theme(legend.position = "bottom") + 
  ylim(c(50,370))

write.csv(df.all,  file=file.path(dir.out, 'Q1_extrapolation_estimation.csv'))
