library(tmvtnorm)
library(truncnorm)

################################################################
# total sample size: 10000
################################################################


#--------------------tail proportation 20%----------------------
n <- 500
p <- 0.98
#scheme 1 
set.seed(1234)
sp.tail <- rexp(n*(1-p))

set.seed(1234)
sp.bulk <- rnorm(n, mean=10)

u <- quantile(sp.bulk,p)

sp.bulk.1 <- sp.bulk[sp.bulk<u]
sp.tail.1 <- sp.tail + u

plot(density(sp.bulk.1))
plot(density(c(sp.bulk.1, sp.tail.1)), main='scheme 1')


#scheme 2
set.seed(1234)
sp.bulk.2 <- rtmvnorm(n*p, mean=10, sigma=1, upper=u)
set.seed(1234)
sp.bulk.2.1 <- rtruncnorm(n*p,  b=u, mean = 10, sd = 1)
plot(density(c(sp.bulk.2, sp.tail.1)), main='scheme 2')

#plot(density(c(sp.bulk.2.1, sp.tail.1)), main='scheme 2')

qqplot(c(sp.bulk.1, sp.tail.1),c(sp.bulk.2, sp.tail.1), 
       main=paste('Scheme1 against Scheme 2, p=',p,', n=',n, sep=''))
abline(0, 1, col = 'red')
#qqplot(c(sp.bulk.1, sp.tail.1),c(sp.bulk.2.1, sp.tail.1), main='Scheme1 against Scheme 2')


#scheme 3

set.seed(1234)
uni <- runif(n)

set.seed(1234)
sp.1 <- rtmvnorm(n, mean=10, sigma=1, upper=u)

set.seed(1234)
sp.2 <- rexp(n) + u

sp <- c(sp.1[uni<p], sp.2[uni>=p])
p1 <- 0.99
sp1 <- c(sp.1[uni<p1], sp.2[uni>=p1])

plot(density(sp.1[uni<p]))
plot(density(sp), main = 'scheme 3')
plot(density(sp1), main = 'scheme 3.1')


qqplot(c(sp.bulk.1, sp.tail.1),sp, 
       main=paste('Scheme1 against Scheme 3, p=',p,', n=',n, sep=''))
abline(0, 1, col = 'red')

qqplot(c(sp.bulk.1, sp.tail.1),sp1, main=bquote("Scheme1 against Scheme 3," ~ pi == .(p1)))
abline(0, 1, col = 'red')