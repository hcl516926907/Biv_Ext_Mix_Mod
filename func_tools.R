


# find the quantile of the ecdf, i.e. applying F(X) to X
my.ecdf <- function(x){
  #adjust to avoid having F(x)=1, according to p36 in cole's book.
  ecdf <- rank(x,ties.method='max')/(length(x)+1)
  return(ecdf)
}



tail.coef <- function(x,y,u){
  n <- length(x)
  
  x.0 <- which(x<=u)
  x.1 <- which(x>u)
  y.0 <- which(y<=u)
  y.1 <- which(y>u)
  
  x.p.0 <- length(x.0)/n
  x.p.1 <- length(x.1)/n
  y.p.0 <- length(y.0)/n
  y.p.1 <- length(y.1)/n
  
  #conditional distribution
  p.0c0 <- length(intersect(x.0,y.0))/length(y.0)
  p.0c1 <- length(intersect(x.0,y.1))/length(y.1)
  p.1c0 <- length(intersect(x.1,y.0))/length(y.0)
  p.1c1 <- length(intersect(x.1,y.1))/length(y.1)
  
  res <- c(x.p.0, x.p.1, y.p.0, y.p.1,
           p.0c0, p.0c1, p.1c0, p.1c1)
  
  names(res) <- c('x.p.0', 'x.p.1', 'y.p.0', 'y.p.1',
                  'p.0c0', 'p.0c1', 'p.1c0', 'p.1c1')
  return(res)
}