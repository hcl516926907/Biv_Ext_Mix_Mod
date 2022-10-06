library(ggplot2)
library(dplyr)
library(GPfit)
library(rBayesianOptimization)
#example function
f <- function(x) {
  f <- (2 * x - 10)^2 * sin(32 * x - 4)
  return(f)
}


# Create noise-free f for n0 based on 5 points within range of [0,1].
x <- c(0, 1/3, 1/2, 2/3, 1)

eval <- data.frame(x = x, y = f(x)) %>% as.matrix()
eval


#Create a gaussian process with GP_fit() with power exponential correlation function. 
#You can also use Matern correlation function list(type = "matern", nu = 5/2).

fit <- GP_fit(X = eval[ , "x"], 
              Y = eval[ , "y"], 
              corr = list(type = "exponential", power = 1.95))


# After we fitted GP model, we can calculate the expected value μ(x) at 
# each possible value of x and the corresponding uncertainty σ(x). 
# These will be used when computing the acquisition functions over the 
# possible values of x.
x_new <- seq(0, 1, length.out = 100)
pred <- predict.GP(fit, xnew = data.frame(x = x_new))
mu <- pred$Y_hat
sigma <- sqrt(pred$MSE)



ggplot(as.data.frame(eval))+
  geom_line(data = data.frame(x = x_new, y = mu),
            aes(x = x, y = y), color = "red", linetype = "dashed")+
  geom_ribbon(data = data.frame(x = x_new, y_up = mu + sigma, y_low = mu - sigma), 
              aes(x = x_new, ymax = y_up, ymin = y_low), fill = "skyblue", alpha = 0.5) +
  geom_point(aes(x,y), size = 2)+
  theme_minimal() +
  labs(title = "Gaussian Process Posterior of f(x)",
       subtitle = "Blue area indicate the credible intervals",
       y = "f(x)")

y_best <- min(eval[,2])

eps <- 0.01
ei_calc <- function(m, s) {
  if (s == 0) {
    return(0)
  }
  Z <- (m - y_best - eps)/s
  expected_imp <- (m - y_best - eps) * pnorm(Z) + s * dnorm(Z)
  return(expected_imp)
}

expected_improvement <- numeric()
for (i in 1:length(mu)) {
  expected_improvement[i] <- ei_calc(m = mu[i],s =  sigma[i])
}


exp_imp <- data.frame(x = x_new,
                      y = expected_improvement)

exp_best <- exp_imp %>% filter(y == max(y))



ggplot(exp_imp, aes(x, y))+
  geom_line()+
  geom_ribbon(aes(ymin = 0, ymax = y), fill = "skyblue", alpha = 0.5, color = "white")+ 
  geom_vline(xintercept = exp_best$x, linetype = "dashed", color = "red")+
  geom_point(data = exp_best, size = 2)+
  theme_minimal() +
  theme(panel.grid = element_blank())+
  scale_x_continuous(breaks = c(seq(0,1,0.25), round(exp_best$x,2)))+
  labs(title = "Expected Improvement",
       subtitle = "x with the highest expected improvement will be evaluated",
       y = "Expected Improvement")


f.example <- function(x,y) {
  
  res <- (2 * x - 10)^2 * sin(32 * x - 4) + y
  return(list(Score = res, Pred = 0))
}

search_bound <- list(x = c(0,1),
                     y = c(0,1))

set.seed(123)
search_grid <- data.frame(x = runif(20,0,1),
                          y = runif(20,0,1))
head(search_grid)

set.seed(1)
bayes_finance_ei <- BayesianOptimization(FUN = f.example, bounds = search_bound, 
                                         init_grid_dt = search_grid, init_points = 0, 
                                         n_iter = 10, acq = "ucb")




f.example <- function(x) {
  
  res <- (2 * x - 10)^2 * sin(32 * x - 4)
  return(list(Score = res, Pred = 0))
}

which(f.example(1:1000)$Score==max(f.example(1:1000)$Score))

search_bound <- list(x = c(0L,1000L))

set.seed(123)
search_grid <- data.frame(x = sample(1:1000,20))
head(search_grid)

set.seed(1)
bayes_finance_ei <- BayesianOptimization(FUN = f.example, bounds = search_bound, 
                                         init_grid_dt = search_grid, init_points = 0, 
                                         n_iter = 20, acq = "ucb")
