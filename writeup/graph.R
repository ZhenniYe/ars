library(stats)
## test 1: standard normal with -inf, inf
set.seed(3)
set1 <- ars(10000, dnorm, l = -Inf, u = Inf)
x1 <- seq(-4, 4, length=100)
hist(set1, prob=TRUE, xlab = "", ylab = "", main = "Standard Normal Distribution")
lines(x1, dnorm(x1), type="l", col="blue", lwd=2)


## test 2 : standard normal with (-1, 2)
set.seed(0)
set2 <- ars(10000, dnorm, l = -1, u = 2)
x2 <- seq(-1, 2, by = 0.05)
hist(set2, prob = TRUE, xlab = "", ylab = "", main = "Truncated Normal Distribution")
lines(x2, dnorm(x2)/(pnorm(max)-pnorm(-1)), type="l", col = "blue", lwd = 2)



## test 3 : Gamma Distribution
set.seed(0)
gamma_test <- function(x){
  k <- dgamma(x, shape=7.5, scale=1)
  return(k)
}
set3 <- ars(10000, gamma_test, 0.1, 20)
x3 <- seq(0.1, 20, by=0.05)
hist(set3, prob=TRUE, xlab = "", ylab = "", main = "Gamma Distrbution")
lines(x3, gamma_test(x3), type="l", col="blue", lwd = 2)


## test 4 : Beta Distrbution
set.seed(3)
beta_test <- function(x){
  c <- dbeta(x, shape1 = 3, shape2 = 5)
  return(c)
}
set4 <- ars(10000, beta_test, 0.01, 0.99)
x4 <- seq(0.01, 0.99, by=0.05)
hist(set4, prob=TRUE, xlab = "", ylab = "", main = "Beta Distrbution")
lines(x4, beta_test(x4), type="l", col="blue", lwd = 2)



## test 5 : Logistic distribution
set.seed(0)
logistic_test <- function(x){
 l <- dlogis(x)
 return(l)
}
set5 <- ars(10000, logistic_test)
x5 <- seq(-10, 10, by=0.05)
hist(set5, prob=TRUE, xlab = "", ylab = "", main = "Logistic Distribution", ylim = c(0, 0.25))
lines(x5, logistic_test(x5), type="l", col="blue", lwd = 2)

## test 6 : Laplace distribution
set.seed(0)












