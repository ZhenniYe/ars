context('Primary ars() Function Testing')
set.seed(2018)


## functions used for testing
normal_test <- function(x){return(dnorm(x, mean = 50, sd = 10))}
uniform_test <- function(x){return(dunif(x, -1, 1))}
gamma_test <- function(x){return(dgamma(x, 6))}
gamma_test2 <- function(x){return(dgamma(x, 2))}
gamma_test3 <- function(x){return(dgamma(x, 10))}
beta_test <- function(x){return(dbeta(x, shape1=3, shape2=5))}
beta_test2 <- function(x){return(dbeta(x, shape1=1, shape2=2))}
chi_test <- function(x){return(dchisq(x, 1))}
chi_test2 <- function(x){return(dchisq(x, 6))}
chi_test3 <- function(x) {return(dchisq(x, df = 10))}
t_test <- function(x){return(dt(x, 3))}
t_test2 <- function(x){return(dt(x, 1))}
f_test <- function(x){return(df(x, 9, 11))}
f_test2 <- function(x){return(df(x, 1, 2))}
binom_test <- function(x){return(dbinom(x, 50, 0.3))}
other_test <- function(x){return(exp(x^2))}
other_test2 <- function(x){return(exp(x^3))}


test_that("Main function testing",{
  print("test Normal Distribution")
  expect_true(ks.test(ars(1000, dnorm), pnorm)$p.value > 0.05)
  expect_true(ks.test(ars(1000, normal_test, l = 0, u = 100),
                      function(x) {pnorm(x, mean=50, sd=10)})$p.value > 0.05)

  print("test Beta Distribution")
  expect_true(ks.test(ars(1000, beta_test, 0.01, 0.99),
                      function(x) {pbeta(x, 3, 5)})$p.value > 0.05)
  expect_true(ks.test(ars(1000, beta_test2, 0.01, 0.99),
                      function(x) {pbeta(x, 1, 2)})$p.value > 0.05)

  print("test Gamma Distribution")
  expect_true(ks.test(ars(1000, gamma_test, l = 0.1, u = 20),
                      function(x) {pgamma(x, 6)})$p.value > 0.05)
  expect_true(ks.test(ars(1000, gamma_test3, l = 1),
                      function(x) {pgamma(x, 10)})$p.value > 0.05)

  print("test Chi-square Distribution")
  expect_true(ks.test(ars(1000, chi_test2, l = 1),
                      function(x) {pchisq(x, df=6)})$p.value > 0.05)
  expect_true(ks.test(ars(1000, chi_test3, l = 1, u = 100),
                      function(x) {pchisq(x, df=10)})$p.value > 0.05)

  print("test Uniform Distribution")
  expect_true(ks.test(ars(1000, dunif, 0.01, 0.99), punif)$p.value > 0.05)
  expect_true(ks.test(ars(1000, uniform_test, -0.99, 0.99),
                      function(x) {punif(x, -1, 1)})$p.value > 0.05)

  cat('\nMain function tests passed\n')
})



test_that("Output sample points within the boundary",{
  print("test Normal Distribution")
  expect_true(min((ars(1000, dnorm, l = -1, u = 2))) >= -1 &&
                max((ars(1000, dnorm, l = -1, u = 2))) <= 2)

  print("test Beta Distribution")
  expect_true(min((ars(1000, beta_test, l = 0.1, u = 0.9))) >= 0.1 &&
                max((ars(1000, beta_test, l = 0.1, u = 0.9))) <= 0.9)

  print("test Gamma Distribution")
  expect_true(min((ars(1000, gamma_test, l = 1, u = 100))) >= 1 &&
                max((ars(1000, gamma_test, l = 1, u = 100))) <= 100)

  print("test Chi-square Distribution")
  expect_true(min((ars(1000, chi_test2, l = 0.01))) >= 0.01 &&
                max((ars(1000, chi_test2, l = 0.01))) <= Inf)

  print("test Uniform Distribution")
  expect_true(min((ars(1000, dunif, l = 0.1, u = 0.9))) >= 0.1 &&
                max((ars(1000, dunif, l = 0.1, u = 0.9))) <= 0.9)

  cat("\nBoundary tests passed\n")
})



test_that("Throw error message for input error",{
  expect_error(ars('1000', dnorm), "Please provide numeric sample size N")
  expect_error(ars(-100, dnorm), 'Please provide positive integer for sample size')
  expect_error(ars(66.66, dnorm), "Please provide positive integer for sample size")
  expect_error(ars(100, 'dnorm'), 'Please provide fn as an function')
  expect_error(ars(100, 12), 'Please provide fn as an function')
  expect_error(ars(100, dnorm, l='-1', u='1'), 'Please provide the numeric boundary')
  expect_error(ars(100, dnorm, l=-1, u='1'), 'Please provide the numeric boundary')
  expect_error(ars(100, dnorm, l='-1', u=1), 'Please provide the numeric boundary')
  expect_error(ars(100, dnorm, l=1, u=1), 'Please provide the valid boundary')
  expect_warning(ars(100, dnorm, l=5, u=-5), 'The lower bound should be less than the upper bound, swaping values')
  expect_error(ars(100, dnorm, center=1e8), 'Please provide valid boundary instead of using the defualt')

  cat('\nInput tests passed\n')
})



test_that("Throw error message when density is not log-concave", {
  print("test Chi-square Distribution")
  expect_error(ars(1000, chi_test, l=1), 'Please provide the log-concave density')
  expect_error(ars(1000, chi_test, l=3, u=5), 'Please provide the log-concave density')

  print("test Student t Distribution")
  expect_error(ars(100, t_test, l = -5, 5), 'Please provide the log-concave density')
  expect_error(ars(100, t_test2,l = -10, u = 10), 'Please provide the log-concave density')

  print("test f Distribution")
  expect_error(ars(100, f_test, l=1), 'Please provide the log-concave density')
  expect_error(ars(1000, f_test2, l=3), 'Please provide the log-concave density')

  print("test exp(x^2)")
  expect_error(ars(10, other_test, l=-5, u=5), 'Please provide the log-concave density')
  expect_error(ars(100, other_test2, l=-5, u=8), 'Please provide the log-concave density')

  cat('\nLog-concave tests passed\n')
})



test_that("Throw error message when function is not differentiable",{
  print("test f Distribution")
  expect_error(ars(100, f_test, l = -1),
               "Input function is not differentiable in the given boundary. Please check for boundary")

  print("test Chi-square Distribution")
  expect_error(ars(100, chi_test, l = -10),
               "Input function is not differentiable in the given boundary. Please check for boundary")

  print("test Beta Distribution")
  expect_error(ars(100, beta_test, u = 2),
               "Input function is not differentiable in the given boundary. Please check for boundary")

  print("test Gamma Distribution")
  expect_error(ars(100, gamma_test, l = -5),
               "Input function is not differentiable in the given boundary. Please check for boundary")

  print("test Binomial Distribution")
  expect_error(ars(100, binom_test, l = -5),
               "Input function is not differentiable in the given boundary. Please check for boundary")

  cat('\nDifferentiability tests passed\n')
})



test_that("Throw error message when function is not continuous",{
  print("test Beta Distribution")
  expect_error(ars(100, beta_test, l = -1, u = 2),
               "Target function DO NOT EXIST in the given boundary")

  print("test chi-square Distribution")
  expect_error(ars(100, chi_test, l = -5, u = 0),
               "Target function DO NOT EXIST in the given boundary")

  print("test Uniform Distribution")
  expect_error(ars(100, dunif, l = -1, u = 0),
               "Target function DO NOT EXIST in the given boundary")

  print("test f Distribution")
  expect_error(ars(100, f_test, l = -10, u = -1),
               "Target function DO NOT EXIST in the given boundary")

  print("test Gamma Distribution")
  expect_error(ars(100, gamma_test, l = -10, u = 2),
               "Target function DO NOT EXIST in the given boundary")

  cat('\nExistence tests passed\n')
})





