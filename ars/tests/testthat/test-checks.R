context('Check Functions Testing')
set.seed(2018)

# functions used for testing
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


test_that("Check Log-concave", {
  # log-concave densities
  expect_true(Check_logconcave(dnorm, c(-10:10)))
  expect_true(Check_logconcave(gamma_test, c(1:5)))
  expect_true(Check_logconcave(gamma_test2, c(1:10)))
  expect_true(Check_logconcave(beta_test, c(0.25, 0.5, 0.75)))
  expect_true(Check_logconcave(beta_test2, c(0.25, 0.5, 0.75)))

  # too wide boundary (long-tail problem)
  expect_false(Check_logconcave(dnorm, c(-1000:1000)))
  expect_false(Check_logconcave(dexp, c(2:100000)))
  expect_false(Check_logconcave(gamma_test, c(1:10000)))

  # Not log-concave densities
  expect_false(Check_logconcave(chi_test, c(1:10)))
  expect_false(Check_logconcave(chi_test, c(100:1000)))
  expect_false(Check_logconcave(t_test, c(-5:5)))
  expect_false(Check_logconcave(t_test2, c(-1:10)))
  expect_false(Check_logconcave(f_test, c(1:10)))
  expect_false(Check_logconcave(f_test2, c(1:10)))
  expect_false(Check_logconcave(other_test, c(-10:10)))
  expect_false(Check_logconcave(other_test2, c(-100:100)))

  cat('\nIndividual log-concave-check function tests passed\n')
})


test_that("Check Differentiability",{
  # differentiable
  expect_equal(Deriv(0, dnorm, -5, 5), 0)
  expect_equal(Deriv(2, function(x) log(dnorm(x)), -Inf, Inf), -2)
  expect_equal(Deriv(2, function(x) log(gamma_test2(x)), 0, 6), -0.5)
  expect_equal(Deriv(1, function(x) log(other_test(x)), -5, 5), 2)
  expect_equal(Deriv(2, function(x) log(other_test2(x)), -1, 6), 12)

  # NOT differentiable
  expect_error(Deriv(-1, function(x) log(f_test(x)), -5, 5), "Target function is not differentiable in the given boundary. Please check for boundary")
  expect_error(Deriv(-0.5, function(x) log(chi_test(x)), -1, 1), "Target function is not differentiable in the given boundary. Please check for boundary")
  expect_error(Deriv(3, function(x) log(beta_test(x)), 0, 5), "Target function is not differentiable in the given boundary. Please check for boundary")
  expect_error(Deriv(-Inf, function(x) log(gamma_test(x)), -Inf, 5), "Target function is not differentiable in the given boundary. Please check for boundary")
  expect_error(Deriv(-5, function(x) log(binom_test(x)), -5, 5), "Target function is not differentiable in the given boundary. Please check for boundary")

  cat('\nIndividual differentiability-check function tests passed\n')
})


test_that("Check for Existance", {
  expect_error(define_check(-1000, other_test), "Target function DO NOT EXIST in the given boundary")
  expect_error(define_check(1000, other_test2), "Target function DO NOT EXIST in the given boundary")

  cat('\nIndividual existance-check function tests passed\n')
})

