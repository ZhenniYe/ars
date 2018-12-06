context('Primary ars() Function Testing')
set.seed(2018)

test_that("Throw error message for input error",{
  expect_error(ars('1000', dnorm), "Please provide the numeric sample size n_iter")
  expect_error(ars(-100, dnorm), 'Please provide positive integer for sample size')
  expect_error(ars(66.66, dnorm), "Please provide positive integer for sample size")
  expect_error(ars(100, 'dnorm'), 'Please provide fn as an function')
  x <- 6
  expect_error(ars(100, x^2), 'Please provide fn as an function')
  expect_error(ars(100, dnorm, l='-1', u='1'), 'Please provide the numeric boundary')
  expect_error(ars(100, dnorm, l=-1, u='1'), 'Please provide the numeric boundary')
  expect_error(ars(100, dnorm, l='-1', u=1), 'Please provide the numeric boundary')
  expect_error(ars(100, dnorm, l=1, u=1), 'Please provide the valid boundary')
  expect_warning(ars(100, dnorm, l=5, u=-5), 'The lower bound should be less than the upper bound, swaping values')
})


test_that("Throw error message when density is not log-concave", {
  chi_test <- function(x){return(dchisq(x, 2))}
  expect_error(ars(1000, chi_test, l=1), 'Please provide the log-concave density')
  chi_test2 <- function(x){return(dchisq(x, 1))}
  expect_error(ars(100, chi_test2, l=2), 'Please provide the log-concave density')

  t_test <- function(x){return(dt(x, 3))}
  expect_error(ars(10, t_test), 'Please provide the log-concave density')
  t_test2 <- function(x){return(dt(x, 1))}
  expect_error(ars(100, t_test2), 'Please provide the log-concave density')

  f_test <- function(x){return(df(x, 9, 11))}
  expect_error(ars(100, f_test, l=1), 'Please provide the log-concave density')
  f_test2 <- function(x){return(df(x, 1, 2))}
  expect_error(ars(1000, f_test2, l=3), 'Please provide the log-concave density')

  other_test <- function(x){return(exp(x^2))}
  expect_error(ars(10, other_test, l=-5, u=5), 'Please provide the log-concave density')
  other_test2 <- function(x){return(exp(x^3))}
  expect_error(ars(100, other_test2, l=-5, u=8), 'Please provide the log-concave density')
})


test_that("Throw error message when function is not differentiable",{

})




