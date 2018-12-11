context('Primary ars() Function Testing')
set.seed(2018)


test_that("Main function testing",{

  print("test Normal")
  expect_true(ks.test(ars(1000, dnorm), pnorm)$p.value > 0.05)
  expect_true(ks.test(ars(1000, normal_test, l = 0, u = 100),
                      function(x) {pnorm(x, mean=50, sd=10)})$p.value > 0.05)

  print("test Beta")
  expect_true(ks.test(ars(1000, beta_test, 0.01, 0.99),
                      function(x) {pbeta(x, 3, 5)})$p.value > 0.05)
  expect_true(ks.test(ars(1000, beta_test2, 0.01, 0.99),
                      function(x) {pbeta(x, 1, 2)})$p.value > 0.05)

  print("test gamma")
  expect_true(ks.test(ars(1000, gamma_test, l = 0.1, u = 20),
                      function(x) {pgamma(x, 6)})$p.value > 0.05)
  expect_true(ks.test(ars(1000, gamma_test3, l = 1),
                      function(x) {pgamma(x, 10)})$p.value > 0.05)

  print("test Chi square")
  expect_true(ks.test(ars(1000, chi_test2, l = 1),
                      function(x) {pchisq(x, df=6)})$p.value > 0.05)
  expect_true(ks.test(ars(1000, chi_test3, l = 1, u = 100),
                      function(x) {pchisq(x, df=10)})$p.value > 0.05)

  print("test Uniform")
  expect_true(ks.test(ars(1000, dunif, 0.01, 0.99), punif)$p.value > 0.05)
  expect_true(ks.test(ars(1000, uniform_test, -0.99, 0.99),
                      function(x) {punif(x, -1, 1)})$p.value > 0.05)

  cat('\nMain function tests passed\n')
})


test_that("Output sample points within the boundary",{
  print("test Normal")
  expect_true(min((ars(1000, dnorm, l = -1, u = 2))) >= -1 &&
                max((ars(1000, dnorm, l = -1, u = 2))) <= 2)

  print("test Beta")
  expect_true(min((ars(1000, beta_test, l = 0.1, u = 0.9))) >= 0.1 &&
                max((ars(1000, beta_test, l = 0.1, u = 0.9))) <= 0.9)

  print("test Gamma")
  expect_true(min((ars(1000, gamma_test, l = 1, u = 100))) >= 1 &&
                max((ars(1000, gamma_test, l = 1, u = 100))) <= 100)

  print("test Chi square")
  expect_true(min((ars(1000, chi_test2, l = 0.01))) >= 0.01 &&
                max((ars(1000, chi_test2, l = 0.01))) <= Inf)

  print("test Uniform")
  expect_true(min((ars(1000, dunif, l = 0.1, u = 0.9))) >= 0.1 &&
                max((ars(1000, dunif, l = 0.1, u = 0.9))) <= 0.9)

  cat("\nBoundary tests passed\n")
})





test_that("Throw error message for input error",{
  expect_error(ars('1000', dnorm), "Please provide the numeric sample size N")
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

  cat('\nInput tests passed\n')
})


test_that("Throw error message when density is not log-concave", {
  expect_error(ars(1000, chi_test, l=1), 'Please provide the log-concave density')
  expect_error(ars(1000, chi_test, l=3, u=5), 'Please provide the log-concave density')

  expect_error(ars(100, t_test, l = -5, 5), 'Please provide the log-concave density')
  expect_error(ars(100, t_test2,l = -10, u = 10), 'Please provide the log-concave density')

  expect_error(ars(100, f_test, l=1), 'Please provide the log-concave density')
  expect_error(ars(1000, f_test2, l=3), 'Please provide the log-concave density')

  expect_error(ars(10, other_test, l=-5, u=5), 'Please provide the log-concave density')
  expect_error(ars(100, other_test2, l=-5, u=8), 'Please provide the log-concave density')

  cat('\nLog-concave tests passed\n')
})


test_that("Throw error message when function is not differentiable",{
  expect_error(ars(100, f_test, l = -1), "Target function is not differentiable in the given boundary. Please check for boundary")
  expect_error(ars(100, chi_test, l = -10), "Target function is not differentiable in the given boundary. Please check for boundary")
  expect_error(ars(100, beta_test, u = 2), "Target function is not differentiable in the given boundary. Please check for boundary")
  expect_error(ars(100, gamma_test, l = -5), "Target function is not differentiable in the given boundary. Please check for boundary")
  expect_error(ars(100, binom_test, l = -5), "Target function is not differentiable in the given boundary. Please check for boundary")

  cat('\nDifferentiability tests passed\n')
})


test_that("Throw error message when function is not continuous",{
  expect_error(ars(100, beta_test, l = -1, u = 2), "Target function DO NOT EXIST in the given boundary")
  expect_error(ars(100, chi_test, l = -5, u = 0), "Target function DO NOT EXIST in the given boundary")
  expect_error(ars(100, dunif, l = -1, u = 0), "Target function DO NOT EXIST in the given boundary")
  expect_error(ars(100, f_test, l = -10, u = -1), "Target function DO NOT EXIST in the given boundary")
  expect_error(ars(100, gamma_test, l = -10, u = 2), "Target function DO NOT EXIST in the given boundary")

  cat('\nExistence tests passed\n')
})





