context('Auxiliary Functions Testing')
set.seed(2018)


## functions used for testing
square_test <- function(x){return(x^2)}


test_that("Test `calcDeriv()` function",{
  expect_equal(calcDeriv(dnorm, 0), 0)
  expect_equal(round(calcDeriv(square_test,1)), 2)

  cat('\n`calcDeriv()` function tests passed`\n')
})


test_that("Test `lowerPDF()` function",{
  expect_equal(lowerPDF(dnorm, 100, c(-10, 10)), 0)
  expect_equal(lowerPDF(dnorm, -11, c(-10, 10)), 0)
  expect_is(lowerPDF(dnorm, 2, c(-10, 10)), "numeric")
  expect_true(lowerPDF(dnorm, 2, c(-10, 10)) != 0)
  expect_is(lowerPDF(dexp, 5, c(0, 10)), "numeric")
  expect_true(lowerPDF(dexp, 5, c(0, 10)) != 0)

  cat('\n`lowerPDF()` function tests passed`\n')
})


test_that("Test `expPDF()` function",{
  expect_is(expPDF(1, c(-1, 0, 2, 5), c(-0.5, 1.25, 3), c(3, 2, -2, -1), c(5, 10, 8, 6)), "numeric")
  expect_is(expPDF(c(0, 1, 2, 3, 4), c(-1, 0, 2, 5), c(-0.5, 1.25, 3), c(3, 2, -2, -1), c(5, 10, 8, 6)), "numeric")
  expect_true(is.vector(expPDF(c(0, 1, 2, 3, 4), c(-1, 0, 2, 5), c(-0.5, 1.25, 3), c(3, 2, -2, -1), c(5, 10, 8, 6))))
  expect_is(expPDF(10, c(-10, 0, 5, 30), c(-5, -1, 6, 22), c(5, 3, 0, -9), c(10, 30, 50, 9)), "numeric")

  cat('\n`expPDF()` function tests passed`\n')
})


test_that("Test `expCDF()` function",{
  expect_is(expCDF(c(0,1,2,3,4), c(-1, 0, 2, 5), c(-0.5, 1.25, 3), c(3, 2, -2, -1), c(5, 10, 8, 6)), "list")
  expect_is(expCDF(c(0,1,2,3,4), c(-1, 0, 2, 5), c(-0.5, 1.25, 3), c(3, 2, -2, -1), c(5, 10, 8, 6))$cdf, "numeric")
  expect_true(is.vector(expCDF(c(0,1,2,3,4), c(-1, 0, 2, 5), c(-0.5, 1.25, 3), c(3, 2, -2, -1), c(5, 10, 8, 6))$cdf))
  expect_equal(expCDF(c(0,1,2,3,4), c(-1, 0, 2, 5), c(-0.5, 1.25, 3), c(3, 2, -2, -1), c(5, 10, 8, 6))$cdf[1], 0)
  expect_equal(tail(expCDF(c(0,1,2,3,4), c(-1, 0, 2, 5), c(-0.5, 1.25, 3), c(3, 2, -2, -1), c(5, 10, 8, 6))$cdf, n=1), 1)
  expect_true(is.vector(expCDF(c(0,1,2,3,4), c(-1, 0, 2, 5), c(-0.5, 1.25, 3), c(3, 2, -2, -1), c(5, 10, 8, 6))$difs))
  expect_is(expCDF(c(0,1,2,3,4), c(-1, 0, 2, 5), c(-0.5, 1.25, 3), c(3, 2, -2, -1), c(5, 10, 8, 6))$difs, "numeric")
  expect_is(expCDF(c(0,1,2,3,4), c(-1, 0, 2, 5), c(-0.5, 1.25, 3), c(3, 2, -2, -1), c(5, 10, 8, 6))$norm.const, "numeric")
  expect_is(expCDF(c(0,1,2,3,4), c(-1, 0, 2, 5), c(-0.5, 1.25, 3), c(3, 2, -2, -1), c(5, 10, 8, 6))$adj, "numeric")

  expect_equal(expCDF(c(-2,-1), c(-2,0), -1, c(2,0), c(-2.9, -0.9))$cdf[1], 0)
  expect_equal(expCDF(c(-2,-1), c(-2,0), -1, c(2,0), c(-2.9, -0.9))$cdf[2], 1)
  expect_equal(round(expCDF(c(-2,-1), c(-2,0), -1, c(2,0), c(-2.9, -0.9))$norm.const, 3), 0.176)
  expect_true(expCDF(1, c(0.09,0.99), 1, c(-1,-100), c(0.6, -4))$cdf == "NaN")

  cat('\n`expCDF()` function tests passed`\n')
})


test_that("Test `invCDF()` function",{
  expect_is(invCDF(1, c(1,2,3), c(0, -2), c(-0.919, 1.08),  0.571, c(0, -0.600), 0), "numeric")
  expect_is(invCDF(0.2, 1, c(0, -2), c(-0.919, 1.08),  0.571, c(0, -0.600), 0), "numeric")
  expect_equal(round(invCDF(1, 1, c(0, -2), c(-0.919, 1.08),  0.571, c(0, -0.600), 0), 3), 1.964)
  expect_equal(round(invCDF(0.2, 1, c(0, -2), c(-0.919, 1.08),  0.571, c(0, -0.600), 0), 3), 0.286)
  expect_equal(round(invCDF(1, -1, c(2, 0), c(-2.919, -0.919), 0.571, c(0, -0.598), 0.027)), 0)
  expect_equal(round(invCDF(-0.14, 0, c(1, -1), c(-0.419, -0.419),  0.832, c(0, -1.315), 0.242), 3), -1.656)

  cat('\n`invCDF()` function tests passed`\n')
})


test_that("Test `setParams()` function",{
  expect_is(setParams(dnorm, -10, 10, c(-5, -1, 0, 1, 5)), "list")
  expect_is(setParams(dnorm, -10, 10, c(-5, -1, 0, 1, 5))$m_p, "numeric")
  expect_true(is.vector(setParams(dnorm, -10, 10, c(-5, -1, 0, 1, 5))$m_p))
  expect_is(setParams(dnorm, -10, 10, c(-5, -1, 0, 1, 5))$lf_p, "numeric")
  expect_true(is.vector(setParams(dnorm, -10, 10, c(-5, -1, 0, 1, 5))$lf_p))
  expect_is(setParams(dnorm, -10, 10, c(-5, -1, 0, 1, 5))$int_x, "numeric")
  expect_true(is.vector(setParams(dnorm, -10, 10, c(-5, -1, 0, 1, 5))$int_x))
  expect_is(setParams(dnorm, -10, 10, c(-5, -1, 0, 1, 5))$b_vec, "numeric")
  expect_true(is.vector(setParams(dnorm, -10, 10, c(-5, -1, 0, 1, 5))$b_vec))
  expect_is(setParams(dnorm, -10, 10, c(-5, -1, 0, 1, 5))$shift, "numeric")
  expect_true(is.vector(setParams(dnorm, -10, 10, c(-5, -1, 0, 1, 5))$shift))
  expect_is(setParams(dnorm, -10, 10, c(-5, -1, 0, 1, 5))$nc, "numeric")
  expect_is(setParams(dnorm, -10, 10, c(-5, -1, 0, 1, 5))$adj, "numeric")
  expect_is(setParams(dnorm, -2, 0, c(-2,0))$m_p[1], "numeric")
  
  expect_equal(round(setParams(dnorm, -2, 0, c(-2,0))$m_p[1]), 2)
  expect_equal(setParams(dnorm, -2, 0, c(-2,0))$m_p[2], 0)
  expect_equal(round(setParams(dnorm, -2, 0, c(-2,0))$int_x), -1)
  expect_equal(setParams(dnorm, -2, 0, c(-2,0))$shift[1], 0)

  cat('\n`setParams()` function tests passed`\n')
})
