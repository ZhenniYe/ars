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


  cat('\n`expPDF()` function tests passed`\n')
})


