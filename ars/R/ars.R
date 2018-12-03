#' Adapative Rejection Sampling
#'
#' `ars()` generates random numbers from log-concave functions by adapative rejection sampling.
#'
#' This function takes a numeric argument n and a function argument f
#' to generate random numbers from the given distribution in the specified sample size,
#' based on adapative rejection sampling method. For more detials, see Gilks and Wild (1992).
#'
#' @param fn function; log-concave function f(x).
#' @param min_bound numeric; the lower bound of the domain as the first starting point. The default if -Inf.
#' @param max_bound numeric; the upper bound of the domain as the second starting point. The default is Inf.
#' @param p numeric vector; starting fixed points for evaluation.
#' @param n_iter numeric; number of points to sample. The default is 1000.
#'
#' @return \code{ars} returns a vector of length n_iter containing sampled values.
#'
#' @export
#' @examples
#' ##  Sample from Standard Normal Distribution
#' ars(dnorm, min_bound = -5, max_bound = 5, p = c(-1, 1))
#'
#' ## Sample from Gamma Distribution
#' gamma_test <- function(x){
#' k <- dgamma(x, shape=7.5, scale=1)
#' return(k)
#' }
#' ars(gamma_test, min_bound = 0.1, max_bound = 20, p = c(3,10), n_iter = 10000)
#'
#' @author
#' Pan, Yanting; Ye, Zhenni; Myers, Vince; based on Gilks and Wild (1992).
#'
#' @references
#' W. R. Gilks, P. Wild. (1992), "Adaptive Rejection Sampling for Gibbs Sampling," Applied Statistics 41:337–348.
#'




ars <- function(fn, min_bound = -Inf, max_bound = Inf, p, n_iter = 1000){
  ### Checks ###
  # Check for function
  if (class(fn) != "function") stop('Please provide fn as an function', call. = FALSE)

  # Check for domain
  if (!is.numeric(min_bound) || !is.numeric(max_bound)) stop('Please provide the numeric domain', call. = FALSE)
  if (min_bound == max_bound) stop('Please provide the valid domain', call. = FALSE)
  if (min_bound > max_bound) {
    warning('The lower bound should be less the upper bound, swaping values', call. = FALSE, immediate. = TRUE)
    tmp <- min_bound
    min_bound <- max_bound
    max_bound <- tmp
  }

  # Check for fixed points
  if (!is.numeric(p)) stop('Please provide the numeric starting points', call. = FALSE)
  if (length(p) < 2) stop('Please provide at least two starting points for evaluation', call. = FALSE)
  if ( min_bound > min(p) || max_bound < max(p) ) stop('Please provide the fixed points within the domain',
                                                       call. = FALSE)

  # Check if the domain is unbounded (shink)
  if (min_bound == -Inf) {min_bound <- min(p) - 1000}
  if (max_bound == Inf) {max_bound <- max(p) + 1000}

  # Check if the input function is log-concave
  if (check_logconcave(fn, p) == FALSE) stop('Please provide the log-concave density', call. = FALSE)

  # Check for sample size
  if(!is.numeric(n_iter)) stop('Please provide the numeric sample size', call. = FALSE)
  if(n_iter%%1 != 0 || n_iter < 1) {
    warning('The sample size should be the positive integer, re-assigning as the default value',
            call. = FALSE, immediate. = TRUE)
    n_iter = 1000
  }





  ### Main ###
  p <- sort(p)
  i <- 1
  set <- rep(0, n_iter)
  par <- setParams(fn, min_bound, max_bound, p) # calculate parameters for initial set of fixed points

  # loop through iterations
  while (i <= n_iter){
    y <- runif(1)
    unif <- runif(1)
    x_val <- invCDF(y, par$int_x, par$m_p, par$b_vec, par$nc, par$shift, par$adj)
    exp_pdf <- expPDF(x_val, p, par$int_x, par$m_p, par$lf_p)
    squeeze <- lowerPDF(fn, x_val, p)
    # if squeeze test is passed, accept and iterate
    if ((squeeze / exp_pdf) > unif){
      set[i] <- x_val
      i <- i + 1
    } else { # if squeeze test failed, perform rejection test
      # if envelope test passed, accept and iterate
      if ((fn(x_val) / exp_pdf) > unif){
        set[i] <- x_val
        i <- i + 1
      }
      # whether rejection test fails or not, add candidate point to fixed points
      p <- sort(c(p, x_val))
      par <- setParams(fn, min_bound, max_bound, p)

    }
  }
  return(set)
}
