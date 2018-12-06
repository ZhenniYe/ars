#' Adapative Rejection Sampling
#'
#' `ars()` generates random numbers from log-concave functions by adapative rejection sampling.
#'
#' This function takes a numeric argument n and a function argument f
#' to generate random numbers from the given distribution in the specified sample size,
#' based on adapative rejection sampling method. For more detials, see Gilks and Wild (1992).
#'
#' @param n_iter numeric: number of points to sample.
#' @param fn - function: log-concave density function to draw sample from.
#' @param l - numeric: the lower bound of the density function. The default is -Inf.
#' @param u - numeric: the upper bound of the density function. The default is Inf.
#' @param center - numeric: the estimated center of of log-scale fn. The default is 0.
#' @param step - numeric: value of bandwidth around center used to initialize function. The default is 0.5.
#'
#' @return \code{ars} returns a vector of length n_iter containing sampled values.
#'
#' @export
#' @examples
#' ##  Sample from Standard Normal Distribution
#' ars(n_iter=10000, fn=dnorm, l = -Inf, u= Inf, center=0, step=0.5)
#' ars(10000, dnorm) # equivalent call
#'
#' ## Sample from Gamma Distribution
#' gamma_test <- function(x){
#' k <- dgamma(x, shape=7.5, scale=1)
#' return(k)
#' }
#' ars(n_iter = 10000, fn= gamma_test, l= 0.1, u=20)
#'
#' @author
#' Pan, Yanting; Ye, Zhenni; Myers, Vince; based on Gilks and Wild (1992).
#'
#' @references
#' W. R. Gilks, P. Wild. (1992), "Adaptive Rejection Sampling for Gibbs Sampling," Applied Statistics 41:337–348.
#'





ars <- function(n_iter, fn, l = -Inf, u = Inf, center = 0, step = 0.5){
  ### Initiatlization ###

  ## check for n_iter
  if(is.numeric(n_iter) == FALSE) {stop('Please provide the numeric sample size n_iter', call. = FALSE)}
  if (n_iter%%1 != 0 || n_iter < 1) {stop("Please provide positive integer for sample size", call. = FALSE)}
  ## check for fn
  if(is.function(fn) == FALSE) {stop('Please provide fn as an function', call. = FALSE)}
  ## check the bounds
  if (!is.numeric(l) || !is.numeric(u)) {stop('Please provide the numeric boundary', call. = FALSE)}
  if (l == u) {stop('Please provide the valid boundary', call. = FALSE)}
  if (l > u) {
    warning('The lower bound should be less the upper bound, swaping values', call. = FALSE, immediate. = TRUE)
    tmp <- l
    l <- u
    u <- tmp
  }
  

  # Take log scale of the density function
  FUN <- function(x,fun=fn){
    return(log(fun(x)))
  }

  # create the starting points
  count  = 1
  ## case 1: user finite input as starting points
  if (l != -Inf && u != Inf){
    ### check whether the point is defined on the function
    define_check(l, FUN)
    define_check(u, FUN)
    inif <- c(l, u)

    ## for uniform case:
    ### test the derivative
    test1 <- Deriv(l, FUN, l, u)
    test2 <- Deriv(u, FUN, l, u)
    if ((abs(test1) < 1e-6)&&(abs(test1)< 1e-6)){
      return (set = runif(n = n_iter, l, u))
      }
    }

  ## case 2: user enter infinite lower bound and finite upper bound
  ## define the lower starting if lower bound is inf
  if (l == -Inf && u != Inf){

    if (center > u) { ### u < 0
      center <- u - step  ### find the lower bound begining with u-0.5
    }
    ll <- center
    test <- Deriv(ll, FUN, l, u)
    ### push the samller starting abscissae left unitl find the first one that the diff is postive
    while (-Inf < test && test <= 0 && count <=100){
      ll <- ll - step
      test <- Deriv(ll, FUN, l, u)
      count = count + 1
    }
    ### check define
    define_check(ll, FUN)
    define_check(u, FUN)

    inif <- c(ll, u)
    }

  ## case 3: user input finite lower bound but infinit upper bound
  ## find the upper starting
  if (l != -Inf && u == Inf){

    if (center < l) {  ### l > 0
      center = l + step ### find the upper bound begining with l+0.5
    }
    uu <- center ### uu= 0 if l<0, uu= l+0.5 if l >0
    test <- Deriv(uu, FUN, l, u)
    ### push the larger starting abscissae right unitl find the first one that the diff is negative
    while (0 <=  test && test < Inf && count <=100){
      uu <- uu + step
      test <- Deriv(uu, FUN, l, u)
      count = count + 1
    }
    ### check define
    define_check(l, FUN)
    define_check(uu, FUN)
    inif <- c(l, uu)
    }

  ## case 4: the default (-inf, inf)
  if (l == -Inf && u == Inf){

    ll <- center - step
    uu <- center + step
    test1 <- Deriv(ll, FUN, l, u)
    test2 <- Deriv(uu, FUN, l, u)

    ### push the samller starting abscissae left unitl find the first one that the diff is postive
    while (-Inf < test1 && test1 <= 0 && count <= 100 ){
      ll <- ll - step
      test1 <- Deriv(ll)
      count = count + 1
    }
    ### push the larger starting abscissae right unitl find the first one that the diff is negative
    while (0 <=  test2 && test2 < Inf && count <= 100){
      uu <- uu + step
      test <- Deriv(uu)
      count = count + 1
    }
    define_check(ll, FUN)
    define_check(uu, FUN)
    inif <- c(ll,uu)
    }

  if (count >= 100) {stop ("Initial points cannot be found under given conditions.
                           Please check the boundary or(and) center of the log-concave density.",
                           .call = FALSE)}


  ### Main ###
  p <- sort(inif)
  # check for log-concave
  if (Check_logconcave(fn, p) == FALSE) stop('Please provide the log-concave density', call. = FALSE)

  i <- 1
  set <- rep(0, n_iter)
  min_bound <- ifelse(l== -Inf, -1e8, l)
  max_bound <- ifelse(u== Inf, 1e8, u)
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
        # if rejection test passed, accept and iterate
        if ((fn(x_val) / exp_pdf) > unif){
          set[i] <- x_val
          i <- i + 1
          }
        # whether the rejection test fails or not, add candidate point to fixed points
        p <- sort(c(p, x_val))
        # check for log-concave
        if (Check_logconcave(fn, p) == FALSE) stop('Please provide the log-concave density', call. = FALSE)
        par <- setParams(fn, min_bound, max_bound, p)
        }
    }

  return(set)
  }

