#' Adapative Rejection Sampling
#'
#' `ars()` generates random numbers from log-concave functions by adapative rejection sampling.
#'
#' This function takes a numeric argument n and a function argument f
#' to generate random numbers from the given distribution in the specified sample size,
#' based on adapative rejection sampling method. For more detials, see Gilks and Wild (1992).
#'
#' @param n_iter numeric; number of points to sample 
#' @param fn - function: density function users interested in 
#' @param l - numeric: the optional lower bound of the density function users interested in 
#' @param u - numeric: the optional upper bound of the density function users interested in
#' @param mode - numeric: the estimated mode of of log-scale fn, optional 
#' @param u - number: value of bandwidth around mode used to initialize function, optional
#'
#' @return \code{ars} returns a vector of length n_iter containing sampled values.
#'
#' @export
#' @examples
#' ##  Sample from Standard Normal Distribution
#' ars(n_iter=10000, fn=dnorm, l = -Inf, u= Inf, mode=0, step=0.5)
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




ars <- function(n_iter, fn, l = -Inf, u = Inf, mode = 0, step = 0.5){ 
  
  ## check whether the input is correct (fn, DFUN and bounds)
  ## check numeric n_iter 
  if(is.numeric(n_iter) == FALSE){
    stop('input n_iter is not numeric', call. = FALSE)
  }
  ## check postive n_iter
  if (n_iter%%1 != 0 | n_iter < 1) stop("sample size is not positive integer")
  ## check whether the fn input is correct or not
  if(is.function(fn) == FALSE){
    stop('input fn is not a function', call. = FALSE)
  }
  ## check the bounds
  if(l >= u){
    stop('please provide resonalbe lower bound and upper bound', call. = FALSE)
  }
  
  
  ## take log scale of the density function
  FUN <- function(x,fun = fn){
    return(log(fun(x)))
  }  
  
  
  ## design the dervative function
  Deriv <- function(x, FUN, l, u){
    eps = 1e-8
    if (x==l) {
      return ((FUN(x + eps)-FUN(x))/eps)
    }
    if (x==u) {
      return ((FUN(x)-FUN(x - eps))/eps)
    }
    if (l <= x && x <= u) {
      return((FUN(x + eps)-FUN(x - eps))/2*eps)
    }
  }
  
  
  ## check  whether the point is defined on the function
  define_check <- function (point, FUN){
    if(!is.finite(FUN(point))) {
      stop('point is not defined on fn', call. = FALSE)}
  }
  
  
  ## create the starting points
  count  = 1
  ## case 1: user input finite lower and upper bound 
  
  if (l != -Inf && u != Inf){
    ## check whether the point is defined on the function
    define_check(l, FUN)
    define_check(u, FUN)
    
    ## test the derivative 
    test1 <- Deriv(l, FUN, l, u)
    test2 <- Deriv(u, FUN, l, u)
    if(test1 < 0 | test2 > 0){
      stop('the function is not log-convexity on the provided bound')
    } else{
      p <- c(l, u)}
  }
  ## case 2: user enter infinite lower bound and finite upper bound
  ## define the lower starting if lower bound is inf 
  if (l == -Inf && u != Inf){
    
    if (mode > u) { # u < 0
      mode <- u - step 
    }
    ll <- mode
    test <- Deriv(ll, FUN, l, u)
    ## push the samller starting abscissae left unitl find the first one that the diff is postive
    while (-Inf < test && test <= 0 && count <=50){  
      ll <- ll - step
      test <- Deriv(ll, FUN, l, u)
      count = count + 1
    }
    ## check define
    define_check(ll, FUN)
    define_check(u, FUN)
    
    p <- c(ll, u)
  }
  ## case 3: user input finite lower bound but infinit upper bound
  if (l != -Inf && u == Inf){
    
    if (mode < l) {  # l > 0
      mode = l + step
    }
    uu <- mode # uu= 0 if a<0, uu= a+0.5 if a >0
    test <- Deriv(uu, FUN, l, u)
    ## push the larger starting abscissae right unitl find the first one that the diff is negative
    while (0 <=  test && test < Inf && count <=50){  
      uu <- uu + step 
      test <- Deriv(uu, FUN, l, u) 
      count = count + 1
    }
    ## check define
    define_check(l, FUN)
    define_check(uu, FUN)
    p <- c(l, uu)
  }
  ## case 4: the default (-inf, inf)
  if (l == -Inf && u == Inf){
    
    ll <- mode - step 
    uu <- mode + step
    test1 <- Deriv(ll, FUN, l, u)
    test2 <- Deriv(uu, FUN, l, u)
    
    ## push the samller starting abscissae left unitl find the first one that the diff is postive
    while (-Inf < test1 && test1 <= 0 && count <= 50 ){
      ll <- ll - step
      test1 <- Deriv(ll)
      count = count + 1
    }
    ## push the larger starting abscissae right unitl find the first one that the diff is negative
    while (0 <=  test2 && test2 < Inf && count <= 50){
      uu <- uu + step
      test <- Deriv(uu)
      count = count + 1
    }
    define_check(ll, FUN)
    define_check(uu, FUN)
    p <- c(ll,uu)
  }
  
  
  ## main function
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
    } else { # if squeeze test failed, perform envelope test
      # if envelope test passed, accept and iterate
      if ((fn(x_val) / exp_pdf) > unif){
        set[i] <- x_val
        i <- i + 1
      } else { # if envelope test fails, reject and add candidate point to fixed points
        p <- sort(c(p, x_val))
        par <- setParams(fn, min_bound, max_bound, p)
      }
    }
  }
  return(set)
}


