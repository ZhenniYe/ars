## the whole function used to test


# function to calculate slope of tangent lines
calcDeriv <- function(fn, p){
  dif <- 1e-10
  
  p1 <- fn(p)
  p2 <- fn(p + dif)
  
  slope <- (p2 - p1) / dif
  return(slope)
}

# function to calculate lower squeeze function
lowerPDF <- function(fn, x, p){
  
  # squeeze function only applies if x is between the min and max fixed points
  if (x <= min(p) | x >= max(p)){
    return(0)
    
  } else {
    
    w <- findInterval(x, c(-Inf, p, Inf)) # find which section of squeeze line applies to x
    y <- ((p[w] - x)*log(fn(p[w-1])) + (x - p[w-1])*log(fn(p[w])))/(p[w] - p[w-1]) # value of squeeze at x
    return(exp(y))
  }
}

# function to compute piecewise exponential PDF
expPDF <- function(x, p, int_x, m_p, lf_p){
  
  b_vec <- lf_p - p*m_p # calculate b in y=mx+b format for each tangent line
  w <- findInterval(x, c(-Inf, int_x, Inf)) # find section of enveloping function that applies to x
  y <- m_p[w]*x + b_vec[w] # calculate y=log(env(x))
  return(exp(y))
  
}

# function to compute piecewise exponential CDF
expCDF <- function(x, p, int_x, m_p, lf_p){
  exp_cdf <- rep(0, length(x))
  b_vec <- lf_p - p*m_p # calculate b in y=mx+b format for each tangent line
  which_section <- findInterval(x, c(-Inf, int_x, Inf)) # break x vector into sections denoting related enveloping line
  nsections <- max(which_section)
  
  # calculate starting point of each segment of piecewise CDF
  piece_difs <- rep(0, length(int_x)+1)
  for (j in 1:length(int_x)){
    left <- (1/m_p[j])*exp(m_p[j]*int_x[j] + b_vec[j])
    right <- (1/m_p[j+1])*exp(m_p[j+1]*int_x[j] + b_vec[j+1])
    piece_difs[j+1] <- right - left
  }
  piece_difs_cumul <- cumsum(piece_difs) # vertical adjustment to pieces of CDF to create continuity
  
  for (k in 1:nsections){ # loop through pieces of piecewise function, computing F(x) for each x
    x_section <- x[which_section == k]
    y_section <- (1/m_p[k])*exp(m_p[k]*x_section + b_vec[k])
    exp_cdf[which_section == k] <- y_section - piece_difs_cumul[k]
  }
  
  adj <- min(exp_cdf) # set minimum point of CDF to zero
  nc <- max(exp_cdf - adj) # normalizing constant
  norm_exp_cdf <- (exp_cdf - adj) / nc # normalized piecewise CDF
  return(list(cdf=norm_exp_cdf, difs=piece_difs_cumul, norm.const = nc, adj = adj))
}

# function to take invert the CDF at a given point for F(x)
invCDF <- function(y, int_x, m_p, b_vec, nc, shift, adj){
  n <- length(m_p) - 1
  breaks <- ((1/m_p[1:n])*exp(m_p[1:n]*int_x + b_vec[1:n]) - shift[1:n])/nc
  loc <- findInterval(y, c(-Inf, breaks, Inf))
  x_output <- (log(m_p[loc]*(nc*y + shift[loc] + adj)) - b_vec[loc])/m_p[loc]
  return(x_output)
}

# calculate the parameters for the enveloping function for a given set of fixed points p
setParams <- function(fn, min_bound, max_bound, p){
  
  x <- sort(c(min_bound, max_bound, p)) # critical x-values
  m_p <- calcDeriv(function(x) log(fn(x)), p) # calculate slope of tangent lines
  lf_p <- log(fn(p)) # calculate value of log function at eval point
  
  # compute intersection points
  n <- length(p)
  int_x <- (lf_p[2:n] - lf_p[1:n-1] - p[2:n]*m_p[2:n] + p[1:n-1]*m_p[1:n-1]) / (m_p[1:n-1] - m_p[2:n])
  
  # calculate parameters for inverse CDF
  exp_cdf <- expCDF(x, p, int_x, m_p, lf_p) # y-output for piecewise linear PDF
  nc <- exp_cdf$norm.const
  shift <- exp_cdf$difs
  adj <- exp_cdf$adj
  b_vec <- lf_p - p*m_p
  
  # return the parameters needed by the invCDF function
  return(list(m_p = m_p, lf_p = lf_p, int_x = int_x, b_vec = b_vec, nc=nc, shift=shift, adj = adj))
}


# Design the dervative function
Deriv <- function(x, FUN, l, u){
  eps = 1e-8
  if (x==l) {return ((FUN(x + eps)-FUN(x))/eps)}
  if (x==u) {return ((FUN(x)-FUN(x - eps))/eps)}
  if (l <= x && x <= u) {return((FUN(x + eps)-FUN(x - eps))/2*eps)}
}


# Check  whether the point is defined on the function
define_check <- function (point, FUN){
  if(!is.finite(FUN(point))) {
    stop('point is not defined on fn', call. = FALSE)}
}


# Check log-concave
Check_logconcave <- function(fn, p){
  x <- sort(p)
  hp <- calcDeriv(function(x) log(fn(x)), p)
  
  result <- all(diff(hp) <= 0)
  return(result)
}



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




## test 1: standard normal with -inf, inf
set.seed(0)
set1 <- ars (N=10000, dnorm, l = -Inf, u = Inf)
x1 <- seq(-4, 4, length=100)
z1 <- dnorm(x)
plot(density(set1))
lines(x1, z1, type="l", col="blue")

## test2 : uniform distribution (0,1)
set.seed(0)
set2 <- ars (N=10000, dunif, l = 0, u = 1)
plot(density(set2))


## test3: gamma with 0.1, inf
set.sed(0)
gamma_test <- function(x){
  k <- dgamma(x, shape=7.5, scale=1)
  return(k)
}
z3 <- seq(0.1, 20, by=0.05)
set3 <- ars(n_iter = 10000, fn= gamma_test, l= 0.1, u=Inf)
plot(density(set3))
lines(z3, gamma_test(z3), type="l", col="blue")


## test 4: truncate t (-3, -2)
set.seed(0)
set4 <- ars (n_iter=10000, dnorm, l = -3, u = -2)
plot(density(set4))
