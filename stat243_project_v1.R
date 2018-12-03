### The two main functions are ars() and setParams(). The calcDeriv() and expCDF() functions are used
### by setParams(); the setParams(), lowerPDF(), expPDF(), and invCDF() functions are used by ars().
### setParams() calculates the parameters of the enveloping function for a set of fixed points; those
### parameters are used in ars() for all of the other function calls.

### KEY:
### p = fixed points; intially provided by user, and updated after every rejection
### min_bound / max_bound = the boundaries on the function; all fixed points should be within min_bound and max_bound
### m_p = slope of lines in log of enveloping function at each fixed point p (uses y = mx + b form)
### b_vec = b in (y = mx + b)
### lf_p = value of log of user-provided function at fixed points p
### int_x = intersection point for lines in log of enveloping function
### shift = vertical shift for each segment of enveloping function to create a continuous CDF
### adj = vertical shift for the entire enveloping CDF so that F(x) = 0 at min_bound
### nc = normalizing constant so that max point of enveloping CDF is 1


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

# main ARS sampling function
# Requires user to input min, max, and fixed points p
# I would suggest having tests to make sure the user provides at least 2 fixed points, and 
# making sure all fixed points are between the min and the max
ars <- function(fn, min_bound, max_bound, p, n_iter = 1000){ 
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


#### test 1: standard normal distribution

min <- -5
max <- 5
p <- c(-3, 2)
set1 <- ars(dnorm, min, max, p, n_iter = 10000)
z <- seq(min, max, by=0.05)
hist(set1, prob=TRUE)
lines(z, dnorm(z), type="l", col="blue")


### test 2: truncated standard normal distribution

min <- -1
max <- 2
p <- c(-0.8, 1.3)
set2 <- ars(dnorm, min, max, p, n_iter = 10000)
z <- seq(min, max, by=0.05)
hist(set2, prob=TRUE)
lines(z, dnorm(z)/(pnorm(max) - pnorm(min)), type="l", col="blue")


### test 3: gamma distribution

gamma_test <- function(x){
  k <- dgamma(x, shape=7.5, scale=1)
  return(k)
}

min <- 0.1
max <- 20
p <- c(3, 10)
set3 <- ars(gamma_test, min, max, p, n_iter = 10000)
z <- seq(min, max, by=0.05)
hist(set3, prob=TRUE)
lines(z, gamma_test(z), type="l", col="blue")



testfn <- function(x){
  k <- dgamma(x, shape=7.5, scale=1)
  return(k)
}
