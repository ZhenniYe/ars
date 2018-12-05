### Functions needed for in-function tests ###

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

