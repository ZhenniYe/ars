# Functions needed for in-function tests

## function to take the derivative
Deriv <- function(x, FUN, l, u){
  eps = 1e-8
  if (x==l) {res <- ((FUN(x + eps)-FUN(x))/eps)}
  if (x==u) {res <- ((FUN(x)-FUN(x - eps))/eps)}
  if (l <= x && x <= u) {res <- ((FUN(x + eps)-FUN(x - eps))/(2*eps))}

  if (res == 'NaN' || !is.finite(res)){
    stop('Input function is not differentiable in the given boundary. Please check for boundary',
                                            call. = FALSE)}
  return(res)
}


## function to check  whether the point is defined on the function
define_check <- function (point, FUN){
  if(!is.finite(FUN(point))) {
    stop('Target function DO NOT EXIST in the given boundary', call. = FALSE)}
}


## function to check log-concave
Check_logconcave <- function(fn, p){
  x <- sort(p)
  hp <- calcDeriv(function(x) log(fn(x)), p)
  result <- all(diff(hp) < 0)
  return(result)
}


