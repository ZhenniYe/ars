
# Assign mode
Mode <- function(fn, l, u){
  intv <- seq(l, u, length.out = 1000)
  if (l<0 && 0<u) {intv <- c(intv, 0)}
  results <- fn(intv)
  mode <- intv[results == max(results)]
  if (length(mode) > 1) {
    diff <- abs(mode)
    mode <- mode[diff == min(diff)]
  }
  return(mode)
}


# Check log-concave
Check_logconcave <- function(fn, p){
  x <- sort(p)
  hp <- calcDeriv(function(x) log(fn(x)), p)

  result <- all(diff(hp) < 0)
  return(result)
}






