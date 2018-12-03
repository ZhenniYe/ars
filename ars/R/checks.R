
# Check log-concave
check_logconcave <- function(fn, p){
  x <- sort(c(min_bound, max_bound, p)) # critical x-values
  hp <- calcDeriv(function(x) log(fn(x)), p)

  result <- all(diff(hp) < 0)
  return(result)
}




