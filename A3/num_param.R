num.param <- function(x, span, m) {
  n <- length(x)
  that <- NULL
  for (i in 1:m) {
    # generate rv V to minimize variance of the trace estimate
    v <- NULL
    for (j in 1:n){
      u <- runif(1, 0,1)
      if (u < 0.5) {
        v <- c(v, 1)
      }
      else {
        v <- c(v, -1)
      }
    }
    # estimate vhat generated from loess function
    r <- loess(v~x, span=span)
    vhat <- r$fitted
    that <- c(that, sum(v*vhat))
  }
  # compute the average (expectation) of v^T*S*v
  numparam = mean(that)
  # estimate the standard error
  sderror=sd(that)/sqrt(n)
  result <- list(numparam = numparam, sderror = sderror)
}

x <- rnorm(10000, 0, 1)
r <- num.param(x, span = 0.7, 100)
r
