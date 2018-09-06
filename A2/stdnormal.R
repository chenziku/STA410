#Question 1 Rejection Sampling
stdnormal <- function(n) {
  u <- NULL
  2
  v <- NULL
  rejections <- 0
  uplus <- 1
  vminus <- -sqrt(2/exp(1))
  vplus <- sqrt(2/exp(1))
  for (i in 1:n) {
    reject <- T
    while(reject) {
      ustar <- runif(1, 0, uplus)
      vstar <- runif(1, vminus, vplus)
      #check if ustar and vstar are in Ch
      if (ustar <= exp((-vstar^2)/(4*ustar^2))) {
        u <- c(u, ustar)
        v <- c(v, vstar)
        reject <- F
      }
      else rejections <- rejections + 1
    }
  }
  x <- v/u
  # calculate proposal acceptance rate
  accept.rate <- n/(n+rejections)
  r <- list(x=x,accept.rate=accept.rate)
  r
}
y <- stdnormal(1000)
#Acceptence Rate
y$accept.rate
#Normal Q-Q Plot
qqnorm(y$x)
qqline(y$x)
#Histogram
hist(y$x, main="Histogram of X", xlab = "X", xlim = c(-4,4))
