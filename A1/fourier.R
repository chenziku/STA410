#fast fourier transform

lambda = 7
epsilon = 10E-5
# PDF for X
p <- function(x) {
  p <- 0
  if (x >= 0 & x <=6){
    p <- (4-abs(3-x))/16
  }
  return(p)}
# PGF for X
phi <- function(s) {
  phi = 0
  for (x in 0:6) {
    phi <- phi + p(x) * s^x
  }
  return(phi)}
# function for M
M <- function(s) {
  (-log(epsilon) - lambda*(1-phi(s)))/log(s)
}
# fit a curve for M and look for the minimum, which turns out to be around 71.5
curve(M, from=1.3, to=1.5, xname="s", ylab="M", main = "Plot for Bound M")
M = round(optimize(M, interval=c(1, 2))$objective)

# Fast Fourier Transformation with Bound M
x=c(0:M-1)
px <- sapply(x, p)
pjhat <- fft(px)
gj <- exp(-lambda * (1-pjhat))
py <- Re(fft(gj, inv=T)/M)

# Plot the Distribution of Y
plot(x, py, main = "Probability Distribution of Y", xlab = "y", ylab = "P(Y = y)")