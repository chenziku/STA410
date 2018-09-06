# Question 2 G-S Method
seidel <- function(y,lambda,theta,max.iter=50,eps=1.e-6) {
  n <- length(y)
  # Define initial estimates if unspecified
  if (missing(theta)) theta <- rep(mean(y),n)
  # Compute objective function for initial estimates
  obj <- sum((y-theta)^2)+lambda*sum(diff(theta)^2)
  # The function diff(theta) computes first differences of the vector theta
  no.conv <- T
  # Do Gauss-Seidel iteration until convergence or until max.iter iterations
  iter <- 0
  theta.old <- theta
  while(no.conv) {
    theta[1] <- (y[1] + lambda*theta[2])/(1+lambda)
    # Update theta[2],..., theta[n-1]
    for (j in 2:(n-1)) {
      theta[j] <- (y[j]+lambda*(theta[j+1]+theta[j-1]))/(1+2*lambda)
    }
    theta[n] <- (y[n]+lambda*theta[n-1])/(1+lambda)
    # Compute new objective function for current estimates
    obj.new <- sum((y-theta)^2)+lambda*sum(diff(theta)^2)
    iter <- iter + 1
    # Now set no.conv to F if either convergence or iter=max.iter
    # and update the value of the objective function variable obj
    if (obj==obj.new) no.conv <- F
    if (iter==max.iter) no.conv <- F
    obj <- obj.new
    theta.old <- theta
  }
  r <- list(y=y,theta=theta,obj=obj,niters=iter)
  r
}

y <- c(rep(0,250),rep(1,250),rep(0,50),rep(1,450)) + rnorm(1000,0,0.1)

#plot y and its estimate theta with different values of lambda
lambda <- 1
plot(c(1:1000), y, col="blue")
points(c(1:1000), seidel(y,lambda)$theta, col="red")

lambda <- 10
plot(c(1:1000), y, col="blue")
points(c(1:1000), seidel(y,lambda)$theta, col="red")

lambda <- 100
plot(c(1:1000), y, col="blue")
points(c(1:1000), seidel(y,lambda)$theta, col="red")


