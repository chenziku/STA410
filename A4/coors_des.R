coor.des <- function(y,lambda,theta1,phis,eps=1.e-6) {
  n <- length(y)
  lambda2 <- lambda/2
  phis <- c(NA, phis)
  # initial objective funciton value
  term2 <- 0
  for (i in 2:n) {
    term2 <- term2 + (y[i]-theta1-sum(phis[2:i]))^2
  }
  new.obj <- (y[1]-theta1)^2 + term2 + lambda*sum(abs(phis[2:n]))
  no.conv <- T
  iter <- 0
  while (no.conv) {
    obj <- new.obj
    # update theta1
    theta1 <- y[1]
    for (i in 2:n) {
      theta1 <- theta1 + y[i]-sum(phis[2:i])
    }
    theta1 <- theta1/n
    # update phi's
    term2 <- 0
    for (j in 2:n) {
      sumj <- 0
      for (i in j:n) {
        sumj <- sumj + y[i]-theta1-sum(phis[2:i][-j])
      }
      if (abs(sumj)<=lambda2) {
        phis[j] <- 0
      }
      else if (abs(sumj)>lambda2){
        phis[j] <- (sumj-lambda2)/(n-j+1)
      }
      else {
        phis[j] <- (sumj+lambda2)/(n-j+1)
      }
      term2 <- term2 + (y[j]-theta1-sum(phis[2:j]))^2
    }
    #print(phis)
    iter <- iter + 1
    # compute the new objective function value
    new.obj <- (y[1]-theta1)^2 + term2 + lambda*sum(abs(phis[2:n]))
    
    print(new.obj)
    
    if (abs(obj-new.obj)<eps) no.conv <- F
  }
  r <- list(theta1=theta1, phis=phis, iter=iter)
  r
}

lambda <- 0
y <- c(rep(0,250),rep(1,250),rep(0,50),rep(1,450)) + rnorm(1000,0,0.1)
# using the result of the seidal function from A2 for initial estimates
thetas <- seidel(y, lambda)$theta
theta1 <-thetas[1]
n <- length(thetas)
phis <- thetas[2:n] - thetas[1:n-1]

coor.des(y,lambda,theta1,phis)

