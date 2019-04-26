normalmixture <- function(x,k,mu,sigma,lambda,eps=1e-6,max.iter=5000) {
  n <- length(x)
  x <- sort(x)
  vars <- sigma^2
  means <- mu
  lam <- lambda/sum(lambda)  # guarantee that lambdas sum to 1
  delta <- matrix(rep(0,n*k),ncol=k) 
  # initial deltas
  for (i in 1:n) {
    xi <- x[i]
    for (j in 1:k) {
      mj <- means[j]
      varj <- vars[j]
      denom <- 0
      for (u in 1:k) {
        mu <- means[u]
        varu <- vars[u]
        denom <- denom + lam[u]*dnorm(xi,mu,sqrt(varu))
      }
      delta[i,j] <- lam[j]*dnorm(xi,mj,sqrt(varj))/denom
    }
  }
  # initial log likelihood value
  new.loglik <- 0
  for (i in 1:n) {
    xi <- x[i]
    for (j in 1:k) {
      mj <- means[j]
      varj <- vars[j]
      fkxi <- dnorm(xi,mj,sqrt(varj))
      new.loglik <- new.loglik + log(fkxi)*delta[i,j] + delta[i,j]*log(lam[j])
    }
  }
  
  iter <- 1
  no.conv <- T
  while (no.conv && iter <= max.iter) {
    loglik <- new.loglik
    # compute updates of deltas 
    for (i in 1:n) {
      xi <- x[i]
      for (j in 1:k) {
        mj <- means[j]
        varj <- vars[j]
        denom <- 0
        for (u in 1:k) {
          mu <- means[u]
          varu <- vars[u]
          denom <- denom + lam[u]*dnorm(xi,mu,sqrt(varu))
        }
        delta[i,j] <- lam[j]*dnorm(xi,mj,sqrt(varj))/denom
      }
    }
    # compute updated estimates of means, variances, and probabilities - the 
    # function weighted.mean may be useful here for computing the estimates of
    # the means and variances.
    for (j in 1:k) {
      deltaj <- as.vector(delta[,j])
      sum_dj <- sum(deltaj)
      lambda[j] <- sum_dj/n
      means[j] <- sum(x*deltaj)/sum_dj
      vars[j] <- sum((x-means[j])^2*deltaj)/sum_dj
    }
    lam <- lambda/sum(lambda) 
    iter <- iter + 1
    # Log-likelihood computation
    new.loglik <- 0
    for (i in 1:n) {
      xi <- x[i]
      for (j in 1:k) {
        mj <- means[j]
        varj <- vars[j]
        fkxi <- dnorm(xi,mj,sqrt(varj))
        # rule out the cases when the probability is too small 
        # to evaluate (negligible to the sum)
        if (fkxi != 0) {
          new.loglik <- new.loglik + log(fkxi)*delta[i,j] + delta[i,j]*log(lam[j])
        }
      }
    }
    print(new.loglik)
    if (abs(new.loglik-loglik)<eps) no.conv <- F
  }
  r <- list(mu=means,var=vars,delta=delta,lambda=lam,loglik=loglik,iter=iter)
  r
}

# getting data from txt file
setwd("/Users/zikunchen/Desktop/GitHub/STA410/A4")
stamp <- read.table("stamp.txt", fill = TRUE)
stamp <- unname(unlist(stamp))
stamp <- stamp[!is.na(stamp)]

# 5 modes
k <- 5
# initial estimate based on graph
plot(density(stamp,bw=0.0026))
lambda <- rep(1/k, k)
mu <- c(0.079, 0.09, 0.1, 0.11, 0.12)
sigma <- c(0.1414, 0.1, 0.1, 0.1, 0.12)
r5 <- normalmixture(stamp,k,mu,sigma,lambda)

# 6 modes
k <- 6
# initial estimate based on graph
plot(density(stamp,bw=0.0024))
lambda <- rep(1/k, k)
mu <- c(0.079, 0.09, 0.1, 0.11, 0.12, 0.13)
sigma <- c(0.1414, 0.1, 0.1, 0.1, 0.1, 0.1)
r6 <- normalmixture(stamp,k,mu,sigma,lambda)

# 7 modes
k <- 7
# initial estimate based on graph
plot(density(stamp,bw=0.0015))
lambda <- rep(1/k, k)
mu <- c(0.071, 0.08, 0.09, 0.1, 0.11, 0.12, 0.124)
sigma <- c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1)
r7 <- normalmixture(stamp,k,mu,sigma,lambda)

# compare models using log-likelihood test


