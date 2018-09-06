gamma.nr <- function(x,alpha,lambda,eps=1.e-8,max.iter=100) {
  n <- length(x)
  m1 <- mean(x)
  var <- sum((x-m1)^2)/n
  
  # use MoM to estimate alpha and lambda if missing
  if (missing(alpha)) {
    lambda <- m1/var
    alpha <- lambda*m1
  }
  
  lambda1 <- log(lambda)
  alpha1 <- log(alpha)
  theta1 <- c(alpha1,lambda1)
  
  # compute the scores based on the initial estimates
  sum_x <- m1*n
  sum_lnx <- sum(log(x))
  score1 <- n*lambda1 + exp(alpha1)*sum_lnx - n*digamma(exp(alpha1))*exp(alpha1)/gamma(exp(alpha1))
  score2 <- n*alpha1 + exp(lambda1)*sum_x
  score <- c(score1,score2)
  iter <- 1
  while (max(abs(score))>eps && iter<=max.iter) {
    # compute observed Fisher information
    info.11 <- n*((trigamma(exp(alpha1))*exp(2*alpha1)+exp(alpha1)*digamma(exp(alpha1)))-(digamma(exp(alpha1))^2*exp(2*alpha1)))/gamma(exp(alpha1))^2 
      - exp(alpha1)*sum_lnx
    info.12 <- -n
    info.22 <- -exp(lambda1)*sum_x
    info <- matrix(c(info.11,info.12,info.12,info.22),ncol=2)
    # Newton-Raphson iteration
    theta1 <- theta1 + solve(info,score)
    
    alpha1 <- theta1[1]
    lambda1 <- theta1[2]
    iter <- iter + 1
    # update score
    score1 <- n*lambda1 + exp(alpha1)*sum_lnx - n*digamma(exp(alpha1))*exp(alpha1)/gamma(exp(alpha1))
    score2 <- n*alpha1 + exp(lambda1)*sum_x
    score <- c(score1,score2)
  }
  if (max(abs(score))>eps) print("No convergence")
  else {
    alpha <- exp(alpha1)
    lambda <- exp(lambda1)
    print(paste("Number of iterations =",iter-1))
    loglik <- n*alpha*log(lambda) - n*log(gamma(alpha)) + (alpha-1)*sum_lnx - lambda*sum_x
    info.11 <- n*(trigamma(alpha)*gamma(alpha) - digamma(alpha)^2)/gamma(alpha)^2
    info.12 <- -n/lambda
    info.22 <- n*alpha/lambda^2
    info <- matrix(c(info.11,info.12,info.12,info.22),ncol=2)
    r <- list(alpha=alpha,lambda=lambda,loglik=loglik,info=info)
    r
  }
}

alpha <- 1
lambda <- 1
x <- rgamma(10000, shape = alpha, rate = lambda)
r <- gamma.nr(x)
r
# estimated variance-covariance matrix
solve(r$info)
