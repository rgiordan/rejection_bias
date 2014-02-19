
library(mvtnorm)


SimulateMVNorm <- function(n, beta, sigma.chol) {
  # I do this myself to save the time calculating the cholesky decomposition each
  # time.
  k <- length(beta)
  stopifnot(k == ncol(sigma.chol))
  x <- matrix(rnorm(n * k), nrow=n)
  return(x %*% sigma.chol + beta)
}

SimulatRejectionProb <- function(beta, x.sim, x.lower, x.upper, n=1e4) {
  AllRejected <- function(x.row) {
    all((x.row + beta) < x.lower | (x.row + beta) > x.upper)
  }
  stopifnot(length(x.lower) == length(beta) && length(x.upper) == length(beta))
  stopifnot(all(x.lower <= x.upper))
  
  # Dig that you don't actually need to resimulate x each time -- just recenter it.
  z <- apply(x.sim, MARGIN=1, AllRejected) 
  return(sum(z) / length(z))
}



ConditionalLogLik <- function(x, x.sim, x.lower, x.upper, beta, sigma) {
  if (any(x > x.lower & x < x.upper)) {
    return(-Inf)
  }
  l <- dmvnorm(x=x, mean=beta, sigma=sigma, log=TRUE)
  rejection.prob <- SimulatRejectionProb(beta=beta, x.sim=x.sim,
                                         x.lower=x.lower, x.upper=x.upper)
  log.rejection.prob <- log(rejection.prob + 1e-9)
  
  log.lik <- l - log.rejection.prob
  #print(log.lik)
  return(log.lik)
}




PackageAreaCalculation <- function(beta, sigma, x.lower, x.upper, t=0.5) {
  
  # Just need to guarantee that this is within the bounds.
  x <- t * x.upper + (1 - t) * x.lower
  base.log.lik <- dtmvnorm(x=x, mean=beta, sigma=sigma, log=TRUE)
  
  k <- length(beta)
  GetLogArea <- function(i) {
    lower <- rep(-Inf, k)
    upper <- rep(Inf, k)
    lower[i] <- x.lower[i]
    upper[i] <- x.upper[i]
    return(base.log.lik - dtmvnorm(x=x, mean=beta, sigma=sigma,
                                   lower=lower, upper=upper, log=TRUE))
  }
  areas <- exp(sapply(1:k, GetLogArea))
  #print(areas)
  #print("--------")
  joint.area <-  exp(base.log.lik - dtmvnorm(x=x, mean=beta, sigma=sigma,
                                             lower=x.lower, upper=x.upper, log=TRUE))
  #print(joint.area)
  #print("--------")
  
  total.area <- sum(areas) - (k - 1) * joint.area
  return(total.area)
}



PackageConditionalLogLik <- function(x, x.lower, x.upper, beta, sigma) {
  if (any(x > x.lower & x < x.upper)) {
    return(-Inf)
  }
  l <- dmvnorm(x=x, mean=beta, sigma=sigma, log=TRUE)
  rejection.prob <- 1 - PackageAreaCalculation(beta=beta, sigma=sigma,
                                               x.lower=x.lower, x.upper=x.upper)
  if (rejection.prob < 0) {
    print(rejection.prob)
  }
  stopifnot(rejection.prob > 0)
  
  log.rejection.prob <- log(rejection.prob + 1e-9)
  
  log.lik <- l - log.rejection.prob
  #print(log.lik)
  return(log.lik)
}


