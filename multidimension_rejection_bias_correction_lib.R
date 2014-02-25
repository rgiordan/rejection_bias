
library(mvtnorm)


SimulatRejectionProb <- function(beta, x.sim, x.lower, x.upper) {
  # For a given beta, estimate the rejection probability. with simulated data.
  #
  # Args:
  #   beta:    A k-vector of the centers of the distribution.
  #   x.sim:   An n by k matrix of simulations from the desired distribution centered
  #            at zero.  Here, n should be large to get good simulated results.
  #   x.lower: A k-vector of lower bounds for rejecting the kth column of <x.sim>.
  #   x.upper: A k-vector of upper bounds for rejecting the kth column of <x.sim>.
  #
  # Returns:
  #   The proportion of rows of <x.sim> that are outside the bounds specified
  #   by <x.lower> and <x.upper> when it is recentered to <beta>.  This is an
  #   estimate of the probability of rejecting data whose covariance structure
  #   is approximated by <x.sim> is the data were to be centered at <beta>.
  
  AllRejected <- function(x.row) {
    # A boolean whether or not all the columns of x.row were rejected according
    # to <x.lower> and <x.upper>, when re-centered at <beta>.
    return(all((x.row + beta) < x.lower | (x.row + beta) > x.upper))
  }
  
  # Sanity check the inputs.
  stopifnot(length(x.lower) == length(beta) && length(x.upper) == length(beta))
  stopifnot(all(x.lower <= x.upper))

  z <- apply(x.sim, MARGIN=1, AllRejected) 
  return(sum(z) / length(z))
}



SimulatedConditionalLogLik <- function(x, x.sim, x.lower, x.upper, beta, sigma) {
  # Use simulated data to return log likelihood conditional on rejection.
  #
  # Args:
  #   x:       The actual observed data (e.g. your regression estimates)
  #   x.sim:   An n by k matrix of n draws from the covariance structure of 
  #            <x>, but centered at zero.  This will be used to estimate the
  #            probability of rejection.  See SimulatRejectionProb for details.
  #   x.lower: A k-vector of lower bounds for rejecting the kth element of <x>.
  #   x.upper: A k-vector of upper bounds for rejecting the kth element of <x>.
  #   beta:    The center of the distribution at which the likelihood is to
  #            be evaluated (e.g. the candidate true regression parameters).
  #   sigma:   The covariance matrix of the observations.  Note that <x.sim>
  #            is expected to be drawn from this covariance structure.
  #
  #  Returns:
  #    An estimated log likelihood of the observations <x>, given that they were drawn
  #    from a multivariate normal distribution centered at <beta> and with
  #    covariance <sigma>, conditional on each <x> being outside the bounds
  #    <x.lower> to <x.upper>.
  #
  #   This uses <x.sim> to estimate the probability of being outside the bounds.
  
  # The likelihood is zero conditional on an observation not being rejected.
  if (any(x > x.lower & x < x.upper)) {
    return(-Inf)
  }
  
  l <- dmvnorm(x=x, mean=beta, sigma=sigma, log=TRUE)
  rejection.prob <- SimulatRejectionProb(beta=beta, x.sim=x.sim,
                                         x.lower=x.lower, x.upper=x.upper)
  
  # Add 1e-9 to prevent numerical issues.
  log.rejection.prob <- log(rejection.prob + 1e-9)
  log.lik <- l - log.rejection.prob
  
  return(log.lik)
}


NumericAcceptanceProbability <- function(beta, sigma, x.lower, x.upper, t=0.5, force=FALSE) {
  # The probability of at least one multivariate normal draw being within bounds.
  #
  # Args:
  #   beta:    The k-dimensional center of the multivariate normal.
  #   sigma:   The k by k covariance matrix of the multivariate normal.
  #   x.lower: A k-vector of lower bounds for rejecting the kth element.
  #   x.upper: A k-vector of upper bounds for rejecting the kth element.
  #   t:       Used to select a point within the bounds.  For debugging only.
  #
  # Returns:
  #   A numeric approximation of the probability that a draw from the multivariate
  #   normal centered at <beta> with covariance <sigma> has at least one component
  #   between its corresponding <x.lower> and <x.upper>.
  #
  #   This function backs this probability out of the full likelihood returned by
  #   dtmvnorm, which is rather a hack and requires an extra likelihood evaluation.
  #   It could probably be done more effeciently.
  
  # <x> is the point at which we will evaluate the likelihood.  All that matters
  # is that it be within the bounds.  As long as this is the case, the answer
  # will not depend on the particular x you choose.
  stopifnot(t > 0 && t < 1)
  x <- t * x.upper + (1 - t) * x.lower
  
  # The way I've handled the areas only works for d=2
  if (length(beta) > 2 && !force) {
    stop("Sorry, this doesn't work yet in dimension > 2.")
  }

  k <- length(beta)
  GetLogArea <- function(i) {
    # Return the marginal log probabilty of the i_th component lying within
    # the bounds (x.lower[i], x.upper[i]).
    base.log.lik <- dtmvnorm(x=x[i], mean=beta[i], sigma=sigma[i, i], log=TRUE)
    return(base.log.lik - dtmvnorm(x=x[i], mean=beta[i], sigma=sigma[i, i],
                                   lower=x.lower[i], upper=x.upper[i], log=TRUE))
  }
  
  areas <- exp(sapply(1:k, GetLogArea))
  
  # This is the area of the intersection of all the hyper-rectangular strips
  # of the individual variables.  This area is counted <k>-1 extra times in
  # <areas>.
  base.log.lik <- dtmvnorm(x=x, mean=beta, sigma=sigma, log=TRUE)
  joint.area <-  exp(base.log.lik - dtmvnorm(x=x, mean=beta, sigma=sigma,
                                             lower=x.lower, upper=x.upper, log=TRUE))
  total.area <- sum(areas) - (k - 1) * joint.area
  return(total.area)
}


NumericConditionalLogLik <- function(x, x.lower, x.upper, beta, sigma) {
  # Use numeric area probabilitiess to return log likelihood conditional on rejection.
  #
  # Args:
  #   x:       The actual observed data (e.g. your regression estimates)
  #   x.lower: A k-vector of lower bounds for rejecting the kth element of <x>.
  #   x.upper: A k-vector of upper bounds for rejecting the kth element of <x>.
  #   beta:    The center of the distribution at which the likelihood is to
  #            be evaluated (e.g. the candidate true regression parameters).
  #   sigma:   The covariance matrix of the observations.  Note that <x.sim>
  #            is expected to be drawn from this covariance structure.
  #
  #  Returns:
  #    An estimated log likelihood of the observations <x>, given that they were drawn
  #    from a multivariate normal distribution centered at <beta> and with
  #    covariance <sigma>, conditional on each <x> being outside the bounds
  #    <x.lower> to <x.upper>.
  #
  #   This uses NumericAcceptanceProbability to estimate the probability of being
  #   outside the bounds.
  
  # The likelihood of any observation being within the bounds is zero.
  if (any(x > x.lower & x < x.upper)) {
    return(-Inf)
  }
  
  # A bit rude of them not to check this in dmvnorm I think.
  sigma <- as.matrix(sigma)
  
  l <- dmvnorm(x=x, mean=beta, sigma=sigma, log=TRUE)
  rejection.prob <- 1 - NumericAcceptanceProbability(beta=beta, sigma=sigma,
                                                     x.lower=x.lower, x.upper=x.upper)
  
  # This can happen if NumericAcceptanceProbability has numeric issues.
  if (rejection.prob < 0 || rejection.prob > 1) {
    stop(sprintf("Rejection probability = %f.  This indicates numeric problems.",
                 rejection.prob))
  }
  
  # Add 1e-9 to avoid numeric issues.
  log.rejection.prob <- log(rejection.prob + 1e-9)
  
  log.lik <- l - log.rejection.prob
  return(log.lik)
}



SimulatedEstimateTruncatedMean <- function(x, x.lower, x.upper, sigma,
                                           verbose=F, n=1e4, ...) {
  # Use simulation and optim to find the MLE of the center of the truncated distribution.
  #
  # Args:
  #     n:  The number of data points to simulate.  The higher this number is,
  #         the more accurate the estimated rejection probabilities will be,
  #         but the slower the estimation will be.
  #     For other args, see SimulatedConditionalLogLik.
  #
  # Returns:
  #     The output of optim applied to SimulatedConditionalLogLik.
  
  x.sim <- rmvnorm(n=n, mean=0, sigma=sigma)
  SimulatedConditionalLogLikForOptim <- function(beta) {
    if (verbose) {
      print(beta)
    }
    return(-SimulatedConditionalLogLik(x=x, x.sim=x.sim, beta=beta,
                                       x.lower=x.lower, x.upper=x.upper,
                                       sigma=sigma))
  }
  
  results <- optim(par=0.9 * x, fn=SimulatedConditionalLogLikForOptim,
                   method="L-BFGS-B", lower=rep(0, length(x)), upper=x, ...)
  return(results) 
}


NumericEstimateTruncatedMean <- function(x, x.lower, x.upper, sigma,
                                         verbose=F, ...) {
  # Use optim to find the MLE of the center of the truncated distribution.
  #
  # Args:
  #     See NumericConditionalLogLik.
  #
  # Returns:
  #     The output of optim applied to NumericConditionalLogLik.
  
  NumericConditionalLogLikForOptim <- function(beta) {
    if (verbose) {
      print(beta)
    }
    return(-NumericConditionalLogLik(x=x, beta=beta,
                                     x.lower=x.lower, x.upper=x.upper,
                                     sigma=sigma))
  }
  
  results <- optim(par=0.9 * x, fn=NumericConditionalLogLikForOptim,
                   method="L-BFGS-B", lower=rep(0, length(x)), upper=x, ...)
  return(results) 
}