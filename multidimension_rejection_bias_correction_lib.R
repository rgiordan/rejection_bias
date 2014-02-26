
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
  #
  # TODO: This could surely benefit greatly from importance sampling.
  
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


NumericRejectionProbability <- function(beta, sigma, x.lower, x.upper) {
  # The probability of all the observations being out of bounds.
  #
  # Args:
  #   beta:    The k-dimensional center of the multivariate normal.
  #   sigma:   The k by k covariance matrix of the multivariate normal.
  #   x.lower: A k-vector of lower bounds for rejecting the kth element.
  #   x.upper: A k-vector of upper bounds for rejecting the kth element.
  #
  # Returns:
  #   A numeric approximation of the probability that a draw from the multivariate
  #   normal centered at <beta> with covariance <sigma> has at least all components
  #   between their corresponding <x.lower> and <x.upper>.
  #
  #   This function backs this probability out of the full likelihood returned by
  #   dtmvnorm, which is rather a hack and requires an extra likelihood evaluation.
  #   It could probably be done more effeciently.

  k <- length(beta)
  
  if (k > 8) {
    stop(paste("You probably don't want to run this with dimension > 8.",
               "The algorithm requires 2^d operations to handle a d-dimensional vector.",
               "For that, you should consider switching to simulations.",
               sep="\n"))
  }
  
  RecursiveGetRejectionArea <- function(current.lower, current.upper) {
    # A recursive function that ultimately returns the full rejection area.
    #
    # Args:
    #   current.lower: A vector of lower bounds, possibly incomplete.
    #   current.upper: A vector of upper bounds, possibly incomplete.
    #
    # Other variable usage:
    #   This function uses x.lower, x.upper, beta, and sigma from the containing
    #   namespace.
    #
    # Returns:
    #   The area of the rejection region corresponding to the probability that
    #   all the observations are outside the bounds specified by x.lower and
    #   x.upper in the containing namespace.
    #
    # The rejection region consists of 2^k infinite hyper-rectangles over which
    # we want to integrate the likelihood.  This function recursively calls
    # itself to evaluate the areas from (-Inf, x.lower[i]) and (x.upper[i], Inf)
    # for each combination of i in 1:k.  There will, of course, be 2^k such
    # area evaluations, so this is prohibitive for even moderate k.  (8 seems
    # to be a reasonable limit.)
    #
    # For large dimensions, simulations are probably better.
    # If you have better ideas for the analytic soultion I'd love to hear them
    # (Ryan). 
    
    k <- length(current.lower)
    stopifnot(length(current.lower) == length(current.upper))
    # If we are at the leaf of the tree, return the probability
    if (k == length(x.lower)) {
      # This particular point only needs to be in the desired region, otherwise,
      # its value doesn't matter.  Note that, for each i, either
      # current.lower[i] = -Inf or current.upper[i] = Inf.
      delta <- 0.1 * sqrt(diag(sigma))
      x <- rep(0, k)
      finite.lower.bound <- is.finite(current.lower)
      stopifnot(!any(finite.lower.bound & is.finite(current.upper)))
      x[finite.lower.bound] <- current.lower[finite.lower.bound] + delta[finite.lower.bound]
      x[!finite.lower.bound] <- current.upper[!finite.lower.bound] - delta[!finite.lower.bound]
      
      base.log.lik <- dtmvnorm(x=x, mean=beta, sigma=sigma, log=TRUE)
      return(exp(base.log.lik - dtmvnorm(x=x, mean=beta, sigma=sigma,
                                        lower=current.lower, upper=current.upper,
                                        log=TRUE)))
    } else {
      # Otherwise, proceed down the tree.
      return(RecursiveGetRejectionArea(current.lower=c(current.lower, x.upper[k + 1]),
                                       current.upper=c(current.upper, Inf)) +
             RecursiveGetRejectionArea(current.lower=c(current.lower, -Inf),
                                       current.upper=c(current.upper, x.lower[k + 1])))
    }
  }
  
  return(RecursiveGetRejectionArea(c(), c()))
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
  #   This uses NumericRejectionProbability to estimate the probability of being
  #   outside the bounds.
  
  # The likelihood of any observation being within the bounds is zero.
  if (any(x > x.lower & x < x.upper)) {
    return(-Inf)
  }
  
  # A bit rude of them not to check this in dmvnorm I think.
  sigma <- as.matrix(sigma)
  
  l <- dmvnorm(x=x, mean=beta, sigma=sigma, log=TRUE)
  rejection.prob <- NumericRejectionProbability(beta=beta, sigma=sigma,
                                                x.lower=x.lower, x.upper=x.upper)
  
  # This can happen if NumericRejectionProbability has numeric issues.
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