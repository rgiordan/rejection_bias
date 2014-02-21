
library(ggplot2)
library(tmvtnorm)

setwd("~/Documents/rejection_bias")
source("~/Documents/rejection_bias/multidimension_rejection_bias_correction_lib.R")


##############################
# 1d simulation with no ground truth.

beta.hat <- 2.2
beta.hat.lower <- -2
beta.hat.upper <- 2
sigma <- matrix(0.1)

# Check the area properties. 
x.sim <- rmvnorm(n=1e4, mean=0, sigma=sigma)
1 - NumericAcceptanceProbability(beta=0, sigma=sigma, x.lower=beta.hat.lower, x.upper=beta.hat.upper)
SimulatRejectionProb(beta=0, x.sim=x.sim, x.lower=beta.hat.lower, x.upper=beta.hat.upper)

size <- 0.05
1 - NumericAcceptanceProbability(beta=0, sigma=sigma, x.lower=-size, x.upper=size)
SimulatRejectionProb(beta=0, x.sim=x.sim, x.lower=-size, x.upper=size)

# Note that the simulation method is really slow, even in one dimension,
# and demonstrably pretty random when <n> is small.  I will not use it from
# here on out.
results.simulation <- SimulatedEstimateTruncatedMean(x=beta.hat, sigma=sigma,
                                                     x.lower=beta.hat.lower,
                                                     x.upper=beta.hat.upper,
                                                     n=5e4, verbose=T)

results.numeric <- NumericEstimateTruncatedMean(x=beta.hat, sigma=sigma,
                                                x.lower=beta.hat.lower,
                                                x.upper=beta.hat.upper, verbose=T)

print(results.simulation$par)
print(results.numeric$par)



############
# 2d simulation
beta.hat <- c(2.2, 2.2)
beta.hat.lower <- rep(-2, 2)
beta.hat.upper <- rep(2, 2)

z <- runif(2)
sigma <- z %*% t(z) + diag(0.8, 2)

size <- 0.05
1 - NumericAcceptanceProbability(beta=c(0, 0), sigma=sigma, x.lower=rep(-size, 2), x.upper=rep(size, 2))
SimulatRejectionProb(beta=c(0, 0), x.sim=x.sim, x.lower=rep(-size, 2), x.upper=rep(size, 2))

results <- NumericEstimateTruncatedMean(x=beta.hat, sigma=sigma,
                                        x.lower=beta.hat.lower,
                                        x.upper=beta.hat.upper, verbose=T)
print(results$par)


############
# Arbitrary dimension simulation.  This seems broken for d >= 3.
if (F) {
  k <- 3
  beta.hat <- rep(2.2, k)
  beta.hat.lower <- rep(-2, k)
  beta.hat.upper <- rep(2, k)
  
  z <- runif(k)
  sigma <- 0.1 * z %*% t(z) + diag(1, k)
  
  x.sim <- rmvnorm(n=1e4, mean=rep(0, k), sigma=sigma)
  
  size <- 0.1
  1 - NumericAcceptanceProbability(beta=rep(0, k), sigma=sigma, x.lower=rep(-size, k), x.upper=rep(size, k))
  SimulatRejectionProb(beta=rep(0, k), x.sim=x.sim, x.lower=rep(-size, k), x.upper=rep(size, k))
  
  
  results <- NumericEstimateTruncatedMean(x=beta.hat, sigma=sigma,
                                          x.lower=beta.hat.lower,
                                          x.upper=beta.hat.upper, verbose=T)
  print(results$par)
}

#####################################
# Simulation with ground truth

# This just has to be large enough to get two rejections on average.
k <- 50
true.beta <- rep(2, k)
beta.hat.lower.val <- -2
beta.hat.upper.val <- 2
dimension <- 1

z <- runif(k)
sigma <- z %*% t(z) + diag(1, k)

all.beta.hat <- rmvnorm(n=1, mean=beta, sigma=sigma)
reject <- all.beta.hat < beta.hat.lower.val | all.beta.hat > beta.hat.upper.val
print(k.rej <- sum(reject))

# Choose the first <dimension>, which should be random.
beta.hat <- as.numeric(all.beta.hat[reject])[1:dimension]
beta.hat.sigma <- (sigma[reject, reject])[1:dimension, 1:dimension, drop=F]

# Set up the variables needed for the log likelihood
beta.hat.lower <- rep(beta.hat.lower.val, dimension)
beta.hat.upper <- rep(beta.hat.upper.val, dimension)

results <- NumericEstimateTruncatedMean(x=beta.hat, sigma=beta.hat.sigma,
                                        x.lower=beta.hat.lower,
                                        x.upper=beta.hat.upper, verbose=T)
print(results$par)



###########################
# Playground and plots

# Plot simulated stuff

# You may as well keep this constant.
beta.hat.sigma <- matrix(1)

# Set these and play
beta.hat <- 5
range <- beta.hat - 0.5
beta.hat.lower <- -range
beta.hat.upper <- range

results <- NumericEstimateTruncatedMean(x=beta.hat, sigma=beta.hat.sigma,
                                        x.lower=beta.hat.lower,
                                        x.upper=beta.hat.upper, verbose=T)
print(results$par)

grid.values <- seq(-1, beta.hat * 1.5, length.out=500)
lik.grid <- sapply(grid.values,
                   function(x) {
                     NumericConditionalLogLik(x=beta.hat, beta=x,
                                              x.lower=beta.hat.lower, x.upper=beta.hat.upper,
                                              sigma=beta.hat.sigma)
                   })

ggplot() + geom_line(aes(x=grid.values, y=exp(lik.grid))) +
  geom_vline(aes(xintercept=beta.hat), col="blue") +
  geom_vline(aes(xintercept=beta.hat.upper), color="gray", lwd=2) +
  geom_vline(aes(xintercept=results$par), color="red")
  

