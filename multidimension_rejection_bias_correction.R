
library(ggplot2)
library(tmvtnorm)

setwd("~/Documents/rejection_bias")
source("~/Documents/rejection_bias/multidimentsion_rejection_bias_correction_lib.R")

# Below, k always denotes the dimension and n the data points in the simulation.


#############################
# Sanity check the area calculations.

k <- 10
n <- 1e4

x.lower <- rep(-2, k)
x.upper <- rep(2, k)

beta <- rep(2, k)
z <- runif(k)
sigma <- z %*% t(z) + diag(0.1, k)
sigma.chol <- chol(sigma)

SimulatRejectionProb(beta=rep(0, k), sigma.chol=sigma.chol, x.lower=x.lower, x.upper=x.upper)
SimulatRejectionProb(beta=rep(2, k), sigma.chol=sigma.chol, x.lower=x.lower, x.upper=x.upper)
SimulatRejectionProb(beta=rep(20, k), sigma.chol=sigma.chol, x.lower=x.lower, x.upper=x.upper)






SimulatedConditionalLogLikForOptim <- function(beta) {
  print(beta)
  return(-SimulatedConditionalLogLik(x=beta.hat, beta=beta,
                                     x.lower=beta.hat.lower, x.upper=beta.hat.upper,
                                     sigma=beta.hat.sigma, x.sim=x.sim))
}



##############################
# 1d simulation with no ground truth.

beta.hat <- 2.2
beta.hat.lower <- -2
beta.hat.upper <- 2
sigma <- matrix(0.1)

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



############3
# 2d simulation
beta.hat <- c(2.2, 2.2)
beta.hat.lower <- rep(-1, 2)
beta.hat.upper <- rep(1, 2)

#z <- runif(2)
#sigma <- z %*% t(z) + diag(0.1, 2)
sigma <- diag(0.8^2, 2)
sigma.chol <- chol(sigma)
beta.hat.sigma <- sigma

#x.sim <- SimulateMVNorm(n=1e5, beta=rep(0, 2), sigma.chol=sigma.chol)
x.sim <- rmvnorm(n=1e3, mean=rep(0, 2), sigma=sigma)

size <- 0.05
1 - NumericAcceptanceProbability(beta=c(0, 0), sigma=sigma, x.lower=rep(-size, 2), x.upper=rep(size, 2))
SimulatRejectionProb(beta=c(0, 0), x.sim=x.sim, x.lower=rep(-size, 2), x.upper=rep(size, 2))

results <- optim(par=rep(0, k), fn=SimulatedConditionalLogLikForOptim,
                 method="L-BFGS-B", lower=0, upper=beta.hat)

results.package <- optim(par=rep(0, k), fn=NumericConditionalLogLikForOptim,
                         method="L-BFGS-B", lower=0, upper=beta.hat)


############
# high-d simulation
k <- 20
beta <- rep(2, k)
beta.hat.lower.val <- -2
beta.hat.upper.val <- 2

z <- runif(k)
sigma <- z %*% t(z) + diag(1, k)
beta.hat.sigma <- sigma
sigma.chol <- chol(sigma)

x.sim <- rmvnorm(n=1e4, mean=rep(0, k), sigma=sigma)

size <- 0.05
1 - NumericAcceptanceProbability(beta=rep(0, k), sigma=sigma, x.lower=rep(-size, k), x.upper=rep(size, k))
SimulatRejectionProb(beta=rep(0, k), x.sim=x.sim, x.lower=rep(-size, k), x.upper=rep(size, k))


results <- optim(par=rep(0, k), fn=SimulatedConditionalLogLikForOptim,
                 method="L-BFGS-B", lower=0, upper=beta.hat)

results.package <- optim(par=rep(0, k), fn=NumericConditionalLogLikForOptim,
                         method="L-BFGS-B", lower=0, upper=beta.hat)



###################################
# MV simulation
k <- 50
beta <- rep(2, k)
beta.hat.lower.val <- -2
beta.hat.upper.val <- 2

z <- runif(k)
sigma <- z %*% t(z) + diag(1, k)

all.beta.hat <- SimulateMVNorm(1, beta, chol(sigma))
reject <- all.beta.hat < beta.hat.lower.val | all.beta.hat > beta.hat.upper.val
print(k.rej <- sum(reject))
beta.hat <- as.numeric(all.beta.hat[reject])
beta.hat.sigma <- sigma[reject, reject]
beta.hat.sigma.chol <- chol(beta.hat.sigma)

# Set up the variables needed for the log likelihood
x.sim <- SimulateMVNorm(n=1e3, beta=rep(0, k.rej), sigma.chol=beta.hat.sigma.chol)
beta.hat.lower <- rep(beta.hat.lower.val, k.rej)
beta.hat.upper <- rep(beta.hat.upper.val, k.rej)

# Much better
results.package <- optim(par=rep(0, k.rej), fn=NumericConditionalLogLikForOptim,
                         method="L-BFGS-B", lower=0, upper=beta.hat)




# Super slow
results <- optim(par=rep(0, k.rej), fn=SimulatedConditionalLogLikForOptim,
                 method="L-BFGS-B", lower=rep(0, k.rej), upper=beta.hat)

x.sim <- SimulateMVNorm(n=1e4, beta=rep(0, sum(reject)), sigma.chol=beta.hat.sigma.chol)
results2 <- optim(par=results$par, fn=SimulatedConditionalLogLikForOptim,
                  method="L-BFGS-B", lower=rep(0, k.rej), upper=beta.hat)

