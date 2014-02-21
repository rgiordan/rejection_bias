# This is a preliminary stab at an empirical Bayes solution and is
# not yet ready for public consumption.


# Normal Log Lik
NormLogLik <- function(x, m, sigma) {
  log.lik <- -log(sigma) - ((x - m)^2) / (2 * sigma^2)
  return(sum(log.lik))
}


# Number of equations
K <- 20

# Number of observations per equation
N <- 100

# Predictors.  Each column is an equation.
x <- matrix(runif(N * K), nrow=N, ncol=K)

# True mean and variance of your betas
beta.mean <- 0
beta.sigma <- 0.3
beta <- rnorm(K, mean=beta.mean, sd=beta.sigma)

# Sanity check the log likelihood.  par[1] is the mean and par[2] is the
# log of the standard deviation.
RawBetaLogLik <- function(par) {
  m <- par[1]
  s <- exp(par[2])
  return(-NormLogLik(beta, m, s))
}
beta.fit <- optim(c(0, 0), RawBetaLogLik)$par
print("Fit to prior if we actually observed the betas.")
cat(sprintf("Actual mean: %f, fit mean: %f.  \nActual sd: %f, fit sd %f.",
        beta.mean, beta.fit[1], beta.sigma, exp(beta.fit[2])))

# Outcome variables.  Each column is an equation with its own beta.
sigma <- 0.1
y <- x * rep(beta, each=N) + rnorm(K * N, mean=0, sd=sigma)

beta.hat <- colSums(y * x) / colSums(x * x)
y.hat <- x * rep(beta.hat, each=N)
sigma.hat <- sqrt(colSums((y - y.hat) ^ 2) / (N - 1))
beta.hat.sd <- sigma.hat / sqrt(colSums(x * x))


# Method of moments EB on the prior of beta
print("Fit to prior using method of moments")
cat(sprintf("Actual mean: %f, fit mean: %f.\nActual sd: %f, fit sd %f.",
        beta.mean, mean(beta.hat), beta.sigma, sqrt(var(beta.hat) - sum(beta.hat.sd^2))))



# Fit the same thing with log likelihood
BetaPriorLogLik <- function(par) {
  m <- par[1]
  s <- exp(par[2])
  return(-NormLogLik(beta.hat, m, s + beta.hat.sd))
}
beta.prior.fit <- optim(c(0, 0), BetaPriorLogLik)$par
print("Fit to prior using maximum likelihood")
cat(sprintf("Actual mean: %f, fit mean: %f.\nActual sd: %f, fit sd %f.",
        beta.mean, beta.prior.fit[1], beta.sigma, exp(beta.prior.fit[2])))






