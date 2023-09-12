# function to calculate GP CDF at an integer x
# when delta = 0. This does the same thing as ppois()
gp.cdf <- function(x, mu, delta) {
  # grid of values at which to calculate the mass function
  y <- seq(0, x, by = 1)
  # return CDF value
  return(sum(exp(-mu - delta * y) * (mu * (mu + delta * y)^(y - 1)) / factorial(y)))
}



# n: sample size
rgp <- function(n = n, mu = mu, delta = delta) {
  U <- runif(n)
  X <- rep(0, n)
  if (length(mu) > 1) {
    mu <- mu
  } else {
    (mu <- rep(mu, n))
  }

  # loop through each uniform
  for (i in 1:n)
  {
    # first check if you are in the first interval
    if (U[i] < gp.cdf(0, mu[i], delta)) {
      X[i] <- 0
    } else {
      # while loop to determine which subinterval,I, you are in
      # terminated when B = TRUE
      B <- FALSE
      I <- 0
      while (B == FALSE) {
        # the interval to check
        int <- c(gp.cdf(I, mu[i], delta), gp.cdf(I + 1, mu[i], delta))
        # see if the uniform is in that interval
        if ((U[i] > int[1]) & (U[i] < int[2])) {
          # if so, quit the while loop and store the value
          X[i] <- I + 1
          B <- TRUE
        } else {
          # If not, continue the while loop and increase I by 1
          I <- I + 1
        }
      }
    }
  }
  return(X)
}

get_gen_pois_y_reps <- function(z_jags, max_iters = 100) {
  n_obs <- dim(z_jags$lambda)[1]
  n_strata <- dim(z_jags$delta)[1]

  n_iters <- min(dim(z_jags$delta)[2], max_iters)
  n_chains <- dim(z_jags$delta)[3]

  gpy <- array(NA, dim = c(n_obs, n_iters, n_chains))

  for (this_stratum in 1:n_strata) { # this isn't doing anything!
    for (this_chain in 1:n_chains) {
      for (this_iter in 1:n_iters) { #
        print(sprintf(
          "progress: %s of %s for chain %s of %s",
          this_iter, n_iters, this_chain, n_chains
        ))

        # for(this_obs in 1:n_obs) {
        this_delta <- z_jags$delta[this_stratum, this_iter, this_chain]
        # this_lambda <- z_jags$lambda[this_obs, this_iter, this_chain]

        # gpy[this_obs, this_iter, this_chain] <-
        #   rgp(1, mu = this_lambda, delta = this_delta)
        these_lambdas <- z_jags$lambda[, this_iter, this_chain]
        gpy[, this_iter, this_chain] <-
          rgp(n_obs, mu = these_lambdas, delta = this_delta)
        # vrgp(n_obs, mu = these_lambdas, this_delta)
        # }
      }
    }
  }

  gpy
}

# vrgp <- Vectorize(rgp, 'mu')
