# Moment match for the beta distribution.
params_from_moments <- function(mu, sigma) {
  alpha <- (mu^2 - mu^3 - mu * sigma^2) / sigma^2
  beta <- (mu - 2 * mu^2 + mu^3 - sigma^2 + mu * sigma^2) / sigma^2
  c(alpha = alpha, beta = beta)
}
moments_from_params <- function(alpha, beta) {
  mu <- alpha / (alpha + beta)
  sigma <- sqrt((alpha * beta) / ((alpha + beta)^2 * (alpha + beta + 1)))
  c(mu = mu, sigma = sigma)
}
jitter_guess <- function(guess, w = 0.1) {
  boot::inv.logit(rnorm(1, boot::logit(guess), w))
}
get_inits <- function(p_mean_guess, p_sd_guess, max_tries = 1000) {
  z <- NaN
  n_tries <- 1
  while (is.nan(z) & n_tries <= max_tries) {
    params_guess <-
      params_from_moments(jitter_guess(p_mean_guess), jitter_guess(p_sd_guess))
    guess_info <- paste(round(params_guess, 2), collapse = ", ")
    message(glue("Testing validity of guess {n_tries}: ({guess_info})"))
    z <- rbeta(1, params_guess["alpha"], params_guess["beta"])
    if (n_tries == max_tries) stop("Failed to find a valid init!")
    n_tries <- n_tries + 1
  }
  moments_guess <-
    moments_from_params(params_guess[["alpha"]], params_guess[["beta"]])
  c(mu = moments_guess[["mu"]], sigma = moments_guess[["sigma"]])
}
