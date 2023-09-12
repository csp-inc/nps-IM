library(dplyr)


# ---- posterior predictive loss -----------------------------------------------

get_ppl <- function(mcarray_list, data, likelihood) {
  y_rep_obj <- ifelse(likelihood == "hurdle-ordinal-latent-beta",
    "y.rep.beta", "y.rep"
  )

  y_sim_mean <- summary(mcarray_list[[y_rep_obj]], mean)$stat
  y_sim_sd <- summary(mcarray_list[[y_rep_obj]], sd)$stat

  data_for_ppl <- data %>% get_data_for_ppl(likelihood)
  y_obs <- get_y_from_data_for_ppl(data_for_ppl)

  sum((y_obs - y_sim_mean)^2) + sum(y_sim_sd^2)
}

get_data_for_ppl <- function(data, likelihood) {
  data %>%
    {
      `if`(
        likelihood == "hurdle-ordinal-latent-beta",
        filter(., response > 0), .
      )
    } %>%
    drop_na(matches("response|hits"))
}

get_y_from_data_for_ppl <- function(data) {
  c(data[["hits"]], data[["response"]])
}

# ---- deviance information criterion ------------------------------------------

get_dic <- function(mcarray_list, data, likelihood, jags_data) {
  data_for_ppl <- data %>% get_data_for_ppl(likelihood)
  y_obs <- get_y_from_data_for_ppl(data_for_ppl)

  if (likelihood == "gen-pois" |
    likelihood == "zero-inflated-poisson" |
    likelihood == "zero-inflated-negative-binomial") {
    return(NA)
  }
  if (likelihood == "poisson") {
    mean_mu <- summary(mcarray_list$mu, mean)$stat
    d_hat <- -2 * (sum(dpois(x = y_obs, lambda = mean_mu, log = TRUE)))
  } else if (grepl("negative-binomial", likelihood)) {
    mean_mu <- summary(mcarray_list$mu, mean)$stat
    mean_kappa <- summary(mcarray_list$kappa, mean)$stat
    d_hat <- -2 * (sum(dnbinom(
      x = y_obs,
      mu = mean_mu, size = mean_kappa[jags_data$y.strata],
      log = TRUE
    )))
  } else if (likelihood == "binomial" |
    likelihood == "beta-binomial" |
    likelihood == "zero-inflated-binomial" |
    likelihood == "zero-inflated-beta-binomial") {
    mean_p <- summary(mcarray_list$p, mean)$stat
    d_hat <- -2 * (sum(dbinom(
      x = y_obs,
      prob = mean_p, size = data_for_ppl[["trials"]],
      log = TRUE
    )))
  } else if (likelihood == "beta") {
    mean_alpha <- summary(mcarray_list$alpha, mean)$stat
    mean_beta <- summary(mcarray_list$beta, mean)$stat
    d_hat <- -2 * (sum(dbeta(y_obs, mean_alpha, mean_beta, log = TRUE)))
  } else if (likelihood == "lognormal") {
    mean_log_mu <- summary(log(mcarray_list$mu), mean)$stat
    mean_sigma_log_y <- summary(mcarray_list$sigma.log.y, mean)$stat
    d_hat <- -2 * (sum(dlnorm(
      x = y_obs, mean_log_mu, mean_sigma_log_y,
      log = TRUE
    )))
  } else if (likelihood == "gamma") {
    mean_mu <- summary(mcarray_list$mu, mean)$stat
    mean_sigma_y <- summary(mcarray_list$sigma.y, mean)$stat
    d_hat <- -2 * (sum(dgamma(
      x = y_obs,
      mean_mu^2 / mean_sigma_y^2,
      mean_mu / mean_sigma_y^2,
      log = TRUE
    )))
  }

  if (likelihood == "hurdle-ordinal-latent-beta") {
    DIC <- beta_hurdle_DIC(jags_data, mcarray_list)
  } else if (likelihood == "ordinal-latent-normal") {
    DIC <- oln_DIC(jags_data, mcarray_list)
  } else {
    d_bar <- summary(mcarray_list$deviance, mean)$stat
    pD_DIC <- d_bar - d_hat
    DIC <- d_hat + 2 * pD_DIC
  }

  DIC
}

beta_hurdle_DIC <- function(data, jags_samples) {
  B <- summary(jags_samples$B, mean)$stat
  if ("Beta" %in% names(jags_samples)) {
    Beta <- summary(jags_samples$Beta, mean)$stat
    X <- data$X
  } else {
    Beta <- 0 # a hack for the calculation below to ensure this piece is 0
    X <- matrix(numeric(length(data$y.beta)), ncol = 1)
  }
  sigma <- summary(jags_samples$sigma.ss, mean)$stat
  Dev.Beta.mn <- summary(jags_samples$dev.Beta, mean)$stat # the mean of the deviance
  mu <- a <- b <- deltaF <- numeric(length(data$y.beta))
  lims <- data$lim.dint # I think this should be lim.dint, not lim0!!
  y.site <- data$y.site
  y.strata <- data$y.strata
  x <- data$x
  for (k in 1:length(data$y.beta)) {
    lp_add_cov <- X[k, ] %*% Beta
    mu[k] <- inv_logit(B[y.site[k], 1, y.strata[k]] +
      ifelse(dim(B)[2] == 2, B[y.site[k], 2, y.strata[k]] * x[k], 0) +
      lp_add_cov[, , drop = TRUE])
    a[k] <- (mu[k]^2 - mu[k]^3 - mu[k] * sigma[y.site[k], y.strata[k]]^2) /
      sigma[y.site[k], y.strata[k]]^2
    b[k] <- (mu[k] - 2 * mu[k]^2 + mu[k]^3 - sigma[y.site[k], y.strata[k]]^2 + mu[k] * sigma[y.site[k], y.strata[k]]^2) /
      sigma[y.site[k], y.strata[k]]^2
    deltaF[k] <-
      pbeta(lims[data$y.beta[k] + 1], a[k], b[k]) - pbeta(lims[data$y.beta[k]], a[k], b[k])
  }

  # deltaF: the probability associated with the interval between the lower and the upper cutpoint
  plike <- log(deltaF[deltaF > 0]) # the log of deltaF (just the nonnegative values)
  llike <- sum(plike) # the total log likelihood
  D.thetabar.Beta <- -2 * llike # the first part of the DIC calculation
  pd.Beta <- Dev.Beta.mn - D.thetabar.Beta # the penalty (the effective number of parameters)
  DIC.Beta <- Dev.Beta.mn + 2 * pd.Beta

  # DIC is defined as -2 * (the log of the probability of the data conditional on the model parameters)
  DIC.Beta
}

oln_DIC <- function(data, jags_samples) {
  B <- summary(jags_samples$B, mean)$stat
  if ("Beta" %in% names(jags_samples)) {
    Beta <- summary(jags_samples$Beta, mean)$stat
    X <- data$X
  } else {
    Beta <- 0 # a hack for the calculation below to ensure this piece is 0
    X <- matrix(numeric(length(data$y)), ncol = 1)
  }
  sigma <- summary(jags_samples$sigma.ss, mean)$stat
  deviance.mean <- summary(jags_samples$deviance, mean)$stat # the mean of the deviance
  mu <- a <- b <- deltaF <- numeric(length(data$y))

  lims <- c(-Inf, summary(jags_samples$theta, mean)$stat, Inf)
  y.site <- data$y.site
  y.strata <- data$y.strata
  x <- data$x

  for (k in 1:length(data$y)) {
    lp_add_cov <- X[k, ] %*% Beta
    mu[k] <- B[y.site[k], 1, y.strata[k]] +
      B[y.site[k], 2, y.strata[k]] * x[k] +
      lp_add_cov[, , drop = TRUE]
    deltaF[k] <- pnorm(lims[data$y[k] + 1], mu[k], sigma[y.site[k], y.strata[k]]) -
      pnorm(lims[data$y[k]], mu[k], sigma[y.site[k], y.strata[k]])
  }

  # deltaF: the probability associated with the interval between the lower and the upper cutpoint
  plike <- log(deltaF[deltaF > 0]) # the log of deltaF (just the nonnegative values)
  llike <- sum(plike) # the total log likelihood
  D.thetabar.oln <- -2 * llike # the first part of the DIC calculation
  pd.oln <- deviance.mean - D.thetabar.oln # the penalty (the effective number of parameters)
  DIC.oln <- deviance.mean + 2 * pd.oln

  # DIC is defined as -2 * (the log of the probability of the data conditional on the model parameters)
  DIC.oln
}
