library(digest)
# library(extraDistr)

get_generating_values <- function(n_strata, link,
                                  # P, # TODO: number of covariates in W
                                  params = list(
                                    # Intercept hyperparameters.
                                    mu_mu_B0_l = 2, sigma_mu_B0_l = 1,
                                    mu_sigma_B0_l = 0.125 / 2, sigma_sigma_B0_l = 0.125 / 4,
                                    # Time-slope hyperparameters.
                                    mu_mu_B1_l = 0, sigma_mu_B1_l = 0.25,
                                    mu_sigma_B1_l = 0.125 / 4, sigma_sigma_B1_l = 0.125 / 8,
                                    # Covariate effects.
                                    mu_Beta_l = 0, sigma_Beta_l = 1,
                                    # Random effect parameter correlations.
                                    logit_mu_rho_l = 0, logit_sigma_rho_l = 1.5,
                                    # Error term(s).
                                    mu_sigma_y_l = 1, sigma_sigma_y_l = 0.5,
                                    mu_mu_sigma_y_l = 2, sigma_mu_sigma_y_l = 1,
                                    mu_sigma_sigma_y_l = 1, sigma_sigma_sigma_y_l = 0.5
                                  ),
                                  trials_upper = 100, N = 1,
                                  which_params = NULL,
                                  seed = NULL,
                                  params_override = list(),
                                  var_level = "stratum", var_type = "fixed") {

  # l: subscript for park
  # k: subscript for stratum
  # j: subscript for site (unused)
  # i: subscript for observation (unused)

  set.seed(seed)
  params <- modifyList(params, params_override)
  message(params$mu_mu_B0_l)

  h <- .h(link)

  # B0
  mu_B0_lk <- h(rnorm(n_strata, params$mu_mu_B0_l, params$sigma_mu_B0_l))
  sigma_B0_lk <- rgamma(
    n_strata,
    params$mu_sigma_B0_l^2 / params$sigma_sigma_B0_l^2,
    params$mu_sigma_B0_l / params$sigma_sigma_B0_l^2
  )

  # B1
  mu_B1_lk <- h(rnorm(n_strata, params$mu_mu_B1_l, params$sigma_mu_B1_l))
  sigma_B1_lk <- rgamma(
    n_strata,
    params$mu_sigma_B1_l^2 / params$sigma_sigma_B1_l^2,
    params$mu_sigma_B1_l / params$sigma_sigma_B1_l^2
  )

  # Beta
  mu_Beta_lk <- h(rnorm(1, params$mu_Beta_l, params$sigma_Beta_l)) # TODO: invoke n_strata

  # p.zero
  zi_terms <- if ("p.zero" %in% which_params) {
    list(p.zero = runif(n_strata, 0.2, 0.4))
  } else {
    NULL
  }

  out <- list(
    mu.B0 = mu_B0_lk, sigma.B0 = sigma_B0_lk,
    mu.B1 = mu_B1_lk, sigma.B1 = sigma_B1_lk,
    B1 = mu_B1_lk,
    Beta = mu_Beta_lk,
    rho = boot::inv.logit(
      rnorm(n_strata, params$logit_mu_rho_l, params$logit_sigma_rho_l)
    ) * 2 - 1,
    y.n = rep(trials_upper, N)
  )

  # Error term(s).
  add_terms <- if (var_level == "stratum") {
    sigma_y_lk <- rgamma(
      n_strata,
      params$mu_sigma_y_l^2 / params$sigma_sigma_y_l^2,
      params$mu_sigma_y_l / params$sigma_sigma_y_l^2
    )
    list(sigma = sigma_y_lk)
  } else if (var_level == "site") {
    mu_sigma_y_lk <- rgamma(
      n_strata,
      params$mu_mu_sigma_y_l^2 / params$sigma_mu_sigma_y_l^2,
      params$mu_mu_sigma_y_l / params$sigma_mu_sigma_y_l^2
    )
    sigma_sigma_y_lk <- rgamma(
      n_strata,
      params$mu_sigma_sigma_y_l^2 / params$sigma_sigma_sigma_y_l^2,
      params$mu_sigma_sigma_y_l / params$sigma_sigma_sigma_y_l^2
    )
    list(
      mu.sigma = mu_sigma_y_lk, sigma.sigma = sigma_sigma_y_lk,
      y.q = ifelse(var_type == "hier", 1, 0)
    )
  } else {
    NULL
  }
  out <- c(out, zi_terms, add_terms)

  if (!is.null(which_params)) {
    return(out[names(out) %in% which_params])
  }

  out
} # get_generating_values(2, 'exponential')

.h <- function(link) {
  function(x) {
    if (grepl("lin", link)) exp(x) else x
  }
}

get_hash <- function(args, shorten = TRUE) {
  args_hash <- digest(args, algo = "sha256")
  if (shorten) {
    return(substr(args_hash, 1, 8))
  }
  args_hash
} # get_hash('thing')

get_link <- function(x) {
  str_extract(x, "lin|exp|inv-logit")
} # get_link("assets/count-variables/poisson-exponential-b0/model.jags")

# # Moment match for the beta distribution.
# .mm_beta <- function(mu, sigma) {
#   alpha <- (mu^2 - mu^3 - mu * sigma^2) / sigma^2
#   beta <- (mu - 2 * mu^2 + mu^3 - sigma^2 + mu * sigma^2) / sigma^2
#   c(alpha = alpha, beta = beta)
# }
# mm_beta <- function(mu, sigma) {
#   t(do.call(mapply, c(FUN = .mm_beta, list(mu = mu, sigma = sigma))))
# }

# https://schmidtynotes.com/r/sample%20design/2019/11/19/spsurvey_grts.html
get_design <- function(n_sites_per_year, n_years, n_sites_total, n_strata) {
  n_years_full_cycle <- n_sites_total / n_sites_per_year
  panel_n <- rep(n_sites_per_year, n_years_full_cycle)
  names(panel_n) <- sprintf(
    "panel-%s",
    str_pad(seq_len(n_years_full_cycle), 2, "left", "0")
  )
  stratum_ids <- LETTERS[1:n_strata]
  out <- lapply(stratum_ids, function(x) {
    list(panel = panel_n, seltype = "Equal")
  })
  names(out) <- stratum_ids
  attr(out, "panel_info") <- tibble(
    x = seq(0, n_years - 1),
    panel = rep(names(panel_n), length.out = n_years)
  )
  out
}
get_theme <- function(x, ...) {
  switch(x,
    "sim_theme" = function() sim_theme(...)
  )
}
sim_theme <- function(y_lab_angle = 0) {
  theme_ipsum_rc(
    # base_family = font_an,
    grid = "Y", base_size = 16, axis_title_size = 18, strip_text_size = 18
  ) +
    theme(
      axis.title.y = element_text(angle = y_lab_angle),
      legend.position = "bottom"
    )
}
int_breaks <- function(x, n = 5, x_min = 0) {
  l <- pretty(x, n)
  l[abs(l %% 1) < .Machine$double.eps^0.5]
  l[l >= x_min]
}

source_dir <- function(path, verbose = TRUE, ...) {
  if (verbose) message("Sourcing the following files...")
  for (nm in list.files(path,
    pattern = "\\.[Rr]$", full.names = TRUE,
    recursive = TRUE
  )) {
    if (verbose) cat(nm)
    source(nm, ...) # was: file.path(path, nm)
    if (verbose) cat("\n")
  }
}

get_exit_code <- function(is_converged, is_verifiable) {
  message("Checking exit status...")

  # assign non-zero exit status for a failure to converge (11) or a failure
  # to pass verification (12)

  if (!is_converged) {
    11
  } else if (is_converged & !is_verifiable) {
    12
  } else {
    0
  }
}
