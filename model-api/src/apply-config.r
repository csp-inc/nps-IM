source("model-api/src/stats.r")
source("model-api/src/zone-objects.R")


get_tibble <- function(a_scenario, yaml_info, yaml_branches, output_dir, ...) {
  response_info <- yaml_info$`response info`
  site_loc_info <- yaml_info$`site location info`

  standardize_data(
    file = response_info$file,
    network_code = rev(yaml_branches)[1],
    unit_code = rev(yaml_branches)[2],
    unit_code_col = yaml_info$`unit code column`,
    site_col = yaml_info$`site id column`,
    sample_cols = response_info$`sample id column(s)`,
    response_col = c(a_scenario[["response"]], a_scenario[["hits"]]),
    trials_col = a_scenario[["trials"]],
    sampling_method = response_info$`sampling method`,
    event_date_info = yaml_info$`event date info`,
    # Site location info.
    site_locations_file = site_loc_info$file,
    coord_cols = site_loc_info$`coordinate columns`,
    # Strata info.
    stratum_col = yaml_info$`stratum id column`,
    stratum_area_info = yaml_info$`stratum area`,
    # Covariates.
    covariates_info = yaml_info$`covariate info`,
    covariates = a_scenario[["additional_covariates"]],
    # Censoring and truncation.
    censoring_info = yaml_info$`censoring info`,
    truncation_info = yaml_info$`truncation info`,
    # Parameterization.
    # parameterization_info = yaml_info$`parameterization info`,
    # Logging.
    log_dir = output_dir,
    trend_conditions = yaml_info$`trend conditions`,
    check_response_scale = a_scenario$likelihood == "beta",
    make_finite_pop_corr = isTRUE(yaml_info$`finite population correction`),
    finite_pop_info = yaml_info$`finite population info`,
    eval_drivers = yaml_info$`drivers`,
    me_crossings = yaml_info$`crossings`,
    me_omissions = yaml_info$`omissions`,
    # Offsets (as applicable, for count models)
    exposure_var = yaml_info$`response info`$`state variable`$`offset column`,
    ...
  )
}

default_inits <- function(n_site, n_strata, n_obs, likelihood, var_level,
                          n_chains = 3) {
  B.init <- array(dim = c(max(n_site), 2, n_strata))

  for (k in 1:n_strata) {
    B.init[1:n_site[k], 1, k] <- rep(2, n_site[k]) # -2
    B.init[1:n_site[k], 2, k] <- rep(5, n_site[k]) # -.4
  }

  inits <- list()
  inits <- inits[1:n_chains] # 1 for each chain
  lapply(
    seq_len(n_chains),
    function(x,
             chain_b_init_multiplier = c(1, 1.1, .9),
             chain_mu_B0 = c(3, 2, 1.2),
             chain_mu_B1 = c(.03, .02, .1),
             chain_rho = c(0, -.25, .25),
             chain_sigma_B0 = c(.1, .2, .5),
             chain_sigma_B1 = c(.01, .3, .6),
             chain_sigma = c(.01, .03, .06),
             # chain_beta_sigma = c(-3, -2, -1),
             chain_mu_beta_sigma = c(-5, -4, -3),
             chain_sigma_beta_sigma = c(1, 0.5, 0.25),
             # chain_mu_beta_sigma = c(-5, -4, -3),
             # chain_sigma_beta_sigma = c(.01, .03, .006),
             # chain_sigma_ss = c(.01, .03, .006),  # need Tom's input
             chain_p_zero = c(.6, .1, .8),
             # chain_beta = c(55, 100, 87),  # TODO: ask Tom about this....
             chain_z = ifelse(grepl("zero-inflated-beta-binomial|zero-inflated-binomial", likelihood), 0, 1),
             y_latent = c(.5, .3, .1),
             rng_name = "base::Mersenne-Twister",
             rng_seed = c(1, 245, 22)) {
      sigma_ss <- array(dim = c(max(n_site), n_strata))
      for (k in 1:n_strata) {
        sigma_ss[1:n_site[k], k] <- chain_sigma[x]
      }

      if (var_level == "site") {
        beta_sigma <- array(0, dim = c(max(n_site), n_strata))
        # for(k in 1:n_strata) {
        # beta_sigma[1:n_site[k], k] <- 0
        # }
      } else {
        beta_sigma <- rep(0, n_strata)
      }


      inits[[x]]$B <- B.init * chain_b_init_multiplier[x]
      inits[[x]]$mu.B0 <- rep(chain_mu_B0[x], n_strata)
      inits[[x]]$mu.B1 <- rep(chain_mu_B1[x], n_strata)
      inits[[x]]$rho <- rep(chain_rho[x], n_strata)
      inits[[x]]$sigma.B0 <- rep(chain_sigma_B0[x], n_strata)
      inits[[x]]$sigma.B1 <- rep(chain_sigma_B1[x], n_strata)
      inits[[x]]$sigma <- rep(chain_sigma[x], n_strata)
      inits[[x]]$sigma.ss <- sigma_ss
      inits[[x]]$p.zero <- rep(chain_p_zero[x], n_strata)
      # inits[[x]]$beta <- rep(chain_beta[x], n_strata)
      inits[[x]]$z <- rep(chain_z, n_obs)
      inits[[x]]$y.latent <- rep(y_latent[x], n_obs)
      # inits[[x]]$beta.sigma <- rep(chain_beta_sigma[x], n_strata)
      inits[[x]]$beta.sigma <- beta_sigma
      inits[[x]]$mu.beta.sigma <- rep(chain_mu_beta_sigma[x], n_strata)
      inits[[x]]$sigma.beta.sigma <- rep(chain_sigma_beta_sigma[x], n_strata)
      inits[[x]]$.RNG.name <- rng_name
      inits[[x]]$.RNG.seed <- rng_seed[x]

      inits[[x]]
    }
  )
}

get_stratum_weights <- function(data, yaml_info) {
  stratum_index_lookup <- data %>%
    select(matches("^stratum")) %>%
    distinct()
  d_sw <- yaml_info$`stratum area info` %>%
    as_tibble(.name_repair = "minimal") %>%
    gather(stratum_id, area_sq_m) %>%
    left_join(stratum_index_lookup) %>%
    filter(!is.na(stratum_index)) %>%
    arrange(stratum_index) %>%
    mutate(wt = area_sq_m / sum(area_sq_m))

  out <- d_sw %>% pull(wt)
  attr(out, "park_attr_lookup") <- d_sw %>%
    mutate(is_in_park_attr_yml = TRUE) %>%
    select(stratum_id, stratum_index, is_in_park_attr_yml)
  out
  # browser()
}

get_class_limits <- function(cover_class_info, rescale = TRUE) {
  read_csv(cover_class_info$file) %>%
    {
      `if`(
        rescale,
        mutate_at(
          ., vars(
            cover_class_info$class_midpoint_col,
            cover_class_info$class_break_low,
            cover_class_info$class_break_high
          ),
          function(x) x / 100
        ), .
      )
    } %>%
    # mutate_at(vars(Midpoint, Low, High), function(x) x/100) %>%
    # mutate(Width=High-Low) %>%
    rename(
      response = !!cover_class_info$class_col,
      class_midpoint = !!cover_class_info$class_midpoint_col,
      class_break_low = !!cover_class_info$class_break_low,
      class_break_high = !!cover_class_info$class_break_high
    ) %>%
    arrange(class_break_low) %>%
    select(response, class_midpoint, class_break_low, class_break_high)
}

get_data_list <- function(data, stratum_weights, L,
                          folder = NULL,
                          variances_structure = NULL,
                          cover_class_info = NULL,
                          time_effect = NULL,
                          other_rand_coefs = NULL,
                          finite_pop = NULL,
                          exposure_var = NULL,
                          hat_step = 0.2) { # by default, Xs vary by site

  if (nrow(data) == 0) {
    stop("The data supplied to get_data_list() appears to be empty, please take a look!")
  }

  # ---- objects inherited via attributes ----
  y.n.unobs.site <- attr(data, "y.n.unobs.site")
  X.pred <- attr(data, "X.pred")
  y.zone.index <- attr(data, "y.zone.index")
  eval_mean_for_tv_covariates <- attr(data, "eval_mean_for_tv_covariates")
  in.sample.idx <- attr(data, "in.sample.idx")
  n.x.pred <- attr(data, "n.x.pred")
  x.pred.index <- attr(data, "x.pred.index")
  x.pred.raw <- attr(data, "x.pred.raw")
  j.pred <- attr(data, "j.pred")
  k.pred <- attr(data, "k.pred")
  X.driver <- attr(data, "X.driver")
  which.drivers <- attr(data, "which.drivers")
  # ------------------------------------------


  x_hat_raw <- seq(0, max(data$rel_year), hat_step)
  strat_table <- data %>%
    group_by(stratum_index) %>%
    arrange(stratum_index) %>%
    summarise(n_site = n_distinct(site_id))
  n_site <- strat_table$n_site
  n_strata <- n_distinct(strat_table$stratum_index)
  # attr(data, 'x.hat.raw') <- seq(0, max(data$rel_year), 0.2)
  # data %T>%
  # pipe_assign('x.hat.raw', seq(0, pull(., rel_year) %>% max, .2)) %T>%
  # pipe_assign('strat.table',
  #             group_by(., stratum_index) %>%
  #               arrange(stratum_index) %>%  # new
  #               summarize(n.site=n_distinct(site_id))) %T>%
  # pipe_assign('n.site', strat_table$n.site) %T>% # pull(strat.table, n.site)
  # pipe_assign('n.strata', n_distinct(strat_table$stratum_index)) # pull(strat.table, stratum_index) %>% length

  # hat.site.array <- array(NA, dim=c(max(n.site), length(x_hat_raw), n.strata))

  if (L == "hurdle-ordinal-latent-beta") {
    classes_limits <- get_class_limits(cover_class_info) # sprintf("%.100f",classes_limits$class_break_high[12])

    d_hurdle <- data %>%
      mutate(bern = ifelse(response == 0, 1, 0)) %>%
      left_join(classes_limits)
    gt_0 <- d_hurdle %>% filter(response > 0)

    thresh <- classes_limits %>% pull(class_break_low)
    lim_dint <- c(thresh[2:length(thresh)], 0.99999999)
    thresh_alt <- classes_limits %>% pull(class_break_high)
    names(thresh_alt) <- lead(classes_limits$class_break_low)
    trace_alt <- classes_limits$class_break_low[2]
    lim_dint_alt <- c(trace_alt - ifelse(trace_alt == 1e-08, 1e-09, 1e-08), thresh_alt[2:length(thresh)]) # TODO: validate!
    upper_cutpoints <- classes_limits %>%
      slice(2:(n() - 1)) %>%
      pull(class_break_high)

    nYlevels <- length(upper_cutpoints) + 1
    # hat.site.pr.array <- array(NA, dim=c(dim(hat.site.array), nYlevels))

    L_varying <- list(
      # Data for beta latent process (categories > 0).
      y.beta = gt_0$response,
      # Data for presence absence (0 or 1).
      y.bern = d_hurdle$bern,
      y.bern.site = d_hurdle$site_in_stratum_index,
      y.bern.strata = d_hurdle$stratum_index,
      # All data including 0 category.
      x.raw = gt_0 %>% pull(rel_year),
      x = gt_0 %>% pull(rel_year) %>% scale() %>% as.double(),
      # Categories omitting 0 cover class.
      lim = upper_cutpoints,
      lim.dint = lim_dint_alt, # thresh[2:length(thresh)],
      nYlevels = nYlevels, # number of categories
      # Scale limits from 0 to 1 for computation of likelihood.
      w = d_hurdle %>% pull(rel_year) %>% scale() %>% as.double(),
      y.site = gt_0$site_in_stratum_index,
      y.strata = gt_0$stratum_index # ,
      # hat.site.pr = hat.site.pr.array,
      # hat.site.class.mean = hat.site.array,
      # hat.site.class.new.obs = hat.site.array,
      # hat.site.class0.new.obs = hat.site.array,
      # hat.strat.mean.01 = matrix(NA, nrow=length(x_hat_raw), ncol=n.strata),
      # hat.strat.mean.gt0 = matrix(NA, nrow=length(x_hat_raw), ncol=n.strata)
    )
  } else {
    L_varying <- list(
      y = c(data[["response"]], data[["hits"]]),
      x.raw = data %>% pull(rel_year),
      x = data %>% pull(rel_year) %>% scale() %>% as.double(),
      y.site = data %>% pull(site_in_stratum_index),
      y.strata = data %>% pull(stratum_index)
    )
    if (L == "ordinal-latent-normal") {
      nYlevels <- 6 # the number of soil stability categories
      theta <- rep(NA, nYlevels - 1)
      theta[1] <- 1 + .5
      theta[nYlevels - 1] <- nYlevels - .5
      L_varying$nYlevels <- nYlevels
      L_varying$theta <- theta
    }
    if (L == "binomial" | L == "beta-binomial" |
      L == "zero-inflated-binomial" | L == "zero-inflated-beta-binomial") {
      trials.mean <- ceiling(mean(data$trials))

      trials.hat.array <- array(trials.mean,
        dim = c(max(n_site), length(x_hat_raw), n_strata)
      )

      if (isTRUE(!is.null(finite_pop) & finite_pop)) {
        trials.hat.array <- array(
          trials.mean,
          dim = c(n_distinct(j.pred), length(x_hat_raw), n_strata)
        )
      }

      L_varying$trials.hat <- trials.hat.array
    }
    if (L == "gen-pois") {
      L_varying$Zeros <- rep(0, length(L_varying$y))
    }
    if (L == "lognormal") {
      get_sd_log_y <- function(x) {
        sigma <- sd(x, na.rm = TRUE)
        mu <- mean(x, na.rm = TRUE)
        sqrt(log((sigma^2 + mu^2) / mu^2))
      }

      sd_log_y_by_stratum <- data %>%
        group_by(stratum_id) %>%
        # summarise(sd_log_y = get_sd_log_y(response[!is_censored]))
        summarise(sd_log_y = get_sd_log_y(na.omit(response)))
      L_varying$sd.log.y <- sd_log_y_by_stratum$sd_log_y

      # L_varying$mean.log.y <- mean(log(data$response), na.rm = TRUE)
    }
    if ("is_censored" %in% names(data)) {
      L_varying$is.censored <- as.integer(data$is_censored)
      # L_varying$y.rep <- ifelse(data$is_censored, data$censor_limit_vec, NA)
      L_varying$is.censored.rep <- L_varying$is.censored # or rep(1, nrow(data)), shot in the dark

      L_varying$censor.limit.vec <- data$censor_limit_vec
      L_varying$censor.limit.vec.rep <- data$censor_limit_vec_rep # testing
      L_varying$censor.limit.vec[!data$is_censored] <- Inf
    }
    if ("truncated_below" %in% names(data)) {
      if (all(!is.na(data$truncated_below))) L_varying$trunc.lower <- data$truncated_below
      if (all(!is.na(data$truncated_above))) L_varying$trunc.upper <- data$truncated_above
      # if('is_censored' %in% names(data)) {
      #  L_varying$censor.limit.vec.rep <- cbind(data$truncated_below, data$censor_limit_vec)
      # }
    }
  }

  L_invariant <- list(
    y.n.site = n_site,
    y.n.strata = n_strata,
    x.hat.raw = x_hat_raw,
    x.hat = (x_hat_raw - mean(L_varying$x.raw)) / sd(L_varying$x.raw),
    # hat.site.mean = hat.site.array,
    # hat.site.new.obs = hat.site.array,
    # hat.strat.mean = matrix(NA, nrow=length(x_hat_raw), ncol=n.strata),
    # hat.park.mean = rep(NA, length(x_hat_raw)),
    wt = stratum_weights,
    y.q = ifelse(variances_structure$var_type == "fixed", 0, 1)
  )
  if (variances_structure$var_level == "stratum" |
    L == "poisson" |
    L == "binomial" |
    L == "zero-inflated-binomial") {
    L_invariant <- L_invariant[names(L_invariant) != "y.q"]
  }

  d_for_tv_covariates <- NULL # data for time-varying covariates

  if (eval_mean_for_tv_covariates) {
    # pred.site.array <-
    #   array(NA, dim=c(max(n.site), n.x.pred, n.strata))  # was n_distinct(L_varying$x)

    X.pred.lookup <- array(NA, dim = c(n_distinct(j.pred), n.x.pred, n_distinct(k.pred)))
    idx_vs_X_df <- cbind(j.pred, k.pred, x.pred.index, X.pred) %>% # i.pred
      as_tibble() %>%
      mutate(row.index = 1:n())

    for (k in 1:n_distinct(k.pred)) {
      this_idx_vs_X_df <- filter(idx_vs_X_df, k.pred == k) %>%
        select(j.pred, x.pred.index, row.index) %>%
        pivot_wider(names_from = x.pred.index, values_from = row.index) %>%
        # spread(x.pred.index, row.index) %>%
        select(-j.pred) %>%
        as.matrix()
      X.pred.lookup[1:nrow(this_idx_vs_X_df), 1:ncol(this_idx_vs_X_df), k] <- this_idx_vs_X_df
    }



    d_for_tv_covariates <- list(
      n.x.pred = n.x.pred,
      x.pred.index = x.pred.index,
      x.pred = (x.pred.raw - mean(L_varying$x.raw)) / sd(L_varying$x.raw),
      x.pred.raw = x.pred.raw,
      j.pred = j.pred,
      k.pred = k.pred,
      X.pred = X.pred,
      X.pred.lookup = X.pred.lookup # ,
      # pred.site.mean = pred.site.array,
      # pred.site.new.obs = pred.site.array,
      # pred.strat.mean = matrix(NA, nrow=n.x.pred, ncol=n.strata),
      # pred.park.mean = rep(NA, n.x.pred)
    )
    # d_for_tv_covariates <- c(d_for_tv_covariates, list(X.driver = X.driver))


    # if(L == 'hurdle-ordinal-latent-beta') {
    #   pred.site.pr.array <- array(NA, dim=c(dim(pred.site.array), nYlevels))
    #   d_for_tv_covariates <- c(
    #     d_for_tv_covariates,
    #     list(pred.site.pr = pred.site.pr.array,
    #          pred.site.class.mean = pred.site.array,
    #          pred.site.class.new.obs = pred.site.array,
    #          pred.site.class0.new.obs = pred.site.array,
    #          pred.strat.mean.01 = matrix(NA, nrow=n.x.pred, ncol=n.strata),
    #          pred.strat.mean.gt0 = matrix(NA, nrow=n.x.pred, ncol=n.strata))
    #   )
    # }

    if (L == "binomial" | L == "beta-binomial" |
      L == "zero-inflated-binomial" | L == "zero-inflated-beta-binomial") {
      trials.pred.array <- array(NA, dim = c(max(n_site), n.x.pred, n_strata)) # pred.site.array
      if (isTRUE(!is.null(finite_pop) & finite_pop)) {
        trials.pred.array <- array(NA, dim = c(n_distinct(j.pred), n.x.pred, n_strata))
      }
      trials.pred.array[is.na(trials.pred.array)] <- trials.mean
      d_for_tv_covariates$trials.pred <- trials.pred.array
    }
  }

  deflections_matrix <- NULL
  if (file.exists(file.path(folder, "00-input/deflections-indices.rds"))) {
    deflections_list <- readRDS(file.path(folder, "00-input/deflections-indices.rds"))
    deflections_matrices <- lapply(deflections_list, function(x) {
      matrix(c(length(x), min(x, na.rm = TRUE), max(x, na.rm = TRUE) + 1), nrow = 1)
    })
    deflections_matrix <- list(
      y.n.deflections = length(deflections_list),
      defl.indices = do.call(rbind, deflections_matrices)
    )
  }

  jags_data <- c(
    L_varying[names(L_varying) != "x.raw"],
    L_invariant,
    d_for_tv_covariates,
    list(X.driver = X.driver, which.drivers = which.drivers),
    deflections_matrix
  )


  if (!is.null(data[["hits"]])) {
    jags_data$y.n <- data[["trials"]]
  }
  if (any(grepl("^fe_", names(data)))) {
    jags_data$X <- as.matrix(data %>% select(matches("^fe_")))
    if (L == "hurdle-ordinal-latent-beta") {
      jags_data$X <- as.matrix(gt_0 %>% select(matches("^fe_")))
      jags_data$W <- as.matrix(d_hurdle %>% select(matches("^fe_")))
    }
  }

  if (!is.null(exposure_var)) {
    jags_data$exposure <- data[[exposure_var]]
  }

  if (!is.null(time_effect)) {
    jags_data$B1 <- rep(0, jags_data$y.n.strata)
  }

  if (!is.null(other_rand_coefs)) {
    which_rand_coefs <-
      grepl(paste(other_rand_coefs, collapse = "|"), dimnames(jags_data$X)[[2]])
    jags_data$is.rand.slope <- as.integer(which_rand_coefs)
    n_rand_slopes <- sum(which_rand_coefs) + 1 # this will need to be updated for b0-b1
    jags_data$n.rand.slopes <- n_rand_slopes
    # https://www.mrc-bsu.cam.ac.uk/wp-content/uploads/2019/06/handouts_master.pdf
    # http://mmeredith.net/blog/2020/Correlated_priors.htm
    jags_data$R.wish <- diag(n_rand_slopes) # diag(c(5, 2, 2, 2))
    jags_data$k.wish <- n_rand_slopes + 1 # k = M + 1, where M is the number of variables results in a Uniform(-1, 1) prior for the correlation coefficient
    # if (is.null(time_effect)) browser()
  }

  zone.data <- NULL

  if (isTRUE(finite_pop)) {
    jags_data$tot.site <- jags_data$y.n.site + y.n.unobs.site

    if (!is.null(y.zone.index)) {
      zone.data <- get_zone_objects(
        y.zone.index, x.pred.index, n.x.pred,
        folder, X.pred
      ) # j.pred.zone
    }
    if ("X.pred" %in% names(jags_data)) { # was exists('X.pred')
      jags_data$in.sample.idx <- in.sample.idx
      all.idx <- seq.int(ncol(jags_data$X.pred))
    }
    if ("defl.indices" %in% names(jags_data)) {
      defl.idx <- seq.int(
        jags_data$defl.indices[1, 2],
        jags_data$defl.indices[1, 2] +
          (jags_data$defl.indices[1, 1] - 1) - 1
      )
      jags_data$non.defl.idx <- all.idx[-defl.idx]
    }
  }

  if (!is.null(folder)) {
    save_object(c(jags_data, zone.data), file.path(folder, "00-input"), "jags-data.rds")
  }

  Filter(Negate(is.null), c(jags_data, zone.data))
}

fence <- function(vec, lb = 0.0001, ub = 0.9999) pmax(lb, pmin(vec, ub))

get_jags_inits <- function(type, n_chains, data, likelihood, var_level,
                           output_dir, inits_cache, n_rand_slopes = NULL, ...) {
  inits_cache <- file.path(inits_cache, "inits")

  rngs <- c(
    "base::Wichmann-Hill",
    "base::Marsaglia-Multicarry",
    "base::Super-Duper",
    "base::Mersenne-Twister"
  )

  empty_inits <- lapply(sample(rngs, 3), function(x) list(.RNG.name = x))

  # Inheritance for zero-inflated models.
  chain_z <- ifelse(
    grepl("zero-inflated-beta-binomial|zero-inflated-binomial", likelihood),
    0, 1
  )
  if (chain_z == 0) {
    empty_inits <- lapply(empty_inits, function(x) {
      x$z <- rep(chain_z, nrow(data))
      x
    })
  }

  jags_inits <-
    if (type == "inherited") {
      inits_cache <- file.path(
        inits_cache,
        paste0(digest::digest(output_dir), ".rds")
      )
      if (file.exists(inits_cache)) {
        readRDS(inits_cache)
      } else {
        empty_inits
      }
    } else if (type == "none") {
      empty_inits
    } else if (type == "default") {
      # needs: jags_data, jags_file, likelihood

      di <- get_default_inits(likelihood = likelihood, var_level = var_level, ...)
      if (likelihood == "hurdle-ordinal-latent-beta") {
        source("model-api/src/beta-hurdle-inits.R")
        lapply(di, function(x) {
          tryme <- get_inits(runif(1, 0.1, 0.9), runif(1, 0.1, 0.2)) # list("mu" = tryme['mu'], "sigma" = tryme['sigma'])
          x$mu.B0 <- rep(boot::logit(tryme["mu"]), length(x$mu.B0))
          x$sigma.B0 <- rep(0.1, length(x$mu.B0))
          x$mu.beta.sigma <- rep(boot::logit(tryme["sigma"]), length(x$mu.beta.sigma))
          x$sigma.beta.sigma <- rep(0.1, length(x$mu.beta.sigma))
          x
        })
      } else {
        di
      }
      # empty_inits
    } else if (grepl("*.rds$", type)) {
      readRDS(type) # modifyList(, empty_inits)
    } else if (type == "empirical" & grepl("beta-binomial", likelihood)) { #

      d_inits_raw <- data %>%
        mutate(p_obs = hits / trials) %>%
        group_by(site_in_stratum_index, stratum_index) %>%
        summarise(
          mean_p = fence(mean(p_obs)), sd_p = fence(sd(p_obs)),
          .groups = "drop"
        )
      d_inits_b0 <- d_inits_raw %>%
        pivot_wider(site_in_stratum_index,
          names_from = "stratum_index", values_from = "mean_p"
        ) %>%
        select(-site_in_stratum_index) # %>%
      # mutate(across(everything(), ~ replace_na(.x, 0.0001)))
      d_inits_sigma_ss <- d_inits_raw %>%
        pivot_wider(site_in_stratum_index,
          names_from = "stratum_index", values_from = "sd_p"
        ) %>%
        select(-site_in_stratum_index) # %>%
      # mutate(across(everything(), ~ replace_na(.x, 0.0001)))

      sd_shrink <- 1
      d_inits <- d_inits_raw %>%
        group_by(stratum_index) %>%
        summarise(
          mean_b0 = mean(boot::logit(mean_p)), sd_b0 = sd(boot::logit(mean_p)) * sd_shrink,
          mu_sigma = mean(sd_p), sigma_sigma = sd(sd_p)
        )

      lapply(1:3, function(x) {
        eps <- 0.0001
        mu.B0 <- rnorm(nrow(d_inits), d_inits$mean_b0, eps)
        sigma.B0 <- fence(
          runif(nrow(d_inits), d_inits$sd_b0 - eps, d_inits$sd_b0 + eps)
        )

        sigma.ss <- apply(d_inits_sigma_ss, c(1, 2), function(x, mult = 2) {
          fence(runif(1, x - eps, x + eps))
        })
        mu.sigma <- fence(
          runif(nrow(d_inits), d_inits$mu_sigma - eps, d_inits$mu_sigma + eps)
        )
        sigma.sigma <- fence(
          runif(nrow(d_inits), d_inits$sigma_sigma - eps, d_inits$sigma_sigma + eps)
        )
        B0 <- apply(d_inits_b0, c(1, 2), function(x) {
          rnorm(1, boot::logit(x), eps)
        })

        these_inits <- list(
          mu.B0 = mu.B0, sigma.B0 = sigma.B0,
          # B0 = B0,
          # sigma.ss = sigma.ss,
          mu.sigma = mu.sigma, sigma.sigma = sigma.sigma, sigma = mu.sigma
        )

        if (grepl("zero-inflated", likelihood)) {
          c(these_inits, list(z = rep(0, length(data$hits))))
        } else {
          these_inits
        }
      })
    } else if (type == "empirical" & likelihood == "hurdle-ordinal-latent-beta") {
      # new_inits <- NULL
      message("developing inits using empirical summaries of the data....")
      d_inits_raw <- data %>%
        group_by(site_in_stratum_index, stratum_index) %>%
        summarise(
          p_of_0 = fence(sum(class_midpoint == 0) / n()),
          p_gt_0 = fence(max(0, mean(class_midpoint[class_midpoint > 0]), na.rm = TRUE)),
          sd_gt_0 = fence(max(0, sd(class_midpoint[class_midpoint > 0]), na.rm = TRUE),
            lb = 0.00001, ub = 0.99999
          ),
          .groups = "drop"
        )

      d_inits_b0 <- d_inits_raw %>%
        pivot_wider(site_in_stratum_index,
          names_from = "stratum_index", values_from = "p_gt_0"
        ) %>%
        select(-site_in_stratum_index)
      d_inits_sigma_ss <- d_inits_raw %>%
        pivot_wider(site_in_stratum_index,
          names_from = "stratum_index", values_from = "sd_gt_0"
        ) %>%
        select(-site_in_stratum_index)

      sd_shrink <- 0.5
      d_inits <- d_inits_raw %>%
        group_by(stratum_index) %>%
        summarise(
          mean_g0 = mean(boot::logit(p_of_0)), sd_g0 = sd(boot::logit(p_of_0)) * sd_shrink,
          mean_b0 = mean(boot::logit(p_gt_0)), sd_b0 = sd(boot::logit(p_gt_0)) * sd_shrink,
          # mu_sigma = mean(sd_gt_0), sigma_sigma = sd(sd_gt_0) # coming in far too low!?
          mean_beta_sigma = mean(boot::logit(sd_gt_0)), sd_beta_sigma = sd(boot::logit(sd_gt_0)) * sd_shrink
        )

      lapply(1:3, function(x) {
        eps <- 0.0001
        mu.G0 <- rnorm(nrow(d_inits), d_inits$mean_g0, eps)
        sigma.G0 <- fence(
          runif(nrow(d_inits), d_inits$sd_g0 - eps, d_inits$sd_g0 + eps)
        )
        mu.B0 <- rnorm(nrow(d_inits), d_inits$mean_b0, eps)
        sigma.B0 <- fence(
          runif(nrow(d_inits), d_inits$sd_b0 - eps, d_inits$sd_b0 + eps)
        )
        # sigma.ss <- apply(d_inits_sigma_ss, c(1, 2), function(x, mult = 2) {
        #   fence(runif(1, x * mult - eps, x * mult + eps))
        # })
        # mu.sigma <- fence(
        #   runif(nrow(d_inits), d_inits$mu_sigma - eps, d_inits$mu_sigma + eps)
        # )
        # sigma.sigma <- fence(
        #   runif(nrow(d_inits), d_inits$sigma_sigma - eps, d_inits$sigma_sigma + eps)
        # )
        mu.beta.sigma <- rnorm(nrow(d_inits), d_inits$mean_beta_sigma, eps)
        sigma.beta.sigma <- fence(
          runif(nrow(d_inits), d_inits$sd_beta_sigma - eps, d_inits$sd_beta_sigma + eps)
        )
        beta.sigma <- if (var_level == "stratum") {
          rep(0, nrow(d_inits))
        } else {
          apply(d_inits_sigma_ss, c(1, 2), function(x) {
            0 # rnorm(1, boot::logit(x), eps)
          })
        }
        B0 <- apply(d_inits_b0, c(1, 2), function(x) {
          rnorm(1, boot::logit(x), eps)
        })
        mu.G1 <- rnorm(nrow(d_inits), 0, eps * 10)
        sigma.G1 <- runif(nrow(d_inits), 0.005, 0.01)

        list(
          mu.G0 = mu.G0, sigma.G0 = sigma.G0, mu.B0 = mu.B0, sigma.B0 = sigma.B0,
          B0 = B0,
          # beta.sigma = beta.sigma,
          mu.beta.sigma = mu.beta.sigma, sigma.beta.sigma = sigma.beta.sigma,
          # sigma.ss = sigma.ss,
          # mu.sigma = mu.sigma, sigma.sigma = sigma.sigma, sigma = mu.sigma,
          mu.G1 = mu.G1, sigma.G1 = sigma.G1
        )
      })
    }




  if ("is_censored" %in% names(data)) {
    jags_inits_with_y_inits <- lapply(jags_inits, function(x) {
      y_init <- rep(NA, length(data$response))
      y_init[data$is_censored] <- data$censor_limit_vec[data$is_censored] +
        runif(sum(data$is_censored))
      # y_init_rep <- data$response
      # y_init_rep[data$is_censored] <- data$censor_limit_vec[data$is_censored] -
      #  runif(sum(data$is_censored)) # potentially unsafe
      sd_log_y <- sd(log(data$response), na.rm = TRUE)

      c(x, list(
        y = y_init, y.rep = y_init
        # sigma.log.y = runif(n_distinct(data$stratum_id),
        #                     0.001*sd_log_y, 10*sd_log_y)
      ))
    })
    return(jags_inits_with_y_inits)
  }

  if (likelihood == "gen-pois") {
    return(lapply(jags_inits, function(x) {
      x$delta <- 0
      x
    }))
  }

  if (!is.null(n_rand_slopes)) {
    return(lapply(jags_inits, function(x) {
      K <- n_distinct(data$stratum_id)
      x$xi <- array(NA, dim = c(n_rand_slopes, K))
      x$mu.Beta.raw <- array(NA, dim = c(n_rand_slopes, K))
      x$Tau.Beta.raw <- array(NA, dim = c(n_rand_slopes, n_rand_slopes, K))
      for (k in 1:n_distinct(data$stratum_id)) {
        x$xi[, k] <- runif(n_rand_slopes)
        x$mu.Beta.raw[, k] <- rnorm(n_rand_slopes)
        x$Tau.Beta.raw[, , k] <- MCMCpack::rwish(n_rand_slopes + 1, diag(n_rand_slopes))
      }
      x
    }))
  }


  return(jags_inits)
}

get_default_inits <- function(jags_data, jags_file, likelihood, var_level) {
  jags_inits <- default_inits(
    n_site = jags_data$y.n.site,
    n_strata = jags_data$y.n.strata,
    n_obs = ifelse(likelihood == "hurdle-ordinal-latent-beta",
      length(jags_data$y.beta),
      length(jags_data$y)
    ),
    likelihood = likelihood, var_level = var_level
  )


  params_needed <- system(
    paste("model-api/src/get-unknowns.sh", jags_file),
    intern = TRUE
  )
  params_needed <- grep("^y.*$|^B$|^sigma.ss$|^#", params_needed, value = TRUE, invert = TRUE)

  # inits_needed <- sapply(names(jags_inits[[1]]),
  #                        function(x) {
  #                          pattern <- paste0('\\b', gsub('\\.', '\\\\.', x), '\\b[^\\.]')
  #                          any(grepl(pattern, readLines(jags_file)))
  #                        })
  inits_needed <- names(jags_inits[[1]]) %in% params_needed
  lapply(jags_inits, function(x) x[which(inits_needed)])
}

get_additional_covs_desc <- function(a_scenario) {
  if (is.na(a_scenario$additional_covariates)) {
    "null"
  } else {
    paste(tolower(str_split(a_scenario$additional_covariates, ", ")[[1]]), collapse = "_") %>%
      gsub(" ", "-", .) %>%
      gsub("\\(", "", .) %>%
      gsub("\\)", "", .) %>%
      gsub("\\*", "X", .)
  }
}

invoke_os_command <- function(cmd, capture_stdout = FALSE) {
  if (Sys.info()["sysname"] == "Linux") {
    system(cmd, intern = capture_stdout)
  } else if (Sys.info()["sysname"] == "Windows") {
    tryCatch(
      {
        shell(cmd, '"C:/Program Files/Git/git-bash.exe"', flag = "")
      },
      warning = function(w) {
        git_bash_loc <-
          system('WHERE /R "C:\\Users" git-bash.exe', intern = TRUE)
        shell(cmd, shQuote(git_bash_loc), flag = "")
      },
      error = function(e) {
        git_bash_loc <-
          system('WHERE /R "C:\\Users" git-bash.exe', intern = TRUE)
        shell(cmd, shQuote(git_bash_loc), flag = "")
      }
    )
  }
}

get_truncation_string <- function(data, log = FALSE) {
  paste(ifelse(all(is.na(data$truncated_below)), "",
    ifelse(log, "log(trunc.lower[n])", "trunc.lower[n]")
  ),
  ifelse(all(is.na(data$truncated_above)), "",
    ifelse(log, "log(trunc.upper[n])", "trunc.upper[n]")
  ),
  sep = ","
  )
}

get_jags_model_file <- function(a_scenario, output_dir, verbose = FALSE) {

  # Use the model-builder API to create individual JAGS model files for each
  # scenario defined in the top-level config file.



  build_mod_leave_crumbs <- function(data) {
    invoke_os_command(data$shell_command)

    jags_file <- readLines(file.path(output_dir, "crumbs.txt"))
    file.remove(file.path(output_dir, "crumbs.txt"))
    invoke_os_command(paste("model-api/src/get-watch-list.sh", jags_file))
    vars_to_watch <- readLines(file.path(output_dir, "00-input/vars-to-watch.txt"))

    data %>%
      mutate(jags_file,
        all_vars = list(sub("* $", "", vars_to_watch))
      )
  }

  pipe_message <- function(data, column, verbose) {
    if (verbose) message("Building model with command:\n", pull(data, !!enquo(column)))
    data
  }
  a_scenario %>%
    mutate(shell_command = paste(
      "model-api/src/model-builder/compile-jags-file.sh",
      likelihood,
      deterministic_model,
      group_level_effects,
      get_additional_covs_desc(a_scenario),
      output_dir,
      ifelse(attr(a_scenario, "eval_mean_for_tv_covariates"),
        "keep_pred_switches", "drop_pred_switches"
      ),
      group_level_effects_zeros,
      var_level, var_type, dqs,
      censoring, truncation, parameterization,
      deflections, B_drop, G_drop,
      finite_pop, omit_pred_strat_block,
      get_marg_effects, add_rand_coefs, apply_offset
    )) %>%
    pipe_message(shell_command, verbose = verbose) %>%
    do(build_mod_leave_crumbs(.))
}

wheat_from_chaff <- function(x, jags_data) {
  is_data <- function(x) !all(is.na(x))
  jags_data_is_data <- lapply(jags_data, is_data)

  chaff <- names(jags_data[which(unlist(jags_data_is_data))])
  chaff <- c(
    chaff, "p.site", "site.wt", "i.pred", "Tau.B", "Sigma.B",
    "hat.site.kappa", "hat.site.p.nb",
    "hat.strat.kappa", "hat.strat.p.nb", "p.nb" # "kappa" required for DIC
  )
  wheat <- c(
    x[!x %in% grep("[^theta]", chaff, value = TRUE)],
    "deviance"
  )

  # str_detect(c('#p', 'p'), '^\\#*')
  deflections <- NULL
  if ("defl.indices" %in% names(jags_data)) deflections <- "D"
  c(wheat[!grepl("^\\#[a-zA-Z]*", wheat)], deflections)
}
