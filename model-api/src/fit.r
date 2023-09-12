rowwise_fit <- function(a_scenario, yaml_info, yaml_branches, n_chains = 3,
                        save_full_jags_obj = FALSE,
                        destroy_xy = FALSE, pass_errors = FALSE,
                        max_adapt_attempts = 3) {
  source("model-api/src/get-config.r")
  source("model-api/src/apply-config.r")
  source("model-api/src/theme.r")
  source("model-api/src/serializers.r")
  sapply(list.files("model-api/src", "utils-.*.R", full.names = TRUE), source)
  sapply(list.files("model-api/src", "checking-.*.R", full.names = TRUE), source)
  sapply(list.files("model-api/src", "inference-.*.R", full.names = TRUE), source)
  source("model-api/src/fit-utils.r")
  source("model-api/src/model-eval.R")
  rjags::load.module("dic")
  rjags::load.module("glm")

  print(a_scenario %>% as.list())
  if (destroy_xy) yaml_info$`site location info` <- NULL


  output_dir <- a_scenario %>% get_jags_outputs_dir(yaml_info)
  if (!is.null(yaml_info$rerun)) {
    if (!output_dir %in% yaml_info$rerun) {
      return(a_scenario %>% mutate(status = "skipping"))
    }
  }
  # output_subdirs <-
  #   file.path(
  #     output_dir,
  #     c(
  #       "00-input", #' 01-diagnostics', '01-diagnostics/funnel', #'01-diagnostics/traceplots'
  #       "02-checking", "02-checking/prior", "02-checking/ppc",
  #       "02-checking/ppc/stratum_id", "02-checking/ppc/unit_code",
  #       "02-checking/variography",
  #       "03-inference", "03-inference/strata", "03-inference/park", "03-inference/zone",
  #       "03-inference/site", "03-inference/me", "03-inference/coef",
  #       "03-inference/strata/all", "03-inference/strata/sampled",
  #       "03-inference/park/all", "03-inference/park/sampled",
  #       "03-inference/zone/all", "03-inference/zone/sampled",
  #       "99-misc"
  #     )
  #   )
  # create_tree <- function(x) dir.create(x, showWarnings = FALSE, recursive = TRUE)
  # sapply(output_subdirs, create_tree)

  d <- a_scenario %>% get_tibble(yaml_info, yaml_branches, output_dir)

  jags_info <- a_scenario %>%
    mutate(
      censoring = ifelse("is_censored" %in% names(d), "yes", "no"),
      truncation = ifelse("truncated_below" %in% names(d),
        get_truncation_string(d), # (a_scenario[['likelihood']] == 'lognormal')
        "no"
      ),
      get_marg_effects = ifelse(is.null(yaml_info$drivers) |
        isFALSE(yaml_info$drivers), "no", "yes"),
      apply_offset = ifelse(
        is.null(yaml_info$`response info`$`state variable`$`offset column`),
        "no", "yes"
      )
    ) %>%
    set_attr("eval_mean_for_tv_covariates", attr(d, "eval_mean_for_tv_covariates")) %>%
    get_jags_model_file(output_dir, verbose = TRUE)


  if (dirname(jags_info$jags_file) != output_dir) warning("Check output_dir")

  jags_n_iters <- c("n_adapt" = n_adapt, "n_update" = n_update, "n_iter" = n_iter) %T>%
    save_object(file.path(output_dir, "00-input"), "jags-n-iters.rds")

  resp_vec <- d %>%
    select(matches("^response$|^hits$")) %>%
    pull(1)

  if (any(is.na(resp_vec)) & (!"is_censored" %in% names(d))) {
    # the root of dana's problem!?!?
    warning("NAs exist in response column")
    d %<>% drop_na(any_of(c("response", "hits")))
    # d %<>% na.omit
  }

  this_L <- a_scenario[["likelihood"]]
  if (this_L == "hurdle-ordinal-latent-beta") {
    # Add class midpoint cover for future reference....
    d %<>% left_join(get_class_limits(yaml_info$cover_class_info))
    save_table(d, file.path(output_dir, "00-input"), "state-variable-data.csv")
  }

  trend_step_size <- yaml_info$`trend step size`

  stratum_weights <- get_stratum_weights(d, yaml_info)
  # browser()
  jags_data <- d %>%
    get_data_list(
      stratum_weights = stratum_weights,
      L = jags_info$likelihood, folder = output_dir,
      variances_structure = a_scenario %>% select(matches("^var_*")),
      cover_class_info = yaml_info$cover_class_info,
      time_effect = yaml_info$`time effect`,
      other_rand_coefs = yaml_info$`other random coefficients`,
      finite_pop = yaml_info$`finite population correction`,
      exposure_var = yaml_info$`response info`$`state variable`$`offset column`,
      hat_step = ifelse(is.null(trend_step_size), 0.2, trend_step_size)
    )

  inits_cache <-
    str_extract(output_dir, sprintf(".*(?=/%s)", tail(yaml_branches, 1)))

  jags_inits <- get_jags_inits(
    type = a_scenario$inits_handling,
    n_chains = n_chains,
    data = d,
    likelihood = this_L,
    var_level = a_scenario$var_level,
    output_dir = output_dir,
    inits_cache = inits_cache,
    n_rand_slopes = jags_data$n.rand.slopes,
    # NEW
    jags_data = jags_data, jags_file = jags_info$jags_file
  )
  if (a_scenario$add_rand_coefs == "yes") {
    jags_inits <- lapply(jags_inits, function(x) x[-grep("mu.B0|Beta.grp.tilde", names(x))])
  }
  # if(a_scenario$inits_handling == 'default') {
  # jags_inits <- get_default_inits(jags_data, jags_info$jags_file, this_L)
  # }

  # Copy calling script and metadata doc to analysis results folder.
  sapply(
    c("calling-script.R", "output-metadata.md"),
    function(x) file.copy(file.path("model-api/docs", x), output_dir)
  )
  capture.output(summary(d),
    file = file.path(output_dir, "00-input", "state-variable-data-summary.txt")
  )
  save_object(jags_info, file.path(output_dir, "00-input"), "jags-info.rds")

  # JAGS model object.
  save_object(jags_inits, file.path(output_dir, "00-input"), "jags-inits.rds")

  is_tuned <- FALSE
  attempts <- 0
  start_time <- Sys.time()

  # if (TRUE) {
  #   source('sandbox/richness-identifiability-sim.R')
  #   jags_data$y <- do_sim(jags_data, d, output_dir)
  # }

  if (grepl("Gap", unlist(a_scenario %>% select(matches("response|hits"))))) {
    jags_data$y.t.index <- as.integer(as.factor(jags_data$x))
    jags_data$y.thresh <- 50
    jags_data$y.n.t <- tibble(y.strata = jags_data$y.strata, x = jags_data$x) %>%
      group_by(y.strata) %>%
      summarise(y.n.t = n_distinct(x)) %>%
      arrange(y.strata) %>%
      pull(y.n.t)
    jags_data$y.tot.trans.length <- d %>%
      group_by(stratum_index, rel_year) %>%
      summarise(
        y.n.transects = n_distinct(interaction(site_id, sample_id)),
        y.tot.trans.length = y.n.transects * (5000 - 40)
      ) %>%
      ungroup() %>%
      select(stratum_index, rel_year, y.tot.trans.length) %>%
      spread(rel_year, y.tot.trans.length) %>%
      arrange(stratum_index) %>%
      select(-stratum_index) %>%
      as.matrix()
    ttttt <- lapply(1:jags_data$y.n.strata, function(x) {
      out <- rep(NA, max(jags_data$y.t.index))
      in_years <- sort(unique(jags_data$y.t.index[jags_data$y.strata == x]))
      out[1:length(in_years)] <- in_years
      out
    })
    jags_data$spill <- array(0, dim = c(
      length(jags_data$y),
      max(jags_data$y.t.index),
      jags_data$y.n.strata
    ))
    for (i in 1:length(jags_data$y)) {
      for (t in 1:max(jags_data$y.t.index)) {
        for (k in 1:jags_data$y.n.strata) {
          jags_data$spill[i, t, k] <- ifelse(
            jags_data$y.t.index[i] == t & jags_data$y.strata[i] == k,
            NA, 0
          )
        }
      }
    }
    jags_data$y.rep.gap.gt.thresh <- jags_data$spill
    jags_data$y.years.in.sample <- do.call(rbind, ttttt)
  }


  # n_strata <- jags_data$y.n.strata  # TODO: consider whether something like this might help models converge faster
  # jags_inits=lapply(jags_inits, function(x) {
  #   empirical_B0 <- d %>%
  #     group_by(stratum_index) %>%
  #     summarise(logit_B0 = boot::logit(mean(hits/trials))) %>%
  #     pull(logit_B0)
  #   c(x,
  #     list(mu.B0 = rnorm(n_strata, empirical_B0, 0.01) # sqrt(1 / 0.4444444)
  #          ))  #sigma.B0 = runif(n_strata, 0, 0.01))
  #
  # })

  jags_vars <- jags_info$all_vars[[1]] %>% wheat_from_chaff(jags_data)
  if (this_L == "hurdle-ordinal-latent-beta") {
    jags_vars <- jags_vars[!jags_vars %in% c(
      "mu", "a", "b", "deltaF.0",
      "pllike.0", "tau.gamma", "tau.B0",
      "mean.bern", "sd.bern", "mean.sim.bern", "sd.sim.bern",
      "sd.y", "mean.y", "sd.sim", "mean.sim",
      "hat.site.pr", "pred.site.pr",
      jags_vars[grepl("*.bb$", jags_vars)]
    )]
  }
  if (this_L == "ordinal-latent-normal") {
    jags_vars <- jags_vars[!jags_vars %in% c(
      "pr",
      jags_vars[grepl("tau", jags_vars)],
      jags_vars[grepl("(strat|site).mu$", jags_vars)],
      jags_vars[grepl("*.(strat|site).(pr|pr.oos)$", jags_vars)]
    )]
  }
  if (this_L == "gen-pois") {
    jags_vars <- jags_vars[!jags_vars %in% c(
      "L", paste0("l", 1:4), "eta",
      "C", "Zeros.mean"
    )]
  }
  save_object(jags_vars, file.path(output_dir, "00-input"), "jags-vars.rds")

  # jags_data$X.pred.lookup[is.na(jags_data$X.pred.lookup)] <- 1
  while (!is_tuned & attempts <= max_adapt_attempts) {
    attempts <- attempts + 1
    if (attempts <= (max_adapt_attempts - 1)) {
      # First, try default inits....
      jags_model <- tryCatch(
        {
          # start_time <- Sys.time()
          rjags::jags.model(
            file = jags_info$jags_file,
            data = jags_data,
            inits = jags_inits,
            n.chains = n_chains, n.adapt = n_adapt
          ) # n_adapt
          # end_time <- Sys.time()
          # end_time - start_time
        },
        error = function(e) {
          # browser()
          file.remove(file.path(output_dir, "adapt-phase-errors.txt"))
          save_stdout(e, output_dir, "adapt-phase-errors.txt")
          if (!pass_errors) browser()
          NA
          # sink(file.path(output_dir, 'adapt-phase-errors.txt')); print(e); sink()
        }
      )
    } else {
      # If inits auto-generated by JAGS fail after (max_adapt_attempts - 1) attempts, omit 'em.
      # sink(file.path(output_dir, 'adapt-phase-warnings.txt'))
      # print('Trying old, likely crappy default inits...')
      # sink()
      save_stdout(
        "Trying old, likely crappy default inits...",
        output_dir, "adapt-phase-warnings.txt"
      )
      jags_inits <- get_default_inits(jags_data, jags_info$jags_file, this_L,
        var_level = a_scenario$var_level
      )
      jags_model <- tryCatch(
        {
          rjags::jags.model(
            file = jags_info$jags_file,
            data = jags_data, # c(jags_data, zone.data),
            inits = jags_inits,
            n.chains = n_chains, n.adapt = n_adapt
          )
        },
        error = function(e) {
          save_stdout(e, output_dir, "adapt-phase-errors.txt")
          # sink(file.path(output_dir, 'adapt-phase-errors.txt')); print(e); sink()
          NA
        }
      )
    }
    # do.call(cbind, jags_data[c('y', 'is.censored', 'censor.limit.vec', 'trunc.lower')]) %>% head
    # do.call(cbind, jags_data[c('censor.limit.vec', 'trunc.lower')]) %>% as_tibble() %>%
    #   mutate(test = censor.limit.vec >= trunc.lower) %>%
    #   filter(!test)
    is_tuned <- !is.na(jags_model[1])
    if (attempts == max_adapt_attempts) {
      return(a_scenario %>% mutate(status = "failed, reached max attempts"))
    }
  }
  if (pass_errors & file.exists(file.path(output_dir, "adapt-phase-errors.txt"))) {
    return(
      a_scenario %>% mutate(status = "failed, see adapt-phase-errors.txt")
    )
  }
  save_object(jags_inits, file.path(output_dir, "00-input"), "jags-inits.rds")

  update(object = jags_model, n.iter = n_update)

  # JAGS samples.
  y_cens <- if ("is_censored" %in% names(d)) "y" else NULL

  thin <- 1

  do_thin <- this_L == "hurdle-ordinal-latent-beta" |
    any(grepl("GapSize", c(a_scenario$response, a_scenario$hits))) |
    d$unit_code[1] == "TSRA" | d$unit_code[1] == "GRASS"
  if (do_thin) thin <- 5


  y_rep_pattern <- ifelse(do_thin, "y.rep*", "999*")

  z_jags <- rjags::jags.samples(
    model = jags_model,
    variable.names = c(jags_vars[!grepl(y_rep_pattern, jags_vars)], "lambda", "delta"),
    # 'k.for.j.draw.zone', 'p.site.zone', 'j.draw.zone'),  # temporary
    # variable.names=jags_vars[!grepl(y_rep_pattern, jags_vars)],
    thin = thin,
    n.iter = n_iter
  )


  if (do_thin & jags_info$likelihood != "gen-pois") { # was : this_L=='hurdle-ordinal-latent-beta'
    z_jags_y_rep <- rjags::jags.samples(
      model = jags_model,
      variable.names = c(jags_vars[grepl("y.rep*", jags_vars)], y_cens),
      thin = thin,
      n.iter = n_iter * .5
    )
    z_jags <- c(z_jags, z_jags_y_rep)
  }

  # Coda samples.
  get_coda_vars <- function(x, K = jags_data$y.n.strata) {
    length(dim(x)) <= 3 & dim(x)[1] <= K
  }
  z_jags_keepers <- unlist(lapply(z_jags, get_coda_vars))
  coda_vars <-
    intersect(
      names(z_jags_keepers)[which(z_jags_keepers)],
      c(
        "B", "mu.B0", "mu.B1", "B0", "B1", "sigma.B0", "sigma.B1", "Beta", "Beta.hat",
        "G", "mu.G0", "mu.G1", "G1", "sigma.G0", "sigma.G1", "Gamma", # "G0.fixed",
        "p.zero", "p.mean", "p.sd", "deviance",
        "sigma", "mu.sigma", "sigma.sigma",
        "mu.beta.sigma", "sigma.beta.sigma", "beta.sigma"
      )
    ) # jags_vars[!grepl('^hat.*|^pred.*|^trend.*', jags_vars)]
  coda_vars <- c(coda_vars, if ("Beta" %in% names(z_jags)) "Beta" else NULL)

  if ("Beta.hat" %in% names(z_jags)) {
    coda_vars <- c(grep("^Beta$", coda_vars, value = TRUE, invert = TRUE), "Beta.hat")
  }

  save_object(coda_vars, file.path(output_dir, "00-input"), "coda-vars.rds")

  z_coda <- rjags::coda.samples(
    model = jags_model,
    variable.names = coda_vars,
    thin = thin,
    n.iter = n_iter
  )

  # sink(file.path(output_dir, '03-inference', 'coda-samples-summary.txt'))
  # print(summary(z_coda))
  # sink()
  save_stdout(
    summary(z_coda), file.path(output_dir, "03-inference"),
    "coda-samples-summary.txt"
  )

  coda_samples_quantiles <- as_tibble(summary(z_coda)$quantiles, rownames = "param") %>%
    filter(!param %in% c("p.mean", "p.sd", "deviance"))
  fe_names <- dimnames(jags_data$X)[[2]]
  if (!is.null(fe_names)) {
    # needs to work for Beta.hat[\\d+]!
    fe_lookup <- tibble(covariate = sub("^fe_", "", fe_names)) %>%
      mutate(param = sprintf(
        "%s[%s]",
        ifelse("Beta.hat" %in% names(z_jags),
          "Beta.hat", "Beta"
        ), 1:n()
      ))
    coda_samples_quantiles <- coda_samples_quantiles %>%
      left_join(fe_lookup, by = "param")
  }
  save_table(coda_samples_quantiles, file.path(output_dir, "03-inference"), "coda-samples-quantiles.csv")

  convergence_diagnostic <- coda::gelman.diag(z_coda, multivariate = FALSE)
  # sink(file.path(output_dir, '01-diagnostics', 'convergence-diagnostics.txt'))
  # print(convergence_diagnostic)
  # sink()
  save_stdout(
    convergence_diagnostic,
    file.path(output_dir, "01-diagnostics"),
    "convergence-diagnostics.txt"
  )
  (converged <- all(convergence_diagnostic$psrf[, 1] < 1.1))

  if (this_L == "gen-pois") { # TODO: wrap in function?
    source("model-api/src/rgp-hilbe.r")
    gen_pois_y_reps <- get_gen_pois_y_reps(z_jags)
    class(gen_pois_y_reps) <- class(z_jags$y.hat)
    z_jags$y.rep <- gen_pois_y_reps

    jags_step <- function(y.rep, y, fun = "mean") {
      if (fun == "mean") {
        as.numeric(mean(y.rep[]) - mean(y[]) >= 0)
      } else if (fun == "sd") {
        as.numeric(sd(y.rep[]) - sd(y[]) >= 0)
      }
    }
    gen_pois_p_sd <- apply(gen_pois_y_reps, 2:3, jags_step, y = jags_data$y, fun = "sd")
    # class(gen_pois_p_sd) <- class(z_jags$y.hat)
    gen_pois_p_mean <- apply(gen_pois_y_reps, 2:3, jags_step, y = jags_data$y, fun = "mean")
    # class(gen_pois_p_mean) <- class(z_jags$y.hat)
    z_jags$p.sd <- gen_pois_p_sd
    z_jags$p.mean <- gen_pois_p_mean

    # hist(apply(test[, , ], 2:3, mean)); abline(v = mean(jags_data$y), col = 'red')
    # hist(apply(test[, , ], 2:3, sd)); abline(v = sd(jags_data$y), col = 'red')
  }

  post_pred_loss <- z_jags %>% get_ppl(d, this_L)

  tryCatch(
    {
      DIC <- z_jags %>% get_dic(d, this_L, jags_data)
    },
    error = function(e) {
      browser()
    }
  )

  mod_summary <- z_jags %>%
    get_mod_summary(post_pred_loss, DIC, convergence_diagnostic, output_dir)

  end_time <- Sys.time()
  save_stdout(end_time - start_time, output_dir, "system-time.txt")

  # Write output to disk ------------------------------------------------------

  a_label <- a_scenario$description
  a_deterministic_model <- a_scenario$deterministic_model
  make_dq <- a_scenario$dqs == "all"

  # Save residuals, y reps, and other inferential results....
  source("model-api/src/ordinal-serializers.r")
  source("model-api/src/serializers.r")
  source("model-api/src/serializers-slim.r")
  source("model-api/src/theme.r")
  source("model-api/src/funnel.r")
  source("model-api/src/deflections.r")
  source("model-api/src/misc-plotting.R")
  source("model-api/src/zone-means.R")
  source("model-api/src/get-density.R")
  # source('model-api/src/purge/cnh-results.R')
  source("model-api/src/cdpr-results.R")

  d_raw <- NULL
  d_for_checking <- d
  if (!is.null(yaml_info$`ppc facets`)) {
    d_raw <- a_scenario %>%
      get_tibble(yaml_info, yaml_branches, output_dir, return_raws = TRUE)
    d_for_checking <- d_raw
  } # TODO: remove!
  # yaml_info$`ppc facets` <- c("site_id", "cal_year") #c(yaml_info$`ppc facets`, 'unit_code', 'stratum_id')  #    "unit_code"  "stratum_id"

  sapply(list.files("model-api/src", "utils-.*.R", full.names = TRUE), source)
  sapply(list.files("model-api/src", "inference-.*.R", full.names = TRUE), source)
  source("model-api/src/model-checking.R")
  source("model-api/src/diagnostics.R")

  # TODO: use this with z_jags: attr(d, 'intx_effects_list')!!!!!!!!!
  # print(make_dq) # cli: TRUE; gui: FALSE
  # print(eval_mean_for_tv_covariates) # cli: FALSE, gui: FALSE
  eval_mean_for_tv_covariates <- attr(d, "eval_mean_for_tv_covariates")
  # browser()
  z_jags %T>%
    {
      `if`(
        any(grepl("CDPR/TSRA", output_dir)),
        create_cdpr_summaries(., d, jags_data, d_raw, output_dir, a_label,
          interval = 0.80, method = "eti"
        ), .
      )
    } %T>%
    {
      `if`(
        a_scenario$dqs != "none" & any(grepl("MISC/LAVO|MISC/SEKI", output_dir)),
        get_cnh_maps(., d, output_dir, a_label), .
      )
    } %T>%
    # source('model-api/src/serializers-slim.r')
    # z_jags %T>%
    get_op_plots_new(d, output_dir) %T>%
    # z_jags %T>%
    get_me_inference(output_dir, a_label) %T>% # TODO: reinstate
    get_density(d, a_deterministic_model, this_L, output_dir) %T>%
    {
      `if`(
        d$network_code[1] == "MISC",
        get_misc_contrasts(., d, output_dir, a_label, "mu.D"), .
      )
    } %T>%
    {
      `if`(
        "D" %in% names(z_jags),
        get_deflection_coefs(
          ., d, output_dir, a_label,
          a_deterministic_model
        ), .
      )
    } %T>%
    get_funnel_plot(jags_data, output_dir) %T>%
    # get_op_plots(d, this_L, output_dir, mod_summary)
    get_random_hyperparams(d, output_dir) %T>%
    # z_jags %>%
    get_coef_post_dists(d, output_dir, a_label, a_deterministic_model,
      data_raw = if (!is.null(d_raw)) select_at(d_raw, vars(site_index, yaml_info$`ppc facets`))
    ) %T>%
    {
      `if`(
        grepl("ordinal", this_L) & make_dq,
        get_site_inference_ord(
          ., d, output_dir, a_label,
          jags_data$x.hat.raw, "hat", this_L
        ), .
      )
    } %T>%
    {
      `if`(
        grepl("ordinal", this_L) & eval_mean_for_tv_covariates & make_dq,
        get_site_inference_ord(
          ., d, output_dir, a_label,
          unique(jags_data$x.pred.raw), "pred", this_L
        ), .
      )
    } %T>%
    {
      `if`(
        grepl("hurdle-ordinal-latent-beta", this_L) & make_dq,
        get_stratum_inference_ord(., d, output_dir, a_label,
          jags_data$x.hat.raw, "hat", this_L,
          n_draws = 250
        ), .
      )
    } %T>%
    {
      `if`(
        grepl("hurdle-ordinal-latent-beta", this_L) & eval_mean_for_tv_covariates & make_dq,
        get_stratum_inference_ord(., d, output_dir, a_label,
          unique(jags_data$x.pred.raw), "pred", this_L,
          n_draws = 250
        ), .
      )
    } %T>%
    get_residuals(d, this_L, output_dir) %T>% # was d_for_ppl
    # source('model-api/src/model-checking.R')
    # z_jags %T>%
    get_y_reps(d_for_checking, this_L, output_dir, facet_by = yaml_info$`ppc facets`) %T>%
    {
      `if`(
        make_dq,
        get_trend_inference(., d, output_dir, a_label), .
      )
    } %T>%
    {
      `if`(
        make_dq,
        get_zone_trend_inference(., output_dir, a_label), .
      )
    } %T>%
    {
      `if`(
        make_dq,
        get_park_scale_inference(., d, jags_data, output_dir, a_label,
          n_draws = 500, seed = 123
        ), .
      )
    } %T>%
    {
      `if`(
        make_dq,
        get_site_inference(., d, output_dir, a_label, jags_data$x.hat.raw, "hat", this_L), .
      )
    } %T>%
    # z_jags %T>%
    {
      `if`(
        eval_mean_for_tv_covariates & make_dq,
        get_site_inference(
          ., d, output_dir, a_label,
          unique(jags_data$x.pred.raw), "pred", this_L
        ), .
      )
    } %T>%
    {
      `if`(
        n_distinct(d$stratum_id) >= 1 & make_dq,
        get_stratum_inference(., d, output_dir, a_label,
          jags_data$x.hat.raw, "hat", this_L,
          n_draws = 250
        ), .
      )
    } %T>%
    # z_jags %T>%
    {
      `if`(
        eval_mean_for_tv_covariates & make_dq,
        get_stratum_inference(., d, output_dir, a_label,
          unique(jags_data$x.pred.raw), "pred", this_L,
          n_draws = 250,
          x_hat_raw = jags_data$x.hat.raw
        ), .
      )
    } %T>%
    {
      `if`(
        eval_mean_for_tv_covariates & make_dq & ("hat.zone.mean" %in% names(z_jags)),
        get_zone_inference(., d, output_dir, a_label,
          unique(jags_data$x.pred.raw), this_L,
          n_draws = 250,
          x_hat_raw = jags_data$x.hat.raw
        ), .
      )
    } %T>%
    # New stratum-level results without those pesky CIs!
    {
      `if`(
        n_distinct(d$stratum_id) >= 1 & make_dq,
        get_stratum_inference_no_cis(., d, output_dir, a_label,
          jags_data$x.hat.raw, "hat", this_L,
          n_draws = min(n_iter, 500)
        ), .
      )
    } %T>%
    {
      `if`(
        eval_mean_for_tv_covariates & make_dq,
        get_stratum_inference_no_cis(., d, output_dir, a_label,
          unique(jags_data$x.pred.raw), "pred", this_L,
          n_draws = min(n_iter, 500),
          x_hat_raw = jags_data$x.hat.raw
        ), .
      )
    } %T>%
    # z_jags %T>%
    {
      `if`(
        save_full_jags_obj,
        save_object(., file.path(output_dir, "99-misc"), "z-jags.rds"), .
      )
    }


  # Traceplots.
  tryCatch(
    {
      get_param_trace(z_coda, output_dir)
    },
    error = function(e) {
      browser()
    }
  )

  z_jags %>% cache_inits(n_chains, output_dir,
    chain_z = jags_inits[[1]]$z[1],
    mod_summary$gelman_diag,
    inits_cache
  )

  return(a_scenario %>% mutate(status = "fitted"))
}
