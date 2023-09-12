library(tidyverse)
library(yaml)


get_analysis_scenarios <- function(yaml_info) {
  covariates_vec <- if (!is.null(yaml_info$`additional covariates`)) {
    c(NA, yaml_info$`additional covariates`)
  } else {
    NA
  }

  state_variable_info <- yaml_info$`response info`$`state variable`
  svs <- sapply(state_variable_info, length) # state variable size / length

  if (svs[grep("response|hits column", names(svs))] != svs["description"]) {
    stop("Check that the response column and description entries in your YAML file are 1:1!")
  }
  state_variable_descs <- as_tibble(state_variable_info, .name_repair = "minimal")
  names(state_variable_descs) %<>% str_replace(" column", "")

  responses <- state_variable_info$`response column`
  hits <- trials <- NA
  if (is.null(responses)) {
    responses <- NA
    hits <- state_variable_info$`hits column`
    trials <- state_variable_info$`trials column`
  }

  row_expand <- function(x) {
    gle_zeros <- yaml_info$`group-level effects (zeros)`
    gle_zeros_arg <- if (is.null(gle_zeros)) {
      "none"
    } else {
      gle_zeros
    }
    mod_grid <- expand.grid(list(
      likelihood = yaml_info$likelihood,
      additional_covariates = covariates_vec,
      deterministic_model = yaml_info$`deterministic model`,
      group_level_effects = yaml_info$`group-level effects`,
      group_level_effects_zeros = gle_zeros_arg
    ), stringsAsFactors = FALSE) %>%
      mutate(deterministic_model = ifelse(
        likelihood == "gen-pois",
        "zero-trick", deterministic_model
      )) %>%
      distinct()
    bind_cols(
      x %>% as_tibble(.name_repair = "minimal") %>%
        slice(rep(1:n(), each = nrow(mod_grid))),
      mod_grid
    )
  }

  tibble(response = responses, hits = hits, trials = trials) %>%
    rowwise() %>%
    do(row_expand(.)) %>%
    ungroup() %>%
    {
      `if`(
        all(is.na(.$response)),
        select(., -response), select(., -hits, -trials)
      )
    } %>%
    left_join(state_variable_descs)
}


get_config <- function(yaml_file, rm_dqs = FALSE) {
  yaml_info <- yaml.load_file(yaml_file)

  yaml_tree <- dirname(yaml_file)
  yaml_branches <- NULL # basename(yaml_tree)

  basename(yaml_tree) != "config"

  while (!grepl("config$", basename(yaml_tree)) & basename(yaml_tree) != "_config") {
    this_yaml_file <- list.files(yaml_tree, "^\\_.*[^DEPREC].yml", full.names = TRUE)
    if (!identical(this_yaml_file, character(0))) {
      yaml_info %<>% c(yaml.load_file(this_yaml_file))
      yaml_branches %<>% c(basename(yaml_tree))
    }
    yaml_tree <- dirname(yaml_tree)
    # print(yaml_tree)
  }

  analysis_scenarios <- yaml_info %>%
    get_analysis_scenarios() %>%
    mutate(inits_handling = yaml_info$`initial values`)

  # prepare variance-type data

  resp_df <- analysis_scenarios %>%
    select(matches("^response$|^hits$")) %>%
    distinct()
  get_variances_types <- function(data) {
    if (is.null(data$`response column(s)`)) {
      return(NULL)
    } else {
      vt_df <- data %>%
        as_tibble(.name_repair = "minimal") %>%
        select(level, type, `response column(s)`) %>%
        setNames(c("var_level", "var_type", names(resp_df)))
      return(vt_df)
    }
  }

  variance_type_data <- do.call(bind_rows, lapply(yaml_info$variances, get_variances_types))
  if (nrow(variance_type_data) > 0) {
    analysis_scenarios %<>% left_join(variance_type_data)
  } else {
    # stop('Fix me!')
    analysis_scenarios %<>% mutate(var_level = "stratum", var_type = "fixed")
  }
  parameterization <- ifelse(is.null(yaml_info$`parameterization info`),
    "centered", yaml_info$`parameterization info`
  )

  analysis_scenarios %<>%
    # Substitions for the deterministic model.
    mutate( # deterministic_model = ifelse(likelihood == 'gamma',
      #                            'restricted-linear', deterministic_model),
      deterministic_model = ifelse(likelihood == "lognormal" & deterministic_model == "linear",
        "restricted-linear", deterministic_model
      )
    ) %>%
    mutate(
      var_level = ifelse(likelihood == "poisson" |
        likelihood == "gen-pois" |
        likelihood == "negative-binomial-simple" | # not quite true
        likelihood == "binomial" |
        likelihood == "zero-inflated-binomial",
      "site", var_level
      ),
      var_type = ifelse(likelihood == "poisson" |
        likelihood == "gen-pois" |
        likelihood == "negative-binomial-simple" |
        likelihood == "binomial" |
        likelihood == "zero-inflated-binomial",
      "hier", var_type
      ),
      dqs = ifelse(rm_dqs, "none", "all"),
      dqs = ifelse(!rm_dqs & !is.null(yaml_info$`time effect`), "no-trend", dqs),
      parameterization = parameterization,
      time_effect = TRUE,
      deflections = ifelse(grepl("\\(deflections\\)", additional_covariates),
        "yes", "no"
      )
    ) %>%
    distinct()
  # analysis_scenarios %>% distinct(likelihood, var_level, var_type)
  # analysis_scenarios %>% mutate()
  if (!is.null(yaml_info$`time effect`)) {
    analysis_scenarios %<>%
      filter(
        group_level_effects == "b0",
        grepl("none|g0", group_level_effects_zeros)
      ) %>%
      mutate(time_effect = FALSE)
  }

  # Dropouts for the beta-hurdle model.
  analysis_scenarios %<>%
    mutate(B_drop = "keepit", G_drop = "keepit")

  as_b <- as_g <- as_bg <- NULL
  if (!is.null(yaml_info$dropouts)) {
    d_out_nonzeros_config <- strsplit(yaml_info$`dropouts`, " ")
    d_out_nonzeros_vars <- unlist(lapply(d_out_nonzeros_config, `[[`, 1))
    fix_d_out_nonzeros <- unlist(
      lapply(d_out_nonzeros_config, function(x) grepl("fixed", x[2]))
    )
    as_b <- analysis_scenarios %>%
      select(-B_drop) %>%
      # left_join(tibble(response = yaml_info$dropouts, B_drop = TRUE),
      #           by = "response") %>%  # TODO, what happens if response == 'hits'?
      left_join(
        tibble(
          response = d_out_nonzeros_vars,
          B_drop = TRUE, fixed_B0 = fix_d_out_nonzeros
        ),
        by = "response"
      ) %>%
      # mutate(B_drop = !is.na(B_drop),
      #        B_drop = ifelse(B_drop & group_level_effects %in% c('b0', 'b0-b1'), '-int-only', 'keepit'),  # might be this
      #        group_level_effects = ifelse(B_drop == '-int-only' &
      #                                       group_level_effects == 'b0-b1',
      #                                     'b0', group_level_effects)) %>%
      mutate(
        B_drop = !is.na(B_drop),
        B_drop = ifelse(B_drop & group_level_effects %in% c("b0", "b0-b1"),
          "-int-only", "keepit"
        ),
        B_drop = ifelse(B_drop == "-int-only" & fixed_B0, "-int-only-fixed", B_drop),
        group_level_effects = ifelse(B_drop %in% c("-int-only", "-int-only-fixed") &
          group_level_effects == "b0-b1",
        "b0", group_level_effects
        )
      ) %>%
      select(-fixed_B0) %>%
      distinct()
  }

  if (!is.null(yaml_info$`dropouts (zeros)`)) {
    d_out_zeros_config <- strsplit(yaml_info$`dropouts (zeros)`, " ")
    d_out_zeros_vars <- unlist(lapply(d_out_zeros_config, `[[`, 1))
    fix_d_out_zeros <- unlist(
      lapply(d_out_zeros_config, function(x) grepl("fixed", x[2]))
    )

    get_g_drops <- function(data) {
      if (is.null(data)) {
        return(NULL)
      }
      data %>%
        select(-G_drop) %>%
        left_join(
          tibble(
            response = d_out_zeros_vars,
            G_drop = TRUE, fixed_G0 = fix_d_out_zeros
          ),
          by = "response"
        ) %>%
        mutate(
          G_drop = !is.na(G_drop),
          G_drop = ifelse(G_drop & group_level_effects_zeros %in% c("g0", "g0-g1"),
            "-int-only", "keepit"
          ),
          G_drop = ifelse(G_drop == "-int-only" & fixed_G0, "-int-only-fixed", G_drop),
          group_level_effects_zeros = ifelse(G_drop %in% c("-int-only", "-int-only-fixed") &
            group_level_effects_zeros == "g0-g1",
          "g0", group_level_effects_zeros
          )
        ) %>%
        select(-fixed_G0) %>%
        distinct()
    }

    as_g <- get_g_drops(analysis_scenarios)
    as_bg <- get_g_drops(as_b)
  }

  # bind_rows(analysis_scenarios, as_b, as_g) %>% distinct()
  fpi <- yaml_info$`finite population info`
  omit_pred_strat_block <-
    !is.null(fpi) & isFALSE(fpi$`covariate info includes unsampled sites`)
  analysis_scenarios <- bind_rows(analysis_scenarios, as_b, as_g, as_bg) %>% distinct()
  analysis_scenarios <- analysis_scenarios %>%
    mutate(
      dqs = ifelse(B_drop == "-int-only" & grepl("-int-only", G_drop), # was G_drop == '-int-only'
        "none", dqs
      ),
      additional_covariates = ifelse(B_drop == "-int-only" & grepl("-int-only", G_drop), # was G_drop == '-int-only',
        NA, additional_covariates
      )
    ) %>%
    distinct() %>%
    mutate(
      finite_pop = ifelse(isTRUE(yaml_info$`finite population correction`), "yes", "no"),
      omit_pred_strat_block =
        ifelse(finite_pop == "yes" & omit_pred_strat_block, "yes", "no"),
      add_rand_coefs = ifelse(!is.null(yaml_info$`other random coefficients`), "yes", "no")
    )
  # analysis_scenarios %>% distinct(group_level_effects, group_level_effects_zeros, B_drop, G_drop)


  yaml_info$yaml_file <- yaml_file %>%
    str_replace(".yml", "") %>%
    str_replace("config", "output")

  list(
    yaml_info = yaml_info, yaml_branches = yaml_branches,
    analysis_scenarios = analysis_scenarios
  )
}

get_jags_outputs_dir <- function(a_scenario, yaml_info) {
  param_desc <-
    if (!is.null(yaml_info$`parameterization info`)) "_nc" else NULL
  # return(parameterization)

  rand_Beta <- yaml_info$`other random coefficients`
  rand_Beta_desc <- if (!is.null(rand_Beta)) {
    paste("_rand", paste(rand_Beta, collapse = "-and-"), "effects", sep = "-")
  } else {
    NULL
  }

  a_scenario %>%
    rename(outcome = matches("^response$|^hits$")) %>%
    mutate(
      var_desc = paste(var_type, var_level, sep = "-"), # 'var'
      # Describe the droupouts specification.
      dropouts = sprintf(
        "%s%s%s%s%s%s",
        ifelse(B_drop != "keepit" | G_drop != "keepit", "_", ""),
        ifelse(B_drop == "keepit", "", "b"),
        ifelse(grepl("fixed", B_drop), "fix", ""),
        ifelse(G_drop == "keepit", "", "g"),
        ifelse(grepl("fixed", G_drop), "fix", ""),
        ifelse(B_drop != "keepit" | G_drop != "keepit", "-drop", "")
      ),
      group_level_effects =
        ifelse(group_level_effects_zeros == "none",
          group_level_effects,
          paste(group_level_effects, group_level_effects_zeros, sep = "_")
        ),
      output_dir =
        file.path(
          file.path(yaml_info$yaml_file, outcome), #
          paste0(
            paste(
              abbrev_likelihood(likelihood),
              abbrev_det_fun(deterministic_model),
              group_level_effects, var_desc,
              sep = "_"
            ),
            param_desc, dropouts, rand_Beta_desc
          ),
          get_additional_covs_desc(a_scenario)
        )
    ) %>%
    pull(output_dir)
}

abbrev_likelihood <- function(x) {
  switch(x,
    "ordinal-latent-normal" = "ord-lat-norm",
    "hurdle-ordinal-latent-beta" = "hurd-ord-lat-beta",
    "binomial" = "binom",
    "zero-inflated-binomial" = "zi-binom",
    "beta-binomial" = "od-binom", # overdispersed binomial
    "zero-inflated-beta-binomial" = "zi-od-binom",
    "poisson" = "pois",
    "lognormal" = "lnorm",
    "negative-binomial" = "neg-binom",
    "negative-binomial-simple" = "nb-vanilla",
    x
  )
} # abbrev_likelihood("poisson")

abbrev_det_fun <- function(x) {
  switch(x,
    "linear" = "lin",
    "restricted-linear" = "lin",
    "exponential" = "exp",
    "inverse-logit" = "inv-logit",
    "monomolecular" = "saturating",
    x
  )
}
