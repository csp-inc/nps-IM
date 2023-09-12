library(tidyverse)
library(lazyeval)
library(lubridate)
library(glue)
source("model-api/src/wrangle.r")
source("model-api/src/logging.r")


index_from_category <- function(x, ...) {
  # Creates unique sequential indices from a categorical variable.
  #
  # Args:
  #   x:      The categorical variable for which a numeric index is desired.
  #   ...:    Arguments to pass to factor().
  #
  # Returns:
  #   Ordered, numeric indices.

  as.integer(factor(x, ...))
}

index_call <- function(x) interp(~ index_from_category(a), a = as.name(x))

indices_by_group <- function(data, group, in_col, out_col) {
  data %>%
    group_by_(.dots = group) %>%
    mutate_(.dots = setNames(list(index_call(in_col)), out_col)) %>%
    ungroup()
}

scale_x <- function(x) {
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}

rename_call <- function(x) interp(~a, a = as.name(x))

# rel_year_call <- function(x) interp(~ a - min(a), a = as.name(x))
rel_year_call <- function(x, year_0 = NULL) {

  # interp(~ a - max(c(min(a), year_0)), a = as.name(x))
  interp(~ a - max(c(min(a, na.rm = TRUE), year_0)), a = as.name(x))
}

right_cat <- function(x, r) {
  x %>%
    purrr::map(function(xx) paste0(xx, r)) %>%
    unlist()
}

format_dttm_attr <- function(data, cal_year_col, dttm_format, ref_yr = NULL) {
  # if(!is.null(ref_yr)) browser()

  data %>%
    mutate(
      year_from_datetime =
        year(parse_date_time(get(cal_year_col), dttm_format))
    ) %>%
    mutate_(.dots = setNames(list(rename_call("year_from_datetime")), "cal_year")) %>%
    mutate_(.dots = setNames(list(rel_year_call("year_from_datetime", ref_yr)), "rel_year")) %>%
    select(-year_from_datetime)
}

parse_covariates_info <- function(covs) {
  covariates_sans_parns <- str_replace_all(covs, " \\(.*?\\)", "")
  o_raw <- str_split(covariates_sans_parns, ", ")[[1]]
  # o2 <-
  o <- unique(unlist(str_split(o_raw, "\\*")))
  attr(o, "o_raw") <- o_raw
  o
}

prep_and_join_cov_info <- function(data, cov_info, covs, sample_cols = NULL) {
  cov_time_attr_info <- cov_info$`event date info`

  cov_data <- load_data(cov_info$file) %>%
    select(-matches("X1")) %>%
    {
      `if`(
        !is.null(cov_time_attr_info),
        format_dttm_attr(., cov_time_attr_info$`date-time column`,
          cov_time_attr_info$`date-time format`,
          ref_yr = min(data$cal_year)
        ),
        .
      )
    }

  shared_cols <- intersect(names(data), names(cov_data))
  if (is.null(cov_info$mapping)) {
    cov_info$mapping <- "site"
    sample_cols <- NULL
  }

  cov_data %<>%
    select(all_of(c(shared_cols, covs, sample_cols))) %>% # was: select_(.dots=c(shared_cols, covs)) %T>%
    {
      `if`(
        cov_info$mapping == "sample",
        unite_(., "sample_id", sample_cols, sep = "-", remove = FALSE), .
      )
    } # %T>%
  full_cov_data <- distinct(cov_data)

  suppressWarnings({
    out <- data %>%
      select(-any_of(covs)) %>%
      left_join(cov_data)
  })

  attr(out, "full_cov_data") <- full_cov_data
  out
}

get_design_mat <- function(data, covariates_string, log_dir,
                           center_and_scale_params = NULL, skip_check = FALSE) {
  if (!is.null(center_and_scale_params)) {
    if (nrow(center_and_scale_params) > 0) {
      row_count <- 0
      for (x in center_and_scale_params$fe_column) {
        row_count <- row_count + 1
        this_center <- center_and_scale_params$`scaled:center`[row_count]
        this_scale <- center_and_scale_params$`scaled:scale`[row_count]
        data %<>%
          mutate(this_scaled_x = (get(x) - this_center) / this_scale) %>%
          mutate_(.dots = setNames(list(rename_call("this_scaled_x")), x)) %>%
          select(-this_scaled_x)
      }
    }
  }

  covariates_vector <- parse_covariates_info(covariates_string)
  is.non.indicator.numeric <- function(x) {
    is.numeric(x) & !all(x %in% c(0, 1))
  }

  scale_atts_list <- data %>%
    select(all_of(covariates_vector)) %>%
    map_if(., is.non.indicator.numeric, function(x) {
      attributes(scale(x))[-1] %>%
        as_tibble(.name_repair = "minimal")
    })

  X_raw <- data %>%
    select(all_of(covariates_vector)) %>%
    # pipe_assign('scale_atts_list',
    #             map_if(., is.non.indicator.numeric, function(x) attributes(scale(x))[-1] %>%
    #                      as_tibble(.name_repair = 'minimal'))) %>%
    {
      `if`(
        is.null(center_and_scale_params),
        map_if(., is.non.indicator.numeric, function(x) as.numeric(scale(x))), .
      )
    } %>%
    # map_if(is.numeric, function(x) as.numeric(scale(x))) %>%
    map_if(is.character, function(x) as.factor(x)) %>%
    as_tibble(.name_repair = "minimal") %>%
    set_attr("scale_atts_list", scale_atts_list)
  X_raw_orig <- X_raw

  covariates_vector_raw <- attr(covariates_vector, "o_raw")
  intx_effects_list <- NULL
  if (any(grepl("\\*", covariates_vector_raw))) {
    intx_names <- NULL
    intx_exprs <- grep("\\*", covariates_vector_raw, value = TRUE)
    n_svb <- length(intx_exprs)
    m_k <- grep("\\*", covariates_vector_raw) - 1
    save_object(
      list(n.svb = n_svb, m.k = m_k),
      file.path(log_dir, "00-input"),
      "stratum-varying-effect-info.rds"
    )
    X_raw_list <- list()
    counter <- 0
    intx_effects_list <- list()
    for (intx_expr in intx_exprs) {
      counter <- counter + 1
      intx_name <- sub("\\*", "X", intx_expr)

      intx_terms <- unlist(strsplit(intx_expr, "\\*"))
      intx_term_is_fct <- sapply(intx_terms, function(x) !is.numeric(data[[x]]))

      if (any(intx_term_is_fct)) {
        if (all(intx_term_is_fct)) stop("Unsupported: two categorical variables as interactions.")

        x_cont <- intx_terms[!intx_term_is_fct]
        x_cat <- intx_terms[intx_term_is_fct]
        # d_contXcat <- bind_cols(X_raw[x_cont], data[x_cat])
        m_formula <- as.formula(sprintf("~%s+%s+%s*%s", x_cont, x_cat, x_cont, x_cat))

        # if (intx_expr == 'ndjfm*EcoSite') browser()
        m_contXcat_raw <- # ~deficit.pregr + MDCATY + deficit.pregr*MDCATY
          model.matrix(m_formula, X_raw_orig) # was: X_raw
        m_contXcat <- m_contXcat_raw[, grepl("\\:", dimnames(m_contXcat_raw)[[2]])]
        intx_name <- sub(":", "X", dimnames(m_contXcat)[[2]])
        dimnames(m_contXcat)[[2]] <- intx_name
        X_raw_list[[counter]] <- bind_cols(X_raw_orig %>% select(-all_of(x_cat)), as_tibble(m_contXcat))

        # ---- new! for computing the coefficients of terms involving interactions
        empty_list <- vector(mode = "list", length = length(intx_name) + 1)
        names(empty_list) <- levels(X_raw_orig$EcoSite)
        main_effect <- glue("fe_{x_cont}")
        empty_list[[levels(X_raw_orig$EcoSite)[1]]] <- as.character(main_effect)
        intx_terms <- glue("fe_{intx_name}")
        for (i in seq_along(intx_terms)) {
          empty_list[[i + 1]] <- c(main_effect, intx_terms[i])
        }
        intx_effects_list[[counter]] <- empty_list
        names(intx_effects_list)[counter] <- x_cont
      } else { # the usual interactions case
        X_raw_list[[counter]] <- X_raw_orig %>%
          mutate(!!intx_name := eval(parse(text = intx_expr))) # THIS LINE IS NEW!!
      }
      intx_names <- c(intx_names, intx_name)

      # ================
    }
    # Drop duplicated columns
    rm_duplicated_cols <- function(data) {
      as_tibble(data[!duplicated(as.list(data))])
    }
    X_raw <- rm_duplicated_cols(do.call(cbind, X_raw_list))

    X_raw <- X_raw %>%
      select(c(intersect(covariates_vector, covariates_vector_raw), intx_names))
    covariates_vector <- names(X_raw)
  }


  # Save center and scale values for use later on...
  save_scale_params <- function(x, atts_list) {
    att_list_elem <- atts_list[x] # was: scale_atts_list[x]
    if (is.data.frame(att_list_elem[[1]])) {
      att_list_elem[[1]] %>%
        mutate(fe_column = names(att_list_elem))
    }
  }

  # scale_atts_list <- attr(X_raw, 'scale_atts_list') # this is wrong!

  scale_atts_df <- lapply(names(scale_atts_list), save_scale_params,
    atts_list = scale_atts_list
  ) %>%
    bind_rows()



  # Relevels factors using the reference level info supplied in the analysis
  # config file.
  ref_lev_vec <- get_stuff_in_parens(str_split(covariates_string, ", ")[[1]])
  which_ref_lev <- which(!is.na(ref_lev_vec) & !grepl("deflections", ref_lev_vec))
  if (length(which_ref_lev) > 0) {
    for (i_ref_lev in which_ref_lev) {
      X_raw[[covariates_vector[i_ref_lev]]] <-
        relevel(X_raw[[covariates_vector[i_ref_lev]]],
          ref = ref_lev_vec[i_ref_lev]
        )
    }
  }

  # data %>% select(MgmtZone, MDCATY) %>%
  #   mutate(x1 = ifelse(MgmtZone == 'Upland', 1, 0),
  #          x2 = ifelse(MgmtZone == 'IrrigatedPasture', 1, 0),
  #          x3 = ifelse(MgmtZone %in% c('Upland', 'IrrigatedPasture'), 1, 0),
  #          intX_x1x3 = x1 * x3,  # are you upland or not?
  #          intX_x2x3 = x2 * x3) %>%  # are you irrigated or not?
  #   distinct()

  which_deflect <- which(grepl("deflections", ref_lev_vec))
  deflections_list <- NULL
  if (length(which_deflect) > 0) {
    deflections_list <- list()
    for (i_deflect in which_deflect) {
      n_deflections <- length(unique(X_raw[[covariates_vector[i_deflect]]]))
      contrasts(X_raw[[covariates_vector[i_deflect]]]) <-
        contr.sum(n_deflections, contrasts = TRUE)
      i_deflect_idx <- seq(i_deflect, i_deflect + n_deflections - 1)
      i_deflect_idx[which.max(i_deflect_idx)] <- NA
      names(i_deflect_idx) <- levels(X_raw[[i_deflect]])
      deflections_list[[covariates_vector[i_deflect]]] <- i_deflect_idx
    }

    # covariates_vector[which_deflect]
  }


  fixef_form <- formula(paste("~", paste(covariates_vector, collapse = "+")))
  # if (skip_check) browser()
  if ((nrow(X_raw) != nrow(model.frame(fixef_form, X_raw))) & !skip_check) {
    prime_suspects <-
      lapply(1:ncol(X_raw), function(i) which(is.na(X_raw[, i]))) %>% unlist()
    data %>% slice(unique(prime_suspects)) %T>%
      # select(matches('_id|_index|cal_year')) %T>%
      save_table(file.path(log_dir, "99-misc"), "incomplete-covariates-data.csv")
    # sink(file.path(log_dir, 'data-prep-error.txt'))
    stop('Covariates data probably includes NAs! May have happened after
            the join, in which case some sites are missing from the covariates
            data. Perhaps check that your covariates data includes observations
            for all years for which you have response info... for example, if
            you recently updated the response info to include the latest year
            of data, the covariates will need to be updated accordingly. Check
            "misc/incomplete-covariates-data.csv" for clues.')
    # sink()
  }


  # model.matrix(fixef_form,
  #              model.frame(fixef_form, X_raw)) %>%
  #   as_tibble(.name_repair = NULL) %>%
  #   # select(-`(Intercept)`) %>%
  #   setNames(paste0('fe_', names(.))) %>%
  #   # select(-fe_BotanistJA) %>%
  #   # select(-fe_MgmtZoneIrrigatedPasture, -fe_MgmtZoneRiparian) %>%
  #   distinct

  na_action <- if (skip_check) NULL else na.omit
  mm <- model.matrix(
    fixef_form,
    model.frame(fixef_form, X_raw, na.action = na_action)
  ) %>%
    as_tibble(.name_repair = NULL) %>%
    select(-`(Intercept)`) %>%
    setNames(paste0("fe_", names(.)))

  # For later recall, cache information about how each level of a categorical
  # covariate was encoded.
  for (deflected_covariate in covariates_vector[which_deflect]) {
    defl_lookup <- mm %>%
      bind_cols(X_raw %>% select(deflected_covariate)) %>%
      select(matches(deflected_covariate)) %>%
      distinct()
    save_table(
      defl_lookup,
      file.path(log_dir, "00-input"),
      sprintf(
        "encoding-for-%s.csv",
        deflected_covariate
      )
    )
  }

  list(
    mod_matrix = mm,
    scale_atts_df = scale_atts_df,
    deflections_list = deflections_list,
    intx_effects_list = intx_effects_list
  )
}

standardize_data <- function(file,
                             network_code,
                             unit_code, unit_code_col,
                             site_col, sample_cols, response_col, event_date_info, sampling_method, log_dir,
                             trials_col = NULL,
                             site_locations_file,
                             coord_cols,
                             stratum_col = NULL,
                             stratum_area_info = NULL,
                             agg_samples = FALSE,
                             covariates_info = NULL,
                             covariates = NULL,
                             censoring_info = NULL,
                             truncation_info = NULL,
                             add_cols = NULL,
                             verbose = TRUE,
                             trend_conditions = NULL,
                             check_response_scale = FALSE,
                             make_finite_pop_corr = FALSE,
                             finite_pop_info = NULL,
                             return_raws = FALSE,
                             eval_drivers = NULL,
                             me_crossings = NULL,
                             me_omissions = NULL,
                             exposure_var = NULL) {

  # Loads and prepares data for analysis by subsetting, adding indices for
  # strata and sites, and creating relative year.
  #
  # Args:
  #   file:               The name of the file from which the data are to be
  #                       read. The file path must be relative to the current
  #                       working directory (the project — or repo — root).
  #   unit_code:          The specific unit code. If the data consist of records
  #                       from multiple units, it will be filtered accordingly.
  #   unit_code_col:      Column name for park unit codes. The Park Service calls
  #                       each park a 'unit'.
  #   site_col:           Column name for site. Sites consist of collections of
  #                       samples — typically, transects, quadrats, or even
  #                       plots and are sometimes nested within strata (see
  #                       `is_stratified`, below).
  #   sample_cols:        The column(s) required to differentiate samples (e.g.,
  #                       replicates) within sites.
  #   response_col:       Column name of the response variable.
  #   event_date_info:    Sampling event information (the date-time column and
  #                       format) used to get the calendar year during which a
  #                       site was visited.
  #   coord_cols:         The names of the columns containing the projected
  #                       coordinates for samples. Order matters — X, then Y!
  #   stratum_col:        Column name for strata (if they exist). Only necessary
  #                       if strata is TRUE.
  #   agg_samples:        Logical variable. TRUE if there is a need/desire to
  #                       aggregate (sub)samples at the sample or site level.
  #                       For now, this argument is *only* used to compute the
  #                       number of successes and trials. Defaults to FALSE.
  #   covariates_file:    The name of the file containing space- and/or time-
  #                       varying explanatory variables.
  #   add_cols:           Additional columns to keep. Defaults to NULL (none).
  #   verbose:            Print messages to the console and create a logfile?
  #
  # Returns:
  #   A tibble containing *just* the data needed for analysis (by JAGS).


  response_data <- load_data(file)

  # Filter for strata appearing in _park-level-attributes.yml
  omitted_stata_id <- setdiff(
    unique(response_data[[stratum_col]]), names(stratum_area_info)
  )
  if (length(omitted_stata_id) > 0) {
    message(sprintf("Removing: %s", paste(omitted_stata_id, collapse = ", ")))
    message("Inspect your _park-level-attributes.yml if this is unexpected!")
    response_data <- response_data %>%
      filter(.[[stratum_col]] %in% names(stratum_area_info))
  }
  # browser()

  covariates_switch <- !is.null(covariates_info$file) & !is.na(covariates)
  explanatory_vars_raw <- parse_covariates_info(covariates)

  if (!is.null(site_locations_file)) {
    site_location_data <- load_data(site_locations_file) %>%
      select(c(intersect(names(.), names(response_data)), all_of(coord_cols)))
  }

  is_stratified <- ifelse(!is.null(stratum_area_info), TRUE, FALSE) # TODO, add switch to permit running both ways, even if strata exist

  get_truncation_bounds <- function(x, y) {
    if (is.null(x)) {
      return(NA)
    }
    if (x %in% names(y)) {
      y[[x]]
    } else {
      x
    }
  }

# browser()
  jags_df_raw <- response_data %>%
    # Filter the data using `unit_code`. (Reduces the 'overhead' of the joins.)
    filter_(paste(unit_code_col, "==", shQuote(unit_code))) %>%
    # Rename `unit_code_col`, as necessary.
    mutate_(.dots = setNames(list(rename_call(unit_code_col)), "unit_code")) %>%
    # As necessary, rename `site_col` and create indices for sites.
    mutate_(.dots = setNames(list(rename_call(site_col)), "site_id")) %>%
    mutate(site_index = index_from_category(site_id)) %>%
    # Rename the response variable.
    {
      `if`(
        !is.null(trials_col),
        mutate_(., .dots = setNames(list(rename_call(response_col)), "hits")) %>%
          mutate_(.dots = setNames(list(rename_call(trials_col)), "trials")),
        mutate_(., .dots = setNames(list(rename_call(response_col)), "response"))
      )
    } %>%
    # jags_df %>%
    # Deal with any censoring.
    {
      `if`(
        !is.null(censoring_info),
        left_join(., filter_(., paste(censoring_info$values, collapse = "|")) %>%
          mutate(is_censored = TRUE, censor_limit_vec = .[[censoring_info$limit]])) %>%
          mutate(censor_limit_vec_rep = .[[censoring_info$limit]]) %>% # new!
          mutate(
            is_censored = ifelse(is.na(is_censored), FALSE, is_censored),
            response = ifelse(is_censored, NA, response),
            censor_limit_vec = ifelse(is.na(censor_limit_vec), response, censor_limit_vec)
          ), # was:
        .
      )
    } %>%
    # Deal with any truncation.
    {
      `if`(
        !is.null(truncation_info),
        mutate(.,
          truncated_below = get_truncation_bounds(truncation_info$L, .),
          truncated_above = get_truncation_bounds(truncation_info$U, .)
        ),
        .
      )
    } %>%
    ## Final censoring/trucation data manipulations. # TODO: removed b/c this is not general! talk to Carolyn
    # {`if`(!is.null(truncation_info) & !is.null(censoring_info),
    #      mutate(., truncated_below = ifelse(is_censored, 0, truncated_below)),
    #      .)} %>%
    # Format event date and distinguish between calendar and 'relative' year.
    {
      `if`(
        event_date_info$`date-time column` %in% names(.),
        format_dttm_attr(
          ., event_date_info$`date-time column`,
          event_date_info$`date-time format`
        ),
        .
      )
    } %>%
    # Join the spatial and other time-varying covariates (again, only if
    # specified) to the data.

    {
      `if`(
        covariates_switch,
        prep_and_join_cov_info(., covariates_info, explanatory_vars_raw, sample_cols), .
      )
    } %>%
    # (If necessary, join, and then) rename site-level coordinate columns to
    # 'site_x' and 'site_y'.
    {
      `if`(
        !is.null(site_locations_file),
        left_join(., site_location_data) %>%
          mutate_(.dots = setNames(list(rename_call(coord_cols[1])), "site_x")) %>%
          mutate_(.dots = setNames(list(rename_call(coord_cols[2])), "site_y")),
        .
      )
    } %>%
    # Create indices for strata and sites nested in strata, if strata exist.
    # Defaults to stratum_index = 1 if `is_stratified` == FALSE.
    {
      `if`(
        is_stratified,
        mutate_(., .dots = setNames(
          list(rename_call(stratum_col)),
          "stratum_id"
        )),
        mutate(., stratum_id == "None")
      )
    } %>%
    mutate(stratum_index = index_from_category(stratum_id))
# browser()
  eval_mean_for_tv_covariates <- FALSE
  jags_df <- jags_df_raw %>%
    group_by(stratum_index) %>%
    mutate(site_in_stratum_index = index_from_category(site_id)) %>%
    ungroup() %>%
    # Create unique sequential indices for samples within each site.
    unite_("sample_id", sample_cols, sep = "-", remove = FALSE) %>%
    mutate(sample_index = index_from_category(sample_id)) %>% # was: %T>%
    # set_attr('eval_mean_for_tv_covariates', FALSE) %>%
    set_attr("full_cov_data", attr(jags_df_raw, "full_cov_data"))
  # set_attr('X.driver', NULL) %>%
  # set_attr('which.drivers', NULL)


  # pipe_assign('eval_mean_for_tv_covariates', FALSE) %T>%
  # pipe_assign('X.driver', NULL) %T>%
  # pipe_assign('which.drivers', NULL)

  save_object(
    eval_mean_for_tv_covariates,
    file.path(log_dir, "00-input"),
    "eval-mean-for-tv-covariates.rds"
  )


  if (make_finite_pop_corr & !is.null(finite_pop_info$file)) {
    n_sites_cols <- finite_pop_info$`columns with the number of sampled and unsampled sites`
    y_n_unobs_site <- load_data(finite_pop_info$file) %>%
      left_join(jags_df %>% select(stratum_col, stratum_index) %>% distinct()) %>%
      arrange(stratum_index) %>%
      pull(.[[n_sites_cols[2]]])
    # pipe_assign('y.n.unobs.site', .[[n_sites_cols[2]]])
  }

  if (covariates_switch) {
    covariates_info$mapping == "sample"
    if (is.null(covariates_info$mapping)) {
      covariates_info$mapping <- "site"
    }
    other_selectors <- if (covariates_info$mapping == "sample") {
      c("sample_id", "sample_index")
    } else {
      NULL
    }

    # TODO: deal with scenarios in which we don't have complete observations.
    X_pred_tmp <- attr(jags_df, "full_cov_data") %>%
      {
        `if`(
          "rel_year" %in% names(.),
          filter(., rel_year >= 0, rel_year <= max(jags_df$rel_year)),
          crossing(., rel_year = full_seq(jags_df$rel_year, 1))
        )
      } %>%
      # filter(rel_year >= 0, rel_year <= max(max(jags_df$rel_year))) %>%
      inner_join(jags_df %>%
        select(
          !!stratum_col, !!site_col, site_in_stratum_index, stratum_index, # 11/11/21 added !!stratum_col
          other_selectors
        ) %>%
        distinct()) %>%
      mutate(is_in_sample = TRUE)


    if (make_finite_pop_corr & covariates_switch &
      isTRUE(finite_pop_info$`covariate info includes unsampled sites`)) { # is.null(finite_pop_info)
      unsampled_site_cov_data <- attr(jags_df, "full_cov_data") %>%
        {
          `if`(
            "rel_year" %in% names(.),
            filter(., rel_year >= 0, rel_year <= max(jags_df$rel_year)),
            crossing(., rel_year = full_seq(jags_df$rel_year, 1))
          )
        } %>%
        # filter(rel_year >= 0, rel_year <= max(max(jags_df$rel_year))) %>%
        anti_join(jags_df %>%
          select(!!site_col, site_in_stratum_index, stratum_index) %>%
          distinct()) %>%
        left_join(distinct(jags_df, !!stratum_col := get(stratum_col), stratum_index))

      X_pred_tmp <- X_pred_tmp %>%
        bind_rows(unsampled_site_cov_data) %>%
        group_by(!!stratum_col := get(stratum_col)) %>%
        do(append_unsampled_site_indices(., site_col)) %>%
        ungroup()


      y_n_unobs_site <- unsampled_site_cov_data %>%
        group_by(stratum_index) %>%
        summarise(n_unobs_sites = n_distinct(get(site_col))) %>%
        arrange(stratum_index) %>%
        pull(n_unobs_sites)

      # pipe_assign('y.n.unobs.site', .$n_unobs_sites)
      # X_pred_tmp %>%
      #   group_by(stratum_index) %>%
      #   summarise(n_sites = max(site_in_stratum_index)) %>%
      #   arrange(stratum_index) %T>%
      #   pipe_assign('y.n.unobs.site', .$n_sites)
    }

    X_pred_raw <- X_pred_tmp %>%
      arrange(stratum_index, site_in_stratum_index, rel_year) %T>%
      # pipe_assign('in.sample.idx', which(.$is_in_sample)) %T>%
      # pipe_assign('n.x.pred', .$rel_year %>% n_distinct) %T>%
      # pipe_assign('x.pred.index', group_indices(., rel_year)) %T>%
      # pipe_assign('x.pred.raw', .$rel_year) %T>%
      # # pipe_assign('i.pred', .$sample_index) %T>% # maybe for a future time!
      # pipe_assign('j.pred', .$site_in_stratum_index) %T>%
      # pipe_assign('k.pred', .$stratum_index) %T>%
      # pipe_assign('eval_mean_for_tv_covariates',
      #             (nrow(distinct(., stratum_index, site_in_stratum_index, rel_year)) >= # was just >
      #               nrow(distinct(jags_df, stratum_index, site_in_stratum_index, rel_year))) &
      #               is.null(covariates_info$`destroy pred blocks`)) %T>%  # was nrow(.) > nrow(jags_df)
      save_table(file.path(log_dir, "00-input"), "complete-covariates-data.csv")
    X_pred <- X_pred_raw

    eval_mean_for_tv_covariates <- X_pred_raw %>%
      summarise(
        tmp =
          (nrow(distinct(., stratum_index, site_in_stratum_index, rel_year)) >= # was just >
            nrow(distinct(jags_df, stratum_index, site_in_stratum_index, rel_year))) &
            is.null(covariates_info$`destroy pred blocks`)
      ) %>%
      pull(tmp)
    save_object(
      eval_mean_for_tv_covariates,
      file.path(log_dir, "00-input"),
      "eval-mean-for-tv-covariates.rds"
    )
  }


  # Design matrix for additional covariates.
  if (!is.na(covariates)) {

    # jags_df %>%
    #   select_if(function(x) any(is.na(x))) %>%
    #   summarise_each(list(~sum(is.na(.))))
    if (nrow(jags_df) != nrow(na.omit(select(jags_df, response_col, explanatory_vars_raw)))) {
      warning(sprintf(
        "`jags_df` has %s rows containing NAs. Removing them now (CAUTION)!",
        nrow(jags_df) - nrow(na.omit(jags_df))
      ))
      jags_df <- drop_na(jags_df, response_col, explanatory_vars_raw) %>% # na.omit(jags_df) %>%
        group_by(stratum_index) %>% # THIS IS NEW AS OF Dec 20 2019
        mutate(site_in_stratum_index = index_from_category(site_id)) %>%
        ungroup()
      # rowSums(is.na)
    }

    X_list <- jags_df %>% get_design_mat(covariates, log_dir)
    save_object(
      X_list$scale_atts_df,
      file.path(log_dir, "00-input"),
      "covariate-moments.rds"
    )
    X_scaled <- X_list$mod_matrix
    explanatory_vars <- names(X_scaled)
    jags_df %<>% bind_cols(X_scaled)

    if (!is.null(X_list$deflections_list)) {
      save_object(
        X_list$deflections_list,
        file.path(log_dir, "00-input"),
        "deflections-indices.rds"
      )
    }

    if (!is.null(trend_conditions)) {
      applicable_tc <-
        trend_conditions[names(trend_conditions) %in% explanatory_vars_raw]

      # applicable_tc <- trend_conditions[
      #   grep(paste(names(trend_conditions), collapse = "|"), explanatory_vars)
      # ]

      for (col_name in names(applicable_tc)) {
        if (is.character(X_pred[[col_name]])) {
          unique_levels <- unique(X_pred[[col_name]])
          X_pred[col_name] <- factor(trend_conditions[[col_name]],
            levels = unique_levels
          )
        } else if (all(X_pred[[col_name]] %in% c(0, 1)) &
          applicable_tc[[col_name]] %in% c(0, 1)) {
          X_pred[col_name] <- applicable_tc[[col_name]]
        } else {
          col_mean <- X_list$scale_atts_df %>%
            filter(fe_column == col_name) %>%
            pull(`scaled:center`)
          X_pred[col_name] <- col_mean
        }
      }
    }

    if (!is.null(eval_drivers) & isTRUE(eval_drivers)) {
      which_drivers <- if (!is.null(me_omissions)) {
        grep(paste(paste0("_", me_omissions), collapse = "|"), # was: paste0('_', me_omissions, '$')
          dimnames(X_list$mod_matrix)[[2]],
          invert = TRUE
        )
      } else {
        1:ncol(X_list$mod_matrix)
      }
      driver_names <- dimnames(X_list$mod_matrix)[[2]][which_drivers]
      v_minmax <- lapply(which_drivers, function(v) {
        v_name <- names(X_list$mod_matrix[, v])
        v_name_raw <- sprintf("%s_raw", v_name)
        incr_name <- sprintf("incr_%s", v_name)

        if (all(unique(unlist(X_list$mod_matrix[, v])) %in% c(0, 1))) {

          # If it's categorical handle it one way
          v_increments_vec <- c(zero = 0, pos_one = 1)
          v_increments_df <-
            enframe(sort(v_increments_vec), name = incr_name, value = v_name)
          v_increments_df %>% mutate(!!v_name_raw := get(v_name))
        } else {

          # Do the usual thing.
          v_increments_vec <- c(
            neg_one = -1, zero = 0, pos_one = 1,
            seq(-1.96, 1.96, length.out = 11),
            min_seen = min(X_list$mod_matrix[, v]),
            max_seen = max(X_list$mod_matrix[, v])
          )
          # was: sort(v_increments_vec[-which(duplicated(v_increments_vec))])
          v_increments_df <- enframe(sort(v_increments_vec),
            name = incr_name, value = v_name
          ) %>%
            arrange(get(v_name)) %>%
            mutate(!!incr_name := ifelse(get(incr_name) == "", NA, get(incr_name))) %>%
            distinct()

          # if (v_name == "fe_point_jurisdiction1") browser()
          driver_scale_atts <- X_list$scale_atts_df %>% # X_list$scale_atts_df[v, ]
            filter(grepl(sub("^fe_", "", v_name), fe_column))
          v_scale_raw <- driver_scale_atts$`scaled:scale`
          v_scale <- ifelse(identical(v_scale_raw, character(0)), NA, v_scale_raw)
          v_center_raw <- driver_scale_atts$`scaled:center`
          v_center <- ifelse(identical(v_center_raw, character(0)), NA, v_center_raw)
          # if (!v_name %in% c("fe_Prev3", "fe_tmaxmean")) browser()
          tryCatch(
            expr = {
              v_increments_df %>%
                mutate(!!v_name_raw := get(v_name) * v_scale + v_center) %>%
                filter(
                  get(v_name) >= get(v_name)[which(get(incr_name) == "min_seen")],
                  get(v_name) <= get(v_name)[which(get(incr_name) == "max_seen")]
                )
            },
            error = function(e) {
              browser()
            }
          )
          # v_increments_df %>%
          #   mutate(!!v_name_raw := get(v_name) * v_scale + v_center) %>%
          #   filter(get(v_name) >= get(v_name)[which(get(incr_name) == "min_seen")],
          #          get(v_name) <= get(v_name)[which(get(incr_name) == "max_seen")])
        }
      })

      get_crossings <- function(x) {
        paste(sprintf(".[, 1] == '%s'", x), collapse = "|")
      }

      cond_on_str <- if (is.null(me_crossings)) get_crossings("zero") else get_crossings(me_crossings)

      X_driver <- lapply(seq_along(v_minmax), function(x) {
        # lapply(v_minmax[which(!seq_along(v_minmax) %in% x)], na.omit)
        x_at <- lapply(which(!seq_along(v_minmax) %in% x), function(y) {
          v_minmax[[y]] %>% filter_(cond_on_str)
        })


        expand.grid.df <- function(...) Reduce(function(...) merge(..., by = NULL), list(...))
        this_v_minmax <- if (is.null(me_crossings)) {
          v_minmax[x]
        } else {
          list(v_minmax[[x]] %>%
            drop_na(matches("^incr_")) %>%
            filter(.[[1]] %in% me_crossings)) # TODO: filter for me_crossings!!!
        }
        as_tibble(do.call(expand.grid.df, c(this_v_minmax, x_at))) %>% # was: c(v_minmax[x], x_at)
          mutate(me_idx = x, me_along = driver_names[x]) # was: names(X_list$mod_matrix[, x])
      })


      X_driver_uniq <- distinct(bind_rows(X_driver)) # %T>%
      attr(X_driver_uniq, "X.driver") <- as.matrix(select(X_driver_uniq, all_of(driver_names)))
      attr(X_driver_uniq, "which.drivers") <- which_drivers
      # pipe_assign('X.driver', as.matrix(select(., all_of(driver_names)))) %T>%
      # pipe_assign('which.drivers', which_drivers)
      save_table(X_driver_uniq,
        path = file.path(log_dir, "00-input"),
        filename = "X-driver-lookup.csv"
      )
    }

    X_pred_stuff <- get_design_mat(X_pred, covariates, log_dir,
      center_and_scale_params = X_list$scale_atts_df, skip_check = TRUE
    )

    attr(X_pred_stuff, "X.pred") <- as.matrix(X_pred_stuff$mod_matrix)
    # X_pred_stuff$mod_matrix %T>%
    #   pipe_assign('X.pred', as.matrix(.))
    X_pred_deflections <- X_pred_stuff$deflections_list
    if (!is.null(X_pred_deflections)) {
      attr(X_pred_deflections, "y.zone.index") <-
        as.integer(factor(X_pred[[names(X_pred_deflections)]],
          levels = names(X_pred_deflections[[1]])
        ))
      # factor(X_pred[[names(X_pred_deflections)]],
      #        levels = names(X_pred_deflections[[1]])) %T>%
      #   pipe_assign('y.zone.index', as.integer(.))
    }
  } else {
    explanatory_vars <- NULL
  }

  if (!is.null(site_locations_file)) {
    loc_utils_pattern <- paste0(
      "sample-locations-", tolower(substr(network_code, 1, 4)), ".r"
    )
    tryCatch(
      {
        source(
          list.files("model-api/src", pattern = loc_utils_pattern, full.names = TRUE)
        )
      },
      error = function(e) {
        message("Using non-specific site loc utils, see sample-locations-misc.r
                in model-api/src. Consider contributing your own!")
        source("model-api/src/sample-locations-misc.r")
      }
    )

    # Get point-intercept locations for each transect in each plot.
    # jags_df %>%
    #   group_by(site_id, sample_id) %>%
    #   do(transect_point_intercepts(.)) #%>%
    #   # ungroup %>%
    #   # mutate(pi_geometry = paste0('POINT(',pi_x, ' ',pi_y,')'))
    if (sampling_method == "transect" | sampling_method == "point_intercept") {
      sample_locations <- jags_df %>%
        # rowwise %>%
        group_by(site_id, sample_id) %>%
        do(get_sample_locations(., sampling_method, sample_cols)) %>%
        ungroup() %>%
        distinct() %>%
        select(-contains("fe")) # TODO: fix by removing the need for calling distinct!!!

      # jags_df %>% group_by(site_id, rel_year) %>%
      #   summarise(n_transects=n_distinct(Transect)) %>%
      #   filter(n_transects != 3)

      # jags_df %>% anti_join(sample_locations)
      # bind_cols(jags_df %>% arrange(site_index, rel_year) %>% select(contains('fe')),
      #           sample_locations %>% arrange(site_index, rel_year) %>% select(contains('fe'))) %>%
      #   mutate(fe1 = fe_Prev.Oct.Apr.P == fe_Prev.Oct.Apr.P1,
      #          fe2=fe_Prev.July.Sept.P == fe_Prev.July.Sept.P1) %>%
      #   filter(!fe2)


      jags_df %<>%
        left_join(sample_locations)
    } else if (sampling_method == "plot") {
      jags_df %<>%
        group_by(site_id, sample_id) %>%
        do(get_sample_locations(., sampling_method, sample_cols)) %>% # TODO: is sample_cols necessary?
        ungroup()
    }
  }
  if (return_raws) {
    return(jags_df)
  }

  # xy_coords_to_lat_lon('site', PCS_EPSG) %>%
  # Compute centroids for each plot-based sample.

  relevant_vars <-
    bar_concat(
      "^", c(
        "unit_code",
        right_cat(
          c("site", "sample", "stratum", "site_in_stratum"),
          c("_id", "_index", "_x", "_y")
        ),
        right_cat(c("cal", "rel"), "_year"), add_cols
      ),
      "$"
    )
  rhs_vars <- bar_concat("^", explanatory_vars, "$")

  if (agg_samples) {
    grouped_jags_df <- jags_df %>%
      group_by_at(vars(dplyr::matches(relevant_vars)))
    rhs <- grouped_jags_df %>% # the right-hand side of the equation
      summarise_at(vars(dplyr::matches(rhs_vars)), list(~mean))
    # Calculate the total number of successes and trials to allow a binomial (as
    # opposed to a Bernoulli) likelihood. Vastly increases speed.
    lhs <- grouped_jags_df %>% # the left-hand side
      summarise_at(vars("response"), list(hits = ~sum, trials = ~ n()))
    d_out <- left_join(lhs, rhs) %>% ungroup()
  } else {
    all_vars <- bar_concat(c(
      relevant_vars, rhs_vars,
      "^response$", "^hits$", "^trials$",
      "^is_censored$", "^censor_limit_vec$",
      "^censor_limit_vec_rep$",
      "^truncated_below$", "^truncated_above$"
    ))

    d_out <- jags_df %>%
      select(matches(all_vars), any_of(exposure_var)) %>%
      # select_at(vars(matches(all_vars))) %>%
      mutate(network_code = network_code)

    # any(d_out %>% mutate(sanity_test=hits>trials) %>% pull(sanity_test))
  }


  if (check_response_scale) {
    if (any(d_out$response)) warning("NAs exist in the response column, if not censored, check this!")
    if (max(d_out$response, na.rm = FALSE) > 1) {
      d_out$response <- d_out$response / 100
    }
    # prevents logit(0) = -Inf and prevents logit(1) = Inf
    d_out$response[d_out$response == 0] <- .00001
    d_out$response[d_out$response == 1] <- .99999
  }

  # If `verbose`, print miscellaneous structural information for the data.
  if (verbose) {
    environment(message_for_jags_prep) <- environment()
    message_for_jags_prep()

    logger(d_out, log_dir)
  }

  tc <- function(x) {
    tryCatch(x, error = function(e) NULL)
  }

  attr(d_out, "y.n.unobs.site") <- tc(y_n_unobs_site)
  attr(d_out, "X.pred") <- tc(attr(X_pred_stuff, "X.pred"))
  attr(d_out, "y.zone.index") <- tc(attr(X_pred_deflections, "y.zone.index"))
  attr(d_out, "eval_mean_for_tv_covariates") <- tc(eval_mean_for_tv_covariates)
  attr(d_out, "in.sample.idx") <- tc(which(X_pred_raw$is_in_sample))
  attr(d_out, "n.x.pred") <- tc(X_pred_raw$rel_year %>% n_distinct())
  attr(d_out, "x.pred.index") <- tc(group_indices(X_pred_raw, rel_year))
  attr(d_out, "x.pred.raw") <- tc(X_pred_raw$rel_year)

  attr(d_out, "j.pred") <- tc(X_pred_raw$site_in_stratum_index)
  attr(d_out, "k.pred") <- tc(X_pred_raw$stratum_index)

  attr(d_out, "X.driver") <- tc(attr(X_driver_uniq, "X.driver"))
  attr(d_out, "which.drivers") <- tc(attr(X_driver_uniq, "which.drivers"))

  attr(d_out, "intx_effects_list") <- tc(X_list$intx_effects_list)

  d_out
}

match_value <- function(x, pattern, ignore.case = TRUE) {
  grep(pattern, x, ignore.case, value = TRUE)
}

append_unsampled_site_indices <- function(data, site_col) {
  data$site_in_stratum_index[is.na(data$site_in_stratum_index)] <-
    group_indices(
      filter(data, is.na(data$site_in_stratum_index)),
      get(site_col)
    ) +
    max(data$site_in_stratum_index, na.rm = TRUE)
  data
}
