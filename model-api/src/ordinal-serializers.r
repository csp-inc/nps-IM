get_site_inference_ord <- function(mcarray_list, data, folder, resp_var_desc,
                                   timesteps, jags_obj_prefix, likelihood,
                                   wrap_char_n = 15, n_sites_to_show = 5, ...) {
  print("Inference on sites for ordinal outcomes...")

  jags_obj_suffix <- ifelse(likelihood == "ordinal-latent-normal",
    "site.mean", "site.class.mean"
  )
  # Preliminaries.
  stratum_lookup <- get_stratum_lookup(data)
  site_lookup <- get_site_lookup(data)
  year_lookup <- get_year_lookup(data, timesteps)
  site_means_obj <- paste(jags_obj_prefix, jags_obj_suffix, sep = ".")
  # stratum_means_obj <- paste(jags_obj_prefix, 'strat.mean', sep='.')

  # site_mu_array <-
  #   summary(mcarray_list[[site_means_obj]], mean, na.rm=TRUE)$stat

  get_site_stats_df <- function(site_stat_array, stat) {
    x_index <- 0
    apply(site_stat_array, MARGIN = 3, function(x) {
      x_index <<- x_index + 1
      x %>%
        na.omit() %>%
        as_tibble(.name_repair = NULL) %>%
        mutate(site_in_stratum_index = 1:n(), stratum_index = x_index) %>%
        left_join(stratum_lookup) %>%
        left_join(site_lookup) %>%
        gather(year_col, !!stat, -matches("stratum"), -site_id) %>%
        left_join(year_lookup) %>%
        select(-matches("index|col"))
    }) %>% reduce(rbind)
  }

  site_stats_df <- left_join(
    # Mean and SE.
    lapply(c("mean", "sd"), function(x) {
      get_site_stats_df(summary(mcarray_list[[site_means_obj]], x)$stat, x)
    }) %>% reduce(left_join),
    # Upper and lower bounds (95% CI).
    lapply(c(.025, .975), function(x) {
      get_site_stats_df(
        summary(mcarray_list[[site_means_obj]],
          quantile,
          probs = x, na.rm = TRUE
        )$stat,
        as.character(x)
      )
    }) %>% reduce(left_join)
  )

  # Save plot components for future plotting.
  site_stats_df %T>%
    save_object(
      file.path(folder, "03-inference/site"),
      paste(jags_obj_prefix,
        "site-in-stratum-class-means.rds",
        sep = "-"
      )
    )

  # Save tabular output.
  site_stats_df %>%
    # We used fractional years to get nice smooth lines for plotting. Tabular
    # outputs are probably best summarized by whole years.
    filter(year %% 1 == 0) %>%
    # Reorder columns for easier reading.
    select(stratum_id, site_id, year, mean, sd, `0.025`, `0.975`) %T>%
    save_table(
      file.path(folder, "03-inference/site"),
      paste(jags_obj_prefix, "site-in-stratum-class-means.csv", sep = "-")
    )

  # Plotting.
  y_lev <- "Site-level"
  y_lab <- as.expression(bquote(atop(
    .(y_lev),
    .(tolower(resp_var_desc)) ~ (mu)
  )))
  selected_sites <- data %>%
    distinct(stratum_id, site_id) %>%
    group_by(stratum_id) %>%
    sample_n(n_sites_to_show) %>%
    mutate(selected_site_index = as.numeric(as.factor(site_id))) %>%
    ungroup()
  site_stats_df_select <- site_stats_df %>%
    filter(site_id %in% selected_sites$site_id) %>%
    left_join(selected_sites)
  data_select <- data %>%
    filter(response > 0, site_id %in% selected_sites$site_id) %>%
    left_join(selected_sites)
  p_ss_class_means <- ggplot() +
    facet_grid(selected_site_index ~ stratum_id,
      labeller = label_wrap_gen(wrap_char_n)
    ) +
    geom_ribbon(
      data = site_stats_df_select, aes(x = year, ymin = `0.025`, ymax = `0.975`, group = site_id),
      fill = ima_pal("default")[2],
      alpha = .2, color = NA
    ) +
    geom_line(
      data = site_stats_df_select, aes(x = year, y = mean, group = site_id),
      color = ima_pal("default")[2]
    ) +
    geom_jitter(
      data = data_select,
      aes(x = cal_year, y = response, fill = response), alpha = .8,
      width = .2, height = .1, pch = 21, size = 2
    ) +
    scale_fill_viridis_c("Cover class") +
    geom_label(
      data = data_select %>%
        distinct(selected_site_index, stratum_id, site_id),
      aes(x = -Inf, y = Inf, label = site_id),
      vjust = "inward", hjust = "inward"
    ) +
    labs(x = "Year", y = y_lab) +
    scale_x_continuous(breaks = scales::pretty_breaks()) +
    theme_ima("default")
  save_figure(
    plot = p_ss_class_means,
    path = file.path(folder, "03-inference/site"),
    filename = paste(jags_obj_prefix,
      "site-in-stratum-class-means.png",
      sep = "-"
    ),
    p_width = 1.5, p_height = 1.5, device = NULL
  )
}

get_stratum_inference_ord <- function(mcarray_list, data, folder, resp_var_desc,
                                      timesteps, jags_obj_prefix, likelihood,
                                      wrap_char_n = 15, ...) {
  print("Inference on strata for ordinal outcomes...")

  # strat_mean_suffixes <- c('', '.oos')
  oos_entry <- ifelse(paste(jags_obj_prefix, "strat", "mean", "oos", sep = ".") %in% names(mcarray_list),
    ".oos", NA
  )
  strat_mean_suffixes <- na.omit(c("", oos_entry))
  lapply(strat_mean_suffixes, function(strat_mean_suffix) {
    sms_dash <- sub("\\.", "-", strat_mean_suffix)
    sm_sub_dir <- ifelse(grepl("oos", strat_mean_suffix), "all", "sampled")
    lapply(c("phi", "mu"), function(component) {
      # stratum_means_obj <- paste(jags_obj_prefix, component, sep='.')

      stratum_means_obj <- paste0(
        paste(component, jags_obj_prefix, "strat", sep = "."), strat_mean_suffix
      )
      # stratum_means_obj <- paste0(
      #   paste(jags_obj_prefix, component, sep='.'), strat_mean_suffix)
      stratum_lookup <- get_stratum_lookup(data)
      stratum_mu_full <- concat_chains(mcarray_list[[stratum_means_obj]], axis = 4)
      year_vec <- data %>%
        get_year_lookup(timesteps) %>%
        pull(year)

      x_index <- 0

      stratum_mu_draws <- apply(stratum_mu_full %>% thin_mcmc_output(axis = 3, ...),
        MARGIN = 2, function(x) {
          x_index <<- x_index + 1
          x %>%
            as_tibble(.name_repair = NULL) %>%
            mutate(year = year_vec, stratum_index = x_index) %>%
            left_join(stratum_lookup) %>%
            gather(k, stratum_mu, -year, -matches("stratum"))
        }
      ) %>% reduce(rbind)

      stratum_mean_mu <- summary(mcarray_list[[stratum_means_obj]], mean)$stat %>%
        as_tibble(.name_repair = NULL) %>%
        mutate(year = year_vec) %>%
        gather(stratum_col, stratum_mean_mu, -year) %>%
        left_join(stratum_lookup)

      x_index <- 0
      stratum_mean_cis <-
        apply(summary(mcarray_list[[stratum_means_obj]], quantile, c(.025, .975))$stat,
          MARGIN = 3, function(x) {
            x_index <<- x_index + 1
            x %>%
              as_tibble(.name_repair = NULL) %>%
              mutate(ci = c("0.025", "0.975"), stratum_index = x_index) %>%
              gather(year_col, ci_value, -ci, -stratum_index) %>%
              left_join(get_year_lookup(data, timesteps)) %>%
              left_join(stratum_lookup) %>%
              select(-matches("index|col"))
          }
        ) %>% reduce(rbind)

      # Save plot components for future plotting.
      component_description <- ifelse(component == "strat.mean.01",
        "prob-bare-ground", "cover-when-present"
      )
      list(
        stratum_mean_mu = stratum_mean_mu, stratum_mu_draws = stratum_mu_draws,
        stratum_mean_cis = stratum_mean_cis
      ) %T>%
        save_object(
          file.path(folder, "03-inference/strata", sm_sub_dir),
          paste(jags_obj_prefix,
            paste0(
              "stratum-means-", component_description,
              sms_dash, ".rds"
            ),
            sep = "-"
          )
        )

      # Save tabular output.
      temp <- stratum_mean_cis %>%
        select(-rel_year, -cal_year) %>%
        spread(ci, ci_value)
      stratum_mu_draws %>%
        # We used fractional years to get nice smooth lines for plotting. Tabular
        # outputs are probably best summarized by whole years.
        filter(year %% 1 == 0) %>%
        rename(mu = stratum_mu) %>%
        # For each year, compute the mean and standard deviation.
        group_by(year, stratum_id) %>%
        summarise_at("mu", list(~ mean(.), ~ sd(.))) %>%
        ungroup() %>%
        left_join(temp, by = c("year", "stratum_id")) %T>%
        save_table(
          file.path(folder, "03-inference/strata", sm_sub_dir),
          paste(jags_obj_prefix,
            paste0(
              "stratum-means-", component_description,
              sms_dash, ".rds"
            ),
            sep = "-"
          )
        )


      # fncols <- function(data, cname) {
      #   # https://stackoverflow.com/questions/45857787/adding-column-if-it-does-not-exist
      #   add <-cname[!cname%in%names(data)]
      #
      #   if(length(add)!=0) data[add] <- 1
      #   data
      # }
      #
      # stratum_data <- data %>%
      #   select(year=cal_year, stratum_id, outcome=matches("^response$|^hits$"),
      #          normalizer = matches('trials')) %>%
      #   fncols('normalizer') %>%
      #   mutate(outcome = outcome / normalizer)

      # Plotting.
      y_lev <- "Stratum-level"
      get_y_lab <- function(showing_y_new) {
        if (showing_y_new) {
          as.expression(
            bquote(atop(.(y_lev), .(tolower(resp_var_desc)) ~
              (mu * "," ~ italic(y)^new)))
          )
        } else {
          as.expression(
            bquote(atop(.(y_lev), .(tolower(resp_var_desc)) ~ (mu)))
          )
        }
      }

      site_new_obs_obj <- if (component == "phi") {
        paste(jags_obj_prefix, "site.class0.new.obs", sep = ".")
      } else {
        paste(jags_obj_prefix, "site.gt0.new.obs", sep = ".")
      }
      obj_year_lookup <- get_year_lookup(data, timesteps) # new!
      # obj_year_lookup <- if(jags_obj_prefix == 'pred') {
      #  get_year_lookup(data, timesteps)
      # } else {
      #  get_year_lookup(data, x.hat.raw)
      # }
      # if(site_new_obs_obj != 'phi') browser()
      if (!is.null(mcarray_list[[site_new_obs_obj]])) {
        n_strata <- dim(mcarray_list[[site_new_obs_obj]])[3]
        if (n_strata == 1) {
          site_new_obs_tmp <- concat_chains(mcarray_list[[site_new_obs_obj]][, , 1, , ], axis = 4)
          site_new_obs <- array(NA, dim = c(
            nrow(site_new_obs_tmp), ncol(site_new_obs_tmp),
            n_strata, tail(dim(site_new_obs_tmp), 1)
          ))
          site_new_obs[, , 1, ] <- site_new_obs_tmp
        } else {
          site_new_obs <- concat_chains(mcarray_list[[site_new_obs_obj]], axis = 5)
        }

        new_obs_df <- lapply(seq_len(dim(site_new_obs)[3]), function(x) {
          stratum_site_new_obs <- site_new_obs[, , x, ]

          new_obs_bcis <- apply(stratum_site_new_obs,
            MARGIN = 2, FUN = quantile,
            probs = c(.025, .975), na.rm = TRUE
          )
          new_obs_bcis %>%
            as_tibble(.name_repair = NULL) %>%
            mutate(stat = c("lwr", "upr")) %>%
            gather(year_col, stat_val, -stat) %>%
            left_join(obj_year_lookup) %>%
            mutate(stratum_col = paste0("V", x))
        }) %>%
          bind_rows() %>%
          left_join(get_stratum_lookup(data))

        p_comp <- ggplot() +
          facet_wrap(~stratum_id, labeller = label_wrap_gen(wrap_char_n)) +
          geom_line(
            data = new_obs_df, aes(x = year, y = stat_val, group = stat),
            linetype = "twodash", alpha = .5
          ) +
          geom_line(
            data = stratum_mean_mu, aes(x = year, y = stratum_mean_mu),
            color = ima_post_col("default")[2], size = 1.01
          ) +
          geom_line(
            data = stratum_mean_cis, aes(x = year, y = ci_value, group = ci),
            color = ima_post_col("default")[2], linetype = "twodash", size = 1.01
          ) +
          labs(x = "Year", y = get_y_lab(TRUE)) +
          scale_x_continuous(breaks = scales::pretty_breaks()) +
          theme_ima("default")
        save_figure(
          plot = p_comp,
          path = file.path(folder, "03-inference/strata", sm_sub_dir),
          filename = paste0(
            jags_obj_prefix,
            paste0("-stratum-means-", component_description, sms_dash), #' -stratum-new-obs',
            ".png"
          ),
          p_width = 4, p_height = 5, device = NULL
        )
        # stratum_mean_mu
      }
    })
  })
}
