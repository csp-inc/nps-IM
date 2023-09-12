library(glue)

get_stratum_inference_no_cis <- function(mcarray_list, data, folder, resp_var_desc,
                                         timesteps, jags_obj_prefix, likelihood,
                                         wrap_char_n = 15, x_hat_raw = timesteps,
                                         ...) {
  print("Inference on strata (no cis)...")
  if (!sprintf("%s.strat.mean", jags_obj_prefix) %in% names(mcarray_list)) {
    return(NULL)
  }

  show_preds_against_data <- show_data_in_output(mcarray_list)
  if (likelihood == "hurdle-ordinal-latent-beta") show_preds_against_data <- FALSE
  # Preliminaries.
  # strat_mean_suffixes <- if(jags_obj_prefix == 'hat') {
  #   c('', '.oos')
  # } else {
  #   ''
  # }
  oos_entry <- ifelse(any(grepl(
    sprintf("%s.strat.*oos$", jags_obj_prefix),
    names(mcarray_list)
  )), ".oos", NA)
  strat_mean_suffixes <- na.omit(c("", oos_entry))
  in_and_out <- lapply(strat_mean_suffixes, function(strat_mean_suffix) {
    sms_dash <- sub("\\.", "-", strat_mean_suffix)
    sm_sub_dir <- ifelse(grepl("oos", strat_mean_suffix), "all", "sampled")
    is_beta_binomial <- grepl("beta-binomial", likelihood)
    jags_obj_suffix <- "mean" # ifelse(is_beta_binomial, 'p', 'mean')
    stratum_means_obj <- paste0(
      paste(jags_obj_prefix, "strat", jags_obj_suffix, sep = "."), strat_mean_suffix
    )
    stratum_lookup <- get_stratum_lookup(data)
    stratum_mu_full <- concat_chains(mcarray_list[[stratum_means_obj]], axis = 4)

    year_vec <- data %>%
      get_year_lookup(timesteps) %>%
      pull(year)

    new_obs_obj <- paste(jags_obj_prefix, "strat.new.obs", sep = ".")

    show_preds_against_data <-
      if (new_obs_obj %in% names(mcarray_list) & !grepl("oos", strat_mean_suffix)) {
        show_preds_against_data
      } else {
        FALSE
      }

    x_index <- 0
    n_strata <- dim(mcarray_list[[stratum_means_obj]])[2]
    stratum_mu_subset <- stratum_mu_full %>% thin_mcmc_output(axis = 3, ...)
    stratum_mu_draws <- lapply(seq_len(n_strata), function(x) {
      x_index <<- x_index + 1

      this_stratum_mu_subset <- if (n_strata == 1) {
        stratum_mu_subset
      } else {
        stratum_mu_subset[, x, ]
      }
      this_stratum_mu_subset %>%
        as_tibble(.name_repair = v_digit) %>%
        mutate(year = year_vec, stratum_index = x_index) %>%
        left_join(stratum_lookup) %>%
        gather(k, stratum_mu, -year, -matches("stratum"))
    }) %>% bind_rows()

    stratum_mean_mu <- summary(mcarray_list[[stratum_means_obj]], mean)$stat %>%
      as_tibble(.name_repair = NULL) %>%
      mutate(year = year_vec) %>%
      gather(stratum_col, stratum_mean_mu, -year) %>%
      left_join(stratum_lookup)

    fncols <- function(data, cname) {
      # https://stackoverflow.com/questions/45857787/adding-column-if-it-does-not-exist
      add <- cname[!cname %in% names(data)]

      if (length(add) != 0) data[add] <- 1
      data
    }

    stratum_data <- data %>%
      select(
        year = cal_year, stratum_id, outcome = matches("^response$|^hits$"),
        normalizer = matches("trials")
      ) %>%
      fncols("normalizer") %>%
      mutate(outcome = outcome / normalizer)

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
    lapply(show_preds_against_data, function(with_data) {
      p_strat_means <- ggplot() +
        facet_grid(. ~ stratum_id, labeller = label_wrap_gen(wrap_char_n)) +
        geom_line(
          data = stratum_mu_draws, aes(x = year, y = stratum_mu, group = k),
          color = "blue", alpha = .01
        ) +
        geom_line(
          data = stratum_mean_mu, aes(x = year, y = stratum_mean_mu),
          color = ima_post_col("default")[2], size = 1.01
        ) +
        # Add layer with actual data.
        {
          if (with_data) {
            geom_point(
              data = stratum_data, aes(x = year, y = outcome),
              color = ima_post_col("default")[2], alpha = .25,
              position = position_jitter(
                width = diff(range(stratum_data$year)) * .01,
                height = diff(range(stratum_data$outcome, na.rm = TRUE)) * .01
              )
            ) # width = 0.1, height = 0.1
          } else {
            NULL
          }
        } +
        labs(x = "Year", y = get_y_lab(with_data)) +
        scale_x_continuous(breaks = scales::pretty_breaks()) +
        theme_ima("default")
      save_figure(
        plot = p_strat_means,
        path = file.path(folder, "03-inference/strata", sm_sub_dir),
        filename = paste0(
          jags_obj_prefix, "-stratum-means",
          ifelse(with_data, "-with-data", ""),
          sms_dash, "-no-cis.png"
        ),
        p_width = 4, p_height = 5, device = NULL
      )
    })
  })
}


y_rep_obj_name <- function(likelihood) {
  "y.rep"
}

rnd <- function(x, digits = n_sig_figs) {
  format(round(x, digits), nsmall = digits)
}

get_op_plots_new <- function(mcarray_list, data, folder,
                             interval = 0.95, method = c("hdi", "eti")) {
  if (!any(grepl("^y.rep", names(mcarray_list)))) {
    return(NULL)
  }
  # output_dir <- file.path(folder, '02-checking/op')
  # dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  y_rep_objs <- grep("^y.rep", names(mcarray_list), value = TRUE)

  lapply(y_rep_objs, function(y_rep_obj) {
    print(
      glue::glue(
        "Observed vs. predicted plots for the `{y_rep_obj}` JAGS object...."
      )
    )
    y_rep_name <- gsub("\\.", "-", y_rep_obj)

    y_rep_posts <- concat_chains(mcarray_list[[y_rep_obj]], 3)
    d_y_obs <- data %>% select(outcome = matches("hits|response"))

    if (y_rep_obj == "y.rep") {
      y_obs <- d_y_obs %>% pull(1) # the usual case
    } else if (y_rep_obj == "y.rep.bern") { # hurdle case
      y_obs <-
        d_y_obs %>%
        mutate(outcome = ifelse(outcome > 0, 0, 1)) %>%
        pull(1)
      y_rep_stats <- apply(y_rep_posts, 2, function(y_rep, y_obs) {
        iter_xtab <- table(y_obs = y_obs, y_rep = y_rep)
        iter_xtab / sum(iter_xtab)
      }, y_obs = y_obs)
      xtab_summary_raw <- apply(
        y_rep_stats, 1, get_interval2,
        interval = interval, method = method[1], include_median = TRUE
      )
      xtab_summary <- rnd(xtab_summary_raw, 3)
      xtab_vec <- glue("{xtab_summary['median', ]} ({xtab_summary['lower', ]}, {xtab_summary['upper', ]})")
      xtab_mat <- matrix(xtab_vec, ncol = 2)
      dimnames(xtab_mat) <- list(y_obs = c(0, 1), y_rep = c(0, 1))
      save_stdout(
        xtab_mat, file.path(folder, "02-checking/op"),
        glue("err-type-xtab-{y_rep_name}.txt")
      )
      return(NULL)
    } else if (y_rep_obj == "y.rep.beta") { # non-zero 0-1 outcomes
      y_obs <- d_y_obs %>%
        filter(outcome > 0) %>%
        pull(1)
    } else if (y_rep_obj == "y.rep.cond") {
      return(NULL)
    }

    y_rep_stats <- apply(
      y_rep_posts, 1, get_interval2,
      interval = interval, method = method[1], include_median = TRUE
    )

    d_op <- t(y_rep_stats) %>%
      as_tibble() %>%
      mutate(y = y_obs) %>%
      drop_na(y) # handles censoring
    if (nrow(d_op) > 1000) {
      message("Subsetting observations in OP plot to avoid overplotting")
      d_op <- d_op %>% sample_n(1000)
    }


    rmse <- sqrt(sum(d_op$y - d_op$median) / nrow(d_op))

    jitter_wh <- ifelse(all(d_op$y %% 1 == 0), 0.2, 0)
    pj <- position_jitter(width = jitter_wh, height = jitter_wh)

    p_op <- ggplot(d_op) +
      geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
      geom_segment(aes(x = lower, xend = upper, y = y, yend = y),
        position = pj, alpha = 0.25
      ) +
      geom_point(aes(x = median, y = y),
        pch = 21, size = 2, fill = "white",
        position = pj, alpha = 0.5
      ) +
      geom_label(
        x = Inf, y = -Inf, hjust = "inward", vjust = "inward",
        label = sprintf("RMSE: %s", round(rmse, 3))
      ) +
      labs(x = "Predicted", y = "Observed") +
      theme_ima("default", rotate_y_lab = FALSE)

    save_figure(
      plot = p_op,
      path = file.path(folder, "02-checking/op"),
      filename = glue("op-{y_rep_name}.jpg"),
      p_width = 4, p_height = 3
    )
  })
}
