library(ggridges)
source("model-api/src/fit-utils.r")
source("model-api/src/model-eval.R")

check_packages <- function(packages) {
  lapply(
    packages,
    FUN = function(x) {
      if (!x %in% rownames(installed.packages())) {
        install.packages(x, dependencies = TRUE)
        # library(x, character.only = TRUE)
      }
    }
  )
}

vcat <- function(names) {
  sprintf("V%s", seq.int(length(names)))
}
increment <- function(x, scale = 10) {
  10^floor(log10(abs(x))) / scale * ifelse(x < 0, -1, 1)
}

fig_dims <- function(which_dim, ratio = 2 / 3, # e.g., 4/6, 5/7
                     n_facet_x = 1, n_facet_y = 1,
                     in_per_x_facet = 3, in_per_y_facet = 3) {
  if (which_dim == "x") {
    w <- in_per_x_facet * n_facet_x
    return(c(width = w, height = w / ratio))
  } else if (which_dim == "y") {
    h <- in_per_y_facet * n_facet_y
    return(c(width = ratio * h, height = h))
  }
}

whip_y_sim_into_shape <- function(mcarray_list, data, x_hat_raw,
                                  obj = "hat.site.new.obs",
                                  ...) {
  index_shift <- ifelse(grepl("site", obj), 0, -1)


  n_strata <- dim(mcarray_list[[obj]])[2]
  if (n_strata == 1) {
    y_sim_tmp <- concat_chains(mcarray_list[[obj]][, 1, , ], axis = 4 + index_shift)
    y_sim_array <- array(NA, dim = c(nrow(y_sim_tmp), n_strata, ncol(y_sim_tmp)))
    y_sim_array[, 1, ] <- y_sim_tmp
  } else {
    y_sim_array <- concat_chains(mcarray_list[[obj]], axis = 5 + index_shift)
  }


  y_sim_stratum_cis_array <- apply(y_sim_array, MARGIN = c(2, 3) + index_shift, function(i) {
    quantile(i, probs = c(.025, .975), na.rm = TRUE)
  })
  x_index <- 0
  y_sim_stratum_cis_df <- apply(y_sim_stratum_cis_array, MARGIN = 1, function(x) {
    x_index <<- x_index + 1
    x %>%
      as_tibble(.name_repair = NULL) %>%
      mutate(rel_year = x_hat_raw, quantile_index = x_index) %>%
      gather(stratum_col, quantile, -rel_year, -quantile_index) %>%
      left_join(get_stratum_lookup(data)) %>%
      left_join(get_year_lookup(data, x_hat_raw)) %>%
      select(-matches("col"))
  }) %>%
    reduce(rbind) %>%
    left_join(tibble(quantile_index = c(1, 2), quantile_lab = c("0.025", "0.975")))
  # ggplot(y_sim_stratum_cis_df, aes(x=year, y=quantile, group=quantile_lab)) +
  #   facet_wrap(~stratum_id) +
  #   geom_line()
  y_sim_stratum_cis_df
}

get_site_inference <- function(mcarray_list, data, folder, resp_var_desc,
                               timesteps, jags_obj_prefix, likelihood,
                               wrap_char_n = 15, ...) {
  print("Inference on sites...")

  # Preliminaries.
  stratum_lookup <- get_stratum_lookup(data)
  site_lookup <- get_site_lookup(data)
  year_lookup <- get_year_lookup(data, timesteps)
  is_beta_binomial <- grepl("beta-binomial", likelihood)
  jags_obj_suffix <- "mean" # ifelse(is_beta_binomial, 'p', 'mean')
  site_means_obj <- paste(jags_obj_prefix, "site", jags_obj_suffix, sep = ".")
  stratum_means_obj <- paste(jags_obj_prefix, "strat", jags_obj_suffix, sep = ".")

  site_mu_array <-
    summary(mcarray_list[[site_means_obj]], mean, na.rm = TRUE)$stat

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

  stratum_mean <- if (stratum_means_obj %in% names(mcarray_list)) {
    summary(mcarray_list[[stratum_means_obj]], mean, na.rm = TRUE)$stat %>%
      as_tibble(.name_repair = NULL) %>%
      mutate(year = year_lookup %>% pull(year)) %>%
      gather(stratum_col, stratum_mean, -year) %>%
      left_join(stratum_lookup) %>%
      select(-matches("index|col"))
  } else {
    NULL
  }

  # Save plot components for future plotting.
  list(site_stats_df = site_stats_df, stratum_mean = stratum_mean) %T>%
    save_object(
      file.path(folder, "03-inference/site"),
      paste(jags_obj_prefix, "site-in-stratum-means.rds", sep = "-")
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
      paste(jags_obj_prefix,
        "site-in-stratum-means.csv",
        sep = "-"
      )
    )

  # Plotting.
  y_lev <- "Site-level"
  y_lab <- as.expression(bquote(atop(
    .(y_lev),
    .(tolower(resp_var_desc)) ~ (mu)
  )))
  # png_dims <- fig_dims('x', ratio = (n_distinct(data$stratum_id)/3 * 16)/9,
  #                      n_facet_x = n_distinct(data$stratum_id),
  #                      n_facet_y = 1,
  #                      in_per_x_facet = 4)

  p <- ggplot() +
    facet_grid(. ~ stratum_id, labeller = label_wrap_gen(wrap_char_n)) + # , scales = 'free'
    # facet_wrap(~stratum_id, labeller = label_wrap_gen(wrap_char_n)) +
    geom_line(
      data = site_stats_df, aes(x = year, y = mean, group = site_id),
      color = ima_pal("default")[2]
    ) +
    labs(x = "Year", y = y_lab) +
    scale_x_continuous(breaks = scales::pretty_breaks()) +
    theme_ima("default")
  if (!is.null(stratum_mean)) {
    p <- p +
      geom_line(
        data = stratum_mean, aes(x = year, y = stratum_mean),
        color = ima_pal("default")[1], size = 1.01
      )
  }
  save_figure(
    plot = p,
    path = file.path(folder, "03-inference/site"),
    filename = paste(jags_obj_prefix,
      "site-in-stratum-means.png",
      sep = "-"
    ),
    p_width = 4, p_height = 5, device = NULL
  )

  if (likelihood != "hurdle-ordinal-latent-beta") {
    lapply(unique(site_stats_df$stratum_id), function(s) {
      this_site_stats_df <- site_stats_df %>%
        filter(stratum_id == s)
      # sites_hl <- this_site_stats_df %>%
      #   filter(rel_year==0) %>%
      #   arrange(desc(mean)) %>%
      #   pull(site_id)
      this_data <- data %>%
        filter(stratum_id == s) %>%
        {
          `if`("midpoint" %in% names(.), mutate(., response = midpoint), .)
        } %>%
        {
          `if`("hits" %in% names(.), mutate(., response = hits / trials), .)
        }
      jitter_frac <- .01
      response_minmax_diff <- diff(range(this_data$response, na.rm = TRUE))
      jitter_scale <- ifelse(response_minmax_diff == 0, .001, response_minmax_diff * jitter_frac)
      jitter_scale <- ifelse(grepl("binomial|beta", likelihood), 0, jitter_scale)
      resp_col <- ifelse("class_midpoint" %in% names(this_data),
        "class_midpoint", "response"
      )

      n_wrap_rows <- ceiling(n_distinct(this_site_stats_df$site_id) / 3)

      p_site_means <- ggplot(this_site_stats_df) +
        facet_wrap(~site_id, scales = "free", ncol = 3) +
        geom_ribbon(
          data = this_site_stats_df, aes(x = year, ymin = `0.025`, ymax = `0.975`),
          fill = ima_pal("default")[1], alpha = .5
        ) +
        geom_line(data = this_site_stats_df, aes(x = year, y = mean)) +
        geom_point(
          data = this_data, aes_string(x = "cal_year", y = resp_col),
          alpha = .25,
          position = position_jitter(height = jitter_scale, width = .1)
        ) +
        labs(x = "Year", y = y_lab) +
        theme_ima("default")
      save_figure(
        plot = p_site_means,
        path = file.path(folder, "03-inference/site"),
        filename = paste0(
          jags_obj_prefix,
          "-site-means-(", s, ").png"
        ),
        p_width = 2, p_height = 1.5, device = NULL
      )
    })
  }


  # summary(mcarray_list$B, mean, na.rm=TRUE)$stat[, , 1]
}

cache_inits <- function(mcarray_list, n_chains, folder, chain_z,
                        gelman_diag, inits_cache) {
  inits_cache <- file.path(inits_cache, "inits")

  dir.create(inits_cache, showWarnings = FALSE, recursive = TRUE) # probs not needed
  inits_csv <- file.path(inits_cache, "inits-file-lookup.csv")
  if (!file.exists(inits_csv)) {
    write_csv(
      tibble(model = NA, inits_file = NA, gelman_diag = NA) %>% slice(-1),
      inits_csv
    )
  }
  # Sys.sleep(sample.int(10, 1))  # temporary
  invoke_os_command(paste(
    "./model-api/src/model-builder/get-inits-needed.sh",
    folder
  ))
  inits_names_raw <- names(mcarray_list)[
    names(mcarray_list) %in%
      str_trim(readLines(file.path(folder, "inits-needed.txt")))
  ]
  inits_names <- grep("^[^y|^hat|^pred|^j.draw]", inits_names_raw, value = TRUE)
  param_quants <- lapply(inits_names, function(this_init) {
    this_init_mat <- summary(mcarray_list[[this_init]], quantile,
      probs = c(0.025, .5, .975), na.rm = TRUE
    )$stat
    dimnames(this_init_mat)[[1]] <- NULL
    if (any(is.na(this_init_mat))) {
      return(NULL)
    } else {
      this_init_list <- lapply(1:nrow(this_init_mat), function(x) {
        if (length(dim(this_init_mat)) > 2) {
          return(NULL)
        }
        this_init_mat[x, ]
      })
      names(this_init_list) <- rep(this_init, length(this_init_list))
      return(this_init_list)
    }
  })
  nonnull_param_quants <- param_quants[lengths(param_quants) != 0]
  inits <- list()
  for (i in 1:n_chains) {
    inits[[i]] <- list()
    for (j in 1:length(nonnull_param_quants)) {
      inits[[i]] <- c(inits[[i]], nonnull_param_quants[[j]][i])
    }
  }

  if (any(grepl("^z$", inits_names_raw))) {
    inits <-
      lapply(inits, function(x) {
        x$z <- rep(chain_z, length(x$z))
        x
      })
  }
  inits <-
    lapply(inits, function(x) {
      x[sapply(x, function(i) !is.null(i))]
    })

  if (!is.na(gelman_diag) & gelman_diag < 20) {
    inits_lookup_file <- file.path(inits_cache, "inits-file-lookup.csv")
    last_gelman_diag <-
      read_csv(inits_lookup_file) %>%
      filter(model == folder) %>%
      pull(gelman_diag)
    convergence_benchmark <- ifelse(length(last_gelman_diag) == 0, 20, last_gelman_diag)

    lookup_entry <- tibble(
      model = folder,
      inits_file = paste0(digest::digest(folder), ".rds"),
      gelman_diag = gelman_diag
    )

    if (gelman_diag <= convergence_benchmark) {
      inits_lookup <- lookup_entry %>%
        bind_rows(read_csv(inits_lookup_file, col_types = "ccd")) %>%
        group_by(model, inits_file) %>%
        filter(gelman_diag == min(gelman_diag)) %>%
        distinct() %>%
        ungroup() %T>%
        write_csv(file.path(inits_cache, "inits-file-lookup.csv"))
      saveRDS(inits, file.path(
        inits_cache,
        paste0(digest::digest(folder), ".rds")
      ))
    }

    # inits_lookup %>% filter(model == folder) %>% pull(gelman_diag)
  }


  inits
}

get_coef_post_dists <- function(mcarray_list, data, folder, resp_var_desc,
                                deterministic_model, data_raw = NULL) {
  print("Inference on coefficients...")

  # grep('^B0$|^B1$|^Beta$|^mu.B0$', names(mcarray_list), value = TRUE)
  transform_coef <- deterministic_model %in% c("exponential", "inverse-logit")
  x_label <- "Coefficient estimate"

  tform_scenarios <- if (transform_coef) {
    "untransformed" # c('transformed', 'untransformed')
  } else {
    "untransformed"
  }
  lapply(tform_scenarios, function(x) {
    if (x == "transformed") x_label <- paste(x_label, "(exponentiated)")


    coef_posts_to_df <- function(arr, data, stat_val, stat_lab) {
      arr %>%
        as_tibble(.name_repair = v_digit) %>%
        mutate(site_index = 1:n()) %>%
        gather(stratum_col, !!stat_val, -site_index) %>%
        left_join(get_stratum_lookup(data)) %>%
        mutate(
          stat_lab = stat_lab,
          stratum_id = gsub(" ", "_", stratum_id)
        ) %>%
        select(-stratum_col)
    }


    coef_arrays <- grep("^B$|^G$", names(mcarray_list), value = TRUE)
    coef_labs <- c("beta", "gamma") # str_to_title(coef_labs)
    lapply(seq_along(coef_arrays), function(i_coef) {
      coef_posts <- mcarray_list[[coef_arrays[i_coef]]]
      if (x == "transformed") coef_posts <- exp(coef_posts)
      coef_stats <- summary(coef_posts, quantile,
        probs = c(.025, .5, .975), na.rm = TRUE
      )$stat
      coef_stats_df <- lapply(seq_len(3), function(x) {
        # B0: B[j, 1, k] and B1: B0: B[j, 2, k]
        b0 <- coef_posts_to_df(
          coef_stats[x, , 1, ], data,
          paste0(coef_labs[i_coef], "[0]"),
          dimnames(coef_stats)[[1]][x]
        )
        if (dim(coef_stats)[3] == 2) {
          b1 <- coef_posts_to_df(
            coef_stats[x, , 2, ], data,
            paste0(coef_labs[i_coef], "[1]"),
            dimnames(coef_stats)[[1]][x]
          )
          left_join(b0, b1)
        } else {
          b0
        }
      }) %>%
        bind_rows() %>%
        gather(coef, stat_val, -site_index, -stratum_id, -stratum_index, -stat_lab) %>%
        spread(stat_lab, stat_val)

      coef_regex <- paste(paste0("^mu.", coef_arrays[i_coef], c(0, 1), "$"),
        collapse = "|"
      )
      if (any(grepl(coef_regex, names(mcarray_list)))) {
        these_hyperparams <-
          grep(coef_regex, names(mcarray_list), value = TRUE)
        hyperparams_stats <- lapply(
          seq_along(these_hyperparams),
          function(i_param) {
            hyperdist_posts <- mcarray_list[[these_hyperparams[i_param]]]
            if (x == "transformed") hyperdist_posts <- exp(hyperdist_posts)
            hyperdist_stats <-
              summary(hyperdist_posts,
                quantile,
                probs = c(.025, .5, .975)
              )$stat

            hyperdist_stats %>%
              as_tibble(.name_repair = NULL) %>%
              mutate(stat_lab = dimnames(hyperdist_stats)[[1]]) %>%
              gather(stratum_col, stat_val, -stat_lab) %>%
              left_join(get_stratum_lookup(data)) %>%
              mutate(
                stratum_id = gsub(" ", "_", stratum_id),
                l_type = ifelse(stat_lab == "50%", "solid", "dashed"),
                coef = paste0(coef_labs[i_coef], "[", i_param - 1, "]")
              )
          }
        ) %>% bind_rows()
      }



      # if (FALSE) {
      #   data_tmp <- data %>%
      #     unite('pj', matches("point_jurisdiction")) %>%
      #     select(site_index, stratum_index, pj) %>%
      #     distinct()
      # }
      # coef_stats_df %>% left_join(zzz)
      # d_cov <- read_csv('data/ROMN/GRKO/GRKO_covariates_20200402.csv')
      # d_cov %>% group_by(MDCATY) %>%
      #   summarise(n_riparian = sum())
      # d_cov_erin <- d_cov %>%
      #   group_by(site_id = SiteName) %>%
      #   summarise(MgmtZone = names(which.max(table(MgmtZone)))) %>%
      #   ungroup()
      # tmp <- coef_stats_df %>% left_join(
      #   data %>%
      #     distinct(site_index, site_id, stratum_id, site_in_stratum_index) %>%
      #     left_join(d_cov_erin) %>%
      #     select(site_index = site_in_stratum_index, stratum_id, MgmtZone)
      # )

      p <- ggplot(coef_stats_df) +
        # ggplot(coef_stats_df %>% left_join(data_tmp)) +
        facet_grid(stratum_id ~ coef, scales = "free", labeller = label_parsed) +
        geom_vline(xintercept = ifelse(x == "transformed", 1, 0), color = "red") +
        scale_y_reverse() +
        geom_errorbarh(aes(xmin = `2.5%`, xmax = `97.5%`, y = site_index)) + # , color = MgmtZone
        geom_point(aes(x = `50%`, y = site_index)) +
        theme_ima("default") +
        scale_linetype_identity() +
        labs(x = x_label, y = "Site index")

      if (exists("hyperparams_stats")) {
        p <- p +
          geom_vline(
            data = hyperparams_stats, aes(xintercept = stat_val, linetype = l_type),
            color = "#7cb5ec"
          )
      }
      save_figure(
        plot = p,
        path = file.path(folder, "03-inference/coef"),
        filename = paste(coef_labs[i_coef], "coef-estimates",
          paste0(x, ".png"),
          sep = "-"
        ),
        p_width = 4, p_height = 5, device = NULL
      )
      # mcarray_list$B
      # coef_posts <- concat_chains(mcarray_list$B, axis=5)

      coef_tc <- str_to_title(coef_labs[i_coef]) # title case
      if (any(grepl("^Beta.hat$", names(mcarray_list)))) {
        coef_tc <- paste(coef_tc, "hat", sep = ".")
      }
      if (any(grepl(paste0("^", coef_tc, "$"), names(mcarray_list)))) {
        # dim(mcarray_list[['Beta']])[1]>1
        Beta_posts <- mcarray_list[[coef_tc]]
        if (x == "transformed") Beta_posts <- exp(Beta_posts)
        Beta_posts_mat <- concat_chains(Beta_posts, axis = 3, drop = FALSE)[, , 1]

        Beta_posts_df <- if (dim(Beta_posts)[1] > 1) {
          as_tibble(t(Beta_posts_mat), .name_repair = NULL)
        } else {
          as_tibble(Beta_posts_mat, .name_repair = "minimal")
        }

        b_name <- paste0(
          coef_labs[i_coef], "[",
          sub("fe_", "", grep("^fe_*", names(data), value = TRUE)),
          "]"
        )
        names(Beta_posts_df) <- b_name
        p_add_coefs <- ggplot(Beta_posts_df %>% gather(key, val)) +
          # facet_grid(key~., scales = 'free', labeller=label_parsed) +
          facet_wrap(~key, scales = "free", labeller = label_parsed) +
          geom_histogram(aes(x = val, y = ..density..),
            color = "black", fill = "white"
          ) +
          geom_vline(xintercept = ifelse(x == "transformed", 1, 0), color = "orange") +
          theme_ima("default") +
          labs(x = x_label, y = "Density")
        save_figure(
          plot = p_add_coefs,
          path = file.path(folder, "03-inference/coef"),
          filename = paste0("additional-coef-estimates-", coef_tc, "-", x, ".png"),
          p_width = 4 * 0.8, p_height = 3 * 0.8, device = NULL
        )

        intx_effects_list <- attr(data, "intx_effects_list")
        if (!is_empty(intx_effects_list)) { # !is.null(intx_effects_list)

          names(Beta_posts_df) <- grep("^fe_*", names(data), value = TRUE)

          tmp <- lapply(intx_effects_list, function(z_outer) {
            lapply(z_outer, function(z_inner) {
              rowSums(Beta_posts_df %>% select(any_of(z_inner)))
            })
          })

          p_intx_effects <- ggplot(
            enframe(tmp) %>%
              unnest_longer(value) %>% unnest_longer(value)
          ) +
            facet_grid(value_id ~ name, scales = "free") +
            geom_histogram(aes(x = value, y = ..density..),
              color = "black", fill = "white"
            ) +
            geom_vline(xintercept = ifelse(x == "transformed", 1, 0), color = "orange") +
            theme_ima("default") +
            labs(x = x_label, y = "Density")

          save_figure(
            plot = p_intx_effects,
            path = file.path(folder, "03-inference/coef"),
            filename = paste0("derived-intx-coefs", coef_tc, "-", x, ".png"),
            p_width = 4 * 0.8, p_height = 3 * 0.8, device = NULL
          )
        }
      }
    })
  })
}

get_random_hyperparams <- function(mcarray_list, data, folder,
                                   ...) {
  print("Inference on random effects hyperparameters...")
  hypers <- grep("mu.B0$|mu.B1$|mu.G0$|mu.G1$", names(mcarray_list), value = TRUE)

  lapply(sapply(hypers, function(x) strsplit(x, "\\.")[[1]][2]), function(hyper) {
    get_mat <- function(par_name) {
      out <- concat_chains(mcarray_list[[par_name]], axis = 3, drop = FALSE)[, , 1]
      if (dim(mcarray_list[[par_name]])[1] > 1) {
        t(out)
      } else {
        matrix(out, ncol = 1)
      }
    }

    mu_name <- paste("mu", hyper, sep = ".")
    mu_hyper <- as_tibble(get_mat(mu_name), .name_repair = vcat) %>%
      mutate(iteration = 1:n()) %>%
      gather(stratum_col, !!mu_name, -iteration) %>%
      left_join(get_stratum_lookup(data))

    sigma_name <- paste("sigma", hyper, sep = ".")
    sigma_hyper <- if (sigma_name %in% names(mcarray_list)) {
      as_tibble(get_mat(sigma_name), .name_repair = vcat) %>%
        mutate(iteration = 1:n()) %>%
        gather(stratum_col, !!sigma_name, -iteration) %>%
        left_join(get_stratum_lookup(data))
    } else { # triggered if the intercept term for the hurdle is fixed
      get_stratum_lookup(data) %>% mutate(!!sigma_name := NA)
    }

    # p_pairs <- ggplot(left_join(sigma_hyper, mu_hyper)) +
    #   facet_wrap(~stratum_id) +
    #   geom_point(aes_string(x = mu_name, y = sigma_name), alpha = .1)
    # save_figure(
    #   plot = p_pairs, path = file.path(folder, "03-inference/coef"),
    #   filename = sprintf("pairs-%s.png", hyper),
    #   p_width = 4, p_height = 5
    # )

    both_params <- left_join(sigma_hyper, mu_hyper) %>%
      select(iteration, stratum_id, all_of(sigma_name), all_of(mu_name)) %>%
      gather(parameter, value, -iteration, -stratum_id)
    p_hypers_hist <- ggplot(both_params) +
      facet_grid(stratum_id ~ parameter, scales = "free") +
      geom_histogram(aes(x = value, y = ..density..), fill = "white", color = "black") +
      theme_ima("default")
    save_figure(
      plot = p_hypers_hist, path = file.path(folder, "03-inference/coef"),
      filename = sprintf("hypers-hist-%s.png", hyper),
      p_width = 8, p_height = 3
    )
  })
}
