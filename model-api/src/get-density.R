get_density <- function(mcarray_list, data, deterministic_model, likelihood,
                        folder, thresh = 100) {
  if (!likelihood %in% c("lognormal", "gamma") | data$network_code[1] == "MISC") {
    return(NULL)
  }

  h <- function(det_fun) {
    function(x) {
      if (grepl("exponential", det_fun)) {
        exp(x)
      } else if (grepl("restricted-linear", det_fun)) {
        pmax(x, .00001)
        # sapply(x, function(z) max(0.00001, z), simplify = FALSE)
      } else {
        x
      }
    }
  }

  intercept_samples <- if ("B.tilde" %in% names(mcarray_list)) {
    mcarray_list$B.tilde
  } else {
    mcarray_list$B0.tilde
  }
  # mu <- summary(h(deterministic_model)(intercept_samples),
  #               mean, na.rm=TRUE)$stat
  # if (!is.null(dim(mu))) mu <- mu[1, ]

  variance_samples <- if ("sigma.tilde" %in% names(mcarray_list)) {
    mcarray_list$sigma.tilde
  } else {
    mcarray_list$sigma
  }
  # sigma <- summary(variance_samples, mean, na.rm=TRUE)$stat

  stratum_lookup <- get_stratum_lookup(data)


  # max samples thing
  y_rep_maxima <- apply(mcarray_list$y.rep, c(2, 3), function(w) {
    # note: ignoring years I guess....
    tibble(y_rep = w, stratum_id = data$stratum_id, stratum_index = data$stratum_index) %>%
      group_by(stratum_id, stratum_index) %>%
      summarise(max_y_rep = max(y_rep)) %>%
      ungroup()
  })
  save_table(bind_rows(y_rep_maxima), file.path(folder, "99-misc"), "y-rep-maxima.csv")

  y_moments_fun <- function(this_jags_obj, desc = NA) {
    apply(this_jags_obj, c(2, 3), function(w) {
      tibble(
        stat = w,
        stratum_id = unique(data$stratum_id), stratum_index = unique(data$stratum_index)
      ) %>%
        mutate(moment = desc)
    }) %>% bind_rows()
  }

  y_moments <- bind_rows(
    y_moments_fun(h(deterministic_model)(intercept_samples), desc = "mu"),
    y_moments_fun(variance_samples, desc = "sigma")
  )
  save_table(y_moments, file.path(folder, "99-misc"), "y-moments.csv")
  # h(deterministic_model)(intercept_samples)

  moments_array <- array(NA, dim = c(2, dim(variance_samples)))
  moments_array[1, , , ] <- h(deterministic_model)(intercept_samples)
  moments_array[2, , , ] <- variance_samples

  these_den_results <- apply(moments_array, MARGIN = c(3, 4), function(w) {
    d_tmp <- mapply(function(mu_k, sigma_k, k) {
      if (likelihood == "lognormal") {
        mu_k_log_y <- log(mu_k) - 1 / 2 * log((sigma_k^2 + mu_k^2) / mu_k^2)
        sigma_k_log_y <- sqrt(log((sigma_k^2 + mu_k^2) / mu_k^2))
        x_range <- qlnorm(c(0.005, 0.995), mu_k_log_y, sigma_k_log_y)

        if (!is.finite(x_range[1])) browser()
        x <- seq(floor(x_range[1]), ceiling(x_range[2]))
        d <- dlnorm(x, mu_k_log_y, sigma_k_log_y)
        p_gt_thresh <- 1 - plnorm(q = thresh, mu_k_log_y, sigma_k_log_y)
      } else if (likelihood == "gamma") {
        alpha_k <- mu_k^2 / sigma_k^2
        beta_k <- mu_k / sigma_k^2
        x_range <- qgamma(c(0.005, 0.995), alpha_k, beta_k)
        x <- seq(floor(x_range[1]), ceiling(x_range[2]))
        d <- dgamma(x, alpha_k, beta_k)
        p_gt_thresh <- 1 - pgamma(q = thresh, alpha_k, beta_k)
      }

      list(
        tibble(x = x, d = d, p_gt_thresh) %>%
          mutate(mu = mu_k, sigma = sigma_k, stratum_index = k) %>%
          left_join(stratum_lookup %>% select(stratum_id, stratum_index))
      )
    }, mu_k = w[1, ], sigma_k = w[2, ], k = seq_along(w[1, ]))
    bind_rows(d_tmp)
  })
  # names(these_den_results) <- 1:length(these_den_results)
  testing123 <- bind_rows(these_den_results, .id = "m")

  # d_tmp_comb <- these_den_results[[1]]  # bind_rows(d_tmp)
  # d_tmp_labs <- d_tmp_comb %>%
  #   distinct(mu, sigma, p_gt_thresh, stratum_id) %>%
  #   mutate(moments = sprintf('mu==%s*","~sigma == %s',
  #                            round(mu, 2), round(sigma, 2)),
  #          p_gt_thresh = sprintf('P(italic(y)>%s)==%s',
  #                                thresh, round(p_gt_thresh, 2)))

  # ggplot(testing123) +
  #   facet_wrap(~stratum_id) +
  #   # geom_label(data = d_tmp_labs, aes(x = Inf, y = Inf, label = moments),
  #   #            hjust = 'inward', vjust = 'inward', parse = TRUE) +
  #   # geom_ribbon(aes(x = x, ymin = 0, ymax = d, group = m),
  #   #             fill = 'gray', color = 'black', alpha = 0.2) +
  #   geom_line(aes(x = x, y = d, group = m),
  #               color = 'black', alpha = 0.01) +
  #   # geom_label(data = d_tmp_labs, aes(x = Inf, y = -Inf, label = p_gt_thresh),
  #   #            hjust = 'inward', vjust = 'inward', parse = TRUE, fill = 'white') +
  #   ggthemes::theme_hc()
  # ggsave(sprintf('%s/99-misc/basic-density.png', folder),
  #        width = 6 * 1.5, height = 4 * 1)

  u_testing123 <- testing123 %>% distinct(m, p_gt_thresh, mu, sigma, stratum_index, stratum_id)
  save_table(testing123, file.path(folder, "99-misc"), "basic-density.csv")
}
