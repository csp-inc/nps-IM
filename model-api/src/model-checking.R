get_y_reps <- function(mcarray_list, data, likelihood, folder,
                       facet_by = NULL, ...) {
  message("Model checking...")

  jags_obj <- grep("y.rep", names(mcarray_list), value = TRUE)

  lapply(jags_obj, function(this_jags_obj) {
    hurdle_filter <- ifelse(this_jags_obj == "y.rep.beta",
      "outcome > 0", "outcome == outcome"
    )
    hurdle_mutate <- ifelse(this_jags_obj == "y.rep.bern",
      "ifelse(outcome > 0, 0, 1)", "outcome"
    )

    if ("censor_limit_vec" %in% names(data)) {
      data <- data %>%
        rename(censored_response = response, response = censor_limit_vec)
    }

    d_obs <- data %>%
      select(outcome = matches("^response$|^hits$")) %>%
      mutate(draw = "0", row_idx = 1:n()) %>%
      mutate_(.dots = setNames(hurdle_mutate, "outcome"))
    out <- concat_chains(mcarray_list[[this_jags_obj]], axis = 3) %>%
      as_tibble() %>%
      gather(draw, outcome) %>%
      bind_rows(d_obs) %>%
      mutate(grp = ifelse(draw == "0", "Actual", "Replicate")) %>%
      filter_(hurdle_filter)

    if (this_jags_obj == "y.rep.bern") {
      x_axis <- scale_x_discrete(limit = c(0, 1), labels = c("P", "A"))
      x_axis_bw <- 0.5
    } else {
      x_axis <- scale_x_continuous(breaks = scales::pretty_breaks(n = 3))
      x_axis_bw <- NULL
    }

    default_facets <- c("unit_code", "stratum_id")
    facet_by <- unique(c(facet_by, default_facets))

    keepers <- d_obs %>%
      filter_(hurdle_filter) %>%
      pull(row_idx)
    d_facet <- data %>%
      select(one_of(facet_by)) %>%
      slice(keepers) %>%
      mutate(obs_index = 1:n())

    y_reps_draws_with_fct <- out %>%
      group_by(draw, grp) %>%
      mutate(obs_index = 1:n()) %>%
      ungroup() %>%
      left_join(d_facet, by = "obs_index")

    lapply(facet_by, function(this_fct) {
      message(this_fct)
      # if(this_jags_obj=='y.rep.beta' & this_fct == 'site_id') browser()
      # if (this_jags_obj == 'y.rep.cond') browser()
      y_reps_stats_by_fct <- y_reps_draws_with_fct %>%
        group_by_at(vars("draw", "grp", this_fct)) %>%
        summarise(mu = mean(outcome), sigma = sd(outcome)) %>%
        ungroup() %>%
        gather(stat, value, -one_of("draw", "grp", this_fct))

      bayesian_p_by_fct <- y_reps_stats_by_fct %>%
        filter(grp == "Replicate") %>%
        left_join(y_reps_stats_by_fct %>%
          filter(grp == "Actual") %>%
          select(-draw, -grp) %>%
          rename(act_val = value)) %>%
        group_by_at(vars(this_fct, "stat")) %>%
        summarise(p_val = sum(value > act_val) / n()) %>%
        ungroup() %>%
        spread(stat, p_val)
      save_table(
        bayesian_p_by_fct, file.path(folder, "02-checking/ppc", this_fct),
        sprintf(
          "%s-bayes-p-by-%s.csv",
          gsub("\\.", "-", this_jags_obj), this_fct
        )
      )

      draw_idx_for_rep <- y_reps_draws_with_fct %>%
        filter(grp == "Replicate") %>%
        distinct(draw)
      this_p3 <- ggplot(y_reps_draws_with_fct %>%
        filter(draw %in% c("0", sample(draw_idx_for_rep$draw, 9)))) +
        facet_grid(.data[[this_fct]] ~ draw) +
        geom_histogram(aes(x = outcome, y = ..density.., fill = grp),
          color = "black", binwidth = x_axis_bw
        ) +
        ggthemes::theme_hc(base_size = 16) +
        scale_fill_manual("Group",
          labels = c(
            expression(italic(y)),
            expression(italic(y)^rep)
          ),
          values = c("#7cb5ec", "#434348")
        ) +
        x_axis +
        theme(
          strip.background.x = element_blank(),
          strip.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_text(angle = 0, vjust = 0.5)
        ) +
        labs(x = expression(italic(y)), y = expression("[" * italic(y) * "]"))


      save_figure(
        plot = this_p3,
        path = file.path(folder, "02-checking/ppc", this_fct),
        filename = sprintf(
          "%s-by-%s-9-draws.jpg",
          gsub("\\.", "-", this_jags_obj), this_fct
        ),
        p_width = 1.25, p_height = (3 / 4) * 1.25, device = NULL
      )

      vlines_data <- y_reps_stats_by_fct %>%
        filter(grp == "Actual") %>%
        left_join(bayesian_p_by_fct %>% gather(stat, bayesian_p, -matches(this_fct)))
      this_p4 <- ggplot(y_reps_stats_by_fct %>% filter(grp == "Replicate")) +
        facet_grid(.data[[this_fct]] ~ stat,
          scales = "free",
          labeller = labeller(stat = label_parsed)
        ) +
        geom_histogram(aes(x = value, y = ..density..),
          fill = "#434348",
          color = "black", size = .1
        ) +
        geom_vline(
          data = vlines_data,
          aes(xintercept = value), color = "#7cb5ec", size = 1.1
        ) +
        geom_label(
          data = vlines_data,
          aes(x = value, y = 0, label = format(round(bayesian_p, 2),
            nsmall = 2
          )),
          vjust = "inward"
        ) +
        labs(x = "", y = "Density") +
        ggthemes::theme_hc(base_size = 16) +
        scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
        theme(panel.spacing = unit(1.5, "lines"))

      save_figure(
        plot = this_p4,
        path = file.path(folder, "02-checking/ppc", this_fct),
        filename = sprintf(
          "%s-by-%s-stats.jpg",
          gsub("\\.", "-", this_jags_obj), this_fct
        ),
        p_width = 2.5, p_height = (3 / 4) * 2.5, device = NULL
      )
      # if (this_fct == "cal_year") browser()
      if (this_fct %in% default_facets) {
        n_draws <- n_distinct(y_reps_draws_with_fct$draw) - 1
        p_ridge_yrep <- ggplot(y_reps_draws_with_fct %>%
          filter(draw %in% c("0", sample(
            unique(draw),
            min(1000, n_draws)
          )))) +
          geom_density_ridges(aes(x = outcome, y = .data[[this_fct]], fill = grp),
            alpha = 0.75, scale = 0.95, stat = "binline", bins = 50
          ) +
          ggthemes::theme_hc(base_size = 16) +
          x_axis +
          labs(x = expression(italic(y)^rep), y = "") +
          ggthemes::scale_fill_hc(name = expression(italic(y)))
        save_figure(
          plot = p_ridge_yrep,
          path = file.path(folder, "02-checking/ppc", this_fct),
          filename = sprintf(
            "%s-by-%s-binline.jpg",
            gsub("\\.", "-", this_jags_obj), this_fct
          ),
          p_width = 7, p_height = 5, device = NULL
        )

        p_ridge_yrep_stats <- ggplot(y_reps_stats_by_fct %>% filter(stat == "mu", grp == "Replicate")) +
          geom_density_ridges(aes(x = value, y = NA, fill = .data[[this_fct]]),
            alpha = 0.75, scale = 0.95, stat = "binline", bins = 50
          ) +
          geom_vline(
            data = y_reps_stats_by_fct %>% filter(stat == "mu", grp == "Actual"),
            aes(xintercept = value, color = .data[[this_fct]]), linetype = "dashed",
            size = 0.8
          ) +
          ggthemes::theme_hc(base_size = 16) +
          x_axis +
          labs(x = expression(mean(italic(y)^rep) ~ vs ~ mean(italic(y))), y = "") +
          ggthemes::scale_fill_hc(name = "") +
          ggthemes::scale_color_hc(name = "") +
          theme(axis.text.y = element_blank())
        save_figure(
          plot = p_ridge_yrep_stats,
          path = file.path(folder, "02-checking/ppc", this_fct),
          filename = sprintf(
            "%s-by-%s-binline-stats.jpg",
            gsub("\\.", "-", this_jags_obj), this_fct
          ),
          p_width = 7, p_height = 5, device = NULL
        )
      }


      if (grepl("beta|bern", this_jags_obj)) {
        p_cat <- classwise_ppc(y_reps_draws_with_fct, this_fct, folder)
        save_figure(
          plot = p_cat,
          path = file.path(folder, "02-checking/ppc", this_fct),
          filename = sprintf(
            "%s-by-%s-kruschke.jpg",
            gsub("\\.", "-", this_jags_obj), this_fct
          ),
          p_width = 2.5, p_height = (3 / 4) * 2.5, device = NULL
        )
      }
    })

    NULL
  })
}

classwise_ppc <- function(x, facet, folder) {
  n_obs_by_cat_and_facet <- x %>%
    group_by_at(vars("draw", "grp", facet, "outcome")) %>%
    tally() %>%
    ungroup()

  stats_by_cat_and_facet <- n_obs_by_cat_and_facet %>%
    filter(grp == "Replicate") %>%
    group_by_at(vars(facet, "outcome")) %>%
    summarise(lwr = hdi(n)[1], med = median(n), upp = hdi(n)[2])

  ggplot(stats_by_cat_and_facet) +
    facet_wrap(as.formula(paste("~", facet))) +
    geom_bar(
      data = n_obs_by_cat_and_facet %>% filter(grp == "Actual"),
      aes(x = outcome, weight = n), color = "black", fill = "#7cb5ec"
    ) +
    geom_errorbar(aes(x = outcome, ymin = lwr, ymax = upp),
      width = 0, size = 1.5
    ) +
    geom_point(aes(x = outcome, y = med),
      pch = 21, fill = "white", size = 3, stroke = 1.5
    ) +
    ggthemes::theme_hc()
}
