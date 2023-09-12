library(glue)

get_funnel_plot <- function(mcarray_list, jags_data, output_dir) {
  if ("sigma.B0" %in% names(mcarray_list)) {
    lapply(seq_len(jags_data$y.n.strata), function(this_stratum) {
      this_chain <- 1
      n_groups_to_visualize <- 8

      n_groups <- jags_data$y.n.site[this_stratum]

      these_groups <- sample(seq_len(n_groups), min(n_groups, n_groups_to_visualize))


      sigma_b0 <- tibble(value = mcarray_list$sigma.B0[this_stratum, , this_chain]) %>%
        mutate(iteration = 1:n(), name = "sigma[beta[0]]", group = "NA")
      beta_js <- mcarray_list$B[these_groups, 1, this_stratum, , this_chain]
      dimnames(beta_js) <- list(these_groups, NULL)

      combined_samples <- as_tibble(t(beta_js)) %>%
        mutate(iteration = 1:n(), name = "beta[j]") %>%
        gather(group, value, -iteration, -name) %>%
        bind_rows(sigma_b0)

      p_funnel_trace <- ggplot(combined_samples) +
        facet_wrap(~name, nrow = 2, scales = "free", labeller = label_parsed) +
        geom_line(aes(x = iteration, y = value, color = group), alpha = 0.5) +
        scale_color_brewer(name = "Site", type = "qual", palette = "Set1") +
        labs(x = "Iteration", y = "Value") +
        ggthemes::theme_hc()

      save_figure(
        plot = p_funnel_trace,
        path = file.path(output_dir, "01-diagnostics/funnel"),
        filename = glue("funnel-traces-stratum-{this_stratum}.jpg"),
        p_width = 7, p_height = 3.5
      )

      selected_group <- sample(these_groups, 3)
      funnel_plot_data <- lapply(selected_group, function(this_group) {
        combined_samples %>%
          filter(group == this_group | group == "NA") %>%
          select(-group) %>%
          spread(name, value) %>%
          mutate(group = this_group)
      }) %>% bind_rows()

      p_funnel <- ggplot(funnel_plot_data) +
        facet_wrap(~group, ncol = 3) +
        geom_point(aes(x = `beta[j]`, y = `sigma[beta[0]]`),
          alpha = min(1000 / n_iter, 0.5)
        ) +
        labs(
          x = substitute(beta[j == jj], list(jj = selected_group)),
          y = expression(sigma[beta[0]])
        ) +
        ggthemes::theme_hc()

      save_figure(
        plot = p_funnel,
        path = file.path(output_dir, "01-diagnostics/funnel"),
        filename = glue("funnel-plot-stratum-{this_stratum}.jpg"),
        p_width = 2.5, p_height = 2.5
      )
    })
  } else {
    NULL
  }
}
