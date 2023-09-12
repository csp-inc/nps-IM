get_zone_inference <- function(mcarray_list, data, folder, resp_var_desc,
                               timesteps, likelihood,
                               wrap_char_n = 15, x_hat_raw = timesteps,
                               ...) {
  jags_obj_prefix <- str_extract(
    grep("zone.mean$", names(mcarray_list), value = TRUE),
    ".*(?=\\.zone\\.mean)"
  )
  lapply(jags_obj_prefix, function(x) {
    jags_stuff <- concat_chains(mcarray_list[[sprintf("%s.zone.mean", x)]], axis = 4)
    zone_summary <- lapply(1:ncol(jags_stuff), function(l) {
      means <- matrix(apply(jags_stuff[, l, ], 1, mean), nrow = 1)
      dimnames(means) <- list("mean", NULL)
      l_summary <-
        rbind(means, apply(jags_stuff[, l, ], 1, HDInterval::hdi, credMass = 0.95))

      # l_summary <- apply(jags_stuff[, l, ], 1, quantile, probs = c(0.025, 0.5, 0.975))
      as_tibble(l_summary) %>%
        mutate(stat = dimnames(l_summary)[[1]]) %>%
        gather(year_col, val, -stat) %>%
        spread(stat, val) %>%
        left_join(get_year_lookup(data, timesteps)) %>%
        mutate(zone = filter(get_zone_lookup(folder), zone_idx == l)$zone)
    }) %>% bind_rows()
    p_zone_means <- ggplot(zone_summary) +
      facet_wrap(~zone) +
      geom_line(aes(x = cal_year, y = mean)) +
      geom_ribbon(aes(x = cal_year, ymin = lower, ymax = upper),
        fill = "steelblue", alpha = 0.2
      ) +
      labs(x = "Year", y = "Zone mean") +
      theme_ima("default")
    save_figure(
      plot = p_zone_means,
      path = file.path(folder, "03-inference/zone"),
      filename = sprintf("%s-zone-means.png", x),
      p_width = 3, p_height = 5, device = NULL
    )
    save_table(zone_summary, file.path(folder, "03-inference/zone"), sprintf("%s-zone-means.csv", x))
    save_object(
      jags_stuff,
      file.path(folder, "03-inference/zone"),
      sprintf("%s-jags-zone-means.rds", x)
    )
  })
}
