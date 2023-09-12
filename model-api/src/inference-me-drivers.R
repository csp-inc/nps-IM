get_me_inference <- function(mcarray_list, folder, resp_var_desc,
                             interval = 0.95, method = c("hdi", "eti")) {
  if (is.null(mcarray_list$mu.driver.strat)) {
    return(NULL)
  }
  message("Inference on the marginal effects of drivers...")

  tag_as_scaled <- function(x) {
    paste(x, "scaled", sep = "_")
  }

  X_driver_lookup <- read_csv(file.path(folder, "00-input/X-driver-lookup.csv"))
  X_driver_lookup <- X_driver_lookup %>%
    rename_with(
      tag_as_scaled,
      any_of(sub("_raw", "", grep("fe_.*_raw$", names(X_driver_lookup),
        value = TRUE
      )))
    )

  mu_driver_strat <- concat_chains(mcarray_list$mu.driver.strat[, 1, , ], axis = 3)
  dimnames(mu_driver_strat) <-
    list(NULL, sprintf("iter_%s", 1:ncol(mu_driver_strat)))

  # draws <- sample(ncol(mu_driver_strat), min(250, ncol(mu_driver_strat)))
  tmp <- X_driver_lookup %>%
    bind_cols(as_tibble(mu_driver_strat)) %>% # [, draws]
    pivot_longer(matches("iter_"),
      names_to = "iter_id", values_to = "iter_val"
    ) %>%
    select(matches("_scaled|_raw|incr_|iter_|me_"))



  out_group <- function(x, y) {
    z <- x[!x %in% y]
    if (identical(z, character(0))) {
      NA
    } else {
      z
    }
  }

  all_along <- unique(tmp$me_along)
  me_grps <- X_driver_lookup %>% distinct(me_idx, me_along)
  tmp2_precursor <- lapply(1:max(X_driver_lookup$me_idx), function(i) {
    along <- me_grps %>%
      filter(me_idx == i) %>%
      pull(me_along)
    cond_on <- out_group(all_along, along)

    tmp %>%
      filter(me_idx == i) %>%
      mutate(
        incr_time = "zero", isNA = is.na(get(sprintf("%s_raw", along))),
        x_val = ifelse(isNA, get(sprintf("%s_scaled", along)),
          get(sprintf("%s_raw", along))
        )
      ) %>%
      unite("me_at", any_of(sprintf("incr_%s", c(cond_on, "time"))),
        sep = ":", remove = FALSE
      ) %>%
      mutate(
        me_at_key = paste(c(cond_on, "time"), collapse = ":"),
        intx = interaction(me_at, me_at_key, iter_id)
      ) %>%
      ungroup() # %>% distinct(me_at, me_at_key, intx)
  })
  tmp2 <- bind_rows(tmp2_precursor) %>%
    group_by(me_along) %>%
    mutate(is01 = all(x_val %in% c(0, 1))) %>%
    ungroup() %>%
    mutate(
      col = "black",
      col = ifelse(is01, "orange", col), col = ifelse(isNA, "steelblue", col)
    ) # interaction(is01, isNA)

  tmp3 <- tmp2 %>%
    group_by(me_along, x_val, me_at, me_at_key) %>%
    summarise(
      lwr = get_interval(iter_val, interval, method[1])[1],
      mid = median(iter_val),
      upr = get_interval(iter_val, interval, method[1])[2],
      .groups = "drop"
    )
  # tmp3 %>% distinct(me_along, me_at, me_at_key) %>% View()

  save_object(tmp3, file.path(folder, "03-inference/me"), "me-X.rds")
  tmp4 <- tmp3 %>%
    pivot_longer(all_of(c("lwr", "mid", "upr")),
      names_to = "stat", values_to = "val"
    ) %>%
    filter(!grepl("seen", me_at)) %>%
    mutate(lty = ifelse(stat == "mid", "solid", "dashed"))

  p <- ggplot(
    tmp2 %>%
      filter(
        !grepl("seen", me_at),
        iter_id %in% sample(
          unique(iter_id),
          min(250, ncol(mu_driver_strat))
        )
      ),
    aes(x = x_val, y = iter_val, group = intx, color = col)
  ) +
    facet_wrap(~me_along, scales = "free") +
    geom_line(alpha = 0.05) +
    ggthemes::theme_hc(base_size = 16) +
    labs(
      x = "Covariate value", y = resp_var_desc,
      title = "Marginal effects of covariates",
      subtitle = "Continuous (black) and indicator (orange) variables (raw values), and interactions (blue, scaled only)"
    ) +
    geom_line(
      data = tmp4,
      aes(x = x_val, y = val, group = stat, linetype = lty),
      color = "black"
    ) +
    scale_linetype_identity() +
    theme(legend.position = "none") +
    scale_color_identity()
  # scale_color_brewer('', type = 'qual', palette = 'Set1')
  save_figure(
    plot = p,
    path = file.path(folder, "03-inference/me"),
    filename = "me-X.jpg",
    p_width = 3 * 1.5, p_height = 2 * 1.5,
    device = NULL
  )



  # tmp3
  # lapply(all_along, function(i) {
  #   ggplot(tmp3 %>% filter(me_along == i)) +
  #     facet_grid(~me_along, scales = 'free') +
  #     geom_line(aes(x = x_val, y = mid,
  #                   color = interaction(me_at_key, me_at, sep = '\n'))) +
  #     ggthemes::theme_hc(base_size = 16) +
  #     scale_color_brewer('', type = 'qual', palette = 'Set1')
  # })
  # tmp3 %>% filter(me_along == i) %>% separate(me_at, into = NA, sep = ':')

  # scale_color_brewer('', guide = F)
}
