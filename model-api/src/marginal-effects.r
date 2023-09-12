library(tidyverse)
library(cowplot)


# ---- upper panel -----

vis_vars <- c("AMJP", "JASP") # L_AMJP, AMJP, NDJFMP

weather_wide <- read_csv("data/SCPN/CHCU/AllClimCov.csv") %>%
  select(keys, vis_vars)
weather_long <- weather_wide %>%
  gather(var, val, -keys)
weather_long_labs <- weather_long %>%
  filter(var == "AMJP") %>%
  summarise(
    min = min(val), max = max(val),
    x_min = Year[which.min(val)],
    x_max = Year[which.max(val)]
  ) %>%
  ungroup() %>%
  gather(empirical_stat, stat, -matches("x_")) %>%
  mutate(
    x_start = ifelse(empirical_stat == "min", x_min, x_max),
    x_end = Inf, lab = sprintf("%s mm", format(stat, nsmall = 1))
  )

p <- ggplot(weather_long) +
  # facet_wrap(~var) +
  geom_line(aes(x = Year, y = val, group = interaction(Plot, var), color = var),
    alpha = 0.8
  ) +
  ggthemes::theme_hc(base_size = 16) +
  labs(x = "Year", y = "Accumulated rainfall (mm)") +
  scale_color_grey("Period",
    labels = c("Spring (AMJP)", "Monsoon (JASP)"),
    start = 0.8, end = 0.2
  ) +
  theme(legend.position = c(0.05, 0.9))
p1 <- p +
  geom_segment(
    data = weather_long_labs,
    aes(x = x_start, xend = 2019.25, y = stat, yend = stat),
    linetype = "dotted"
  ) +
  geom_point(
    data = weather_long_labs,
    aes(x = x_start, y = stat), size = 2
  ) +
  geom_label(
    data = weather_long_labs, aes(x = 2019.25, y = stat, label = lab),
    hjust = "inward", vjust = "inward", size = 5
  ) +
  geom_text(aes(x = Inf, y = Inf, label = sprintf("(%s)", "a")),
    size = 5,
    vjust = "inward", hjust = "inward"
  ) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5))
p1

# ---- lower panel -----

X_driver_lookup <- read_csv(file.path(output_dir, "00-input/X-driver-lookup.csv"))
# mutate(incr_fe_AMJP = ifelse(is.na(incr_fe_AMJP), fe_AMJP, incr_fe_AMJP))

mu_driver_strat <- concat_chains(z_jags$mu.driver.strat[, 1, , ], axis = 3)
dimnames(mu_driver_strat) <-
  list(NULL, sprintf("iter_%s", 1:ncol(mu_driver_strat)))

me_of <- sprintf("fe_%s_raw", "JASP") # marginal effects of which variable?
conditional_on <- sprintf("incr_fe_%s", "AMJP") # conditional on another?
bookends <- c("min_seen", "max_seen") # c('-1.96', '1.96')

me_test <- X_driver_lookup %>%
  bind_cols(as_tibble(mu_driver_strat)) %>%
  gather(iter, val, -one_of(names(X_driver_lookup))) %>% # overall_mean
  filter(get(conditional_on) %in% bookends) %>%
  mutate(line_grp = interaction(get(conditional_on), iter))

cred_mass <- 0.5
me_test_summary <- me_test %>%
  group_by(!!me_of := get(me_of), !!conditional_on := get(conditional_on)) %>%
  summarise( # lwr = HDInterval::hdi(val, credMass = cred_mass)[1],  #
    lwr = quantile(val, probs = c(0.25)),
    mid = median(val),
    # upr = HDInterval::hdi(val, credMass = cred_mass)[2],
    upr = quantile(val, probs = c(0.75))
  ) %>% #
  ungroup() %>%
  gather(this_stat, line_val, -me_of, -!!conditional_on) %>%
  mutate(
    line_grp = interaction(get(conditional_on), this_stat),
    lty = ifelse(this_stat == "mid", "solid", "dashed")
  )

p2 <- ggplot(me_test %>% filter(iter %in% sample(unique(iter), 200))) +
  # geom_line(aes_string(x = me_of, y = 'val',
  #               group = 'line_grp', color = conditional_on),
  #           alpha = 0.1) +
  geom_line(
    data = me_test_summary,
    aes_string(
      x = me_of, y = "line_val", group = "line_grp",
      color = conditional_on, linetype = "lty"
    ),
    size = 0.8
  ) +
  ggthemes::theme_hc(base_size = 16) +
  scale_linetype_identity() +
  # labs(x = me_of, y = 'Mean cover') +
  scale_color_grey() +
  labs(x = "Monsoonal rainfall (JASP, mm)", y = "Mean cover") +
  guides(color = FALSE) +
  geom_text(aes(x = Inf, y = Inf, label = sprintf("(%s)", "b")),
    size = 5,
    vjust = "inward", hjust = "inward"
  )

plot_grid(p1, p2, ncol = 1)
ggsave("sandbox/example-2-drivers.png", height = 3 * 3, width = 2 * 3)
