library(tidyverse)
library(lubridate)
library(hrbrthemes)

# Helpers.
degrees_to_radians <- function(degrees) {
  degrees * pi / 180
}
to_slope_percent <- function(degrees) {
  tan(degrees_to_radians(degrees)) * 100
} # to_slope_percent(45)
get_effects <- function(data) {
  data %>%
    mutate(
      # sample_year_squared = sample_year^2,
      date = date_of_first_sample + years(sample_year),
      date_aligned = min(date_of_first_sample) + years(sample_year),
      treat_by_time = ifelse(plot_type == "B", 1 * sample_year, 0),
      # treat_by_time_squared = treat_by_time^2,
      z0 = ifelse(plot_type == "B", 1, 0),
      z1 = ifelse(plot_type == "B" & sample_year > 0, 1, 0),
      z2 = ifelse(plot_type == "B" & sample_year == 2, 1, 0),
      # z3 = ifelse(plot_type == "B" & sample_year == 3, 1, 0),
      z3 = ifelse(plot_type == "B" & sample_year > 1 & sample_year <= 3,
        sample_year - 1, 0
      ),
      roi = "TSRA"
    )
}

files <- list.files("assets/uplands-data/CASP", pattern = "*0to25.csv$", full.names = TRUE)

plot_info <- read_csv("assets/uplands-data/CASP/metadata/plot_info.csv") %>%
  mutate(
    staggering_raw = date_of_first_sample - min(date_of_first_sample),
    staggering = as.numeric(staggering_raw)
  ) %>%
  distinct(plot_name, plot_type, park_name, dominant_spp, staggering)

park_lookup <- plot_info %>% distinct(plot_name, park_name)

# The effect of aspect should involve an interaction with slope:
# https://www.fs.fed.us/rm/pubs_journals/1976/rmrs_1976_stage_a001.pdf
plot_covariates <- read_csv("assets/uplands-data/CASP/metadata/plot_biophysical_data.csv") %>%
  select(plot_name, plot_type, elevation, northness,
    slope_angle = slope,
    matches("tpi")
  ) %>%
  mutate(
    slope_percent = to_slope_percent(slope_angle),
    s = slope_percent / 100
  ) %>%
  # s_x_northness = s * northness) %>%
  select(-slope_angle, -slope_percent, -matches("tpi")) %>% # was -northness,
  left_join(park_lookup, by = "plot_name") %>%
  left_join(plot_info) %>%
  mutate(dominant_genus = substr(dominant_spp, 1, 2))
write_csv(plot_covariates, "assets/uplands-data/CASP/modified/plot-level-covariates.csv")

s_deg <- seq(0, 50)
s_perc <- to_slope_percent(s_deg)
azimuth <- seq(0, 360)
# plot(azimuth, cos(degrees_to_radians(azimuth))); abline(v = seq(0, 360, 90), col = 'red')

d_terrain <- tibble(azimuth, northness = cos(degrees_to_radians(azimuth))) %>%
  crossing(s_perc = seq(0, 500, 25)) %>%
  arrange(s_perc) %>%
  mutate(
    s_perc_x_northness = s_perc * northness,
    terrain_effect_naive = 0.08 * s_perc_x_northness,
    terrain_effect_stage = 0.08 * s_perc_x_northness + (0.13) * s_perc
  ) %>%
  pivot_longer(matches("terrain_effect"), names_to = "approach", values_to = "terrain_effect")
ggplot(d_terrain) +
  facet_wrap(~approach) +
  geom_line(aes(
    x = azimuth, y = terrain_effect, group = s_perc,
    color = factor(s_perc)
  ))
# plot(s_perc~s_deg)

plot_loc_info <- read_csv("assets/uplands-data/CASP/metadata/plot_info.csv") %>%
  distinct(plot_name, park_name, easting, northing)
write_csv(plot_loc_info, "assets/uplands-data/CASP/modified/plot-locs.csv")

sampling_lookup <- read_csv("assets/uplands-data/CASP/metadata/plot_info.csv") %>%
  distinct(plot_name, park_name, date_of_first_sample)

library(tidyverse)

d_z <- crossing(
  b0 = 10, b1 = -0.5, t = seq(0, 5), p = c("Burn", "Control"),
  b_z0 = c(0, -1), b_z1 = c(0, -1.5), b_z2 = c(0, 1.5),
  b_z3 = c(0, 1.5 / 2),
  b_tbt = c(0, -0.5)
) %>%
  mutate(
    z0 = as.integer(p == "Burn"), # control is reference level
    z1 = ifelse(p == "Burn" & t > 0, 1, 0),
    z2 = ifelse(p == "Burn" & t == 2, 1, 0),
    z3 = ifelse(p == "Burn" & t > 1 & t <= 3, t - 1, 0),
    treat_by_time = z0 * t, # was: ifelse(p == "Burn", 1 * t, 0),
    lp = b0 + b1 * t + b_z0 * z0 + b_z1 * z1 + b_z2 * z2 + b_z3 * z3 +
      b_tbt * treat_by_time
  )
d_z_scenario <- d_z %>%
  filter(b_z0 == -1, b_z2 == 0, b_tbt == -0.5)
# filter(b_z0 == -1, b_z3 == 0, b_tbt == -0.5)
ggplot(d_z_scenario) + # filter(b_tbt == 0)
  facet_grid(b_z1 ~ b_z3, labeller = label_both) +
  # facet_grid(b_z1~b_z2, labeller = label_both) +
  geom_line(aes(x = t, y = lp, color = p), size = 0.8) +
  labs(x = "Year", y = "Mean") +
  scale_color_brewer("Plot type", type = "qual", palette = "Set1") +
  theme_ipsum_rc(
    grid = "Y", base_size = 16, axis_title_size = 18, strip_text_size = 18
  )

d_z_scenario <- d_z %>%
  # filter(b_z0 %in% c(0, -1), b_z1 == 0, b_z2 == 0, b_z3 == 0, b_tbt == 0) # z0-concept.jpg
  # filter(b_z0 %in% c(-1), b_z2 == 0, b_z3 == 0, b_tbt == 0) # z1-concept.jpg
  # filter(b_z0 %in% c(-1), b_z1 == -1.5, b_z3 == 0, b_tbt == 0) # z2-concept.jpg
  # filter(b_z0 %in% c(-1), b_z1 == -1.5, b_z2 == 0, b_tbt == 0) # z3-concept.jpg
  filter(b_z0 %in% c(-1), b_z1 == 0, b_z2 == 0, b_z3 == 0) %>%
  rename(b_z0Xrel_year = b_tbt) # z0-rel_year.jpg
p_z_scenario <- ggplot(d_z_scenario) + # filter(b_tbt == 0)
  # facet_grid(.~b_z0, labeller = label_both) +
  facet_grid(b_z0 ~ b_z0Xrel_year,
    labeller = function(labs) {
      label_both(labs, multi_line = TRUE)
    }
  ) +
  # geom_vline(xintercept = 3, linetype = 'dashed', alpha = 0.5, size = 0.5) +
  # geom_hline(yintercept = seq(5, 10, 1), linetype = 'dashed', color = 'orange', alpha = 0.5, size = 0.5) +
  geom_line(aes(x = t, y = lp, color = p), size = 1, alpha = 0.75) +
  labs(x = "Relative year", y = "Linear predictor") +
  scale_color_brewer("", type = "qual", palette = "Set1") +
  theme_ipsum_rc(
    grid = "Y", base_size = 16, axis_title_size = 18, strip_text_size = 18
  ) +
  theme(
    legend.position = "none", # c(0.875, 0.95), #'bottom', #legend.justification="left",
    legend.margin = margin(0, 0, 0, 0),
    # legend.box.margin=margin(-10,-10,-10,-10),
    plot.margin = margin(0, 0, 0, 0),
    # axis.text.x = element_blank(), axis.text.y = element_blank(),
    axis.title.x = element_blank(), axis.title.y = element_blank(),
    panel.spacing = unit(0.5, "lines")
  ) +
  scale_y_continuous(
    breaks = scales::pretty_breaks(n = 3),
    limits = c(min(d_z_scenario$lp) * 0.8, max(d_z_scenario$lp) * 1.2)
  )
save_figure(
  plot = p_z_scenario, path = "tmp",
  filename = "z0Xrel_year-concept.jpg", p_width = 2, p_height = 1.5
)

ok <- lapply(files, function(x) {
  out <- read_csv(x) %>%
    filter(
      # FABCO1D10-15 was set up in 1992 as a burn plot but never burned.
      # FABCO1D10-45 was within a prescribed burn perimeter but the severity
      # reading is 5, which means the burn didn’t actually go through the plot.
      !(plot_name == "FABCO1D10-45" & sample_year == 1)
    ) %>%
    left_join(sampling_lookup) %>%
    drop_na(date_of_first_sample) %>%
    mutate_at(vars(matches("_per_ha")), function(x) x * 0.1) %>%
    rename_at(vars(matches("_per_ha")), function(x) sub("_per_ha", "_per_plot", x)) %>%
    get_effects() %>%
    mutate(time_x_treat = interaction(sample_year, plot_type, sep = "/"))
  if (grepl("forest_structure", x)) {
    sq_m_per_overstory_plot <- 20 * 50
    sq_m_per_acre <- 4046.85642
    out <- out %>%
      mutate(n_acres_per_plot = sq_m_per_overstory_plot / sq_m_per_acre)
  }
  write_csv(
    out %>% select(-treat_by_time, -z1, -z2),
    sprintf("assets/uplands-data/CASP/modified/%s", basename(x))
  )
  out
})
ok[[1]]
# ok[[1]] %>% distinct(sample_year, year_fct)
# levels(ok[[1]]$year_fct)

d_community <- ok[[4]] %>%
  group_by(plot_name, plot_type) %>%
  filter(sample_year == min(sample_year)) %>%
  summarise( # ratio_white_fir = white_fir_trees_per_plot / pine_trees_per_plot,
    prop_white_fir = white_fir_trees_per_plot /
      (pine_trees_per_plot + white_fir_trees_per_plot), .groups = "drop"
  )

d_burn_severity_raw <- read_csv("assets/uplands-data/CASP/metadata/burn_severity_(combined_data_sources).csv")
d_burn_severity <- d_burn_severity_raw %>%
  group_by(plot_name) %>%
  mutate(combined_severity = mean(c(stanton_severity, mandeno_severity), na.rm = TRUE)) %>%
  ungroup()
# ok[[1]] %>% group_by(plot_name, plot_type) %>% tally()
plot_covariates_full <- ok[[1]] %>%
  distinct(plot_name, plot_type, sample_year, roi, date_of_first_sample) %>%
  complete(
    sample_year = seq(min(sample_year), max(sample_year)),
    nesting(plot_name, plot_type, roi, date_of_first_sample)
  ) %>%
  get_effects() %>%
  left_join(plot_covariates) %>%
  left_join(d_community) %>%
  left_join(d_burn_severity %>% select(plot_name, combined_severity)) %>%
  mutate(
    combined_severity = ifelse(plot_type == "C", 5, combined_severity),
    combined_severity = ifelse(is.na(combined_severity), 5, combined_severity),
    # Flip the severity scale so it is more intuitive (i.e., so that larger
    # numbers correspond to higher burn severities).
    burn_severity = 5 - combined_severity
  ) %>%
  select(-combined_severity)

# d_tmp = plot_covariates_full %>% left_join(d_burn_severity_raw) %>%
#   select(plot_type, plot_name, matches("severity")) %>%
#   distinct() %>% filter(plot_type == 'B')
# write_csv(d_tmp, 'sandbox/joined-burn-severity.csv')
# d_tmp %>% filter(is.na(stanton_severity) & is.na(mandeno_severity))

# filter(plot_name %in% c('FABCO1D10-09', 'FABCO1D10-10')) %>% arrange(plot_name)
plot_covariates_full <- plot_covariates_full %>%
  mutate(
    event = year(date_aligned) - min(year(date_aligned)),
    event_cat = sprintf("yr_%s", event),
    cal_year_is_gte_2010 = as.integer(year(date) >= 2010)
  )
write_csv(plot_covariates_full, "assets/uplands-data/CASP/modified/plot-level-covariates-complete.csv")
# plot_covariates_full %>% group_by(plot_name) %>% tally() %>% as.data.frame()

# TODO: Luke start here!

p_hist <- ggplot(plot_covariates_full) +
  geom_histogram(aes(x = s * northness, y = ..density..), # burn_severity, prop_white_fir, s
    color = "black", fill = "orange", bins = 10
  ) +
  theme_ipsum_rc(
    grid = "Y", base_size = 16, axis_title_size = 18, strip_text_size = 18
  ) +
  labs(x = "North facing-ness", y = "Density")
save_figure(
  plot = p_hist, path = "tmp",
  filename = "sXnorthness-hist.jpg", p_width = 4, p_height = 1.5
)

ok_interim1 <- ok[[1]] %>%
  group_by(plot_type) %>%
  mutate(
    plot_idx = as.integer(as.factor(plot_name)),
    rel_year = year(date_aligned) - min(year(date_aligned))
  )

p_rel_year <- ggplot(ok_interim1) +
  facet_wrap(~plot_type) +
  geom_point(aes(x = date, y = plot_idx, fill = rel_year), # burn_severity, prop_white_fir, s
    color = "black", pch = 21, size = 2
  ) +
  theme_ipsum_rc(
    grid = "Y", base_size = 16, axis_title_size = 18, strip_text_size = 18
  ) +
  labs(x = "Relative year", y = "Plot index") +
  scale_fill_viridis_c("Relative year", direction = -1, option = "inferno") +
  theme(legend.position = "bottom", legend.key.width = unit(1, "cm"), plot.margin = margin(0, 0, 0, 0))
p_rel_year
save_figure(
  plot = p_rel_year, path = "tmp",
  filename = "rel_year-sched.jpg", p_width = 2, p_height = 1.5
)

ok_interim2 <- ok[[1]] %>%
  left_join(
    plot_covariates_full %>% select(plot_name, date_of_first_sample, date, date_aligned, staggering)
  ) %>%
  group_by(plot_type) %>%
  mutate(plot_idx = as.integer(as.factor(plot_name)))

p_pt <- ggplot(ok_interim2) +
  facet_wrap(~plot_type) +
  geom_segment(aes(
    x = date_of_first_sample, xend = date_of_first_sample - days(staggering),
    y = plot_idx, yend = plot_idx, color = staggering
  )) +
  geom_point(aes(x = date_of_first_sample, y = plot_idx), # burn_severity, prop_white_fir, s
    color = "black", fill = "white", pch = 21, size = 2
  ) +
  theme_ipsum_rc(
    grid = "Y", base_size = 16, axis_title_size = 18, strip_text_size = 18
  ) +
  labs(x = "Date", y = "Plot index") +
  scale_color_viridis_c("Staggering", direction = -1, option = "inferno") +
  theme(legend.position = "bottom", legend.key.width = unit(1, "cm"), plot.margin = margin(0, 0, 0, 0))
p_pt
save_figure(
  plot = p_pt, path = "tmp",
  filename = "staggering-sched.jpg", p_width = 2, p_height = 1.5
)
# save_figure(plot = p_hist, path = 'tmp',
#             filename = 'sXnorthness-hist.jpg', p_width = 4, p_height = 1.5)
