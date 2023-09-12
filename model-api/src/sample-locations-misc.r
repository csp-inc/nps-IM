
get_sample_locations <- function(data, type, ...) {
  # This is creally just a pass-through.
  data %>%
    rename(sample_x = site_x, sample_y = site_y)
}
