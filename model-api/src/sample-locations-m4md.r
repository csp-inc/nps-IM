get_sample_locations <- function(data, type, ...) {
  data %>%
    rename(sample_x = x_coord, sample_y = y_coord)
}
