get_sample_locations <- function(data, type, ...) {
  data %>%
    rename(sample_x = easting, sample_y = northing)
}
