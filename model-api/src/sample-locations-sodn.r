# TODO: specify libraries/dependencies


POINT_M <- seq(.5, 20, .5)
PCS_EPSG <- 26912 # NAD83 / UTM zone 12N

transect_azimuth_and_centroid <- function(x, y) {
  tibble(
    "transect_azimuth" = atan2(y[2] - y[1], x[2] - x[1]) * 180 / pi,
    "transect_centroid_x" = mean(x),
    "transect_centroid_y" = mean(y),
    "pi_order" = ifelse(y[2] > y[1], 1, -1)
  )
}

get_transect_params <- function(data) {
  data %>%
    bind_cols(data %>% st_coordinates(geometry) %>% as_tibble(.name_repair = NULL)) %>%
    mutate(transect_position = factor(transect_position, levels = c("origin", "end"))) %>%
    arrange(transect_position) %>% # ordered from origins to ends
    group_by(Transect) %>%
    do(transect_azimuth_and_centroid(.$X, .$Y)) %>%
    ungroup()
}

point_along_transect <- function(data, d = seq(.25, 10, .5)) {
  browser()
  m <- tan(rad(data[["transect_azimuth"]]))
  origin <- data %>%
    select(transect_centroid_x, transect_centroid_y) %>%
    as.matrix()
  origin_mat <- matrix(rep(origin, each = length(d)),
    nrow = length(d),
    dimnames = list(NULL, c("X", "Y"))
  )
  d_mat <- matrix(d, ncol = 1)
  pis_mat <- cbind(d_mat * (1 / sqrt(1 + m^2)), d_mat * (m / sqrt(1 + m^2)))
  if (data[["pi_order"]] == 1) {
    pi_m <- POINT_M
  } else {
    pi_m <- rev(POINT_M)
  }
  rbind(origin_mat + pis_mat, origin_mat - pis_mat) %>%
    as_tibble(.name_repair = NULL) %>%
    rename(pi_x = X, pi_y = Y) %>%
    arrange(pi_y) %>%
    mutate(point = pi_m)
}

transect_point_intercepts <- function(data) {
  transect_params <- data %>% get_transect_params()
  transect_params %>%
    group_by(Transect) %>%
    do(point_along_transect(.)) %>%
    ungroup() %>%
    mutate(pi_geometry = paste0("POINT(", pi_x, " ", pi_y, ")"))
}

get_sample_locations <- function(data, type, ...) {
  # The 'precursors' already finds the center of each transect, as such this is
  # really just a pass-through.
  data %>%
    rename(sample_x = UTM_x, sample_y = UTM_y)
  # rename(sample_x=transect_centroid_x, sample_y=transect_centroid_y)
}

macroplot_boundary <- function(data) {
  macroplot_geometry <- data %>%
    st_union() %>%
    st_convex_hull() %>%
    st_buffer(dist = 5.25, nQuadSegs = -2)
  tibble(
    plot_id = data[["plot_id"]][1],
    macroplot_bounds = st_as_text(macroplot_geometry)
  )
}

get_most_recent_hits_sf <- function(data, filter_string) {
  df <- data %>%
    group_by(plot_id) %>%
    filter(year == max(year)) %>%
    ungroup() %>%
    filter_(filter_string) %>%
    na.omit()
  st_as_sf(df, coords = c("pi_lon", "pi_lat"), crs = GCS_EPSG)
}

get_most_recent_plot_cover_sf <- function(data, macroplots_sf) {
  data %>%
    group_by(plot_id) %>%
    filter(year == max(year)) %>%
    ungroup() %>%
    group_by(year, plot_id) %>%
    summarise(
      c3_cover = sum(any_is_live_c3_npg) / 300 * 100,
      c4_cover = sum(any_is_live_c4_npg) / 300 * 100
    ) %>%
    right_join(macroplots_sf) %>%
    na.omit() %>%
    mutate(geometry = st_as_text(macroplot_bounds)) %>%
    st_as_sf(wkt = "geometry", crs = GCS_EPSG)
}
