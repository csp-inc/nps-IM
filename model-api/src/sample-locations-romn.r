source("model-api/src/utils-math-stats.R")
source("model-api/src/wrangle.r")


FIA_ANNULAR_RADIUS <- 36.6
DIST_TO_ABC_PLOTS <- 5
DIST_TO_D_PLOT <- 7.32
ANGLE <- 30
PI_METER_MARKS <- seq(1, 36, .5)
PCS_EPSG <- 26913 # NAD83 / UTM zone 13N


centroid_math <- function(component, distance) {
  if (component == 1 | component == "A") {
    c(0, distance)
  } else if (component == 2 | component == "B") {
    c(cos(rad(ANGLE)) * distance, -sin(rad(ANGLE)) * distance)
  } else if (component == 3 | component == "C") {
    c(-cos(rad(ANGLE)) * distance, -sin(rad(ANGLE)) * distance)
  } else if (component == "D") {
    c(cos(rad(ANGLE)) * distance, sin(rad(ANGLE)) * distance)
  }
}

get_sample_locations <- function(data, type, ...) {
  # Obtain sample locations.
  switch(type,
    plot = get_plot_locs(data, ...),
    transect = get_transect_locs(data, ...) # ,
    # point_intercept = get_point_intercept_locs(data, ...)  # TODO: deprecate?
  )
}

get_plot_locs <- function(data, ...) {

  # Pull plot and transect info for the record.
  transect_id_col <- unlist(...) %>% match_value("transect")
  transect_id <- data %>%
    pull(transect_id_col) %>%
    unique()
  plot_id_col <- unlist(...) %>% match_value("plot")
  plot_id <- data %>%
    pull(plot_id_col) %>%
    unique()

  # Calculate the shifts necessary to map individual samples.
  shift_vec_for_monument_center <-
    centroid_math(transect_id, FIA_ANNULAR_RADIUS)
  dist_to_plot <- ifelse(plot_id == "D", DIST_TO_D_PLOT, DIST_TO_ABC_PLOTS)
  shift_vec_for_plot <- centroid_math(plot_id, dist_to_plot)

  data %>%
    mutate(
      transect_terminus_x = site_x + shift_vec_for_monument_center[1],
      transect_terminus_y = site_y + shift_vec_for_monument_center[2],
      sample_x = ifelse(plot_id == "D", site_x, transect_terminus_x) +
        shift_vec_for_plot[1],
      sample_y = ifelse(plot_id == "D", site_y, transect_terminus_y) +
        shift_vec_for_plot[2]
    ) %>%
    select(-matches("terminus"))
}

get_point_intercept_locs <- function(data, ...) {
  transect_id_col <- unlist(...) %>% match_value("transect")

  rel_angles <- c(90, -30, -150)
  transect_pis <- lapply(seq_len(3), function(x) {
    m <- tan(rad(rel_angles[x]))
    dir <- ifelse(rel_angles[x] < -90, -1, +1)
    sample_x <- data$site_x[1] + (PI_METER_MARKS / sqrt(1 + m^2)) * dir
    sample_y <- data$site_y[1] + (PI_METER_MARKS * m / sqrt(1 + m^2)) * dir
    cbind(sample_x, sample_y, transect = x, meter_mark = PI_METER_MARKS)
  })
  do.call("rbind", transect_pis) %>%
    as_tibble(.name_repair = NULL) %>%
    mutate(site_id = data$site_id[1]) %>%
    rename_(.dots = setNames(list(rename_call("transect")), transect_id_col))
}

get_transect_locs <- function(data, ...) {
  get_point_intercept_locs(data, ...) %>%
    group_by(Transect, site_id) %>%
    summarise_at(c("sample_x", "sample_y"), mean) %>%
    ungroup()
}

get_most_recent_hits_sf <- function(data, grouping_vars, filter_string, coords) {
  df <- data %>%
    group_by_(.dots = grouping_vars) %>%
    filter(year == max(year)) %>%
    ungroup() %>%
    filter_(filter_string) %>%
    na.omit()
  st_as_sf(df, coords = coords, crs = PCS_EPSG)
}

make_plot_bboxes <- function(plot_data) {
  bb <- st_bbox(plot_data$geometry)
  bbox <- st_linestring(rbind(
    c(bb[1], bb[2]), c(bb[3], bb[2]),
    c(bb[3], bb[4]), c(bb[1], bb[4]), c(bb[1], bb[2])
  ))
  bbox <- st_polygonize(bbox)
  plot_data %>%
    as_tibble(.name_repair = NULL) %>%
    mutate(geometry = st_as_text(bbox))
}
