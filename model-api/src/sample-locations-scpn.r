source("model-api/src/utils-math-stats.R")


QUADRAT_CENTERS_AC <- c(2.5, 12.5, 22.5, 32.5, 42.5)
QUADRAT_CENTERS_B <- c(7.5, 17.5, 27.5, 37.5, 47.5)
PCS_EPSG <- 26912 # NAD83 / UTM zone 12N


line_angle_degrees <- function(x, y) {
  atan2(y[2] - y[1], x[2] - x[1]) * 180 / pi
}

offset_coords <- function(from_coords, angle, d, direction) {
  xy_offset <- c(cos(rad(angle)) * d, sin(rad(angle)) * d) * direction
  from_coords + xy_offset
}

point_along_transect <- function(origin, angle, d, direction, o) {
  if (o %in% c("NW", "SE")) direction <- direction * -1
  m <- tan(rad(angle))
  origin_mat <- matrix(rep(origin, each = length(d)),
    nrow = length(d),
    dimnames = list(NULL, c("X", "Y"))
  )
  d_mat <- matrix(d, ncol = 1)
  origin_mat + cbind(
    d_mat * (1 / sqrt(1 + m^2)), d_mat * (m / sqrt(1 + m^2))
  ) * direction
}

quadrant <- function(x) {
  if (x > 0 & x < 90) {
    "SE"
  } else if (x > 90 & x < 180) {
    "NE"
  } else if (x < 0 & x > -90) {
    "SW"
  } else {
    "NW"
  }
}

direction <- function(o) {
  if (o == "NE") {
    c(-1, -1)
  } else if (o == "SE") {
    c(-1, -1) # c(-1, +1)
  } else if (o == "SW") {
    c(+1, +1)
  } else {
    c(+1, +1) # c(+1, -1)
  }
}

macroplot_boundary <- function(data) {
  macroplot_geometry <- st_sf(data) %>%
    filter(!stringr::str_detect(line_field, "B$")) %>%
    st_union() %>%
    st_convex_hull() %>%
    st_buffer(dist = 10.5, nQuadSegs = -2) # not sure why -2 works, but it does!
  data %>%
    select(site_id) %>%
    distinct() %>%
    mutate(macroplot_bounds = st_as_text(macroplot_geometry))
}

quadrat_centroids <- function(data) {
  data %<>%
    mutate(origin = grepl("[ABC]S", Location)) %>%
    arrange(desc(origin))

  coords <- data %>%
    select(Easting, Northing) %>%
    as.matrix()
  angle <- line_angle_degrees(coords[, 1], coords[, 2])
  orientation <- quadrant(angle)

  shift_dir <- direction(orientation)
  layout_dir <- diff(coords[, 1]) / abs(diff(coords[, 1]))
  layout_dir <- ifelse(orientation %in% c("SE", "NW"), layout_dir * -1, layout_dir)
  shift_angle <- angle + 90 * layout_dir

  offset_origin <- offset_coords(coords[1, , drop = FALSE], shift_angle, 1, shift_dir)

  transect_id <- stringr::str_sub(data$Location[1], 1, 1)
  quadrat_centroids <-
    if (transect_id == "B") {
      point_along_transect(offset_origin, angle, QUADRAT_CENTERS_B, layout_dir, orientation)
    } else {
      point_along_transect(offset_origin, angle, QUADRAT_CENTERS_AC, layout_dir, orientation)
    }


  quadrat_centroids %>%
    as_tibble(.name_repair = NULL) %>%
    rename(Quadrat_X = X, Quadrat_Y = Y) %>%
    mutate(Quadrat = 1:5)

  # mutate(plot=paste(data$park[1], data$site_id[1], sep='_'),
  #        transect=paste(plot, transect_id, sep='_'),
  #        sample_id=paste(transect, 1:5, sep='_'))
}

get_sample_locations <- function(data, type, ...) {
  # The 'precursors' already finds the center of each quadrat, as such this is
  # really just a pass-through.
  data %>%
    rename(sample_x = site_x, sample_y = site_y) # %>%
  # mutate(site_x=NA, site_y=NA)
}
