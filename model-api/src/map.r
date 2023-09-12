library(leaflet)
library(leaflet.extras)
library(htmltools)


GCS_EPSG <- 4326 # GCS_EPSG
MAPBOX_ACCESS_TOKEN <- # default public token
  "pk.eyJ1IjoibHphY2htYW5uIiwiYSI6ImNpcW1oODczZTAwcjBnc2pmaGRhYjVudHIifQ.LeGAHvHXv36-vorTmuNtSg"

style_url <- function(style, access_token = MAPBOX_ACCESS_TOKEN) {
  paste0(
    "https://api.mapbox.com/styles/v1/mapbox/",
    style,
    "/tiles/256/{z}/{x}/{y}?access_token=",
    access_token
  )
}

mapbox_streets <- style_url("streets-v10")
mapbox_outdoors <- style_url("outdoors-v10")
mapbox_satellite <- style_url("satellite-v9")

create_leaflet_map <- function() {
  # Map widget plus Mapbox basemaps.
  # options = providerTileOptions(minZoom = 4, maxZoom = 20)
  leaflet() %>%
    addTiles(urlTemplate = mapbox_streets, group = "Streets") %>%
    addTiles(urlTemplate = mapbox_outdoors, group = "Terrain") %>%
    addTiles(urlTemplate = mapbox_satellite, group = "Satellite") %>%
    addFullscreenControl()
}

xy_coords_to_lat_lon <- function(data, coord_cols_prefix, pcs_epsg) {
  new_coord_cols <- c("lon", "lat") %>%
    map_chr(function(x) paste(coord_cols_prefix, x, sep = "_"))
  lat_lon <- data %>%
    select(dplyr::matches(paste0(coord_cols_prefix, "_[xy]$"))) %>%
    st_as_sf(coords = names(.)) %>%
    st_set_crs(pcs_epsg) %>%
    st_transform(GCS_EPSG) %>%
    st_coordinates() %>%
    as_tibble(.name_repair = NULL) %>%
    rename_(.dots = setNames(c("X", "Y"), new_coord_cols))
  data %>% bind_cols(lat_lon)
}
