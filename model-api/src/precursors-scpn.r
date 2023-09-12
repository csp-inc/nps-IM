library(lubridate)
source("model-api/src/wrangle.r")
source("model-api/src/sample-locations-scpn.r")


dir.create("data/SCPN/modified", showWarnings = FALSE)
site_locs <- load_data("data/SCPN/PefoPlots.xlsx") %>%
  bind_rows(
    tibble(
      Park = "PEFO",
      SiteIdLong = c("PEFOS24", "PEFOC21", "PEFOC21"),
      Location = c("CS", "CS", "CF"),
      Easting = c(612567, 614410, 614397),
      Northing = c(3860105, 3869150, 3869097)
    )
  ) %>%
  mutate(
    Plot = str_sub(SiteIdLong, 5, 7),
    Transect = str_sub(Location, 1, 1),
    EcoSite = ifelse(grepl("^C", Plot), "PEFO_C", "PEFO_S")
  ) %>%
  filter(!grepl("C$", Location))


# plots_missing_monuments <-
#   site_locs %>% group_by(Plot) %>% tally %>% filter(n!=6) %>% pull(Plot)
# sf::st_as_sf(site_locs %>% filter(Plot %in% plots_missing_monuments),
#              coords=c('Easting', 'Northing'))


quadrat_locs <- site_locs %>%
  # filter(!Plot %in% plots_missing_monuments) %>%  # stopgap
  group_by(Park, EcoSite, Plot, Transect) %>%
  # Get the centroids of each sample (i.e., quadrat).
  do(quadrat_centroids(.)) %>%
  ungroup() %T>%
  write_csv("data/SCPN/modified/PefoPlots.csv")
# select(matches('Quadrat')) %>% plot


read_csv("data/SCPN/PEFOFunctionalGroup08312018.csv") %>%
  mutate(Plot = str_sub(Plot, 6, 8)) %>%
  filter(Park == "PEFO", Plot %in% unique(quadrat_locs$Plot)) %>%
  left_join(quadrat_locs %>% distinct(Park, Plot, EcoSite)) %T>%
  write_csv("data/SCPN/modified/filtered-PEFOFunctionalGroup08312018.csv")

# read_csv('data/SCPN/FunctionalGroupDataCSP.csv') %>%
#   filter(Park=='PEFO', Plot %in% unique(quadrat_locs$Plot)) %>%
#   left_join(quadrat_locs %>% distinct(Park, Plot, EcoSite)) %T>%
#   write_csv('data/SCPN/modified/filtered-FunctionalGroupDataCSP.csv')

readRDS("data/SCPN/PEFOspeciesRich08312018.rds") %>% # SpeciesRichPEFOQuad.csv
  mutate(
    Plot = str_sub(TransectQuadratID, 6, 8),
    Transect = str_sub(TransectQuadratID, -3, -3),
    Quadrat = as.numeric(str_sub(TransectQuadratID, -1, -1))
  ) %>%
  left_join(quadrat_locs %>% distinct(Park, Plot, EcoSite)) %T>%
  # select(-X1) %>% filter(Plot != 'C03') %T>%
  saveRDS(., "data/SCPN/modified/PEFOspeciesRich08312018.rds")

# read_csv('data/SCPN/SCPN_PEFO_cov.csv') %>%
#   filter(Plot!='S23')  %T>%
#   write_csv('data/SCPN/modified/SCPN_PEFO_cov.csv')


# Get and format gap size information.
prep_gap_size_data <- function(data, gap_size_col, gap_info = "proportion",
                               gap_start_col = "Gap_Start", gap_end_col = "Gap_End") {
  data %<>%
    group_by(Park, Plot, TransectID, TransectLetter, EventYear) %>%
    summarise(gap_proportion = sum(get(gap_size_col)) / 5000) %>%
    ungroup()

  names(data)[names(data) == "gap_proportion"] <- sub("Length_Segment_cm", "Prop", gap_size_col)

  data
}


pefo_basal_gaps <- load_data("data/SCPN/PEFOBasal.csv") # %>%
# group_by(Park, Plot, EventYear) %>%
# mutate(n_transects=n_distinct(TransectLetter)) %>%
# ungroup %>%
# filter(n_transects != 2)

# Develop basal and canopy gap proportion data.
prep_gap_size_data(pefo_basal_gaps, "BasalGapLength_Segment_cm", "proportion") %T>%
  write_csv(paste0("data/SCPN/", "PEFO_gaps_props.csv"))
