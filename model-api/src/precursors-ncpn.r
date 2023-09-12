source("model-api/src/wrangle.r")
source("model-api/src/patch.r") # tidyjson
source("model-api/src/sample-locations-ncpn.r")
library(magrittr)
library(lubridate)


plants_table <- readRDS("C:/R/NCPN_Master_Species_Database/NCPNplants.rds") %>%
  mutate(
    plant = "plant",
    is_plant = (plant == "plant")
  ) %>%
  select(
    plant_code = Master_PLANT_Code,
    matches("is_plant")
  )


# TODO: add biological soil crust, bare soil, and
plants_table <- load_data("data/NCPN/tlu_NCPN_Plants.xlsx") %>%
  # Binary attributes for lifeform, duration, and nativity.
  mutate(
    plant = "plant",
    is_plant = (plant == "plant"),
    is_graminoid = (Lifeform == "Graminoid"),
    is_perennial = (Duration == "Perennial"),
    is_native = (Nativity == "Native")
  ) %>%
  # Is a species native, perennial, and graminoid? If so, is it C3 or C4?
  rowwise() %>%
  mutate(
    is_npg = all(is_native, is_perennial, is_graminoid),
    is_c3_npg = all(is_npg, Photosynthetic_Pathway == "C3"),
    is_c4_npg = all(is_npg, Photosynthetic_Pathway == "C4"),
    # Perennial grass cover.
    is_pg = all(is_perennial, is_graminoid),
    # Non-native plant cover.
    is_non_native = Nativity == "NonNative",
    # Shrub cover.
    is_shrub = Lifeform == "Shrub"
  ) %>%
  ungroup() %>%
  # Select just the grass-related attributes we need for analysis.
  select(
    plant_code = Master_PLANT_Code,
    matches("is_c\\d_npg|is_pg|is_non_native|is_shrub|is_plant")
  )

care_pi_data <- load_data("data/NCPN/CARE_PI.rds") %>%
  # Coerce `Plot_ID`, if character, to numeric data type.
  mutate(Plot_ID = as.numeric(Plot_ID)) %>%
  # Filter out plots removed from the long-term sampling design.
  filter(!Plot_ID %in% c(78, 207))

# The name of columns that will *not* be gathered below. In other words, these
# are our identifying (or grouping) variables.
design_vars <- c(
  "Unit_Code", "Master_Stratification", "Start_Date", "Plot_ID",
  "Transect", "Point"
)

# BSC and bare soil summaries.
pi_surface <- care_pi_data %>%
  select(one_of(design_vars), Surface) %>%
  mutate(
    is_BSC = (Surface == "CY" | Surface == "LC" | Surface == "M"),
    is_baresoil = (Surface == "S")
  )

# Gather hits.
spp_code_fields <- # species code at each canopy level (from top to surface)
  add_regex_boundaries(c("Top", "LCS\\d+", "Surface")) %>%
  paste(collapse = "|")
pi_spp_long <- care_pi_data %>%
  select(one_of(design_vars), matches(spp_code_fields)) %>%
  gather(canopy_level, plant_code, -one_of(design_vars))

# Gather status info (i.e., was the individual encountered dead or alive?).
status_fields <- # status (alive or dead) for each observation
  add_regex_boundaries(c("Alive", "LCA\\d+", "Surface_Alive")) %>%
  paste(collapse = "|")
pi_is_alive_long <- care_pi_data %>%
  select(one_of(design_vars), matches(status_fields)) %>%
  gather(status, is_alive, -one_of(design_vars))

# Join the long-format hits and status data using a lookup table for their
# respective key fields (`field_lookup`) and, finally, join in NCPN plants info.
field_lookup <- tibble( # the info necessary to crosswalk codes and status
  canopy_level = grep(spp_code_fields, names(care_pi_data), value = TRUE),
  status = grep(status_fields, names(care_pi_data), value = TRUE)
)
pi_data_long <- pi_spp_long %>%
  left_join(field_lookup) %>%
  left_join(pi_is_alive_long) %>%
  left_join(plants_table)


# Construct the queries necessary to develop the final response info.
final_pi_data <- pi_data_long %>% # note: fairly slow!
  # Asks whether a given records meets the following two criteria: namely, that
  # a 'hit' is 1) a native, perennial C3 or C4 grass, and 2) alive.
  rowwise() %>%
  mutate(
    is_live_plant = all(is_plant, is_alive),
    is_live_c3_npg = all(is_c3_npg, is_alive),
    is_live_c4_npg = all(is_c4_npg, is_alive),
    is_live_pg = all(is_pg, is_alive),
    is_live_shrub = all(is_shrub, is_alive),
    is_live_non_native = all(is_non_native, is_alive)
  ) %>%
  ungroup() %>%
  # Asks whether anything encountered (at any canopy level!) at a given point-
  # intercept meets the criteria above (1 = 'yes', 0 = 'no').
  group_by_(.dots = design_vars) %>%
  summarise(
    any_is_live = bin_any(is_live_plant),
    any_is_live_c3_npg = bin_any(is_live_c3_npg),
    any_is_live_c4_npg = bin_any(is_live_c4_npg),
    any_is_live_pg = bin_any(is_live_pg),
    any_is_live_shrub = bin_any(is_live_shrub),
    any_is_live_non_native = bin_any(is_live_non_native)
  ) %>%
  ungroup() %>%
  left_join(pi_surface) %>%
  # Summarise presence/absence data at each point intercept as hits/trials
  # along each transect.
  group_by_(.dots = design_vars[design_vars != "Point"]) %>%
  summarise(
    live_hits = sum(any_is_live),
    c3_grass_hits = sum(any_is_live_c3_npg),
    c4_grass_hits = sum(any_is_live_c4_npg),
    pg_hits = sum(any_is_live_pg),
    shrub_hits = sum(any_is_live_shrub),
    non_native_hits = sum(any_is_live_non_native),
    BSC_hits = sum(is_BSC),
    baresoil_hits = sum(is_baresoil),
    trials = 100
  ) %>%
  ungroup() %T>%
  save_data("data/NCPN/modified/Point_Intercept_CARE_summaries.csv")

care_loc_info <- "data/NCPN/CARE_spatial_data/CARE_transect_ends.shp"
care_transect_ends_sf <- load_data(care_loc_info, stringsAsFactors = FALSE) %>%
  filter(
    # Remove '2nd rebar' records --- they're not necessary for our purposes.
    str_detect(Comment, "2nd rebar") %in% c(NA, FALSE),
    # As above, filter out plots removed from the long-term sampling design.
    !PlotID %in% c(78, 207)
  ) %>%
  # Create transect number and position fields (`?<=`` is a positive lookbehind).
  mutate(
    transect_lab = Transect,
    Transect = as.numeric(str_extract(transect_lab, "(?<=^T)\\d")),
    transect_position = str_extract(transect_lab, regex("(?<=^T\\d_).*"))
  ) %>%
  # Make any patches, as necessary.
  apply_patch(care_loc_info, get_patch_info()) %>%
  select(
    Plot_ID = PlotID, Transect, # careful! 'PlotID' vs 'Plot_ID'...
    transect_lab, transect_position, geometry
  )

# Get point-intercept locations for each transect in each plot.
sample_locs <- care_transect_ends_sf %>%
  group_by(Plot_ID) %>%
  do(get_transect_params(.)) %>%
  ungroup() %T>%
  write_csv("data/NCPN/modified/CARE_transect_centers.csv")

load_data("data/NCPN/CARE.P.cov.rds") %>%
  mutate(
    Unit_Code = str_sub(Unit_Plot, 1, 4),
    Plot_ID = as.numeric(str_sub(Unit_Plot, 5, 8))
  ) %T>%
  write_csv("data/NCPN/modified/CARE.P.cov.csv")

# Prepare the exoctic frequency data.
visitation_dates <- load_data("data/NCPN/UplandLocations.rds") %>%
  filter(Unit_Code == "CARE", !Plot_ID %in% c(78, 207)) %>%
  mutate(Plot_ID = as.numeric(Plot_ID), Year = year(Start_Date)) %>%
  group_by(Unit_Code, Master_Stratification, Plot_ID, Start_Date, Year) %>%
  expand(Transect = seq(1, 3)) %>%
  ungroup()
exotic_freq_data_raw <- load_data("data/NCPN/CARE_EF.rds")
exotic_species <- unique(exotic_freq_data_raw$Species)
exotic_freq_data <- exotic_freq_data_raw %>%
  # pipe_assign('exotic_species', unique(.$Species)) %>%
  mutate(Plot_ID = as.numeric(Plot_ID), Year = as.numeric(Year))
efd_by_spp_long <- exotic_freq_data %>%
  mutate(exotic_hits = rowSums(.[grep("M[0-9]{1,2}", names(.))]))
efd_by_spp_wide <- efd_by_spp_long %>%
  select(one_of(design_vars[design_vars != "Point"]), Year, exotic_hits, Species) %>%
  spread(Species, exotic_hits)
efd_across_spp_wide <- efd_by_spp_long %>%
  group_by_(.dots = design_vars[design_vars != "Point"]) %>%
  summarise_at(vars(matches("M[0-9]{1,2}")), list(~any)) %>%
  ungroup() %>%
  mutate(any_exotics = rowSums(.[grep("M[0-9]{1,2}", names(.))])) %>%
  select(-matches("M[0-9]{1,2}"))
efd_combined <- left_join(efd_by_spp_wide, efd_across_spp_wide)

# any_exotics=rowSums(.[names(.) %in% exotic_species] > 0),
zero_for_na <- function(x) ifelse(is.na(x), 0, x)
exotic_freq_data_complete <- (visitation_dates %>% filter(Year >= 2011)) %>%
  left_join(efd_combined) %>%
  mutate_at(vars(c(exotic_species, "any_exotics")), list(~zero_for_na)) %>%
  mutate(exotics_trials = 10) %>%
  filter(!Plot_ID %in% c(78, 207)) %T>%
  write_csv("data/NCPN/modified/CARE_exotic_frequency.csv")

# exotic_freq_data_complete %>%
#   filter(any_exotics>7)

# Shrub density by species
shrub_species_data_raw <- load_data("data/NCPN/CARE_ShrubSpecies.rds")
shrub_species <- unique(shrub_species_data_raw$Species)
# pipe_assign('shrub_species', unique(.$Species)) %>%
shrub_species_data <- shrub_species_data_raw %>%
  mutate(Plot_ID = as.numeric(Plot_ID), Year = as.numeric(Visit_Year))
shrub_by_spp_wide <- shrub_species_data %>%
  select(one_of(design_vars[design_vars != "Point"]), Year, Total, Species) %>%
  spread(Species, Total)

zero_for_na <- function(x) ifelse(is.na(x), 0, x)
shrub_spp_data_complete <- (visitation_dates %>% filter(Year >= 2009)) %>%
  left_join(shrub_by_spp_wide) %>%
  filter(!Plot_ID %in% c(78, 207))
shrub_spp_data_complete[is.na(shrub_spp_data_complete)] <- 0
CARE_shrub_spp_counts <- shrub_spp_data_complete %>%
  write_csv("data/NCPN/modified/CARE_shrub_spp_counts.csv")



# # Get and format gap size information.
# prep_gap_size_data <- function(data, gap_size_col, gap_info='proportion',
#                                gap_start_col='Gap_Start', gap_end_col='Gap_End') {
#   gap_size_data <- as_tibble(data, .name_repair = NULL) %>%
#     # Remove non-targets and gaps at either transect end.
#     filter(Master_Stratification != 'nontarget') %>%
#     # Update fields as necessary, and calculate gap proportion per transect
#     mutate(Plot_ID=as.numeric(Plot_ID), year=lubridate::year(Start_Date))
#
#   if(gap_info == 'proportion') {

#     gap_size_data %<>%
#       group_by(Unit_Code, Master_Stratification,
#                Plot_ID, Transect, year) %>%
#       summarise(gap_proportion=sum(get(gap_size_col))/5000) %>%
#       ungroup %>%
#       filter(gap_proportion <= 1)  # TODO: remove. This should only be temporary!
#
#     names(gap_size_data)[names(gap_size_data)=='gap_proportion'] <-
#       sub('Size', 'Prop', gap_size_col)
#
#   } else if (gap_info == 'size') {

#     gap_size_data %<>%
#       # Update fields as necessary, and calculate gap midpoints.
#       mutate(Gap_Midpoint=get(gap_end_col) -
#                (get(gap_end_col)-get(gap_start_col))/2,
#              Intersects_Transect_Ends=get(gap_start_col) %in% c(0, 5000) |
#                get(gap_end_col) %in% c(0, 5000),
#              GapSizeClass=cut(get(gap_size_col),
#                               breaks=quantile(get(gap_size_col), probs = seq(0, 1, 0.1)),
#                               include.lowest = TRUE))
#     # print(gap_size_data %>%
#     #         group_by(GapSizeClass) %>%
#     #         summarise(num_tossed=sum(Intersects_Transect_Ends)))
#     gap_size_data %<>%
#       select(-GapSizeClass)  #-get(gap_start_col), -get(gap_end_col)
#
#   }
#
#
#   gap_size_data
#   # write_csv(gap_size_data,
#   #           paste0('data/NCPN/', deparse(substitute(data)), '.csv'))
#
# }
#
#
# load('data/NCPN/CARE_bgaps_raw.rda')
# load('data/NCPN/CARE_cgaps_raw.rda')
# # Develop basal and canopy gap proportion data.
# prep_gap_size_data(CARE_bgaps_raw, 'BasalGapSize', 'proportion') %>%
#   full_join(prep_gap_size_data(CARE_cgaps_raw, 'CanopyGapSize', 'proportion')) %T>%
#   write_csv(paste0('data/NCPN/', 'CARE_gaps_props.csv'))
#
# # Develop basal and canopy gap size data.
# prep_gap_size_data(CARE_bgaps_raw, 'BasalGapSize', 'size') %T>%
#   write_csv(paste0('data/NCPN/', 'CARE_bgap.csv'))
# prep_gap_size_data(CARE_cgaps_raw, 'CanopyGapSize', 'size', 'Start') %T>%
#   write_csv(paste0('data/NCPN/', 'CARE_cgap.csv'))
#
# # oof <-
# prep_gap_size_data(CARE_cgaps_raw, 'CanopyGapSize', 'size', 'Start') %>%
#   group_by(Unit_Code, Master_Stratification, Plot_ID, Transect, year) %>%
#   mutate(CanopySize=lead(Start)-Gap_End) %>%
#   ungroup %>%
#   # filter(CanopySize==-7)
#   # filter(Plot_ID==213, year==2015, Transect==2, Gap_ID %in% c(75, 77)) %>% as.data.frame %>% select(Start, Gap_End)
#   filter(CanopySize>1) %>%
#   pull(CanopySize) %T>%
#   saveRDS('sandbox/vector-of-empirical-canopy-sizes.rds')
#
#
# tmp <- prep_gap_size_data(CARE_cgaps_raw, 'CanopyGapSize', 'size', 'Start') %>%
#   group_by(Unit_Code, Master_Stratification, Plot_ID, Transect, year) %>%
#   summarise(n=n(), mean=mean(CanopyGapSize), median=median(CanopyGapSize))
# # mod <- nls(mean ~ exp(a + b * n), data = tmp, start = list(a = 0, b = 0))
# # n_vec <- seq(min(tmp$n), max(tmp$n), .1)
#
# ggplot() +
#   geom_point(data=tmp, aes(x=n, y=mean), alpha=.25) #+
#   # geom_line(data=tibble(n=n_vec, mean_hat=predict(mod, list(n = n_vec))),
#   #           aes(x=n, y=mean_hat), col='red')

CARE_SS <- load_data("data/NCPN/CARE_SS.rds") %>%
  mutate(Plot_ID = as.numeric(Plot_ID), Year = year(Start_Date)) %>%
  filter(!Plot_ID %in% c(78, 207))
design_vars <- c("Unit_Code", "Master_Stratification", "Start_Date", "Year", "Plot_ID")
ss_long <- CARE_SS %>%
  select(one_of(design_vars), contains("_Pos"), contains("_Rating"), contains("_Veg")) %>%
  gather(variable, value, -one_of(design_vars)) %>%
  mutate(Transect = as.numeric(substring(variable, 2, 2)), var = substring(variable, 7, 12), ID = substring(variable, 2, 5)) %>%
  select(-variable)
ss_wide <- ss_long %>%
  spread(var, value) %>%
  mutate(
    TransectLocation = as.numeric(Pos) / 100,
    SurfaceRating = Rating, VegCoverType = Veg,
    SurfaceRating = ifelse(SurfaceRating == 0, NA, SurfaceRating)
  ) %>%
  select(-ID, -Pos, -Rating, -Veg) %T>%
  write_csv("data/NCPN/modified/CARE_SS.csv")
