library(tidyverse)
library(lubridate)

# d_fine_fuels <- read_csv('assets/uplands-data/CASP/fine_fuels_0to25.csv')
# d_pts_raw <- read_csv('assets/uplands-data/CASP/raw/Cover - Points (metric)_XPT 25 YR.csv')
d_rich_raw <- read_csv("assets/uplands-data/CASP/raw/Cover - Species Composition (metric)_XPT 25 YR.csv")
# after 2010 it's 2.5m belt just WITHIN the plot, BUT after 2010 they also tallied all spp appearing elsewhere in the plot

join_pts_data <- FALSE

spp_pts_raw <- d_pts_raw$`Species Symbol` %>%
  table() %>%
  names()
spp_rich_raw <- d_rich_raw$`Species Symbol` %>%
  table() %>%
  names()
spp_pts_raw[!spp_pts_raw %in% spp_rich_raw]

non_spp_codes <- c("BARE1", "BOLE1", "LITT1", "ROCK1", "WOOD1")

format_spp_data <- function(data) {

  # sq_m_per_overstory_plot <- 20 * 50
  sq_m_per_acre <- 4046.85642

  data %>%
    mutate(
      date = parse_date_time(Date, orders = c("%m/%d/%Y %H:%M:%S Op", "%m/%d/%y %H:%M")),
      plot_name = str_extract(`MacroPlot Name`, "(?<=(B:|C |C:)).*"),
      plot_name = str_replace(plot_name, ":", "-"),
      plot_name = str_replace(
        plot_name,
        "(?<=-)\\d+",
        str_pad(str_extract(plot_name, "(?<=-)\\d+"), 2, "left", "0")
      ),
      sample_year = ifelse(`Monitoring Status` == "00Pre", "01YR00", `Monitoring Status`),
      sample_year = as.integer(str_extract(sample_year, "(?<=YR)\\d+")),
      plot_type = str_extract(`MacroPlot Name`, "[BC](?=(:| ))"),
      # the total area sampled prior to 2010 was 10X50, while after 2010 it would have been 20X50
      sq_m_per_belt = ifelse(year(date) >= 2010, 20 * 50, 10 * 50), # TODO: investigate 'diminishing returns' of search area
      n_acres_per_belt = sq_m_per_belt / sq_m_per_acre,
      cal_year_is_gte_2010 = as.integer(year(date) >= 2010)
    )
  # group_by(`Monitoring Status`, `MacroPlot Name`, plot_type, plot_name, sample_year, sq_m_per_belt) %>%
  # summarise(n_spp = n_distinct(`Species Symbol`, na.rm = TRUE),
  #           .groups = "drop")
}
# test_cases <- c('FPIJE1D09-3', 'FPIJE1D09-13')



d_rich <- format_spp_data(d_rich_raw)

if (join_pts_data) {
  d_rich <- d_rich %>%
    bind_rows(
      d_pts_raw %>%
        filter(!`Species Symbol` %in% non_spp_codes) %>%
        format_spp_data()
    )
}

d_rich <- d_rich %>%
  group_by(
    `Monitoring Status`, `MacroPlot Name`, plot_type, plot_name, sample_year,
    sq_m_per_belt, n_acres_per_belt, cal_year_is_gte_2010
  ) %>%
  summarise(
    n_spp = n_distinct(`Species Symbol`, na.rm = TRUE),
    .groups = "drop"
  )


# tmp = d_rich_raw %>%
#   mutate(date = parse_date_time(Date, '%m/%d/%Y %H:%M:%S Op')) %>%
#   mutate(cal_year = year(date)) %>%
#   distinct(cal_year, `Monitoring Status`) %>%
#   arrange(`Monitoring Status`)
# View(tmp %>% filter(cal_year >= 2010)) #  %>% filter(`Monitoring Status` == "YR10")

# d_rich <- d_rich_raw %>%
#   mutate(
#     date = parse_date_time(Date, '%m/%d/%Y %H:%M:%S Op'),
#     plot_name = str_extract(`MacroPlot Name`, "(?<=(B:|C |C:)).*"),
#     plot_name = str_replace(plot_name, ':', '-'),
#     sample_year = ifelse(`Monitoring Status` == '00Pre', '01YR00', `Monitoring Status`),
#     sample_year = as.integer(str_extract(sample_year, "(?<=YR)\\d+")),
#     plot_type = str_extract(`MacroPlot Name`, "[BC](?=(:| ))"),
#     # the total area sampled prior to 2010 was 10X50, while after 2010 it would have been 20X50
#     sq_m_per_belt = ifelse(year(date) >= 2010, 20*50, 10*50)
#   ) %>%
#   # distinct(`MacroPlot Name`, `Monitoring Status`, plot_type, plot_name, sample_year) %>% View()
#   group_by(`Monitoring Status`, `MacroPlot Name`, plot_type, plot_name, sample_year, sq_m_per_belt) %>%
#   summarise(n_spp = n_distinct(`Species Symbol`, na.rm = TRUE),
#             .groups = "drop")

ggplot(d_rich, aes(x = sample_year, y = n_spp / sq_m_per_belt, group = plot_name)) +
  facet_wrap(~plot_type) +
  geom_line(
    alpha = 0.5
  ) +
  geom_point(
    pch = 21, fill = "white"
  ) +
  labs(x = "Sample year", y = "Number of species")

d_plot_info <- read_csv("assets/uplands-data/CASP/metadata/plot_info.csv")
tmp <- d_rich %>%
  left_join(
    d_plot_info %>%
      # filter(grepl('PIP', plot_name))
      select(park_name, plot_name, plot_type2 = plot_type)
  )
# LJZ: ask Dan about this stuff.

# These are plots for which we have other forest structural info, but not
# species composition info. There were likely plots that had no species in
# the understory, so they didn't include that form. So some are possibly
# true zeros (no richness, no species in the understory).
d_plot_info$plot_name[!d_plot_info$plot_name %in% d_rich$plot_name] %>% unique()
# These are plots in the species comp data that do not appear in the
# forest structure for Tahoe. Some of these are probably Plumas and Malakoff.
d_rich$plot_name[!d_rich$plot_name %in% d_plot_info$plot_name] %>% unique()
# tmp %>% mutate(matches = ifelse(plot_type == plot_type2, T, F)) %>% View()
d_plot_info %>% filter(grepl("1D08", plot_name))
write_csv(d_rich, "assets/uplands-data/CASP/spp_comp_0to25.csv")

range(d_rich$n_spp)
hist(d_rich$n_spp)
