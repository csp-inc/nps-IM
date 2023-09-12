library(tidyverse)
tmp <- sessionInfo()

computer <- "macbook" # 'macbook' or 'monsoon'

otherPkgs <- lapply(tmp$otherPkgs, function(x) {
  tibble(
    package = x$Package, version = x$Version, device = computer,
    component = "otherPkgs"
  )
}) %>% bind_rows()

loadedOnly <- lapply(tmp$loadedOnly, function(x) {
  tibble(
    package = x$Package, version = x$Version, device = computer,
    component = "loadedOnly"
  )
}) %>% bind_rows()

tmp_df <- bind_rows(otherPkgs, loadedOnly) %>%
  mutate(platform = tmp$platform, r_version = tmp$R.version$version.string)

write_csv(tmp_df, paste0(
  "output/",
  computer, "-package-versions-latest.csv"
))


# ====

if (computer == "macbook") {
  bind_rows(tmp_df, read_csv("output/monsoon-package-versions.csv")) %>%
    # distinct(r_version, platform)
    select(-r_version, -platform) %>%
    # mutate(entry = 1:n()) %>%
    spread(device, version) %>%
    mutate(is_same = macbook == monsoon) %>%
    filter(!is_same)
}
