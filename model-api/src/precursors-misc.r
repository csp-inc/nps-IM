library(tidyverse)

d <- read_csv("data/MISC/leatherback_data_5-21-18.csv") %>%
  mutate(stratum = "none", scope = "leatherback-turtles") %>%
  select(-X5, -X6) %T>%
  write_csv("data/MISC/modified/leatherback_data_5-21-18.csv")

d %>%
  filter(stck != "Fla") %T>%
  write_csv("data/MISC/modified/leatherback_data_sans_FL.csv")
