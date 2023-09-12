library(tidyverse)
library(ggthemes)

d <- read.csv("output/every-model-summary.csv", stringsAsFactors = FALSE) %>%
  mutate( # response = paste(hits, trials, sep = '/'),
    response = ifelse(is.na(response), paste(hits, trials, sep = "/"), response) # ,
    # gelman_diag_lte5 = gelman_diag <= 5
  )
# View(d)
# d %>% filter(park == "CARE", response == "non_native_hits/trials") %>% View
d %>%
  group_by(network, park, response) %>%
  summarise(prop_na = sum(is.na(p_mean)) / n()) %>%
  filter(prop_na == 1)

# d <- summaries

ggplot(d, aes(x = group_level_effects, y = dic)) +
  facet_grid(likelihood ~ var_type) +
  geom_boxplot() +
  geom_point(aes(color = gelman_diag_lte5), position = position_jitter(width = .1))

d %>%
  filter(gelman_diag < 1.5, p_sd < .8, p_sd > .2, p_mean < .8, p_mean > .2) %>%
  # slice(1:2) %>% View
  mutate(
    source = paste("/scratch/lz62/uplands-ci", sub("./", "", full_path), sep = "/"),
    # destination = dirname(sub('./', '', full_path)),
    destination = file.path("/projects/lci/public_data", sub("./", "", full_path))
  ) %>%
  # mutate(cmd = paste('mkdir -p', dirname(destination))) %>%
  mutate(cmd = paste("mkdir -p", dirname(destination), "&&", "cp -R", source, dirname(destination))) %>% # paste0(, '/')
  pull(cmd) %>%
  paste(collapse = " && ") %>%
  write_lines(., "sandbox/cmd-1.txt")



# d %>% filter(network == 'SCPN') %>% pull(response) %>% table
d %>% filter(response == "SurfaceRating")


d %>%
  group_by(response) %>%
  summarise(perc_na = sum(is.na(dic)) / n() * 100)
d %>%
  filter(gelman_diag < 1.5, p_sd < .8, p_sd > .2, p_mean < .8, p_mean > .2) %>%
  pull(response) %>%
  table()

d %>% filter(network == "SCPN", response == "SpeciesRichness")


d_in <- d %>%
  filter(response == "ShrubCoverClass_10m") %>%
  filter(gelman_diag < 1.5, p_sd < .8, p_sd > .2, p_mean < .8, p_mean > .2) %>%
  select(
    ppl, dic, gelman_diag, likelihood, response,
    matches("^p_|^group_|^var_|covariates")
  ) %>%
  gather(crit, crit_val, -response, -matches("^group_|^var_|^like|covariates")) %>%
  mutate(additional_covariates = gsub(", ", ",\n", additional_covariates)) %>%
  as.data.frame()
d_in

ggplot(d_in, aes(x = additional_covariates, y = crit_val)) +
  facet_grid(crit ~ likelihood, scales = "free") +
  geom_point(aes(fill = var_level),
    color = "black", size = 3,
    alpha = .5, pch = 21
  ) +
  geom_point(aes(shape = var_type), size = 2) +
  geom_point(
    data = d_in %>% filter(group_level_effects == "b0-b1"),
    size = 5, pch = 21
  ) +
  # theme_hc() +
  scale_shape_manual("Variance\ntype", values = c(4, 21)) +
  # scale_color_hc() +
  scale_fill_hc(name = "Variance\nlevel") +
  labs(y = "Value", x = "Additional covariates")
