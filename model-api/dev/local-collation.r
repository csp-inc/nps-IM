library(tidyverse)

file_list <- list.files(
  path = "output-tmp", pattern = "mod-summary.csv",
  full.names = TRUE, recursive = TRUE
)

summaries_list <- lapply(file_list, function(x) {
  dir_metadata <- strsplit(x, "/|//")[[1]]
  read_csv(x) %>%
    mutate(
      network = dir_metadata[2], park = dir_metadata[3],
      results_path = dirname(x)
    )
})
summaries <- bind_rows(summaries_list)
write_csv(summaries, "output-tmp/all-model-summaries.csv")

get_response <- function(x, position = 5) {
  sapply(x, function(z) str_split(z, "/")[[1]][position])
}
d <- read_csv("output-tmp/all-model-summaries.csv") %>%
  mutate(response = get_response(results_path))
d %>%
  group_by(network, park, response) %>%
  summarise(n_na = sum(is.na(p_mean))) %>%
  ungroup()

d %>%
  select(-results_path) %>%
  filter(gelman_diag < 1.5, p_sd < .8, p_sd > .2, p_mean < .8, p_mean > .2) %>%
  group_by(network, park, response) %>%
  tally() %>%
  ungroup() %>%
  as.data.frame()

d %>% filter(response == "soil-stability")


View(summaries)


# ======= one of these, I think, was used, check bash history ------
# THIS ONE! find output -type f | grep -i mod-summary.csv | xargs cp {} --parents -t output-tmp
