cluster_init <- function(cores, ...) {
  library(multidplyr)

  new_cluster(cores) %>%
    cluster_library(c("tidyverse", "magrittr")) %>%
    cluster_copy(c("convertWidth", "convertHeight", unlist(...)))
}
