#!/usr/bin/env Rscript

"
Compile model summaries beneath a given directory.

Usage:
  model-api/src/compile-results.R <path> [--max-rhat=<dbl> --ppp-mean=<dbl> --ppp-sd=<dbl> --filter --distinct --which-rhat=<str>]
  model-api/src/compile-results.R (-h | --help)
  model-api/src/compile-results.R --version

Arguments:
  path                 The relative (project) path to an analysis config file

Options:
  --max-rhat=<dbl>     Maximum allowable Rhat, if filtered [default: 1.1]
  --ppp-mean=<dbl>     Interval for an acceptable posterior predictive p-value of the mean [default: 0.9]
  --ppp-sd=<dbl>       Interval for an acceptable posterior predictive p-value of the std. dev. [default: 0.9]
  --filter             Filter results for model meeting optional criteria directory?
  --distinct           Filter results for the 'best' model for each response variable?
  --which-rhat=<str>   Which Gelman diagnostic to use [default: gelman_diag_pt_est]
" -> doc

# Traverses subdirectories recursively. Consolidates fits in a single file. Note
# that NAs in model checking attributes indicate failures.

source("model-api/src/wrangle.r")

if (!interactive()) {
  args <- docopt::docopt(doc)
  print(args)
  path <- args$path
  max_rhat <- as.numeric(args$`--max-rhat`)
  ppp_mean <- as.numeric(args$`--ppp-mean`)
  ppp_sd <- as.numeric(args$`--ppp-sd`)
  apply_filters <- args$`--filter`
  return_distinct <- args$`--distinct`
  which_rhat <- args$`--which-rhat`
} else {
  path <- "assets/uplands-output/SCPN2/GRASS" # "assets/_output/M4MD"
  max_rhat <- 1.1
  ppp_mean <- 0.9
  ppp_sd <- 0.9
  apply_filters <- FALSE
  return_distinct <- FALSE
  which_rhat <- "gelman_diag_pt_est"
}

mod_summaries_raw <- get_mod_summary_files(path, parse_coda_output = TRUE)
message(sprintf("See %s for more!", file.path(path, "mod-summaries.csv")))

if (apply_filters) {
  tail_area_mean <- 0.5 * (1 - ppp_mean)
  tail_area_sd <- 0.5 * (1 - ppp_sd)
  mod_summaries <- mod_summaries_raw %>%
    filter(
      get(which_rhat) <= max_rhat,
      p_mean >= tail_area_mean, p_mean <= 1 - tail_area_mean,
      p_sd >= tail_area_sd, p_sd <= 1 - tail_area_sd
    )
  write_csv(mod_summaries, file.path(path, "mod-summaries-filtered.csv"))
  message(sprintf(
    "Filtered results appear here: %s",
    file.path(path, "mod-summaries-filtered.csv")
  ))

  if (return_distinct) {
    mod_summaries_uniq <- mod_summaries %>%
      group_by(coalesce(response, hits)) %>%
      filter(ppl == min(ppl)) %>%
      arrange(gelman_diag) %>%
      slice(1) %>%
      ungroup()
    fname <- file.path(path, "mod-summaries-filtered-distinct.csv")
    write_csv(mod_summaries_uniq, fname)
    message(sprintf("Filtered results appear here: %s", fname))
  }
}

# combine_mod_summaries <- function(folder, add_dir_metadata = FALSE,
#                                   write_summaries = FALSE) {
#
#   mod_summaries_files <-
#     list.files(path=folder,
#                pattern="^mod-summaries", recursive = TRUE, full.names = TRUE)
#   mod_summaries_list <- lapply(mod_summaries_files, function(x) {
#     if(add_dir_metadata) {
#       dir_metadata <- strsplit(x, '/|//')[[1]]
#       read_csv(x) %>%
#         mutate(network = dir_metadata[2], park = dir_metadata[3])
#     } else {
#       read_csv(x)
#     }
#   })
#   if(write_summaries) {
#     do.call(bind_rows, mod_summaries_list) %T>%
#       write_csv(file.path(folder, 'combined-mod-summaries.csv'))
#   } else {
#     do.call(bind_rows, mod_summaries_list)
#   }
#
# }

# tmp <- combine_mod_summaries('output', add_dir_metadata = TRUE)
# write_csv(tmp, 'output/every-model-summary.csv')
