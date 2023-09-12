#!/usr/bin/env Rscript

"
Run analyses specified in model config files using JAGS.

Usage:
  model-api/analysis-pipeline.R <file> [--n-adapt=<int> --n-update=<int> --n-iter=<int> --n-cores=<int> --rm-dqs --excl-null --save-mcarray --pass-errors]
  model-api/analysis-pipeline.R (-h | --help)
  model-api/analysis-pipeline.R --version

Arguments:
  file               The relative (project) path to an analysis config file

Options:
  --n-adapt=<int>    The number of iterations for adaptation [default: 1000]
  --n-update=<int>   The number of iterations for updating [default: 5000]
  --n-iter=<int>     The number of iterations to monitor [default: 2500]
  --n-cores=<int>    The number of cores for parallelized analyses [default: 1]
  --rm-dqs           Override (remove) derived quantities?
  --excl-null        Exclude null (time-only) models?
  --save-mcarray     Save mcarray objects to output directory?
  --pass-errors      Pass adapt phase errors and move along to the next model?
" -> doc

if (!interactive()) {
  args <- docopt::docopt(doc)
  yaml_file <- args$file
  n_adapt <- as.numeric(args$`--n-adapt`)
  n_update <- as.numeric(args$`--n-update`)
  n_iter <- as.numeric(args$`--n-iter`)
  n_cores <- as.numeric(args$`--n-cores`)
  override_dqs <- args$`--rm-dqs`
  excl_null_mods <- args$`--excl-null`
  save_mcarray <- args$`--save-mcarray`
  pass_errors <- args$`--pass-errors`
} else {
  yaml_file <- # <<------------------------------------------------- CUSTOMIZE!

  #  './assets/uplands-config/NCPN/DINO/pi-cover-bb-2022.yml'
  #  'assets/uplands-config/NCPN/DINO/soil-stability-2022.yml'
  #  'assets/uplands-config/NCPN/DINO/ex-freq-2022.yml'
  #  'assets/uplands-config/NCPN/DINO/shrub-density-2022.yml'
    'assets/uplands-config/NCPN/DINO/canopy-gaps-2023.yml'
  #  'assets/uplands-config/NCPN/ARCH/pi-cover.yml'
  #  'assets/uplands-config/NCPN/ARCH_GRASS/pi-cover.yml'
    
  n_adapt <- 1000
  n_update <- 5000
  n_iter <- 2500
  n_cores <- 1
  override_dqs <- TRUE
  excl_null_mods <- TRUE
  save_mcarray <- FALSE
  pass_errors <- FALSE
}

# print(yaml_file)
# print(args)

# ==== begin pipeline ===========================================================
source("model-api/src/get-config.r") # TODO: reduce redundancy in this call....
source("model-api/src/apply-config.r")
source("model-api/src/theme.r")
source("model-api/src/serializers.r")
source("model-api/src/fit.r")
source("model-api/demo/utils.R")
source("model-api/src/cluster.r")
source_dir("st-kit")
rjags::load.module("dic")
rjags::load.module("glm")
# check_packages('MCMCpack')

config_list <- get_config(yaml_file, rm_dqs = override_dqs)
analysis_scenarios <- config_list$analysis_scenarios

if (n_cores > 1) cl <- cluster_init(n_cores, ls())
if (interactive()) {
  # warning(paste('You are about to run', nrow(analysis_scenarios),
  #               'analyses. Plan accordingly (or slice)!'))
  warning(
    sprintf(
      "Prior to any added filters below, you are about to run %s jobs.",
      nrow(analysis_scenarios %>%
        {
          `if`(excl_null_mods, filter(., !is.na(additional_covariates)), .)
        })
    )
  )
}

this_slice <- # <<-------------------------------------------------- CUSTOMIZE!
   1:nrow(analysis_scenarios)
analysis_summary <- analysis_scenarios %>%
  mutate(scenario = 1:n()) %>%
  {
    `if`(excl_null_mods, filter(., !is.na(additional_covariates)), .)
  } %>%
  {
    `if`(interactive(), filter(., scenario %in% this_slice), .)
  } %>%
  # If n_cores > 1, then multidplyr handles each call to `rowwise_fit`, sending
  # each analysis scenario (row) to a different core (CPU), otherwise it just
  # takes each row one at a time (serially) using `rowwise`.
  group_by(scenario) %>%
  {
    `if`(
      n_cores > 1,
      partition(., cluster = cl),
      .
    )
  } %>%
  do(rowwise_fit(.,
    yaml_info = config_list$yaml_info,
    yaml_branches = config_list$yaml_branches,
    save_full_jags_obj = save_mcarray,
    destroy_xy = FALSE,
    pass_errors = pass_errors
  )) %>%
  {
    `if`(n_cores > 1, collect(.), ungroup(.))
  }
