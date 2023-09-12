#!/usr/bin/env Rscript

"
Simulate data from a JAGS model.

Usage:
  model-api/demo/sim.R <file> [--n-strata=<int> --n-sites-total=<int> --n-sites-per-year=<int> --n-years=<int> --n-reps-per-site=<int> --seed=<int> --n-possible-sites=<int> --n-morph=<int> --gen-vals=<list>]
  model-api/demo/sim.R (-h | --help)
  model-api/demo/sim.R --version

Arguments:
  file

Options:
  -h --help                 Show CLI help.
  --version                 Show program version.
  --n-strata=<int>          Number of strata.
                            [default: 2]
  --n-sites-total=<int>     Total number of sites allocated to each stratum over
                            all years of sampling.
                            [default: 20]
  --n-sites-per-year=<int>  Number of sites to sample in each year in a stratum.
                            Sites are allocated to years following a 'rotating'
                            schedule, with `n-sites-per-year` sites sampled per
                            year until a total of `n-sites-total` sites is
                            reached. Once `n-sites-total` is met, sites visted
                            earlier are revisted, again on a rotating basis. If
                            `n-sites-total / n-sites-per-year > n-years` there
                            will be some sites that haven't been sampled yet!
                            [default: 4]
  --n-years=<int>           Number of years over which sampling will occur.
                            [default: 8]
  --n-reps-per-site=<int>   Number of subsamples (replicates) per site.
                            [default: 4]
  --seed=<int>              RNG seed.
                            [default: 1337]
  --n-possible-sites=<int>  Number of possible sites in the larger sampling
                            frame (e.g., park unit). This is the sum of the
                            number of sites in all strata.
                            [default: 400]
  --n-morph=<int>           Number of times to apply the operator used to create
                            a random field.
                            [default: 3]
  --gen-vals=<list>         Generating values for the simulation, supplied as a
                            string (e.g., 'list(mu_mu_B0_l=5)'), no spaces!
                            [default: list()]
  --n-attempts=<int>        Number of attempts to simulate values using
                            different RNG seeds.
                            [default: 10]
" -> doc


library(tidyverse)
library(sf)
# devtools::install_version(
#   "spsurvey", version = "4.1.4", repos = "http://cran.us.r-project.org"
# )
library(spsurvey) # WARNING:the newer version is way different!
library(R2jags)
library(runjags) # install.packages('runjags')
library(bayesplot) # install.packages('bayesplot')
library(ggridges)
library(hrbrthemes)
library(glue)

source("model-api/demo/utils.R")
# source('st-kit/eda/histograms.R')
source_dir("st-kit")

args <- docopt::docopt(doc, version = "sim.R v0.1")
print(args)

# extrafont::font_import(prompt = FALSE)
# extrafont::loadfonts()
# hrbrthemes::import_roboto_condensed()
# system("fc-cache -f -v") # install font

print(args$file)

# Fictitious labels.
park_unit <- "ELDO"
start_year <- 1992

# ==== variables ====

if (interactive()) {
  file <- "tmp/count-vars/negative-binomial-linear-b0-b1/model.jags"
  n_strata <- 2
  n_sites_total <- 20
  n_sites_per_year <- 4
  n_years <- 8
  n_reps_per_site <- 4
  seed <- 123
  n_possible_sites <- 400
  n_morph <- 3
  gen_vals <- list()
  n_attempts <- 10
} else {
  file <- args$file
  n_strata <- as.integer(args$`--n-strata`)
  n_sites_total <- as.integer(args$`--n-sites-total`)
  n_sites_per_year <- as.integer(args$`--n-sites-per-year`)
  n_years <- as.integer(args$`--n-years`)
  n_reps_per_site <- as.integer(args$`--n-reps-per-site`)
  seed <- as.integer(args$`--seed`)
  n_possible_sites <- as.integer(args$`--n-possible-sites`)
  n_morph <- as.integer(args$`--n-morph`)
  gen_vals <- eval(parse(text = args$`--gen-vals`))
  n_attempts <- as.integer(args$`--n-attempts`)
}

dir.create(file.path(dirname(file)), showWarnings = FALSE, recursive = TRUE)
print("file dirname")
message(file.path(dirname(file)))
print("working dir")
message(getwd())

message("Here are all of the supplied args")
print(args)

# ==== the program ====

# print(list(n_morph, n_strata, n_sites_per_year, n_sites_total, n_years,
#            n_possible_sites, seed))

set.seed(seed) # for reproducibility

# ---- create a spatial (static) covariate(s), W ----
n_pix_per_side <- floor(sqrt(n_possible_sites)) # set dim of square sample frame
W_raw <- raster::raster(
  xmn = 0, xmx = 1, ymn = 0, ymx = 1,
  ncol = n_pix_per_side, nrow = n_pix_per_side,
  crs = NA
)
raster::values(W_raw) <- rnorm(raster::ncell(W_raw))
W <- W_raw
for (i in 1:n_morph) {
  W <- raster::focal(W, w = matrix(1 / 9, nc = 3, nr = 3), na.rm = TRUE, pad = TRUE)
}

# ---- create the sampling frame ("park") and "strata" within park ----
park_feat <- st_as_sfc(st_bbox(W))
strat_nuclei <- st_sample(park_feat, n_strata)
if (n_strata > 1) {
  v <- st_voronoi(do.call(c, strat_nuclei), envelope = park_feat) %>% # was: envelope = park_feat
    st_collection_extract()
  strat_sfc <- st_intersection(v, park_feat)
  if (sum(st_area(strat_sfc)) != 1) stop("Bad topology!")
} else {
  strat_sfc <- park_feat
}
strat_sf <- st_as_sf(strat_sfc) %>%
  rename(geometry = x) %>%
  mutate(stratum = LETTERS[1:n_strata])
# plot(W, asp = 1); raster::contour(W, add = TRUE)
# plot(strat_sfc, add = TRUE)
# plot(park_feat, add = TRUE)

# ---- the "scaffold" for the complete data, including site-level W ----
scaffold_tbl_raw <- raster::as.data.frame(W, xy = TRUE) %>%
  as_tibble() %>% # just a convenience
  mutate(site_index = row_number()) %>%
  rename(w1 = layer, x_coord = x, y_coord = y) %>%
  crossing(sample_index = seq_len(n_reps_per_site), x = seq_len(n_years) - 1) %>%
  st_as_sf(coords = c("x_coord", "y_coord")) %>%
  st_intersection(strat_sf) %>%
  mutate(x_coord = st_coordinates(.)[, 1], y_coord = st_coordinates(.)[, 2]) %>%
  st_drop_geometry()

# ---- the sample ----
design <- get_design(n_sites_per_year, n_years, n_sites_total, n_strata)

grts_sample <- tryCatch(
  {
    # The old grts
    grts(
      design = design,
      DesignID = "site",
      type.frame = "area",
      src.frame = "sf.object",
      sf.object = strat_sf,
      stratum = "stratum",
      shapefile = FALSE
    )
  },
  error = function(e) {
    # TODO: the new grts
    stop("Rework grts call using newer version of the package...")
  }
)
grts_sample_sf_raw <- st_as_sf(grts_sample)

# plot(W, asp = 1)
# plot(strat_sfc, add = TRUE)
# plot(grts_sample_sf_raw[, "stratum"], add = TRUE)

# ra_new <- W
# ra_new[raster::cellFromXY(W, st_coordinates(grts_sample_sf_raw))] <- NA
# plot(ra_new, colNA="gray")
# plot(grts_sample_sf_raw[, "stratum"], add = TRUE)

grts_sample_sf <- grts_sample_sf_raw %>%
  left_join(attr(design, "panel_info"), by = "panel") %>%
  mutate(
    site_index = raster::cellFromXY(W, st_coordinates(.)),
    is_in_sample = TRUE
  )
grts_sample_tbl <- grts_sample_sf %>%
  select(x, panel, site_index, stratum, is_in_sample) %>%
  st_drop_geometry() %>%
  as_tibble()
# grts_sample_sf %>% distinct(site_index, stratum) %>% count(stratum)

scaffold_tbl <- attr(design, "panel_info") %>%
  left_join(scaffold_tbl_raw, by = "x") %>%
  left_join(grts_sample_tbl, by = c("x", "panel", "site_index", "stratum")) %>%
  mutate(stratum_index = as.integer(as.factor(stratum))) %>%
  group_by(stratum_index) %>%
  mutate(site_in_stratum_index = as.integer(as.factor(site_index))) %>%
  ungroup()

p_scaffold <- ggplot(scaffold_tbl %>% filter(!is_in_sample %in% TRUE)) +
  facet_wrap(~x, nrow = 1, labeller = label_both) +
  geom_tile(aes(x = x_coord, y = y_coord, fill = w1)) +
  scale_fill_viridis_c() +
  coord_equal() +
  geom_sf(data = grts_sample_sf, aes(color = stratum)) +
  geom_sf(
    data = crossing(x = seq_len(n_years) - 1, stratum = LETTERS[1:n_strata]) %>%
      left_join(strat_sf, by = "stratum") %>%
      st_as_sf(),
    aes(color = stratum), fill = NA
  ) +
  scale_color_brewer("Stratum", type = "qual", palette = "Set2") +
  labs(x = "X", y = "Y") +
  sim_theme() +
  theme(
    axis.text.x = element_blank(), axis.text.y = element_blank(),
    axis.title.x = element_blank(), axis.title.y = element_blank()
  ) +
  guides(fill = guide_colourbar(barwidth = 10))
save_figure(
  plot = p_scaffold, path = dirname(file),
  filename = "p-scaffold.jpg", p_width = 1.5, p_height = 1.5
)

# list.files()
params_needed <- system(
  paste("model-api/src/get-unknowns.sh", file),
  intern = TRUE
)
params_needed <- grep("^y$|^y.rep$|^B$", params_needed, value = TRUE, invert = TRUE)
is_binomial <- grepl("(binomial|beta-binomial)-inverse-logit", file)
other_quants <- if (is_binomial) {
  "y.n"
} else {
  NULL
}
uses_hier_site_var <- grepl("hier-site", file)
if (uses_hier_site_var) {
  other_quants <- c(other_quants, "y.q")
} else {
  other_quants
}
var_type <- str_extract(file, "fixed|hier(?=-(stratum|site))")
var_type <- ifelse(is.na(var_type), "fixed", var_type)
var_level <- str_extract(file, "(?<=(fixed|hier)-)stratum|site")
var_level <- ifelse(is.na(var_level), "stratum", var_level)

attempt <- 1
jags_samples <- NULL
while (attempt <= n_attempts & is.null(jags_samples)) {
  message(glue("Attempt {attempt} with seed {seed}..."))

  generating_values <- get_generating_values(
    n_strata = n_strata,
    link = get_link(file),
    seed = seed,
    N = nrow(scaffold_tbl),
    which_params = c(params_needed, other_quants),
    params_override = gen_vals,
    var_level = var_level, var_type = var_type
  )
  print(dirname(file))
  print(generating_values[!names(generating_values) %in% "y.n"])
  d_gen_vals <- enframe(generating_values, name = "parameter") %>%
    unnest_longer(value, indices_to = "stratum_index") %>%
    filter(parameter != "y.n")
  save_table(d_gen_vals, dirname(file), "generating-values.csv")

  d_jags <- c(
    # data
    list(
      y.site = scaffold_tbl$site_in_stratum_index,
      y.strata = scaffold_tbl$stratum_index, # TODO: need site in stratum index
      y.n.strata = n_distinct(scaffold_tbl$stratum_index),
      y.n.site = scaffold_tbl %>% distinct(site_index, stratum) %>% count(stratum) %>% pull(n),
      x = scale(scaffold_tbl$x)[, 1], X = scale(scaffold_tbl$w1), # was: matrix(scaffold_tbl$w1, ncol = 1)
      y = rep(NA, nrow(scaffold_tbl))
    ),
    # parameter values (passed as data)
    generating_values
  )

  jags_samples <- tryCatch(
    {
      run.jags(file,
        data = d_jags,
        monitor = c("y"), sample = 1, n.chains = 1, # c("y", "B")
        inits = list(
          .RNG.name = "base::Super-Duper", .RNG.seed = seed
        ),
        summarise = FALSE
      )
    },
    error = function(e) {
      NULL
    }
  )

  # update criteria and seed
  attempt <- attempt + 1
  seed <- sample(1E6, 1)
}

# print(str(coda::as.mcmc(jags_samples)))
coda_samples <- coda::as.mcmc(jags_samples) # reformat the outputs
samples <- as.vector(coda_samples) # hist(samples)
# range(samples); hist(samples)

d <- scaffold_tbl %>% mutate(y = samples)
print(d)

d_in_sample <- d %>%
  filter(is_in_sample %in% TRUE) %>%
  # distinct(x, panel, site_index, stratum) %>%
  group_by(stratum) %>%
  arrange(panel, site_index) %>%
  mutate(seq_site_index = as.integer(factor(site_index, levels = unique(site_index)))) %>%
  ungroup()
p_in_sample <- ggplot(d_in_sample) +
  facet_wrap(~stratum, labeller = label_both) +
  geom_jitter(
    aes(x = x, y = seq_site_index, fill = factor(site_index), size = y),
    width = 0.2, height = 0, pch = 21, alpha = 0.75
  ) +
  labs(x = "Year since sampling began", y = "Site") +
  scale_fill_viridis_d("Site index") +
  sim_theme()
save_figure(
  plot = p_in_sample, path = dirname(file),
  filename = "p-in-sample.jpg", p_width = 5, p_height = 5
)

d_in_sample_summary <- d_in_sample %>%
  group_by(x, site_index, stratum) %>%
  summarise(yhat = mean(y), .groups = "drop")
p_in_sample_summary <- ggplot(d_in_sample_summary) +
  facet_wrap(~stratum, labeller = label_both) +
  geom_jitter(
    data = d_in_sample,
    aes(x = x, y = y, color = factor(site_index))
  ) +
  geom_line(aes(x = x, y = yhat, color = factor(site_index)), alpha = 0.75) +
  labs(x = "Year since sampling began", y = "y") +
  sim_theme() +
  guides(color = guide_legend(title = "Site index"))
save_figure(
  plot = p_in_sample_summary, path = dirname(file),
  filename = "p-in-sample-summary.jpg", p_width = 5, p_height = 5
)

d_summary <- d %>%
  group_by(x, site_index, stratum, x_coord, y_coord) %>%
  summarise(yhat = mean(y), .groups = "drop")
# ggplot(d_summary) +
#   facet_wrap(~x) +
#   geom_tile(aes(x = x_coord, y = y_coord, fill = yhat)) +
#   scale_fill_viridis_c(option = "B") +
#   coord_equal()
#   # geom_sf(data = strat_sf %>% complete(x = seq_len(n_years) - 1),
#   #         aes(color = stratum), fill = NA)

p_summary <- ggplot(d_summary) +
  facet_wrap(~x, nrow = 1, labeller = label_both) +
  geom_tile(aes(x = x_coord, y = y_coord, fill = yhat)) +
  scale_fill_viridis_c(option = "B") +
  coord_equal() +
  geom_sf(
    data = crossing(x = seq_len(n_years) - 1, stratum = LETTERS[1:n_strata]) %>%
      left_join(strat_sf, by = "stratum") %>%
      st_as_sf(),
    aes(color = stratum), fill = NA
  ) +
  scale_color_brewer("Stratum", type = "qual", palette = "Set2") +
  labs(x = "X", y = "Y") +
  sim_theme() +
  theme(
    axis.text.x = element_blank(), axis.text.y = element_blank(),
    axis.title.x = element_blank(), axis.title.y = element_blank()
  ) +
  guides(fill = guide_colourbar(barwidth = 10))
save_figure(
  plot = p_summary, path = dirname(file),
  filename = "p-summary.jpg", p_width = 1.5, p_height = 1.5
)



# why 256? (20*4 + 12*4 )*2
p_ridges <- ggplot(d %>% mutate(grp = is_in_sample %in% TRUE)) +
  facet_wrap(~stratum, labeller = label_both) +
  geom_density_ridges(
    aes(x = y, y = x, group = interaction(x, grp), fill = grp),
    stat = "binline", alpha = 0.75, bins = 50
  ) +
  scale_y_reverse(breaks = function(x) int_breaks(x, n = n_years)) +
  labs(x = "Year since sampling began", y = "") +
  sim_theme() +
  scale_fill_brewer("In sample?", type = "qual", palette = "Set1")
save_figure(
  plot = p_ridges, path = dirname(file),
  filename = "p-ridges.jpg", p_width = 5, p_height = 5
)


if ("y.n" %in% other_quants) {
  trials <- generating_values$y.n[1]
  print(d_in_sample)
  d_in_sample <- d_in_sample %>%
    mutate(trials)
  print(d_in_sample)
}
sample_col <- ifelse(is_binomial, "transect", "plot")
d_demo_resp <- d_in_sample %>%
  mutate(
    park_unit = !!park_unit, event_year = start_year + x,
    !!sample_col := sample_index
  ) %>%
  dplyr::select(-x, -panel, -matches("*_index$|*_coord$|^w\\d+$"), -is_in_sample,
    site = site_index, y_sim = y, any_of("trials")
  )
d_demo_resp
save_table(d_demo_resp, dirname(file), "samples.csv")

d_demo_xy <- d %>%
  mutate(park_unit = !!park_unit, event_year = start_year + x) %>%
  dplyr::select(park_unit, site = site_index, stratum, matches("*_coord$")) %>%
  distinct() %>%
  arrange(site)
d_demo_xy
save_table(d_demo_xy, dirname(file), "sample-locs.csv")

d_demo_W <- d %>%
  mutate(park_unit = !!park_unit, event_year = start_year + x) %>%
  dplyr::select(park_unit, site = site_index, stratum, matches("^w\\d+$")) %>% # event_year
  distinct() %>%
  arrange(site)
d_demo_W
# d_demo_W %>% filter(duplicated(w1))
save_table(d_demo_W, dirname(file), "sample-covariates.csv")

message("success")
# p_hist <- stk_hist(
#   x = d_in_sample, var = 'y',
#   facets = 'stratum',
#   theme = get_theme("sim_theme", y_lab_angle = 90)(),
#   bins = 30
# )
# save_figure(plot = p_hist, path = dirname(file), filename = 'p-hist.jpg',
#             p_width = 3, p_height = 3)
#
# # data
# jags.data_raw <- modifyList(d_jags, list(y = samples))
# # params <- c('mu.B0', 'mu.B1', 'Beta', 'sigma.B0', 'sigma.B1', 'rho')
# jags.data <- jags.data_raw[-which(names(d_jags) %in% params_needed)]
#   #list(y = samples, N = length(samples), x = x)
# str(jags.data)
#
#
# # MCMC settings
# ni <- 15000
# nt <- 1
# nb <- 5000
# nc <- 3
#
# inits <- function(true_vals) {
#   init <- mapply(function(param_val, param_name) {
#     n <- length(param_val)
#     param_init <- if (grepl('^mu.B\\d|^Beta|^B\\d', param_name)) {
#       rnorm(n)
#     } else if (grepl('^sigma', param_name)) {
#       runif(n, 0, 0.1)
#     } else if (grepl('^rho', param_name)) {
#       runif(n, -0.5, 0.5)
#     } else {
#       browser()
#     }
#     list(param_init)
#   }, true_vals, names(true_vals))
#   function() init
# }
# inits(d_jags[names(d_jags) %in% params_needed])()
# #
# # inits <- function(){
# #   # list(
# #   #   alpha = rnorm(1), beta = rnorm(1), sigma = runif(1,0,10)
# #   # )
# # }
# #
# # call JAGS from R
# just_params <- grep('y.n', params_needed, value = TRUE, invert = TRUE)
# res <- jags(
#   data = jags.data,
#   inits = inits(d_jags[names(d_jags) %in% just_params]), #d_jags[names(d_jags) %in% params_needed]
#   parameters.to.save = just_params, #names(d_jags)[names(d_jags) %in% params_needed],
#   model.file = file,
#   n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb
# )
# print(res, digits = 3)
# ya <- coda::as.mcmc(res)
# #
# d_gen_vals_new <- d_gen_vals %>%
#   group_by(parameter) %>%
#   mutate(
#     parameter = sprintf('%s%s',
#                         parameter,
#                         ifelse(rep(n() > 1, n()), paste0('[', stratum_index, ']'), ''))
#   ) %>%
#   ungroup()
# color_scheme_set(scheme = "brewer-Spectral")
# p_trace_raw <- mcmc_trace(ya,  regex_pars = unique(d_gen_vals$parameter)) # mcmc_combo
# p_trace <- p_trace_raw +
#   geom_hline(
#     data = d_gen_vals_new,
#     aes(yintercept = value), linetype = 'dashed', alpha = 0.8
#   )
# save_figure(plot = p_trace, path = dirname(file), filename = 'p-trace.jpg',
#             p_width = 3, p_height = 1.5)
#
# # Assign exit status for automated checks.
# d_verify <- d_gen_vals_new %>%
#   left_join(as_tibble(res$BUGSoutput$summary, rownames = "parameter"),
#             by = "parameter") %>%
#   mutate(value_is_within_post = value >= `2.5%` & value <= `97.5%`)
# exit_code <- get_exit_code(
#   is_converged = all(res$BUGSoutput$summary[, 'Rhat'] < 1.1),
#   is_verifiable = all(d_verify$value_is_within_post)
#   )
# q(save = "no", status = exit_code)
