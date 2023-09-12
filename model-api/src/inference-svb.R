# get_site_inference(., d, output_dir, a_label)

get_svbs <- function(mcarray_list, data, folder, resp_var_desc,
                     interval = 0.95, method = c("hdi", "eti")) {
  folder <- ex_path
  mcarray_list <- readRDS(file.path(ex_path, "99-misc", "z-jags.rds"))
  data <- read_csv(file.path(ex_path, "00-input", "state-variable-data.csv"))
  jags_data <- readRDS(file.path(ex_path, "00-input", "jags-data.rds"))

  svb_info_path <- file.path(ex_path, "00-input", "stratum-varying-effect-info.rds")
  if (!file.exists(svb_info_path)) {
    return(NULL)
  }

  svb_outputs_path <- file.path(folder, "03-inference", "svb")
  svb_info <- readRDS(svb_info_path)
  s_lookup <- read_csv(file.path(folder, "00-input", "stratum-ids-and-indices.csv"))

  Beta <- t(concat_chains(mcarray_list$Beta, 3))
  # kappa <- matrix(NA, nrow = svb_info$n.svb, ncol = max(s_lookup$stratum_index))

  Beta_labs <- dimnames(jags_data$X)[[2]][svb_info$m.k]
  Kappa_raw <- list()
  for (i in 1:svb_info$n.svb) {

    # n.svb: the number of stratum-varying effects
    # m.k: the index for the ith stratum-varying effect

    kappa <- tibble(!!s_lookup$stratum_id[1] := Beta[, svb_info$m.k[i]])

    for (j in 2:jags_data$y.n.strata) {
      kappa <- kappa %>%
        mutate(!!s_lookup$stratum_id[j] :=
          Beta[, svb_info$m.k[i]] + Beta[, svb_info$m.k[i] + j - 1])
    }
    Kappa_raw[[i]] <- kappa %>% mutate(variable = Beta_labs[i])
  }
  Kappa <- bind_rows(Kappa_raw) %>%
    pivot_longer(-variable, names_to = "stratum", values_to = "svb")
  p <- ggplot(Kappa) +
    facet_grid(variable ~ stratum, scales = "free") +
    geom_histogram(aes(x = svb, y = ..density..),
      fill = "white", color = "black"
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "orange")
  save_figure(
    plot = p, path = svb_outputs_path,
    filename = "svb-hists.jpg",
    p_width = 2, p_height = 2, device = NULL
  )

  d_x <- data %>%
    mutate(variable_name = Beta_labs[i]) %>%
    rename(variable_val = !!Beta_labs[i])
  p_x <- ggplot(d_x) +
    facet_grid(variable_name ~ stratum_id, scales = "free") +
    geom_histogram(aes(x = variable_val, y = ..density..),
      fill = "white", color = "black"
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "orange")
  d_x %>%
    group_by(variable_name, stratum_id) %>%
    summarise(
      mean = mean(variable_val), sd = sd(variable_val),
      .groups = "drop"
    )

  d_Kappa <- Kappa %>%
    group_by(variable, stratum) %>%
    summarise(
      median = median(svb),
      lower = get_interval(svb, interval, method[1])[1],
      upper = get_interval(svb, interval, method[1])[2],
      .groups = "drop"
    )
  save_table(d_Kappa, svb_outputs_path, "svb-summary.csv")
}
