library(hrbrthemes)
source("model-api/src/utils-mcmc.R")

plot_type_summary <- function(data, interval = 0.95,
                              method = c("hdi", "eti"),
                              by_stratum = TRUE) {
  data %>%
    mutate(stratum_id = ifelse(rep(by_stratum, n()), stratum_id, "All park units")) %>%
    group_by(stratum_id, plot_type, set, timestep, grp_col) %>%
    summarise(
      median = median(outcome, na.rm = TRUE),
      lower = get_interval(outcome, interval, method[1], rm_na = TRUE)[1],
      upper = get_interval(outcome, interval, method[1], rm_na = TRUE)[2],
      .groups = "drop"
    )
}

create_cdpr_summaries <- function(mcarray_list, data, data_list, data_raw,
                                  folder, resp_var_desc, ...) {
  if (!"pred.site.mean" %in% names(mcarray_list)) {
    return(NULL)
  }

  # folder <- 'output/CDPR/TSRA/cwd/cwd-counts/neg-binom_lin_b0-b1_hier-site-var/z1_z2_z0*rel_year_s_s*northness_dominant_genus_staggering'
  # mcarray_list <- readRDS(file.path(folder, '99-misc', 'z-jags.rds'))
  # data <- read_csv(file.path(folder, '00-input', 'state-variable-data.csv'))
  # data_list <- readRDS(file.path(folder, '00-input', 'jags-data.rds'))
  # data_raw <- read_csv('data/CASP/modified/coarse_woody_debris_0to25.csv')

  d_raw_distinct <- if (n_distinct(data$stratum_index) == 1) {
    data_raw %>%
      select(site_id = plot_name, stratum_id = unit_code, plot_type) %>%
      distinct()
  } else {
    data_raw %>% distinct(site_id = plot_name, stratum_id = park_name, plot_type)
  }

  d_select <- data %>%
    select(matches("site|stratum")) %>% # ^cal_year$|^rel_year$
    distinct() %>%
    left_join(d_raw_distinct, by = c("site_id", "stratum_id"))

  pred_site_means_array <- concat_chains(
    mcarray_list[["pred.site.mean"]][, , 1, , ],
    axis = 4
  )
  n_iter <- dim(pred_site_means_array)[3]
  # pred_site_means_array <- # potentially temporary, used to speed things up
  #   pred_site_means_array[, , , sample(n_iter, min(1000, n_iter))]
  pred_site_diffs_array <- apply(pred_site_means_array, 3, function(x) {
    sweep(x, c(1, 2), x[, 1], "-")
    # apply(x, 2, function(z) z - x[, 1])
  })
  dim(pred_site_diffs_array) <-
    c(dim(pred_site_means_array)[1:2], dim(pred_site_diffs_array)[-1])

  arr2tbl <- function(arr, resp_name, summarize = FALSE) {
    dimnames(arr) <- list(
      y.site = 1:dim(arr)[1],
      timestep = 1:dim(arr)[2], # - 1,
      # y.strata = 1:dim(arr)[3],
      iter = 1:dim(arr)[3]
    )

    as.data.frame.table(arr, responseName = resp_name) %>%
      as_tibble() %>%
      mutate(y.strata = 1, timestep = as.integer(timestep) - 1) %>%
      mutate_all(as.double) %>%
      left_join(d_select, by = c(
        "y.site" = "site_in_stratum_index",
        "y.strata" = "stratum_index"
      )) %>%
      filter(!is.na(site_id))
  }

  pred_site_means_df <- arr2tbl(pred_site_means_array, "outcome")
  pred_site_diffs_df <- arr2tbl(pred_site_diffs_array, "outcome")


  sets <- list( # everything = c('Control', 'Burn', 'Burn & Control', 'Burn - Control'),
    control_only = "Control", burn_only = "Burn",
    control_vs_burn = c("Burn & Control", "Burn - Control")
  )
  scenarios <- if (n_distinct(data$stratum_index) == 1) FALSE else c(TRUE, FALSE)
  lapply(scenarios, function(this_grouping_scenario) {
    message(sprintf("By stratum: %s", this_grouping_scenario))
    lapply(c("mean", "difference"), function(which_summary) {
      message(sprintf(".. summary by: %s", which_summary))
      this_df_raw <-
        if (which_summary == "mean") pred_site_means_df else pred_site_diffs_df

      this_df_grp <- this_df_raw %>%
        mutate(stratum_id = ifelse(rep(this_grouping_scenario, n()),
          stratum_id, "All park units"
        ))

      this_df <- bind_rows(
        this_df_grp %>% mutate(set = plot_type),
        this_df_grp %>% mutate(set = "B & C"),
        # This next block computes the difference between burn and controls, for each
        # timestep in each stratum at each iteration....
        this_df_grp %>%
          group_by(iter, stratum_id, plot_type) %>%
          filter(site_id %in% sample(unique(site_id), 1)) %>%
          mutate(y.site = Inf, site_index = NA, site_id = NA) %>%
          ungroup() %>%
          pivot_wider(names_from = plot_type, values_from = outcome) %>%
          mutate(outcome = B - C, set = "B - C", plot_type = set)
      ) %>%
        mutate(
          grp_col = ifelse(plot_type == "C", ipsum_pal()(3)[2], NA),
          grp_col = ifelse(plot_type == "B", ipsum_pal()(3)[1], grp_col),
          grp_col = ifelse(plot_type == "B - C", ipsum_pal()(3)[3], grp_col),
          set = factor(set,
            levels = c("C", "B", "B & C", "B - C"),
            labels = c(
              "Control", "Burn", "Burn & Control",
              "Burn - Control"
            )
          ),
          plot_type = factor(plot_type,
            levels = c("C", "B", "B - C"),
            labels = c(
              "Control", "Burn",
              "Difference\n(Burn - Control)"
            )
          )
        )

      mapply(function(set_name, this_set) {
        message(sprintf(".... set: %s", set_name))
        # print(this_set)
        # if (set_name == "control_vs_burn") browser()
        this_set_df <- this_df %>%
          filter(set %in% this_set)
        # if (this_grouping_scenario == FALSE) browser()
        legend_col_key <- this_set_df %>% distinct(plot_type, grp_col)

        d_plotting <- plot_type_summary(this_set_df,
          by_stratum = this_grouping_scenario,
          ...
        )
        p <- ggplot() +
          facet_grid(set ~ stratum_id, scales = "fixed") +
          geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
          geom_line(
            data = d_plotting,
            aes(
              x = timestep, y = median,
              group = plot_type, color = grp_col
            ), size = 0.8,
            show.legend = set_name == "control_vs_burn"
          ) +
          geom_ribbon(
            data = d_plotting %>% filter(!is.na(median)),
            aes(
              x = timestep, ymin = lower, ymax = upper,
              group = plot_type, color = grp_col, fill = grp_col
            ),
            alpha = 0.8, show.legend = set_name == "control_vs_burn"
          ) +
          theme_ipsum_rc(
            grid = "Y", base_size = 14,
            axis_title_size = 16, strip_text_size = 16
          ) +
          # scale_color_manual("Plot type", values = rev(ipsum_pal()(3))) +
          # scale_fill_manual("Plot type", values = rev(ipsum_pal()(3))) +
          labs(x = "Timestep", y = sprintf("%s (%s)", resp_var_desc, which_summary)) +
          theme(legend.position = "bottom") +
          # guides(fill = set_name == "control_vs_burn",
          #        color = set_name == "control_vs_burn")
          scale_color_identity("", guide = "legend", labels = legend_col_key$plot_type, breaks = legend_col_key$grp_col) +
          scale_fill_identity("", guide = "legend", labels = legend_col_key$plot_type, breaks = legend_col_key$grp_col)
        # p
        p_dir <- file.path(folder, "99-misc", set_name)
        dir.create(p_dir, showWarnings = FALSE, recursive = TRUE)
        p_filename <- sprintf(
          "pred-site-%s-%s.jpg",
          which_summary,
          ifelse(this_grouping_scenario, "by-park", "across-parks")
        )
        save_table(d_plotting, p_dir, sub(".jpg", ".csv", p_filename))
        save_figure(
          plot = p, path = p_dir, filename = p_filename,
          p_width = 3, p_height = 2, device = NULL
        )
      }, names(sets), sets)
    })
  }) # results grouping by park unit (TRUE) or across all parks (FALSE)
}
