get_misc_contrasts <- function(mcarray_list, data, folder, resp_var_desc, jags_obj,
                               crit_thresh = 0.80) {
  deflections_info_file <- file.path(folder, "00-input/deflections-indices.rds")
  # if (!file.exists(deflections_info_file)) return(NULL)
  if (!"mu.D" %in% names(mcarray_list)) {
    return(NULL)
  }

  deflections_info <- readRDS(deflections_info_file)[["point_jurisdiction"]]
  defl_name <- names(deflections_info)

  d_contrasts_mat <- t(concat_chains(mcarray_list[[jags_obj]][1, , 1, , ], axis = 3)) # TODO: hard-coded 1 here in dim 1 is not general, nor is the 1 for stratum in dim 3
  dimnames(d_contrasts_mat) <- list(
    iter = 1:nrow(d_contrasts_mat),
    jurisdiction = defl_name
  )

  d_contrasts_wide <- d_contrasts_mat %>%
    as_tibble()

  contrasts <- combn(defl_name, 2)

  contrast_results <- lapply(1:ncol(contrasts), function(this_contrast) {
    this_diff <- paste(paste0("`", contrasts[, this_contrast], "`"), collapse = " - ")
    diff_label <- paste0("`", paste(contrasts[, this_contrast], collapse = "_minus_"), "`")
    eval(parse(text = sprintf(
      "select(mutate(d_contrasts_wide, %s = %s), %s)",
      diff_label, this_diff, diff_label
    )))
  }) %>%
    reduce(bind_cols) %>%
    gather(contrast, diff)

  p_contr <- ggplot(
    contrast_results %>%
      mutate(contrast = sub("_minus_", " - ", contrast)), #
    aes(x = diff, y = contrast)
  ) +
    geom_density_ridges(fill = "steelblue", alpha = 0.8, stat = "binline", scale = 0.9, bins = 50) +
    geom_vline(xintercept = 0, color = "red") +
    labs(x = sprintf("Difference in mean %s", tolower(resp_var_desc)), y = "Contrast") +
    theme_ima("default")
  save_figure(
    plot = p_contr,
    path = file.path(folder, "03-inference/strata"),
    filename = paste0(jags_obj, "-contrasts.png"),
    p_width = 9, p_height = 6, device = NULL
  )

  outcome <- data %>%
    select(outcome = matches("^response$|^hits$")) %>%
    pull(outcome)
  if (any(grepl("^trials$", names(data)))) outcome <- outcome / data$trials
  cdf_at_raw <- seq(min(outcome), max(outcome), length.out = 11)
  cdf_at <- sort(c(cdf_at_raw * -1, cdf_at_raw * 1))
  cdf_at_lab <- c(
    sprintf("P(z<%s)", sort(cdf_at_raw * -1)),
    sprintf("P(z>%s)", sort(cdf_at_raw * 1))
  )

  split_contrast <- function(x, which) {
    sapply(x, function(z) strsplit(z, "_minus_")[[1]][which])
  }

  # tmp = contrast_results %>% filter(contrast == "USFSNONWILDERNESS_minus_USFSWILDERNESS") %>% pull(diff)
  # hist(tmp)
  d_ecdf <- contrast_results %>%
    # filter(contrast == "USFSNONWILDERNESS_minus_USFSWILDERNESS") %>%
    group_by(contrast) %>%
    summarise(
      !!cdf_at_lab[1] := sum(diff < cdf_at[1]) / n(), # ecdf(diff)(cdf_at[1]),
      !!cdf_at_lab[2] := sum(diff < cdf_at[2]) / n(),
      !!cdf_at_lab[3] := sum(diff < cdf_at[3]) / n(),
      !!cdf_at_lab[4] := sum(diff < cdf_at[4]) / n(),
      !!cdf_at_lab[5] := sum(diff < cdf_at[5]) / n(),
      !!cdf_at_lab[6] := sum(diff < cdf_at[6]) / n(),
      !!cdf_at_lab[7] := sum(diff < cdf_at[7]) / n(),
      !!cdf_at_lab[8] := sum(diff < cdf_at[8]) / n(),
      !!cdf_at_lab[9] := sum(diff < cdf_at[9]) / n(),
      !!cdf_at_lab[10] := sum(diff < cdf_at[10]) / n(),
      !!cdf_at_lab[11] := sum(diff < cdf_at[11]) / n(),
      !!cdf_at_lab[12] := sum(diff > cdf_at[12]) / n(), # 1 - ecdf(diff)(cdf_at[6]),
      !!cdf_at_lab[13] := sum(diff > cdf_at[13]) / n(),
      !!cdf_at_lab[14] := sum(diff > cdf_at[14]) / n(),
      !!cdf_at_lab[15] := sum(diff > cdf_at[15]) / n(),
      !!cdf_at_lab[16] := sum(diff > cdf_at[16]) / n(),
      !!cdf_at_lab[17] := sum(diff > cdf_at[17]) / n(),
      !!cdf_at_lab[18] := sum(diff > cdf_at[18]) / n(),
      !!cdf_at_lab[19] := sum(diff > cdf_at[19]) / n(),
      !!cdf_at_lab[20] := sum(diff > cdf_at[20]) / n(),
      !!cdf_at_lab[21] := sum(diff > cdf_at[21]) / n(),
      !!cdf_at_lab[22] := sum(diff > cdf_at[22]) / n()
    ) %>%
    pivot_longer(-contrast) %>%
    left_join(tibble(thresh = cdf_at, name = cdf_at_lab))
  save_table(
    d_ecdf,
    file.path(folder, "03-inference/strata"),
    paste0(jags_obj, "-contrast-cdf.csv")
  )


  d_ecdf_select <- d_ecdf %>%
    mutate(
      meets_crit_thresh = value >= crit_thresh, # crit_thresh,
      interp = ifelse(thresh < 0 & meets_crit_thresh,
        # sprintf("%s is substantially higher than %s",
        #         split_contrast(contrast, 2), split_contrast(contrast, 1)),
        sprintf(
          "There is a %s%% chance that %s is higher on %s than %s by at least %s%% of its total range", #
          round(value * 100, 2), tolower(resp_var_desc),
          split_contrast(contrast, 2), split_contrast(contrast, 1),
          round(abs(thresh) / max(thresh) * 100, 2)
        ),
        NA
      ),
      interp = ifelse(thresh > 0 & meets_crit_thresh,
        # sprintf("%s is substantially higher than %s",
        #         split_contrast(contrast, 1), split_contrast(contrast, 2)),
        sprintf(
          "There is a %s%% chance that %s is higher on %s than %s by at least %s%% of its total range", #
          round(value * 100, 2), tolower(resp_var_desc),
          split_contrast(contrast, 1), split_contrast(contrast, 2),
          round(abs(thresh) / max(thresh) * 100, 2)
        ),
        interp
      ),
      interp = ifelse(name == "P(z<0)" & meets_crit_thresh,
        sprintf(
          "%s is at least trivially higher than %s (i.e., different from zero)",
          split_contrast(contrast, 2), split_contrast(contrast, 1)
        ),
        interp
      ),
      interp = ifelse(name == "P(z>0)" & meets_crit_thresh,
        sprintf(
          "%s is at least trivially higher than %s (i.e., different from zero)",
          split_contrast(contrast, 1), split_contrast(contrast, 2)
        ),
        interp
      )
    ) %>%
    filter(meets_crit_thresh)
  save_table(
    d_ecdf_select,
    file.path(folder, "03-inference/strata"),
    paste0(jags_obj, "-contrast-cdf-select.csv")
  )


  contrast_summary <- contrast_results %>%
    group_by(contrast) %>%
    summarise(
      median = median(diff),
      lt0 = sum(diff < 0) / n(),
      eq0 = sum(diff == 0) / n(),
      gt0 = sum(diff > 0) / n()
    )
  save_table(
    contrast_summary,
    file.path(folder, "03-inference/strata"),
    paste0(jags_obj, "-contrast-summary.csv")
  )


  d_ecdf_lab <- d_ecdf %>%
    mutate(
      lt0_lab = sprintf("%s is higher", split_contrast(contrast, 2)),
      gt0_lab = sprintf("%s is higher", split_contrast(contrast, 1))
    ) %>%
    mutate(
      contrast = sub("_minus_", " - ", contrast),
      # lab_pos = 0.5, label_col = "P(z=0)",
      lab_pos = ifelse(grepl(">", name), Inf, -Inf),
      label_col = ifelse(grepl(">", name), "P(z > thresh)", "P(z < thresh)")
    )
  line_cols <- RColorBrewer::brewer.pal(3, "Set1")
  p_contr_cdf <- ggplot(
    contrast_results %>%
      mutate(contrast = sub("_minus_", " - ", contrast)), #
    aes(x = diff)
  ) +
    facet_wrap(~ forcats::fct_rev(contrast), ncol = 1) +
    geom_vline(xintercept = cdf_at, linetype = "dashed", color = "orange") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    stat_ecdf(geom = "step", color = line_cols[1]) + # TODO: join and add labels, and save labels df as csv
    geom_step(aes(y = 1 - ..y..), stat = "ecdf", color = line_cols[2]) +
    geom_label(
      data = d_ecdf_lab, aes(
        x = thresh, y = lab_pos,
        label = format(round(value, 3), nsmall = 3),
        color = label_col, fill = value
      ),
      vjust = "inward", alpha = 0.75
    ) +
    scale_fill_gradient("Probability",
      low = "white", high = line_cols[3],
      breaks = seq(0, 1, 0.2), limits = c(0, 1)
    ) +
    theme_ima("default") +
    scale_x_continuous(breaks = cdf_at) +
    labs(x = "Difference (z)", y = "CDF") +
    scale_color_brewer("", type = "qual", palette = "Set1")
  save_figure(
    plot = p_contr_cdf,
    path = file.path(folder, "03-inference/strata"),
    filename = paste0(jags_obj, "-contrasts-cdf.png"),
    p_width = 15, p_height = 8, device = NULL
  )
}
