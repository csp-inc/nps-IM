get_deflection_coefs <- function(mcarray_list, data, folder, resp_var_desc,
                                 deterministic_model) {
  print("Inference on coefficients for deflections...")

  deflections_labels <- readRDS(file.path(folder, "00-input/deflections-indices.rds"))
  deflections_df <-
    tibble(delta = paste0("delta", "[", names(deflections_labels[[1]]), "]"))

  transform_coef <- deterministic_model %in% c("exponential", "inverse-logit")
  x_label <- sub("_", " ", names(deflections_labels)) %>% str_to_sentence()

  tform_scenarios <- if (transform_coef) {
    "untransformed" # c('transformed', 'untransformed')
  } else {
    "untransformed"
  }

  get_defl_plots <- function(coef_posts, defl_type, tform_scenario, ref_line,
                             stratum_id = "", strat_idx = NULL) {
    if (tform_scenario == "transformed") coef_posts <- exp(coef_posts)

    coef_stats <- summary(coef_posts, quantile,
      probs = c(.025, .5, .975), na.rm = TRUE
    )$stat
    if (!is.null(strat_idx)) {
      coef_stats <- coef_stats[, 1, , strat_idx]
    } else {
      coef_stats <- coef_stats[, 1, ]
    }
    coef_stats_df <- lapply(seq_len(3), function(x) {
      delta <- deflections_df %>%
        mutate(
          deflection = coef_stats[x, ],
          stat = rownames(coef_stats)[x]
        )
    }) %>%
      bind_rows() %>%
      select(delta, stat, deflection) %>%
      spread(stat, deflection) %>%
      mutate(deflection_index = 1:n())

    xlab <- ifelse(defl_type == "mu.D", "Mean (intercept + deflection)", "Deflection")
    p_defl <- ggplot(coef_stats_df) +
      geom_vline(xintercept = ref_line, color = "red") +
      geom_errorbarh(aes(xmin = `2.5%`, xmax = `97.5%`, y = deflection_index)) +
      geom_point(aes(x = `50%`, y = deflection_index)) +
      theme_ima("default") +
      labs(x = xlab, y = x_label) +
      scale_y_continuous(
        breaks = coef_stats_df$deflection_index,
        labels = parse(text = coef_stats_df$delta)
      )

    save_figure(
      plot = p_defl,
      path = file.path(folder, "03-inference/coef"),
      filename = sprintf(
        "deflections%s-%s-%s.png",
        stratum_id, defl_type, tform_scenario
      ),
      p_width = 6, p_height = 4, device = NULL
    )

    save_table(
      coef_stats_df %>%
        mutate(response = resp_var_desc, scale = tform_scenario),
      file.path(folder, "03-inference/coef"),
      sprintf(
        "deflections%s-%s-%s.csv",
        stratum_id, defl_type, tform_scenario
      )
    )
  }

  lapply(grep("^D$|^mu.D$", names(mcarray_list), value = TRUE), function(defl_type) {
    lapply(tform_scenarios, function(x) {
      if (x == "transformed") x_label <- paste(x_label, "(exponentiated)")
      if (defl_type == "mu.D") {
        if (x == "untransformed") {
          uniq_strat_ids <- unique(data$stratum_id)
          for (this_strat_idx in 1:length(uniq_strat_ids)) {
            this_strat_id <- uniq_strat_ids[this_strat_idx]
            strat_lab <- ifelse(this_strat_id != "none", paste0("-", this_strat_id), "")
            ref_line <- summary(mcarray_list[["B0.tilde.resp.scale"]], mean)$stat
            get_defl_plots(mcarray_list[[defl_type]], defl_type, x, ref_line, strat_lab, this_strat_idx)
          }
        }
      } else {
        ref_line <- ifelse(x == "transformed", 1, 0)
        get_defl_plots(mcarray_list[[defl_type]], defl_type, x, ref_line)
      }
    })
  })
}
