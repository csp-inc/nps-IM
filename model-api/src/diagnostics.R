library(glue)

get_param_trace <- function(mcmc_list, folder, file_ext = "jpg", ...) {
  message("Diagnostics: trace plot(s)...")

  # output_dir <- file.path(folder, '01-diagnostics/trace')
  # dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  params <- grep("^p.*|deviance", dimnames(mcmc_list[[1]])[[2]],
    value = TRUE, invert = TRUE
  )

  lapply(params, function(param) {
    which_param <- which(dimnames(mcmc_list[[1]])[[2]] == param)
    post_mat <- do.call(cbind, lapply(mcmc_list, function(x) x[, which_param]))
    dimnames(post_mat) <- list(NULL, seq_len(ncol(post_mat)))
    param_mcmc <- as_tibble(post_mat) %>%
      mutate(draw = 1:n()) %>%
      gather(chain, param_val, -draw)

    p <- ggplot(param_mcmc) +
      geom_line(aes(x = draw, y = param_val, color = chain),
        alpha = 0.8, size = 0.5
      ) +
      ggthemes::theme_hc(base_size = 16) +
      scale_color_brewer("Chain", type = "qual", palette = "Dark2") +
      labs(x = "Iteration", y = param)

    save_figure(
      plot = p,
      path = file.path(folder, "01-diagnostics/trace"),
      filename = paste(param, glue("trace.{file_ext}"), sep = "-"),
      p_width = nrow(post_mat) / 1000 * 8, p_height = 3,
      ...
    ) # device = NULL
  })
}
