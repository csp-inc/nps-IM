update_ima_theme <- function(plot_bgcolor, rotate_y_lab = FALSE) {
  {    if (rotate_y_lab) {
    theme(axis.title.y = element_text(angle = 0, vjust = 0.5))
  } else {
    theme(axis.title.y = element_text(angle = 90, vjust = 0))
  }  } +
    ggthemes::theme_hc(style = plot_bgcolor, base_size = 14) +
    theme(panel.spacing = unit(1.5, "lines"))
}

theme_ima <- function(plot_bgcolor, ...) {
  switch(plot_bgcolor,
    default = update_ima_theme("default", ...),
    darkunica = update_ima_theme("darkunica", ...)
  )
}

ima_pal <- function(plot_bgcolor, n = 7) {
  switch(plot_bgcolor,
    default = ggthemes::hc_pal(plot_bgcolor)(n),
    darkunica = ggthemes::hc_pal(plot_bgcolor)(n)
  )
}

ima_post_col <- function(plot_bgcolor) {
  switch(plot_bgcolor,
    default = c("gray", "black"),
    darkunica = c("gray", "white")
  )
}

if (FALSE) {
  show_pal_df <- tibble(
    pal = rep(c("default", "darkunica"), each = 7),
    col = c(ima_pal("default"), ima_pal("darkunica"))
  ) %>%
    group_by(pal) %>%
    mutate(col_index = 1:n()) %>%
    ungroup()

  show_pal <- function(data, plot_bgcolor) {
    ggplot(data %>% filter(pal == plot_bgcolor)) +
      geom_bar(aes(x = col_index, fill = col)) +
      scale_fill_identity() +
      theme_ima(plot_bgcolor, rotate_y_lab = TRUE) +
      labs(x = "col_index", y = "test")
  }
  gridExtra::grid.arrange(
    show_pal_df %>% show_pal("default"),
    show_pal_df %>% show_pal("darkunica"),
    nrow = 2
  )
}
