get_residuals <- function(mcarray_list, data, likelihood, folder) {
  if (likelihood == "hurdle-ordinal-latent-beta") data <- data %>% filter(response > 0)

  # Pulls residuals (used to compute semi-variograms, etc.).
  variography_switch <- any(str_detect(names(data), "sample_x$|sample_y$"))
  epsilon_name <- grep("^epsilon$|^epsilon.beta$", names(mcarray_list),
    value = TRUE
  )
  if (variography_switch) {
    residuals_df <- data %>%
      select(sample_x, sample_y) %>%
      mutate(residual = summary(mcarray_list[[epsilon_name]], mean)$stat) %T>%
      {
        `if`(
          any(is.infinite(.$residual)),
          write_file(
            "Some residuals not finite!",
            file.path(folder, "residuals-warnings.txt")
          ), .
        )
      } %T>%
      {
        `if`(
          any(is.na(.$sample_x)),
          write_file(
            "Some samples are missing location info!",
            file.path(folder, "missing-locs-warnings.txt")
          ), .
        )
      } %T>%
      save_object(file.path(folder, "02-checking/variography"), "residuals-df.rds") %>%
      filter(!is.infinite(residual), !is.na(sample_x))


    sp::coordinates(residuals_df) <- c("sample_x", "sample_y")

    directions_df <- tibble(
      angle = c(0, 45, 90, 135, 999)
    ) %>%
      mutate(angle_lab = factor(angle,
        labels =
          c(
            "Northern", "North-eastern",
            "Eastern", "South-eastern",
            "All"
          )
      ))
    isotropic_variog <- gstat::variogram(residual ~ 1, data = residuals_df) %>%
      mutate(variogram = "Isotropic", dir.hor = 999)
    anisotropic_variog <-
      gstat::variogram(residual ~ 1,
        data = residuals_df,
        alpha = directions_df %>% filter(angle < 360) %>% pull(angle)
      ) %>%
      mutate(variogram = "Anisotropic")

    variography_df <- bind_rows(isotropic_variog, anisotropic_variog) %>%
      left_join(directions_df, by = c("dir.hor" = "angle"))

    p_vario <- ggplot(variography_df, aes(x = dist, y = gamma)) +
      facet_grid(variogram ~ .) +
      geom_line(size = .5, aes(color = angle_lab), alpha = .5) +
      geom_point(size = 3, aes(color = angle_lab)) +
      labs(x = "Distance (m)", y = expression(gamma)) +
      theme_ima("default", rotate_y_lab = TRUE) +
      ggthemes::scale_colour_hc("default") +
      guides(color = guide_legend("Direction"))
    save_figure(
      plot = p_vario,
      path = file.path(folder, "02-checking/variography"),
      filename = "variograms.png",
      p_width = 6, p_height = 2, device = NULL
    )
  } else {
    resids_vec <- summary(mcarray_list[[epsilon_name]], mean)$stat
    save_object(resids_vec, file.path(folder, "02-checking"), "residuals-vec.rds")
  }
}
