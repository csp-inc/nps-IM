.sites <- function(data, in_strata = FALSE) {
  if (in_strata) {
    data %>%
      distinct(
        site_id, stratum_id, site_index, stratum_index,
        site_in_stratum_index
      ) %>%
      arrange(stratum_index, site_in_stratum_index) %>%
      adf(Inf)
  } else {
    data %>%
      distinct(site_id, site_index) %>%
      arrange(site_index) %>%
      adf(Inf)
  }
}

.time <- function(data) {
  data %>%
    distinct(cal_year, rel_year) %>%
    arrange(cal_year) %>%
    adf(Inf)
}

.strata <- function(data) {
  data %>%
    distinct(stratum_id, stratum_index) %>%
    arrange(stratum_index) %>%
    adf(Inf)
}

horiz_line <- function() paste0(paste(rep("-", 78), collapse = ""), "\n")

bt <- function(x) glue::backtick(x)

message_for_jags_prep <- function() {
  # Uses the environment of the parent frame (namely, `prep_data_for_jags`) to
  # create a message containing relevant variable mappings, including R code
  # that can be used at the console to print related (and hopefully useful)
  # summary information.
  #
  # Returns:
  #   A message containing relevant variable mappings. Useful for logging.

  # Wrappers used to make the message 'pretty'.

  hard_remap_msg <- "`->>` indicate hard remappings (e.g., renamed variables)"
  soft_remap_msg <- "`~>` indicate\n# soft remappings (e.g., mutated variables)"
  datetime_col <- event_date_info$`date-time column`
  writeLines(paste0(
    "# Relevant variable mappings for this analysis of data from ", unit_code,
    ". Note that the\n",
    "# ", hard_remap_msg, " while ", soft_remap_msg, ".\n",
    "# ", horiz_line(),
    "# ", bt(response_col), " ->> ", bt("response"), "\n",
    "# ", bt(site_col), " ->> `site_id` ~> ", bt("site_index"),
    " (used to define group-level effects)\n",
    "# ", paste(bt(sample_cols), collapse = ", "), " ->> `sample_id` ~> ",
    bt("sample_index"), "\n",
    "# ", bt(datetime_col), " ->> ", bt("cal_year"), " ~> ", bt("rel_year"), "\n",
    paste0(
      "# ", bt(stratum_col), " ->> `stratum_id` ~> ", bt("stratum_index"),
      "\n",
      "# ", paste(bt(c(site_col, stratum_col)), collapse = ", "), " ~> ",
      bt("site_in_stratum_index"), " (an index for site in stratum)\n"
    ),
    "# ", paste(bt(coord_cols), collapse = ", "),
    " ->> site and sample coordinates `*_x`, `*_y`", "\n\n",
    "# covariates: ", comma_concat(explanatory_vars)
  ))
}

logger <- function(data, log_dir) {

  # Write the data (and relevant indices) to disk for future reference.

  log_dir <- gsub(">", "gt", log_dir)
  dir.create(file.path(log_dir, "00-input"), showWarnings = FALSE, recursive = TRUE)
  write_csv(data, file.path(log_dir, "00-input", "state-variable-data.csv"))
  write_csv(
    data %>% .sites(),
    file.path(log_dir, "00-input", "site-ids-and-indices.csv")
  )
  write_csv(
    data %>% .time(),
    file.path(log_dir, "00-input", "calendar-and-relative-years.csv")
  )
  write_csv(
    data %>% .strata(),
    file.path(log_dir, "00-input", "stratum-ids-and-indices.csv")
  )
  write_csv(
    data %>% .sites(in_strata = TRUE),
    file.path(log_dir, "00-input", "site-in-stratum-info.csv")
  )
  writeLines(paste("\n#", "Success!", "See", bt(log_dir), "for specifics..."))
}
