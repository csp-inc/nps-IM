get_zone_objects <- function(y.zone.index, x.pred.index, n.x.pred,
                             output_dir, X.pred) { # , j.pred.zone


  tot.site.zone <- as_tibble(cbind(y.zone.index, x.pred.index)) %>%
    group_by(y.zone.index, x.pred.index) %>%
    tally() %>%
    ungroup() %>%
    spread(y.zone.index, n) %>%
    arrange(x.pred.index) %>%
    select(-x.pred.index) %>%
    as.matrix()

  X.pred.zone.lookup <-
    array(NA, dim = c(max(tot.site.zone), n.x.pred, max(y.zone.index)))
  # str(X.pred.zone.lookup)

  X_pred <- read_csv(file.path(output_dir, "00-input/complete-covariates-data.csv"))
  j.pred.zone <- X_pred %>%
    group_by(rel_year, MgmtZone) %>%
    mutate(j_pred_zone = index_from_category(SiteName)) %>%
    ungroup() %>%
    pull(j_pred_zone)

  idx_vs_X_df_zone <-
    cbind(j.pred.zone, y.zone.index, x.pred.index, X.pred) %>%
    as_tibble() %>%
    mutate(row.index = 1:n())
  for (l in 1:n_distinct(y.zone.index)) {
    this_idx_vs_X_df_zone <- filter(idx_vs_X_df_zone, y.zone.index == l) %>%
      select(j.pred.zone, x.pred.index, row.index) %>%
      spread(x.pred.index, row.index) %>%
      select(-j.pred.zone) %>%
      as.matrix()
    X.pred.zone.lookup[1:nrow(this_idx_vs_X_df_zone), 1:ncol(this_idx_vs_X_df_zone), l] <-
      this_idx_vs_X_df_zone
  }

  list(
    y.zone.index = y.zone.index,
    tot.site.zone = tot.site.zone,
    X.pred.zone.lookup = X.pred.zone.lookup
  )
}
