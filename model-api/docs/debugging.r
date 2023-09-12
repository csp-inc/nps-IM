# Debugging switch. Change to `TRUE` if filtering `jags_data` for debugging
# purposes and update the call to `filter` in the relevant block below. The
# current arguments to `filter` are for demonstration purposes only.
debugging <- FALSE

# Additional filtering for debugging purposes... Note, this can only be done
# for models without additional covariates (at least for the time being).
if (debugging) {
  # If using this, you'll want to replace `data=jags_data` in the call to
  # `jags.model` with `data=jags_data_filtered`....
  jags_data_filtered <- d %>%
    # Here's where, if desired, you would filter the data to your liking.
    filter(site_in_stratum_index <= 4, response > min(response)) %>%
    # Then you apply our `get_data_list` to recreate the list, this time
    # with the filtered data.
    get_data_list(
      stratum_weights = jags_data$wt,
      L = jags_info$likelihood # likelihood
    )
  jags_inits_filtered <- lapply(jags_inits, function(x) {
    x$B <- x$B[1:max(jags_data_filtered$y.n.site), , ]
    x
  })
}
