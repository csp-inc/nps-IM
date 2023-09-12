library(dplyr)
library(magrittr)
# devtools::install_github("jeremystan/tidyjson");
library(tidyjson)


value_handler <- function(x) {
  tryCatch(
    {
      as.numeric(x)
    },
    warning = function(w) {
      paste0("'", x, "'")
    }
  )
}

json_to_df <- function(tbl_json) {
  # Builds a pipeline to turn 'data/patch.json' into a tidy data frame.

  files_df <- tbl_json %>%
    # Captures the top-level identifier (i.e., file). All patches are specific to
    # a given file --- there can be multiple patches per file, which we'll begin
    # to unpack in the next major series of pipes, below.
    enter_object("file_array") %>%
    gather_array("file_index") %>%
    spread_values(file = jstring("file"))

  patch_attributes_df <- files_df %>%
    # Captures the complete set of attributes (for each file) required to apply
    # a specific patch.
    enter_object("patch_array") %>%
    gather_array("patch_index") %>%
    gather_object("attribute") %>%
    append_values_string("attribute_value") %>%
    select(-document.id)

  filters_df <- files_df %>%
    # Grabs the attributes we need to locate the 'offending' record in each file.
    enter_object("filter_vars") %>%
    gather_array("filter_var_index") %>%
    append_values_string("attribute") %>%
    select(-document.id)

  # Join `patch_attributes_df` and `filters_df` and construct filter strings.
  left_join(patch_attributes_df, filters_df) %>%
    rowwise() %>%
    mutate(filter_string = ifelse(!is.na(filter_var_index),
      paste(
        attribute, "==",
        value_handler(attribute_value)
      ), NA
    )) %>%
    select(-filter_var_index)
}

apply_patch <- function(data, file_path, patch_df) {
  # Create a clone of the 'original', complete, and unpatched data after adding
  # a column called `orig_order` that we can use to reconstruct the original
  # order of the data after all patches have been applied.
  data %<>% mutate(orig_order = seq(1, n()))
  data_clone <- data

  # Subset the complete patches information for the specified file.
  this_patch_df <- patch_df %>% filter(file == file_path)

  # For each patch specified for a given file, pull out the patch parameters.
  patch_indices <- unique(this_patch_df[["patch_index"]])
  for (i in patch_indices) {
    ith_patch_params <- this_patch_df %>% filter(patch_index == i)
    record_locator <-
      paste(na.omit(ith_patch_params[["filter_string"]]), collapse = " & ")
    attr_patches <- ith_patch_params %>% filter(is.na(filter_string))

    # Remove the 'offending' record from the original data object.
    rev_record_locator <- gsub("&", "|", gsub("==", "!=", record_locator))
    data %<>% filter_(rev_record_locator)

    # Create another clone, which we'll use to develop the full ith patch.
    data_for_patch_i <- data_clone
    data_for_patch_i %<>% filter_(record_locator)

    # For each attribute involved in a specific patch, update `data_for_patch_i`.
    for (j in 1:nrow(attr_patches)) {
      attr_j <- attr_patches[["attribute"]][j]
      attr_val_j <- type.convert(attr_patches[["attribute_value"]][j])
      data_for_patch_i[[attr_j]] <- attr_val_j
    }
    data <- rbind(data, data_for_patch_i) %>% arrange(orig_order)
  }
  data %>% select(-orig_order)
}

get_patch_info <- function() {
  read_json("data/patch.json") %>%
    json_to_df()
}
