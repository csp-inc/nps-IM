get_filename <- function(x, ext, ...) {
  sprintf("%s.%s", paste(gsub("\\.", "-", x), ..., sep = "-"), ext)
}
