args <- commandArgs(trailingOnly = TRUE)

.libPaths(c(
  .libPaths(),
  "/packages/R/x86_64-redhat-linux-gnu-library/3.5",
  "/home/lz62/R/3.6.0",
  "/packages/R/3.6.0/lib64/R/library"
))
# "/home/lz62/R/x86_64-redhat-linux-gnu-library/3.6.0",
# "/home/lz62/R/3.6.0",
# "/home/lz62/R/3.5"))
# "/scratch/lz62/R/x86_64-redhat-linux-gnu-library/3.6.0",
# "/scratch/lz62/R/3.6.0"))
print(.libPaths())


load(args[1])
setwd("/scratch/lz62/uplands-ci")
source("src/fit.r")

system.time({
  rowwise_fit(
    a_scenario = a_scenario,
    yaml_info = yaml_info,
    yaml_branches = yaml_branches
  )
})
