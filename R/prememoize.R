# CRAN POLICY: Add precalculated memoization files to the R.cache 
# directory, unless running interactively.  The reason for doing this 
# is solely to make segmentBy[Non]PairedPSCBS examples to run faster
# on R CMD check but not having to create these memoized files.
.prememoize <- function() {
  if (!interactive()) {
    require("R.cache") || throw("Package not loaded: R.cache");
    path <- "PSCBS/segmentByCBS/sbdry"
    pathS <- system.file("misc/_Rcache", path, package="PSCBS");
    pathD <- getCachePath(path);
    copyDirectory(pathS, pathD, recursive=FALSE);
  }
} # .prememoize()

############################################################################
# HISTORY:
# 2012-11-05
# o Created.
############################################################################
