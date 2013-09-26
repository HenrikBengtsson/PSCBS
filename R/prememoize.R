# CRAN POLICY: Add precalculated memoization files to the R.cache
# directory, unless running interactively.  The reason for doing this
# is solely to make segmentBy[Non]PairedPSCBS examples to run faster
# on R CMD check but not having to create these memoized files.
.prememoize <- function() {
  if (!interactive()) {
    # This will load 'R.cache', if not already loaded.
    getCachePath <- R.cache::getCachePath;

    path <- "PSCBS/segmentByCBS/sbdry"
    pathS <- system.file("misc/_Rcache", path, package="PSCBS");
    pathD <- getCachePath(path);
    copyDirectory(pathS, pathD, recursive=FALSE);
  }
} # .prememoize()

############################################################################
# HISTORY:
# 2013-09-26
# o CLEANUP: Now .prememoize() no longer attaches 'R.cache', but only
#   loads its namespace.
# 2012-11-05
# o Created.
############################################################################
