# CRAN POLICY: Add precalculated memoization files to the R.cache
# directory, unless running interactively.  The reason for doing this
# is solely to make segmentBy[Non]PairedPSCBS examples to run faster
# on R CMD check but not having to create these memoized files.
.prememoize <- function() {
  if (!interactive()) {
    ## WORKAROUND: Until 'R.utils' no longer attaches itself in some
    ## of its functions internally, we have to attach it here to avoid
    ## 'R CMD check' example test error:
    ##   > PSCBS:::.prememoize()
    ##   Loading required package: R.utils
    ##   Loading required package: R.oo
    ##   Loading required package: R.methodsS3
    ##   Error in attachNamespace(ns, pos = pos, deps) :
    ##     namespace is already attached
    ## /HB 2013-09-26
    pkg <- "R.utils";
    require(pkg, character.only=TRUE, quietly=TRUE) || throw("Package not attached: ", pkg);

    # This will load 'R.cache', if not already done.
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
