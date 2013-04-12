###########################################################################/**
# @RdocDefault exampleData
#
# @title "Gets an example data set"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{name}{A @character string specifying the name of the data set.}
#   \item{...}{Not used.}
# }
#
# \value{
#   Returns @data.frame.
# }
#
# @author "HB"
#
# @keyword IO
#*/###########################################################################
setMethodS3("exampleData", "default", function(name=c("paired.chr01"), ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'name':
  name <- match.arg(name);

  path <- system.file("data-ex", package="PSCBS", mustWork=TRUE);

  if (name == "paired.chr01") {
    filename <- "PairedPSCBS,exData,chr01.Rbin";
  }

  pathname <- Arguments$getReadablePathname(filename, path=path);
  data <- loadObject(pathname);

  data;
})


############################################################################
# HISTORY:
# 2013-04-11
# o Created.
############################################################################
