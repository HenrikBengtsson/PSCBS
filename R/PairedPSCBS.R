###########################################################################/**
# @RdocClass PairedPSCBS
#
# @title "The PairedPSCBS class"
#
# \description{
#  @classhierarchy
#
#  A PairedPSCBS is an object containing the results from the
#  Paired PSCBS method.
# }
# 
# \usage{PairedPSCBS(fit=list(), ...)}
#
# \arguments{
#   \item{fit}{A @list structure containing the Paired PSCBS results.}
#   \item{...}{Not used.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
# 
# @author
#
# \seealso{
#   The @see "segmentByPairedPSCBS" method returns an object of this class.
# }
#*/###########################################################################
setConstructorS3("PairedPSCBS", function(fit=list(), ...) {
  # Argument 'fit':
  if (!is.list(fit)) {
    throw("Argument 'fit' is not a list: ", class(fit)[1]);
  }

  extend(fit, "PairedPSCBS");
})


setMethodS3("print", "PairedPSCBS", function(x, ...) {
  # To please R CMD check
  fit <- x;

  segs <- as.data.frame(fit, ...);
  print(segs);
}, private=TRUE)


setMethodS3("as.data.frame", "PairedPSCBS", function(x, ...) {
  # To please R CMD check
  fit <- x;

  fit$output;
})



##############################################################################
# HISTORY
# 2011-04-08
# o Added formal constructor for the PairedPSCBS class.
# o Created.
##############################################################################
