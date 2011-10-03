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

  extend(PSCBS(fit=fit, ...), "PairedPSCBS");
})



##############################################################################
# HISTORY
# 2011-10-02
# o CLEANUP: Moved print() and as.data.frame() to PSCBS.
# o Added Rdoc help.
# o Now the constructor of PairedPSCBS calls that of PSCBS.
# 2011-06-28
# o DOCUMENTATION: Added Rd help for as.data.frame() of PairedPSCBS.
# 2011-04-08
# o Added formal constructor for the PairedPSCBS class.
# o Created.
##############################################################################
