###########################################################################/**
# @set "class=data.frame"
# @RdocMethod gapsToSegments
#
# @title "Gets the genomic segments that are complementary to the gaps"
#
# \description{
#  @get "title", with default chromosome boundaries being \code{-Inf}
#  and \code{+Inf}.
# }
#
# @synopsis
#
# \arguments{
#   \item{gaps}{A @data.frame with columns \code{chromosome}, \code{start},
#     and \code{stop}.}
#   \item{resolution}{A non-negative @numeric specifying the minimum
#     length unit, which by default equals one nucleotide/base pair.}
#   \item{...}{Not used.}
# }
#
# \value{
#   Returns @data.frame with columns \code{chromosome} (if that argument
#   is given), \code{start}, \code{stop} and \code{length}.
# }
# 
# @author
#
# \seealso{
#   @see "findLargeGaps".
# }
#
# @keyword IO
# @keyword internal
#*/###########################################################################  
setMethodS3("gapsToSegments", "data.frame", function(gaps, resolution=1L, ...) {
  # Argument 'gaps':
  keys <- colnames(gaps);
  stopifnot(all(is.element(c("chromosome", "start", "end"), keys)));


  keepL <- is.finite(gaps$start);
  gapsL <- gaps[keepL,];
  gapsL$end <- gapsL$start - resolution;
  gapsL$start <- -Inf;

  keepR <- is.finite(gaps$end);
  gapsR <- gaps[keepR,];
  gapsR$start <- gapsR$end + resolution;
  gapsR$end <- +Inf;

  knownSegments <- rbind(gapsL, gaps, gapsR);
  o <- with(knownSegments, order(chromosome, start, end));
  knownSegments <- knownSegments[o,];
  rownames(knownSegments) <- NULL;

  knownSegments$length <- knownSegments$end - knownSegments$start;

  knownSegments;
}) # gapsToSegments()


###############################################################################
# HISTORY:
# 2011-11-22
# o Made gapsToSegments() a method for 'data.frame' class.
# o Renamed gapsToKnownSegments() to gapsToSegments().
# 2011-10-xx
# o Created.
###############################################################################
