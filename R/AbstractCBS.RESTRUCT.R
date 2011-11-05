###########################################################################/**
# @set "class=AbstractCBS"
# @RdocDocumentation "Restructuring AbstractCBS objects"
# @alias RestructPSCBS
#
# \description{
#   This page describes available methods for restructuring an 
#   @see "AbstractCBS" object.
#
#   \itemize{
#     \item @seemethod "append"
#           - Appends one @see "AbstractCBS" to another.
#
#     \item @seemethod "extractChromosomes" /
#           @seemethod "extractChromosome"
#           - Extracts an @see "AbstractCBS" with the specified chromosomes.
#
#     \item @seemethod "extractSegments" /
#           @seemethod "extractSegment"
#           - Extracts an @see "AbstractCBS" with the specified segments.
#
#     \item @seemethod "extractRegions" /
#           @seemethod "extractRegion"
#           - Extracts an @see "AbstractCBS" with the specified regions
#             each of a certain size, where a region is defined as a
#             connected set of segments.
#
#     \item @seemethod "dropRegions" /
#           @seemethod "dropRegion"
#           - Drops specified regions and returns an @see "AbstractCBS"
#             without them.
#
#     \item @seemethod "dropChangePoint" /
#           @seemethod "mergeTwoSegments"
#           - Drops a change point by merging two neighboring segments
#             and recalculating the statistics for the merged segment
#             before returning an @see "AbstractCBS".
#
#     \item @seemethod "mergeThreeSegments"
#           - Merges a segment with its two flanking segments
#             and recalculating the statistics for the merged segment
#             before returning an @see "AbstractCBS".
#   }
#
#   All of the above methods are implemented for @see "CBS" and 
#   @see "PairedPSCBS" objects.
# }
# 
# @author
#
# \seealso{
#   @seeclass
# }
#
# @keyword internal
#*/###########################################################################  


###########################################################################/**
# @set "class=AbstractCBS"
# @RdocMethod append
#
# @title "Appends one segmentation result to another"
#
# \description{
#   @get "title",
#   where both holds segmentation results \emph{of the same sample}.
# }
# 
# @synopsis
#
# \arguments{
#  \item{x, other}{The two @see "AbstractCBS" objects to be combined.}
#  \item{addSplit}{If @TRUE, a "divider" is added between chromosomes.}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a object of the same class as argument \code{x}.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################  
setMethodS3("append", "AbstractCBS", abstract=TRUE);



setMethodS3("extractChromosomes", "AbstractCBS", abstract=TRUE, protected=TRUE);


setMethodS3("extractChromosome", "AbstractCBS", function(x, chromosome, ...) {
  # To please R CMD check
  this <- x;

  # Argument 'chromosome':
  chromosome <- Arguments$getInteger(chromosome, disallow=c("NaN", "Inf"));

  extractChromosomes(this, chromosomes=chromosome, ...);
}, protected=TRUE)



setMethodS3("extractSegments", "AbstractCBS", abstract=TRUE, protected=TRUE);

setMethodS3("extractSegment", "AbstractCBS", function(this, idx, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'region':
  idx <- Arguments$getIndex(idx, max=nbrOfSegments(this, splitters=TRUE));

  extractSegments(this, idxs=idx, ...);
}, private=TRUE) # extractSegment()


setMethodS3("extractRegions", "AbstractCBS", function(this, regions, H=1, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  nbrOfSegments <- nbrOfSegments(this, splitters=TRUE);

  # Argument 'regions':
  regions <- Arguments$getIndices(regions, max=nbrOfSegments);

  # Argument 'H':
  H <- Arguments$getInteger(H, range=c(0,Inf));

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Extract regions of a certain length");

  verbose && cat(verbose, "Left-most segments of regions to be extracted:");
  verbose && str(verbose, regions);
  verbose && cat(verbose, "Number of segments in each region: ", H);


  # Identify segments to keep
  Hs <- seq(length=H);
  regions <- regions - 1L;
  regions <- as.list(regions);
  segments <- lapply(regions, FUN=function(region) region + Hs);
  segments <- unlist(segments, use.names=FALSE);
  segments <- sort(unique(segments));
  verbose && cat(verbose, "Final set of segments to be extracted:");
  verbose && str(verbose, segments);

  res <- extractSegments(this, idxs=segments, ...);

  verbose && exit(verbose);

  res;
}, protected=TRUE) # extractRegions()



setMethodS3("extractRegion", "AbstractCBS", function(this, region, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'region':
  region <- Arguments$getIndex(region, max=nbrOfSegments(this, splitters=TRUE));

  extractRegions(this, regions=region, ...);
}, private=TRUE) # extractRegion()



###########################################################################/**
# @RdocMethod mergeTwoSegments
# @aliasmethod dropChangePoint
#
# @title "Merge two neighboring segments"
#
# \description{
#   @get "title" into one segment, which is done by dropping their 
#   common change point and recalculating the segment statistics.
# }
# 
# @synopsis
#
# \arguments{
#  \item{left}{An @integer specifying the segments (left, left+1)
#    to be merged.}
#  \item{update}{If @TRUE, segment statistics are updated.}
#  \item{verbose}{A @logical or a @see "R.utils::Verbose" object.}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns an @see "AbstractCBS" of the same class with one less segment.
# }
#
# @author
#
# \seealso{
#   To merge a segment and its two flanking segments, see
#   @seemethod "mergeThreeSegments".
#   To drop regions (a connected set of segments) 
#   see @seemethod "dropRegions".
#   @seeclass
# }
#*/###########################################################################  
setMethodS3("mergeTwoSegments", "AbstractCBS", abstract=TRUE, protected=TRUE);


setMethodS3("dropChangePoint", "AbstractCBS", function(fit, idx, ...) {
  # Argument 'idx':
  idx <- Arguments$getIndex(idx, max=nbrOfChangePoints(fit, splitters=TRUE));

  mergeTwoSegments(fit, left=idx, ...);
}, protected=TRUE)



###########################################################################/**
# @RdocMethod mergeThreeSegments
#
# @title "Merge a segment and its two flanking segments"
#
# \description{
#   @get "title" into one segment, and recalculating the segment statistics.
# }
# 
# @synopsis
#
# \arguments{
#  \item{middle}{An @integer specifying the three segments
#    (middle-1, middle, middle+1) to be merged.}
#  \item{...}{Additional arguments passed to @seemethod "mergeTwoSegments".}
# }
#
# \value{
#   Returns an @see "AbstractCBS" of the same class with two less segment.
# }
#
# @author
#
# \seealso{
#   Internally @seemethod "mergeTwoSegments" is used.
#   @seeclass
# }
#*/###########################################################################  
setMethodS3("mergeThreeSegments", "AbstractCBS", function(fit, middle, ...) {
  # Argument 'middle':
  S <- nbrOfSegments(fit, splitters=TRUE);
  middle <- Arguments$getIndex(middle, range=c(2L, S));

  # Assert that the three segments are on the same chromosome
  idxs <- middle + c(-1L, 0, +1L);
  fitT <- extractSegments(fit, idxs);
  chrs <- getChromosomes(fitT);
  if (length(chrs) != 1) {
    throw("Argument 'middle' specifies a segment that is at the very end of a chromosome: ", middle);
  }
  rm(fitT);

  fit <- mergeTwoSegments(fit, left=middle, ...);
  fit <- mergeTwoSegments(fit, left=middle-1L, ...);
  fit;
}) # mergeThreeSegments()



###########################################################################/**
# @RdocMethod dropRegions
# @aliasmethod dropRegion
#
# @title "Drops chromosomal regions (a connected set of segments)"
#
# \description{
#   @get "title" each of a cetain size (number of segments).
#   \emph{None of the statistics are recalculated}.
# }
# 
# @synopsis
#
# \arguments{
#  \item{regions}{An @integer @vector of length R specifying the indices
#    of the left most segment in each of the R regions to be dropped.}
#  \item{H}{A non-negative @integer specifying the size of each region,
#    i.e. the number of segments per region.}
#  \item{...}{Additional arguments passed to @seemethod "extractRegions".}
#  \item{asMissing}{If @TRUE, dropped segments are replaced by missing values,
#    otherwise they are truly dropped.}
#  \item{verbose}{A @logical or a @see "R.utils::Verbose" object.}
# }
#
# \value{
#   Returns an @see "AbstractCBS" object of the same class with (at most)
#   R*H segments dropped.
#   If some regions overlap (share segments), then fewer than R*H segments
#   are dropped.
# }
#
# @author
#
# \seealso{
#   Internally @seemethod "extractRegions" is used.
#   See also @seemethod "dropChangePoint" and @seemethod "mergeTwoSegments".
#   @seeclass
# }
#*/###########################################################################  
setMethodS3("dropRegions", "AbstractCBS", function(this, regions, H=1, ..., asMissing=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  nbrOfSegments <- nbrOfSegments(this, splitters=TRUE);
  # Argument 'regions':
  regions <- Arguments$getIndices(regions, max=nbrOfSegments);

  # Argument 'H':
  H <- Arguments$getInteger(H, range=c(0,Inf));

  # Argument 'asMissing':
  asMissing <- Arguments$getLogical(asMissing);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Dropping regions of a certain length");

  verbose && cat(verbose, "Left-most segments of regions to be dropped:");
  verbose && str(verbose, regions);
  verbose && cat(verbose, "Number of segments in each region: ", H);

  # Nothing to do?
  if (H == 0) {
    verbose && cat(verbose, "Nothing to do. No segments will be dropped.");
    verbose && exit(verbose);
    return(this);
  }

  # Identify segments to drop
  Hs <- seq(length=H);
  regions <- regions - 1L;
  regions <- as.list(regions);
  regions <- lapply(regions, FUN=function(region) region + Hs);
  regions <- unlist(regions, use.names=FALSE);
  regions <- sort(unique(regions));
  verbose && cat(verbose, "Final set of segments to be dropped:");
  verbose && str(verbose, regions);

  # Identify segments to keep
  allRegions <- seq(length=nbrOfSegments);
  keepSegments <- setdiff(allRegions, regions);
  verbose && cat(verbose, "Final set of segments to be kept:");
  verbose && str(verbose, keepSegments);

  dropped <- extractRegions(this, regions=regions, ...);
  res <- this;
  if (length(regions) > 0) {  
    if (asMissing) {
      segs <- getSegments(res, splitters=TRUE);
      pattern <- "(chromosome|id|start|end)$";

      # TODO/AD HOC: Should be class specific /HB 2011-10-17
      pattern <- "(chromosome|id)$";
      excl <- grep(pattern, colnames(segs), ignore.case=TRUE, invert=TRUE);
      segs[regions,excl] <- NA;
      res$output <- segs;

      # TODO/AD HOC: Should be class specific /HB 2011-10-17
      for (ff in grep("segRows", names(res), ignore.case=TRUE, value=TRUE)) {
        res[[ff]][regions,] <- NA;
      }
    } else {
      res <- extractRegions(res, regions=keepSegments, ...);
    }
  }
  res$dropped <- dropped;

  # Sanity check
  if (asMissing) {
    stopifnot(nbrOfSegments(res, splitters=TRUE) == nbrOfSegments(this, splitters=TRUE));
  } else {
    stopifnot(nbrOfSegments(res, splitters=TRUE) + length(regions) == nbrOfSegments(this, splitters=TRUE));
  }

  verbose && exit(verbose);

  res;
}, protected=TRUE)


setMethodS3("dropRegion", "AbstractCBS", function(fit, region, ...) {
  # Argument 'region':
  region <- Arguments$getIndex(region);

  dropRegions(fit, regions=region, ...);
}, protected=TRUE)




# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# DEPRECATED BEGIN
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethodS3("dropByRegions", "AbstractCBS", function(fit, ...) {
  dropRegions(fit, ...);
}, protected=TRUE, deprecated=TRUE)

setMethodS3("dropByRegion", "AbstractCBS", function(fit, ...) {
  dropRegion(fit, ...);
}, protected=TRUE, deprecated=TRUE)

setMethodS3("extractByChromosomes", "AbstractCBS", function(fit, ...) {
  extractChromosomes(fit, ...);
}, protected=TRUE, deprecated=TRUE)

setMethodS3("extractByChromosome", "AbstractCBS", function(fit, ...) {
  extractChromosome(fit, ...);
}, protected=TRUE, deprecated=TRUE)

setMethodS3("extractByRegions", "AbstractCBS", function(fit, ...) {
  extractRegions(fit, ...);
}, protected=TRUE, deprecated=TRUE)

setMethodS3("extractByRegion", "AbstractCBS", function(fit, ...) {
  extractRegion(fit, ...);
}, protected=TRUE, deprecated=TRUE)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# DEPRECATED END
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


############################################################################
# HISTORY:
# 2011-11-04
# o BUG FIX: extractSegment() for AbstractCBS would give an error, because
#   it called itself instead of extractSegments().
# 2011-10-21
# o Added mergeThreeSegments() to AbstractCBS.
# 2011-10-17
# o Added argument 'asMissing' to dropRegions() for AbstractCBS.
# 2011-10-14
# o Added implementation of extractRegions() for AbstractCBS, which
#   utilizes extractSegments().
# o Added abstract extractSegments() and extractSegment() for AbstractCBS.
# 2011-10-10
# o Added extractRegion()/dropRegion() and extractRegions()/dropRegions()
#   for AbstractCBS, where former ones are wrappers for the latter ones.
# o Added dropChangePoint() for AbstractCBS, which is just a
#   "name wrapper" for mergeTwoSegments().
# 2011-10-08
# o Added abstract updateMeans() for AbstractCBS.
# o Added abstract mergeTwoSegments() for AbstractCBS.
# 2011-10-02
# o Created.
############################################################################
