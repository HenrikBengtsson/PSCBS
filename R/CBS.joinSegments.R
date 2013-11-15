###########################################################################/**
# @set "class=CBS"
# @RdocMethod joinSegments
#
# @title "Joins neighboring segments such that there is no gap in between them"
#
# \description{
#  @get "title".
#  For instance, consider two neighboring segments [x1,x2] and [x3,x4]
#  with x1 < x2 < x3 < x4.  After join the segments, they are
#  [x1,x23] and [x23,x4] where x23 = (x2 + x3)/2.
# }
#
# @synopsis
#
# \arguments{
#   \item{range}{(optional) A @numeric @vector of length two.}
#   \item{verbose}{See @see "R.utils::Verbose".}
#   \item{...}{Not used.}
# }
#
# \value{
#   Returns an updated @see "CBS" object.
# }
#
# \details{
#   This function assumes only chromosome exists.
#   If more, an error will be thrown.
# }
#
# @author "HB"
#
# @keyword IO
# @keyword internal
#*/###########################################################################
setMethodS3("joinSegments", "CBS", function(fit, range=NULL, verbose=FALSE, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'range':
  if (!is.null(range)) {
    range <- Arguments$getDoubles(range, length=c(2,2));
    stopifnot(range[2] >= range[1]);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Joining segments");
  segs <- getSegments(fit);
  verbose && cat(verbose, "Segments:");
  verbose && print(verbose, segs);
  verbose && cat(verbose, "Range:");
  verbose && print(verbose, range);

  nbrOfSegs <- nrow(segs);
  if (nbrOfSegs > 1) {
    verbose && enter(verbose, "Centering change points");
    prevSeg <- segs[1L,];
    for (ss in 2:nbrOfSegs) {
      currSeg <- segs[ss,];
      currStart <- currSeg[,"start"];
      prevEnd <- prevSeg[,"end"];

      # Sanity check (will give an error if more than one chromosome)
      stopifnot(all(currStart >= prevEnd, na.rm=TRUE));

      # Center CP
      xMid <- (prevEnd + currStart) / 2;

      # Move previous end and current start to this centered CP
      segs[ss,"start"] <- xMid;
      segs[ss-1L,"end"] <- xMid;

      prevSeg <- currSeg;
    } # for (ss ...)
    verbose && exit(verbose);

    # Sanity checks
    stopifnot(all(segs$start[-1] >= segs$end[-nbrOfSegs], na.rm=TRUE));
    stopifnot(all(diff(segs$start) > 0, na.rm=TRUE));
    stopifnot(all(diff(segs$end) > 0, na.rm=TRUE));

    if (nbrOfSegs > 6) {
      verbose && print(verbose, head(segs));
      verbose && print(verbose, tail(segs));
    } else {
      verbose && print(verbose, segs);
    }
  } # if (nbrOfSegs > 1)


  if (!is.null(range)) {
    verbose && enter(verbose, "Adjust for 'range'");
    verbose && cat(verbose, "Range:");
    verbose && print(verbose, range);
    xMin <- min(range, na.rm=TRUE);
    xMax <- max(range, na.rm=TRUE);
    if (nbrOfSegs > 0) {
      # Sanity check
      stopifnot(xMin <= segs[1L,"start"]);
      segs[1L,"start"] <- xMin;
      # Sanity check
      stopifnot(segs[1L,"end"] <= xMax);
      segs[nbrOfSegs,"end"] <- xMax;

      # Sanity checks
      stopifnot(all(segs$start[-1] >= segs$end[-nbrOfSegs], na.rm=TRUE));
      stopifnot(all(diff(segs$start) > 0, na.rm=TRUE));
      stopifnot(all(diff(segs$end) > 0, na.rm=TRUE));


      if (nbrOfSegs > 6) {
        verbose && print(verbose, head(segs));
        verbose && print(verbose, tail(segs));
      } else {
        verbose && print(verbose, segs);
      }
    } # if (nbrOfSegs > 0)
    verbose && exit(verbose);
  } # if (!is.null(range))

  fit <- setSegments(fit, segs);

  segs <- getSegments(fit, splitters=FALSE);
  if (nbrOfSegs > 6) {
    verbose && print(verbose, head(segs));
    verbose && print(verbose, tail(segs));
  } else {
    verbose && print(verbose, segs);
  }
  verbose && exit(verbose);

  fit;
}, private=TRUE) # joinSegments()


############################################################################
# HISTORY:
# 2013-11-14
# o DOCUMENTATION: Added Rd help for joinSegments().
# o CLEANUP: Removed stray variables.
# 2011-11-17
# o Added more sanity checks to joinSegments().
# 2011-09-04
# o Updated joinSegments() to be aware of new column names in CBS.
# 2011-06-14
# o Updated code to recognize new column names.
# 2010-11-21
# o Extracted from segmentByPairedPSCBS.R
############################################################################
