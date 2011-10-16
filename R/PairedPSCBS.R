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


setMethodS3("getSegmentSizes", "PairedPSCBS", function(fit, by=c("length", "count"), ...) {
  by <- match.arg(by);

  data <- getSegments(fit, ...);
  if (by == "length") {
    res <- data[["tcnEnd"]]-data[["tcnStart"]]+1L;
  } else if (by == "count") {
    res <- data[["tcnNbrOfLoci"]];
  }
  res;
})


setMethodS3("updateMeans", "PairedPSCBS", function(fit, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }
 
  verbose && enter(verbose, "Updating mean level estimates");


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Extract the data and segmentation results
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  data <- getLocusData(fit);

  segs <- getSegments(fit, splitters=TRUE);

  nbrOfSegments <- nrow(segs);
  verbose && cat(verbose, "Number of segments: ", nbrOfSegments);

  chromosome <- data$chromosome;
  x <- data$x;
  CT <- data$CT;
  rho <- data$rho;
  muN <- data$muN;
  isSnp <- !is.na(muN);
  isHet <- isSnp & (muN == 1/2);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Update the TCN segments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  isSplitter <- isSegmentSplitter(fit);
  for (ss in seq(length=nbrOfSegments)[!isSplitter]) {
    verbose && enter(verbose, sprintf("Segment %d of %d", ss, nbrOfSegments));
    seg <- segs[ss,];
    verbose && print(verbose, seg);

    chr <- seg[["chromosome"]];
    chrTag <- sprintf("chr%02d", chr);

    for (what in c("tcn", "dh")) {
      keys <- paste(what, c("Start", "End"), sep="");
      xStart <- seg[[keys[1]]];
      xEnd <- seg[[keys[2]]];
      verbose && printf(verbose, "[xStart,xEnd] = [%.0f,%.0f]\n", xStart, xEnd);
      # Nothing todo?
      if (is.na(xStart) && is.na(xEnd)) {
        next;
      }

      stopifnot(xStart <= xEnd);

      # (b) Identify units
      units <- which(chromosome == chr & xStart <= x & x <= xEnd);

      # (c) Adjust for missing values
      if (what == "tcn") {
        value <- CT;
      } else if (what == "dh") {
        value <- rho;
      }
      keep <- which(!is.na(value[units]));
      units <- units[keep];

      # (d) Update mean
      gamma <- mean(value[units]);

      # Sanity check
      stopifnot(length(units) == 0 || !is.na(gamma));

      # Update the segment boundaries, estimates and counts
      key <- paste(what, "Mean", sep="");
      seg[[key]] <- gamma;
    }

    verbose && print(verbose, seg);

    segs[ss,] <- seg;

    verbose && exit(verbose);
  } # for (ss ...)

  verbose && enter(verbose, "Update (C1,C2) per segment");
  # Append (C1,C2) estimates
  tcn <- segs$tcnMean;
  dh <- segs$dhMean;
  C1 <- 1/2*(1-dh)*tcn;
  C2 <- tcn - C1;
  segs$c1Mean <- C1;
  segs$c2Mean <- C2;
  verbose && exit(verbose);

  # Return results
  res <- fit;
  res$output <- segs;

  verbose && exit(verbose);

  res;
}, private=TRUE) # updateMeans()




##############################################################################
# HISTORY
# 2011-01-16
# o BUG FIX: updateMeans() save to the incorrect column names.
# 2011-01-12
# o Added updateMeans() for PairedPSCBS.
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
