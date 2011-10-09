###########################################################################/**
# @set "class=PairedPSCBS"
# @RdocMethod mergeTwoSegments
#
# @title "Merge two neighboring segments"
#
# \description{
#   @get "title" by recalculating segment statistics.
# }
# 
# @synopsis
#
# \arguments{
#  \item{left}{An @integer specifying the segments (left, left+1)
#    to be merged.}
#  \item{verbose}{A @logical or a @see "R.utils::Verbose" object.}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a @see "PairedPSCBS" with one less segment.
# }
#
# @author
#
# \seealso{
#   To drop regions (a connected set of segments) 
#   see @seemethod "dropByRegions".
#   @seeclass
# }
#*/###########################################################################  
setMethodS3("mergeTwoSegments", "PairedPSCBS", function(this, left, verbose=FALSE, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  nbrOfSegments <- nbrOfSegments(this);
  # Argument 'left':
  left <- Arguments$getIndex(left, max=nbrOfSegments-1L);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Merging to segments");
  verbose && printf(verbose, "Segments to be merged: %s & %s\n", left, left+1);
  verbose && cat(verbose, "Number of segments before merging: ", nbrOfSegments);
  verbose && cat(verbose, "Number of segments after merging: ", nbrOfSegments-1L);

  segs <- this$output;
  tcnSegRows <- this$tcnSegRows;
  dhSegRows <- this$dhSegRows;

  rows <- c(left,left+1);
  segsT <- segs[rows,,drop=FALSE];

  # Sanity check
  chrs <- segsT[["chromosome"]];
  if (chrs[1] != chrs[2]) {
    throw("Cannot merge segments that are on different chromosomes: ", chrs[1], " != ", chrs[2]);
  }

  # Merge segments
  segT <- segsT[1,];
  fields <- colnames(segsT);

  # (chromosome, tcnId, dhId)
  idxsUsed <- 1:3;

  # Starts
  idxs <- grep("Start$", fields);
  segT[,idxs] <- apply(segsT[,idxs,drop=FALSE], MARGIN=2, FUN=min, na.rm=TRUE);
  idxsUsed <- c(idxsUsed, idxs);

  # Ends
  idxs <- grep("End$", fields);
  segT[,idxs] <- apply(segsT[,idxs,drop=FALSE], MARGIN=2, FUN=max, na.rm=TRUE);
  idxsUsed <- c(idxsUsed, idxs);

  # Counts
  idxs <- grep("NbrOf", fields);
  segT[,idxs] <- apply(segsT[,idxs,drop=FALSE], MARGIN=2, FUN=sum);
  idxsUsed <- c(idxsUsed, idxs);

  # "Invalidate" remaining entries
  idxsTodo <- setdiff(seq(along=fields), idxsUsed);
  segT[,idxsTodo] <- NA;

  # Update segment table
  segs[rows[1],] <- segT;
  segs <- segs[-rows[2],];

  # Update 'segRows' tables
  segRows <- tcnSegRows;
  segRows[rows[1],2] <- segRows[rows[2],2];
  segRows <- segRows[-rows[2],];
  tcnSegRows <- segRows;

  segRows <- dhSegRows;
  segRows[rows[1],2] <- segRows[rows[2],2];
  segRows <- segRows[-rows[2],];
  dhSegRows <- segRows;

  # Create results object
  res <- this;
  res$output <- segs;
  res$tcnSegRows <- tcnSegRows;
  res$dhSegRows <- dhSegRows;

  # Update the mean estimates.
  res <- updateMeans(res);

  verbose && exit(verbose);

  res;
}, private=TRUE)



###########################################################################/**
# @set "class=PairedPSCBS"
# @RdocMethod dropByRegions
#
# @title "Drops chromosomal regions (a connected set of segments)"
#
# \description{
#   @get "title" each of a cetain size (number of segments).
# }
# 
# @synopsis
#
# \arguments{
#  \item{regions}{An @integer @vector of length R specifying the indices
#    of the left most segment in each of the R regions to be dropped.}
#  \item{H}{A non-negative @integer specifying the size of each region,
#    i.e. the number of segments per region.}
#  \item{...}{Additional arguments passed to @seemethod "extractByRegions".}
#  \item{verbose}{A @logical or a @see "R.utils::Verbose" object.}
# }
#
# \value{
#   Returns a @see "PairedPSCBS" with (at most) R*H segments dropped.
#   If some regions overlap (share segments), then fewer than R*H segments
#   are dropped.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################  
setMethodS3("dropByRegions", "PairedPSCBS", function(this, regions, H=1, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  nbrOfSegments <- nbrOfSegments(this);
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

  # Identify regions to drop
  Hs <- seq(length=H);
  regions <- regions - 1L;
  regions <- as.list(regions);
  regions <- lapply(regions, FUN=function(region) region + Hs);
  regions <- unlist(regions, use.names=FALSE);
  regions <- sort(unique(regions));
  verbose && cat(verbose, "Final set of segments to be dropped:");
  verbose && str(verbose, regions);

  # Identify regions to keep
  allRegions <- seq(length=nbrOfSegments);
  keepRegions <- setdiff(allRegions, regions);
  verbose && cat(verbose, "Final set of segments to be kept:");
  verbose && str(verbose, keepRegions);

  res <- extractByRegions(this, regions=keepRegions, ...);
  res$dropped <- extractByRegions(this, regions=regions, ...);

  verbose && exit(verbose);

  res;
}, private=TRUE)




setMethodS3("extractLocusLevelC1C2", "PairedPSCBS", function(fit, ...) {
  data <- getLocusData(fit);
  C <- data$CT;
  rho <- data$rho;
  C1 <- 1/2*(1-rho)*C;
  C2 <- C-C1;
  data.frame(C1=C1, C2=C2);
}, private=TRUE) # extractLocusLevelC1C2()


setMethodS3("extractLocusLevelTCN", "PairedPSCBS", function(fit, ...) {
  data <- getLocusData(fit);
  C <- data$CT;
}, private=TRUE) # extractLocusLevelTCN()


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

  segs <- getSegments(fit);
  keep <- is.finite(segs$chromosome);
  segs <- segs[keep,,drop=FALSE];

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
  for (ss in seq(length=nbrOfSegments)) {
    verbose && enter(verbose, sprintf("Segment %d of %d", ss, nbrOfSegments));
    seg <- segs[ss,];

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
      key <- paste(what, "mean", sep=".");
      seg[[key]] <- gamma;
    }

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



############################################################################
# HISTORY:
# 2011-10-8
# o ROBUSTIFICATION: Uses drop=FALSE in mergeTwoSegments() for PairedPSCBS.
# 2011-10-02
# o DOCUMENTATION: Added Rdoc help to mergeTwoSegments() & dropByRegions().
# o Added verbose statements to the above to functions.
# 2011-06-14
# o Updated code to recognize new column names.
# 2011-01-18
# o BUG FIX: Fields 'tcnSegRows' and 'dhSegRows' were not updated by
#   mergeTwoSegments() for PairedPSCBS.
# 2011-01-14
# o Moved extractByRegions() and estimateStdDevForHeterozygousBAF() to
#   psCBS v0.9.36.
# o Now extractByRegions() utilizes the 'segRows' field.
# o Added estimateStdDevForHeterozygousBAF().
# 2011-01-12
# o Added updateMeans() for PairedPSCBS.
# o Added dropByRegions().
# o Added extractByRegions() and extractByRegion().
# o Created.
############################################################################
