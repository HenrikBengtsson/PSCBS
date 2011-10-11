###########################################################################/**
# @set "class=CBS"
# @RdocMethod append
#
# @title "Appends one segmentation result to another"
#
# \description{
#   @get "title".
# }
# 
# @synopsis
#
# \arguments{
#  \item{x, other}{The two @see "CBS" objects to be combined.}
#  \item{other}{A @see "PSCBS" object.}
#  \item{addSplit}{If @TRUE, a "divider" is added between chromosomes.}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a @see "CBS" object of the same class as argument \code{x}.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################  
setMethodS3("append", "CBS", function(x, other, addSplit=TRUE, ...) {
  # To please R CMD check
  this <- x;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'other':
  other <- Arguments$getInstanceOf(other, class(this)[1]);
  for (field in c("data", "output")) {
    thisFF <- this[[field]];
    otherFF <- other[[field]];
    if (ncol(otherFF) != ncol(thisFF)) {
      throw(sprintf("Cannot merge %s objects. Arguments 'other' and 'this' has different number of columns in field '%s': %s != %s", class(this)[1], field, ncol(otherFF), ncol(thisFF)));
    }
  }

  # Argument 'addSplit':
  addSplit <- Arguments$getLogical(addSplit);


  # Allocate results
  res <- this;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Locus data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  data <- getLocusData(this);
  res$data <- rbind(data, getLocusData(other));


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Segmentation data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  indexOffset <- nrow(data);
  fields <- c("output", "segRows");
  for (field in fields[-1]) {
    other[[field]] <- other[[field]] + indexOffset;
  }

  splitter <- if (addSplit) NA else NULL;
  for (field in fields) {
    res[[field]] <- rbind(this[[field]], splitter, other[[field]]);
    rownames(res[[field]]) <- NULL;
  }

  # Sanity check
  ns <- sapply(res[fields], FUN=nrow);
  stopifnot(all(ns == ns[1]));

  res;
}) # append()



setMethodS3("extractRegions", "CBS", function(this, regions, ..., verbose=FALSE) {
  fit <- this;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  updateSegRows <- function(segRows, regions) {
    segRows <- segRows[regions,,drop=FALSE];
    ns <- segRows[,2] - segRows[,1] + 1L;
    from <- c(1L, cumsum(ns)[-length(ns)]);
    to <- from + (ns - 1L);
    segRows[,1] <- from;
    segRows[,2] <- to;
    verbose && str(verbose, segRows);
    segRows;
  } # updateSegRows()

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'regions':
  regions <- Arguments$getIndices(regions, max=nbrOfSegments(fit));

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

 
  verbose && enter(verbose, "Extracting subset by regions");

  verbose && cat(verbose, "Number of regions: ", length(regions));
  verbose && str(verbose, regions);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extract data and estimates
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  data <- getLocusData(fit);
  segRows <- fit$segRows;
  segs <- getSegments(fit);
  params <- fit$params;

  # Sanity checks
  stopifnot(all(!is.na(data$chromosome) & !is.na(data$x)));

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Subset segments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  verbose && enter(verbose, "Update table of segments");
  segsT <- segs[regions,,drop=FALSE];
  verbose && str(verbose, segsT);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Subset data accordingly
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  verbose && enter(verbose, "Update locus data");

  segRows <- segRows[regions,,drop=FALSE];
  from <- segRows[[1]];
  to <- segRows[[2]];
  ok <- (!is.na(from) & !is.na(to));
  from <- from[ok];
  to <- to[ok];
  keep <- logical(nrow(data));
  for (rr in seq(along=from)) {
    keep[from[rr]:to[rr]] <- TRUE;
  }
  keep <- which(keep);
  verbose && printf(verbose, "Identified %d (%.2f%%) of %d data rows:\n", length(keep), 100*length(keep)/nrow(data), nrow(data));
  verbose && str(verbose, keep);

  dataT <- data[keep,,drop=FALSE];
  verbose && str(verbose, dataT);

  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Update 'segRows'
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  verbose && enter(verbose, "Update 'segRows'");

  segRows <- updateSegRows(segRows, regions=regions);
  d <- segRows[regions,] - segRows;
  # Sanity check
  stopifnot(identical(d[,1], d[,2]));
  d <- d[,1];
  verbose && cat(verbose, "Row deltas:");
  verbose && str(verbose, d);

  segRows <- segRows[regions,,drop=FALSE] - d;
  verbose && str(verbose, segRows);
  # Sanity checks
  stopifnot(max(segRows, na.rm=TRUE) <= nrow(dataT));

  verbose && exit(verbose);


  # Create new object
  res <- fit;
  res$data <- dataT;
  res$output <- segsT;
  res$segRows <- segRows;

  verbose && exit(verbose);

  res;
}, protected=TRUE) # extractRegions() 



setMethodS3("mergeTwoSegments", "CBS", function(this, left, verbose=FALSE, ...) {
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

  segs <- getSegments(this);
  segRows <- this$segRows;

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
  idxsUsed <- c();

  # (id) [as in label]
  idxs <- grep("(I|i)d$", fields);
  idxsUsed <- c(idxsUsed, idxs);

  # (chromosome)
  idxs <- grep("chromosome$", fields);
  idxsUsed <- c(idxsUsed, idxs);

  # Starts
  idxs <- grep("(S|s)tart$", fields);
  segT[,idxs] <- apply(segsT[,idxs,drop=FALSE], MARGIN=2, FUN=min, na.rm=TRUE);
  idxsUsed <- c(idxsUsed, idxs);

  # Ends
  idxs <- grep("(E|e)nd$", fields);
  segT[,idxs] <- apply(segsT[,idxs,drop=FALSE], MARGIN=2, FUN=max, na.rm=TRUE);
  idxsUsed <- c(idxsUsed, idxs);

  # Counts
  idxs <- grep("(N|n)brOf", fields);
  segT[,idxs] <- apply(segsT[,idxs,drop=FALSE], MARGIN=2, FUN=sum);
  idxsUsed <- c(idxsUsed, idxs);

  # "Invalidate" remaining entries
  idxsTodo <- setdiff(seq(along=fields), idxsUsed);
  segT[,idxsTodo] <- NA;

  # Update segment table
  segs[rows[1],] <- segT;
  segs <- segs[-rows[2],];

  # Update 'segRows' tables
  segRows[rows[1],2] <- segRows[rows[2],2];
  segRows <- segRows[-rows[2],];
 
  # Create results object
  res <- this;
  res$output <- segs;
  res$segRows <- segRows;

  # Update the mean estimates.
  res <- updateMeans(res);

  verbose && exit(verbose);

  res;
}, protected=TRUE) # mergeTwoSegments()



############################################################################
# HISTORY:
# 2011-10-10
# o Added extractRegions() for CBS.
# 2011-10-08
# o Relabelled column 'id' to 'sampleName' returned by getSegments().
# o BUG FIX: getSegments() for CBS would not set 'id' for "splitter" rows.
# o Added mergeTwoSegments() for CBS.
# o Added updateMeans() for CBS.
# o Added all.equal() for CBS.
# 2011-10-02
# o CLEANUP: Moved getChromosomes(), nbrOfChromosomes(), nbrOfSegments(),
#   nbrOfLoci() and print() to AbstractCBS.
# o Now the CBS class extends the AbstractCBS class.
# 2011-09-04
# o Added writeSegments() for CBS.
# o Added writeLocusData() for CBS.
# o Added getSignalType() for CBS.
# o Added argument 'addCalls' to getLocusData().
# o Added getSampleName() for CBS.
# 2011-09-03
# o Added print() and as.character() for CBS.
# o Added CBS() constructor. Although it rairly will be used
#   it we be a place holder for the documentation.
# 2011-09-02
# o Added nbrOfLoci(), nbrOfSegments(), nbrOfChromosomes() and
#   getChromosomes() for CBS.
# 2010-11-19
# o Added append() for CBS objects.
############################################################################
