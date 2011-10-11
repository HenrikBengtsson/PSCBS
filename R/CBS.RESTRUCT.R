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
