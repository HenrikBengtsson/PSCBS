###########################################################################/**
# @RdocDefault segmentByCBS
#
# @title "Segment genomic signals using the CBS method"
#
# \description{
#  @get "title" of the \pkg{DNAcopy} package.
#  This is a convenient low-level wrapper for the \code{DNAcopy::segment()}
#  method.  It is intended to be applied to one sample and one chromosome
#  at the time.
#  For more details on the Circular Binary Segmentation (CBS) method 
#  see [1,2].
# }
#
# @synopsis
#
# \arguments{
#   \item{y}{A @numeric @vector of J genomic signals to be segmented.}
#   \item{chromosome}{(Optional) An @integer scalar 
#       (or a @vector of length J contain a unique value).
#       Only used for annotation purposes.}
#   \item{x}{Optional @numeric @vector of J genomic locations.
#            If @NULL, index locations \code{1:J} are used.}
#   \item{w}{Optional @numeric @vector in [0,1] of J weights.}
#   \item{undo}{A non-negative @numeric.  If less than +@Inf, then
#       arguments \code{undo.splits="sdundo"} and \code{undo.SD=undo}
#       are passed to \code{DNAcopy::segment()}.}
#   \item{...}{Additional arguments passed to the segmentation function.}
#   \item{joinSegments}{If @TRUE, there are no gaps between neighboring
#     segments.
#     If @FALSE, the boundaries of a segment are defined by the support
#     that the loci in the segments provides, i.e. there exist a locus
#     at each end point of each segment.  This also means that there
#     is a gap between any neighboring segments, unless the change point
#     is in the middle of multiple loci with the same position.
#     The latter is what \code{DNAcopy::segment()} returns.
#   }
#   \item{knownCPs}{Optional @numeric @vector of known 
#     change point locations.}
#   \item{index}{An @integer @vector of length J specifying the 
#     genomewide indices of the loci.}
#   \item{columnNamesFlavor}{A @character string specifying how column names
#     of the returned data frame should be named.}
#   \item{seed}{An (optional) @integer specifying the random seed to be 
#     set before calling the segmentation method.  The random seed is
#     set to its original state when exiting.  If @NULL, it is not set.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#   Returns the fit object.
# }
# 
# \details{
#   Internally @see "DNAcopy::segment" is used to segment the signals.
#   This segmentation method support weighted segmentation.
#
#   The "DNAcopy::segment" implementation of CBS uses approximation
#   through random sampling for some estimates.  Because of this,
#   repeated calls using the same signals may result in slightly 
#   different results, unless the random seed is set/fixed.
# }
#
# \section{Missing and non-finite values}{
#   Signals as well as genomic positions may contain missing
#   values, i.e. @NAs or @NaNs.  If so, they are silently excluded
#   before performing the segmentation.
#
#   None of the input signals may have infinite values,
#   i.e. -@Inf or +@Inf. If so, an informative error is thrown.
# }
#
# @examples "../incl/segmentByCBS.Rex"
#
# @author
#
# \references{
#  [1] A.B. Olshen, E.S. Venkatraman (aka Venkatraman E. Seshan), 
#      R. Lucito and M. Wigler, \emph{Circular binary segmentation for 
#      the analysis of array-based DNA copy number data},
#      Biostatistics, 2004.\cr
#  [2] E.S. Venkatraman and A.B. Olshen, \emph{A faster circular binary
#      segmentation algorithm for the analysis of array CGH data}. 
#      Bioinformatics, 2007.\cr
# }
#
# @keyword IO
#*/########################################################################### 
setMethodS3("segmentByCBS", "default", function(y, chromosome=0, x=NULL, index=seq(along=y), w=NULL, undo=Inf, ..., joinSegments=TRUE, knownCPs=NULL, columnNamesFlavor=c("PSCBS", "DNAcopy"), seed=NULL, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'y':
  disallow <- c("Inf");
  y <- Arguments$getDoubles(y, disallow=disallow);
  nbrOfLoci <- length(y);

  length2 <- rep(nbrOfLoci, times=2);

  # Argument 'chromosome':
  disallow <- c("Inf");
  chromosome <- Arguments$getIntegers(chromosome, range=c(0,Inf), disallow=disallow);
  if (length(chromosome) > 1) {
    chromosome <- Arguments$getIntegers(chromosome, length=length2, disallow=disallow);
##    # If 'chromosome' is a vector of length J, then it must contain
##    # a unique chromosome.
##    chromosomes <- sort(unique(chromosome));
##    if (length(chromosomes) > 1) {
##      throw("Argument 'chromosome' specifies more than one unique chromosome: ", paste(seqToHumanReadable(chromosomes), collapse=", "));
##    }
##    chromosome <- chromosomes;
  }

  # For future usage
  chrom <- rep(chromosome, length.out=nbrOfLoci);

  # Argument 'x':
  if (is.null(x)) {
    x <- seq(length=nbrOfLoci);
  } else {
    disallow <- c("Inf");
    x <- Arguments$getDoubles(x, length=length2, disallow=disallow);
  }

  # Argument 'index':
  index <- Arguments$getIndices(index);

  # Argument 'w':
  hasWeights <- !is.null(w);
  if (hasWeights) {
    disallow <- c("NA", "NaN", "Inf");
    w <- Arguments$getDoubles(w, range=c(0,1), length=length2, disallow=disallow);
  }

  # Argument 'undo':
  undo <- Arguments$getDouble(undo, range=c(0,Inf));

  # Argument 'cpFlavor':
  joinSegments <- Arguments$getLogical(joinSegments);

  # Argument 'knownCPs':
  if (!is.null(knownCPs)) {
    if (is.null(x)) {
      knownCPs <- Arguments$getIndices(knownCPs, max=nbrOfLoci);
    } else {
      knownCPs <- Arguments$getDoubles(knownCPs);
    }
    if (length(knownCPs) != 2) {
      throw("Currently argument 'knownCPs' can be used to specify the boundaries of the region to be segmented: ", length(knownCPs));
      throw("Support for specifying known change points (argument 'knownCPs') is not yet implemented as of 2010-10-02.");
    }
    if (!joinSegments) {
      throw("Argument 'knownCPs' should only be specified if argument 'joinSegments' is TRUE.");
    }
  }

  # Argument 'columnNamesFlavor':
  columnNamesFlavor <- match.arg(columnNamesFlavor);

  # Argument 'seed':
  if (!is.null(seed)) {
    seed <- Arguments$getInteger(seed);
  }
 
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Segmenting by CBS");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Set the random seed
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!is.null(seed)) {
    verbose && enter(verbose, "Setting (temporary) random seed");
    oldRandomSeed <- NULL;
    if (exists(".Random.seed", mode="integer")) {
      oldRandomSeed <- get(".Random.seed", mode="integer");
    }
    on.exit({
      if (!is.null(oldRandomSeed)) {
        .Random.seed <<- oldRandomSeed;
      }
    }, add=TRUE);
    verbose && cat(verbose, "The random seed will be reset to its original state afterward.");
    verbose && cat(verbose, "Seed: ", seed);
    set.seed(seed);
    verbose && exit(verbose);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Setup up data");
  data <- data.frame(chrom=chrom, x=x, y=y, index=index);
  if (hasWeights) {
    verbose && cat(verbose, "Adding locus-specific weights");
    data$w <- w;
  }
  verbose && str(verbose, data);
  rm(chrom, x, index, y, w); # Not needed anymore
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Drop data points without known genomic positions, because that
  # is what DNAcopy::CNA() will do otherwise.  At the end, we will 
  # undo this such that the returned 'data' object is complete.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ok <- (!is.na(data$chrom) & !is.na(data$x));
  if (any(!ok)) {
    verbose && enter(verbose, "Dropping loci with unknown locations");
    verbose && cat(verbose, "Number of loci dropped: ", sum(!ok));
    data <- data[ok,];
    verbose && exit(verbose);
  }
  rm(ok); # Not needed anymore


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Reorder data points along the genome, because that is what
  # DNAcopy::segment() will return.  At the end, we will undo
  # the sort such that the returned 'data' object is always in
  # the same order and number of loci as the input data.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Ordering data along genome");
  o <- order(data$chrom, data$x, decreasing=FALSE, na.last=TRUE);
  # Any change?
  if (any(o != seq(along=o))) {
    data <- data[o,];
  }
  rm(o); # Not needed anymore
  verbose && str(verbose, data);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Multiple chromosomes?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify all chromosomes, excluding missing values
  chromosomes <- sort(unique(data$chrom), na.last=NA);
  nbrOfChromosomes <- length(chromosomes);
  if (nbrOfChromosomes > 1) {
    verbose && enter(verbose, "Segmenting multiple chromosomes");
    verbose && cat(verbose, "Number of chromosomes: ", nbrOfChromosomes);

    fitList <- list();
    for (kk in seq(length=nbrOfChromosomes)) {
      chromosomeKK <- chromosomes[kk];
      chrTag <- sprintf("Chr%02d", chromosomeKK);
      verbose && enter(verbose, sprintf("Chromosome #%d ('%s') of %d", kk, chrTag, nbrOfChromosomes));

      # Extract subset
      dataKK <- subset(data, chrom == chromosomeKK);
      verbose && str(verbose, dataKK);
      fields <- attachLocally(dataKK, fields=c("y", "chrom", "x", "index"));
      rm(dataKK); # Not needed anymore

      fit <- segmentByCBS(y=y,
                chromosome=chrom, x=x,
                undo=undo,
                joinSegments=joinSegments,
                columnNamesFlavor=columnNamesFlavor,
                ..., 
                seed=NULL,
                verbose=verbose);

      # Sanity checks
      stopifnot(nrow(fit$data) == length(y));
      stopifnot(all.equal(fit$data$y, y));

      rm(list=fields); # Not needed anymore

      verbose && print(verbose, head(as.data.frame(fit)));
      verbose && print(verbose, tail(as.data.frame(fit)));
      
      fitList[[chrTag]] <- fit;

      # Not needed anymore
      rm(fit);
      verbose && exit(verbose);
    } # for (kk ...)

    verbose && enter(verbose, "Merging");
    fit <- Reduce(append, fitList);
    # Not needed anymore
    rm(fitList);
    verbose && str(verbose, fit);
    verbose && exit(verbose);

    verbose && print(verbose, head(as.data.frame(fit)));
    verbose && print(verbose, tail(as.data.frame(fit)));
   
    verbose && exit(verbose);

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Return results
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    return(fit);
  } # if (nbrOfChromosomes > 1) 



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Retrieving segmentation function
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Retrieving the fit function");
  pkgName <- "DNAcopy";
  # Assert that package is installed
  isPackageInstalled(pkgName) || throw("Package is not installed: ", pkgName);
  pkg <- packageDescription(pkgName);
  pkgVer <- pkg$Version;
  pkgDetails <- sprintf("%s v%s", pkgName, pkgVer);

  methodName <- "segment";
  verbose && cat(verbose, "Method: ", methodName);
  verbose && cat(verbose, "Package: ", pkgDetails);

  # We need to load package
  require(pkgName, character.only=TRUE) || throw("Package not loaded: ", pkgName);

  # Get the fit function for the segmentation method
#  fitFcn <- getExportedValue(pkgName, methodName);
  fitFcn <- getFromNamespace(methodName, pkgName);
  verbose && str(verbose, "Function: ", fitFcn);
  formals <- formals(fitFcn);
  verbose && cat(verbose, "Formals:");
  verbose && str(verbose, formals);
  verbose && exit(verbose);
 

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setting up arguments to pass to segmentation function
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  verbose && enter(verbose, "Setting up method arguments");

  verbose && enter(verbose, "Setting up ", pkgName, " data structure"); 

  sampleName <- "y";  # This is going to be the name of the data field

  cnData <- DNAcopy::CNA(
    genomdat  = data$y,
    chrom     = data$chrom,
    data.type = "logratio",
    maploc    = data$x,
    sampleid  = sampleName,
    presorted = TRUE
  );
  verbose && str(verbose, cnData);
  names(cnData)[3] <- sampleName;
  verbose && str(verbose, cnData);
  verbose && exit(verbose);
  rm(sampleName);  # Not needed anymore

  # Sanity check
  stopifnot(nrow(cnData) == nrow(data));

  params <- list();
  if (hasWeights) {
    params$weights <- data$w;
  }

  if (undo < Inf) {
    params$undo.splits <- "sdundo";
    params$undo.SD <- undo;
  }

  verbose && cat(verbose, "Segmentation parameters:");
  verbose && str(verbose, params);

  userArgs <- list(...);
  if (length(userArgs) > 0) {
    verbose && cat(verbose, "User arguments:");
    verbose && str(verbose, userArgs);
    # Assign/overwrite by user arguments
    for (ff in names(userArgs)) {
      params[[ff]] <- userArgs[[ff]];
    }
  }

  verbose && cat(verbose, "Segmentation and user parameters:");
  verbose && str(verbose, params);

  # Cleaning out unknown parameters
  keep <- (names(params) %in% names(formals));
  params <- params[keep];

  args <- c(list(cnData), params, verbose=as.logical(verbose));
  verbose && cat(verbose, "Final arguments:");
  verbose && str(verbose, args);

  verbose && exit(verbose);
 

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Calling segmentation function
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, sprintf("Calling %s() of %s", methodName, pkgName));

  # WORKAROUND for the case when there are no data points.
  nbrOfNonMissingLoci <- sum(!is.na(cnData$y));
  if (nbrOfNonMissingLoci == 0) {
    args[[1]] <- DNAcopy::CNA(genomdat=0, chrom=0, maploc=0);
  }

  # In case the method writes to stdout, we capture it
  # Note: DNAcopy::segment() *does* this.
  stdout <- capture.output({
    # Does not work, because some internal function of the fit function
    # may only be accessible from within the namespace
    # How to do this for DNAcopy::segment()? /HB
##    fit <- do.call(fitFcn, args);
    # This works, but requires that one loads the package and that the
    # function is not masked in the search() path.
    t <- system.time({
      fit <- do.call(methodName, args);
    }); 
    # Drop the 'call' (because it will be huge due to the do.call() call)
    fit$call <- NULL;
  });
  attr(fit, "processingTime") <- t;
  attr(fit, "pkgDetails") <- pkgDetails;
  attr(fit, "randomSeed") <- seed;

  # WORKAROUND for the case when there are no data points.
  if (nbrOfNonMissingLoci == 0) {
    # Drop dummy data point...
    fit$data <- cnData; ## fit$data[-1,,drop=FALSE];
    # ...dummy region found
    fit$output <- fit$output[-1,,drop=FALSE];
    fit$segRows <- fit$segRows[-1,,drop=FALSE];
  }

  verbose && cat(verbose, "Captured output that was sent to stdout:");
  stdout <- paste(stdout, collapse="\n");
  verbose && cat(verbose, stdout);

  verbose && cat(verbose, "Fitting time (in seconds):");
  verbose && print(verbose, t);

  verbose && cat(verbose, "Fitting time per 1000 loci (in seconds):");
  verbose && print(verbose, 1000*t/nbrOfLoci);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Restructure
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Restructuring results");

  # Coerce
  fit$output$num.mark <- as.integer(fit$output$num.mark);

  # Coerce 'chrom' to a plain integer
  fit$data$chrom <- unclass(fit$data$chrom);

  # Store genomewide index
  fit$data$index <- data$index;
  rm(data);

  verbose && exit(verbose);



  params <- list(
    joinSegments = joinSegments,
    knownCPs = knownCPs,
    seed = seed
  );

  fit$params <- params;

  class(fit) <- c("CBS", class(fit));

  # Sanity checks
  segRows <- fit$segRows;
  stopifnot(all(segRows[,1] <= segRows[,2], na.rm=TRUE));
  stopifnot(all(segRows[-nrow(segRows),2] < segRows[-1,1], na.rm=TRUE));


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Renaming column names
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (columnNamesFlavor != "DNAcopy") {
    segs <- fit$output;
    names <- colnames(segs);
    if (columnNamesFlavor == "PSCBS") {
      names <- gsub("ID", "id", names, fixed=TRUE);
      names <- gsub("seg.mean", "mean", names, fixed=TRUE);
      names <- gsub("chrom", "chromosome", names, fixed=TRUE);
      names <- gsub("num.mark", "nbrOfLoci", names, fixed=TRUE);
      names <- gsub("loc.", "", names, fixed=TRUE); # loc.start, loc.end
    }
    colnames(segs) <- names;
    fit$output <- segs;
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Join segments?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (joinSegments) {
    fit <- joinSegments(fit, range=knownCPs, verbose=verbose);

    # Sanity checks
    segRows <- fit$segRows;
    stopifnot(all(segRows[,1] <= segRows[,2], na.rm=TRUE));
    stopifnot(all(segRows[-nrow(segRows),2] < segRows[-1,1], na.rm=TRUE));
  }


  verbose && cat(verbose, "Results object:");
  verbose && str(verbose, fit);

  verbose && exit(verbose); 


  verbose && exit(verbose);
  
  fit; 
}) # segmentByCBS()



############################################################################
# HISTORY:
# 2011-09-01
# o GENERALIZATION: Now segmentByCBS() can process multiple chromosomes.
# o Now the random seed is set at the very beginning of the code, which
#   should not make a difference, i.e. it should give identical results.
# 2011-06-14
# o GENERALIZATION: Added argument 'columnNamesFlavor' to segmentByCBS().
# o CONVENTION: Changed the column names of returned data frames. 
#   They now follow the camelCase naming convention and are shorter.
# 2011-05-31
# o Now explicitly using DNAcopy::nnn() to call DNAcopy functions. 
# 2011-04-07
# o ROBUSTNESS: Added 'segRows' field validation in segmentByCBS().
# 2010-12-01
# o Now segmentByCBS() is utilizing 'segRows' from DNAcopy::segment(),
#   which makes it possible to drop all code of trying to infer which
#   loci belong to which segments.
# o Now the 'data' object returned also contains column 'index'.
# 2010-12-01
# o Now the genomewide index is always stored in the 'data' field.
# o Added argument 'index' to segmentByCBS().
# 2010-11-30
# o Now segmentByCBS() returns a field 'lociToExclude'.
# 2010-11-28
# o BUG FIX: The algorithm in segmentByCBS() that infers which loci( of
#   the ones share the same genomic positions) that should be exclude
#   from each segment did not take missing signals into account.
# 2010-11-21
# o Now segmentByCBS(..., joinSegments=TRUE) utilizes joinSegments().
# 2010-11-20
# o Now it is possible to specify the boundaries of the regions to be
#   segmented as known change points via argument 'knownCPs'.
# 2010-11-19
# o Added argument 'joinSegments' to segmentByCBS() in order to specify
#   if neighboring segments should be joined or not.
# o Now segmentByCBS() returns an object of class CBS.
# o Now segmentByCBS() allows for unknown genomic positions.
# o Now segmentByCBS() allows for missing signals.
# o Added argument 'preservOrder' to segmentByCBS().  If TRUE, then
#   the loci in the returned 'data' object are ordered as the input
#   data, otherwise it is ordered along the genome.
# 2010-11-16
# o Now the 'data' object returned by segmentByCBS() contains field
#   'index' if and only if the loci had to be reorder along the genome.
# 2010-11-02
# o Added argument 'undo' to segmentByCBS(), which corresponds to
#   undo.splits="sdundo" and undo.SD=undo, if undo < Inf.
# 2010-10-25
# o Now segmentByCBS() also returns element 'lociNotPartOfSegment',
#   if there are segments that share end points, which can happen if
#   a change point is called in middle of a set of loci that have the
#   same genomic positions.  In such cases, 'lociNotPartOfSegment'
#   specifies which loci are *not* part of which segment.  Then by
#   identifying the loci that are within a segment by their positions
#   and excluding any of the above, one knows exactly which loci
#   CBS included in each segment.
# 2010-10-02
# o Added argument optional 'chromosome'.
# 2010-09-02
# o ROBUSTNESS: Now segmentByCBS() also works if there are no data points.
# 2010-07-09
# o Created from segmentByCBS() for RawGenomicSignals in aroma.core.
#   The latter will eventually call this method.
############################################################################
