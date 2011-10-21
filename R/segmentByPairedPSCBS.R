###########################################################################/**
# @RdocDefault segmentByPairedPSCBS
# @alias segmentByPairedPSCBS.data.frame
# @alias segmentByPairedPSCBS
#
# @title "Segment total copy numbers and allele B fractions using the Paired PSCBS method"
#
# \description{
#  @get "title" [1].
#  This method requires matched normals.
#  This is a low-level segmentation method.
#  It is intended to be applied to one tumor-normal sample at the time.
# }
#
# @synopsis
#
# \arguments{
#   \item{CT}{A @numeric @vector of J tumor total tumor copy number (TCN) ratios in [0,+@Inf) (due to noise, small negative values are also allowed).  The TCN ratios are typically scaled such that copy-neutral diploid loci have a mean of two.}
#   \item{betaT}{A @numeric @vector of J tumor allele B fractions (BAFs) in [0,1] (due to noise, values may be slightly outside as well) or @NA for non-polymorphic loci.}
#   \item{betaN}{A @numeric @vector of J matched normal BAFs in [0,1] (due to noise, values may be slightly outside as well) or @NA for non-polymorphic loci.}
#   \item{muN}{An optional @numeric @vector of J genotype calls in 
#        \{0,1/2,1\} for AA, AB, and BB, respectively, 
#        and @NA for non-polymorphic loci.
#        If not given, they are estimated from the normal BAFs using
#        @see "aroma.light::callNaiveGenotypes" as described in [2].}
#   \item{chromosome}{(Optional) An @integer scalar (or a @vector of length J),
#        which can be used to specify which chromosome each locus belongs to
#        in case multiple chromosomes are segments.
#        This argument is also used for annotation purposes.}
#   \item{x}{Optional @numeric @vector of J genomic locations.
#            If @NULL, index locations \code{1:J} are used.}
#   \item{alphaTCN, alphaDH}{The significance levels for segmenting total
#        copy numbers (TCNs) and decrease-in-heterozygosity signals (DHs),
#        respectively.}
#   \item{undoTCN, undoDH}{Non-negative @numerics.  If less than +@Inf, 
#        then a cleanup of segmentions post segmentation is done.
#        See argument \code{undo} of @see "segmentByCBS" for more
#        details.}
#   \item{...}{Additional arguments passed to @see "segmentByCBS".}
#   \item{flavor}{A @character specifying what type of segmentation and 
#     calling algorithm to be used.}
#   \item{tbn}{If @TRUE, \code{betaT} is normalized before segmentation
#     using the TumorBoost method [2], otherwise not.}
#   \item{joinSegments}{If @TRUE, there are no gaps between neighboring
#     segments.
#     If @FALSE, the boundaries of a segment are defined by the support
#     that the loci in the segments provides, i.e. there exist a locus
#     at each end point of each segment.  This also means that there
#     is a gap between any neighboring segments, unless the change point
#     is in the middle of multiple loci with the same position.
#     The latter is what \code{DNAcopy::segment()} returns.
#   } 
#   \item{knownSegments}{Optional @data.frame specifying 
#     \emph{non-overlapping} known segments.  These segments must
#     not share loci.}
#   \item{seed}{An (optional) @integer specifying the random seed to be 
#     set before calling the segmentation method.  The random seed is
#     set to its original state when exiting.  If @NULL, it is not set.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#   Returns the segmentation results as a @see "PairedPSCBS" object.
# }
# 
# \details{
#   Internally @see "segmentByCBS" is used for segmentation.
#   The Paired PSCBS segmentation method does \emph{not} support weights.
# }
#
# \section{Reproducibility}{
#   The "DNAcopy::segment" implementation of CBS uses approximation
#   through random sampling for some estimates.  Because of this,
#   repeated calls using the same signals may result in slightly 
#   different results, unless the random seed is set/fixed.
# }
#
# \section{Whole-genome segmentation is preferred}{
#   Although it is possible to segment each chromosome independently
#   using Paired PSCBS, we strongly recommend to segment whole-genome
#   (TCN,BAF) data at once.  The reason for this is that downstream
#   CN-state calling methods, such as the AB and the LOH callers,
#   performs much better on whole-genome data.  In fact, they may
#   fail to provide valid calls if done chromsome by chromosome.
# }
#
# \section{Missing and non-finite values}{
#   The total copy number signals as well as any optional positions
#   must not contain missing values, i.e. @NAs or @NaNs.
#   If there are any, an informative error is thrown.
#   Allele B fractions may contain missing values, because such are
#   interpreted as representing non-polymorphic loci.
#
#   None of the input signals may have infinite values, i.e. -@Inf or +@Inf.
#   If so, an informative error is thrown.
# }
#
# @examples "../incl/segmentByPairedPSCBS.Rex"
#
# @author
#
# \references{
#  [1] @include "../incl/OlshenA_etal_2011.Rd" \cr
#  [2] @include "../incl/BengtssonH_etal_2010.Rd" \cr
# }
#
# \seealso{
#   Internally, @see "aroma.light::callNaiveGenotypes" is used to 
#   call naive genotypes, @see "aroma.light::normalizeTumorBoost" is 
#   used for TumorBoost normalization, and @see "segmentByCBS" is used 
#   to segment TCN and DH separately.
#
#   To segment total copy-numbers, or any other unimodal signals,
#   see @see "segmentByCBS".
# }
#
# @keyword IO
#*/########################################################################### 
setMethodS3("segmentByPairedPSCBS", "default", function(CT, betaT, betaN, muN=NULL, chromosome=0, x=NULL, alphaTCN=0.009, alphaDH=0.001, undoTCN=Inf, undoDH=Inf, ..., flavor=c("tcn&dh", "tcn,dh", "sqrt(tcn),dh", "sqrt(tcn)&dh"), tbn=TRUE, joinSegments=TRUE, knownSegments=NULL, seed=NULL, verbose=FALSE) {
  # WORKAROUND: If Hmisc is loaded after R.utils, it provides a buggy
  # capitalize() that overrides the one we want to use. Until PSCBS
  # gets a namespace, we do the following workaround. /HB 2011-07-14
  capitalize <- R.utils::capitalize;

  require("R.utils") || throw("Package not loaded: R.utils");
  require("aroma.light") || throw("Package not loaded: aroma.light");

  # To please R CMD check
  index <- NULL; rm(index);

  # Settings for sanity checks
  tol <- getOption("PSCBS/sanityChecks/tolerance", 0.0005);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'CT':
  disallow <- c("Inf");
  CT <- Arguments$getDoubles(CT, disallow=disallow);
  nbrOfLoci <- length(CT);
  length2 <- rep(nbrOfLoci, times=2);

  # Argument 'betaT':
  betaT <- Arguments$getDoubles(betaT, length=length2, disallow="Inf");

  # Argument 'muN':
  if (!is.null(muN)) {
    muN <- Arguments$getDoubles(muN, length=length2, range=c(0,1), disallow="Inf");
  }

  # Argument 'tbn':
  tbn <- Arguments$getLogical(tbn);

  # Argument 'betaN':
  betaN <- Arguments$getDoubles(betaN, length=length2, disallow="Inf");

  # Argument 'chromosome':
  if (is.null(chromosome)) {
    chromosome <- 0L; 
  } else {
    disallow <- c("Inf");
    chromosome <- Arguments$getIntegers(chromosome, range=c(0,Inf), disallow=disallow);
    if (length(chromosome) > 1) {
      chromosome <- Arguments$getIntegers(chromosome, length=length2, disallow=disallow);
    }
  }

  # Argument 'x':
  if (is.null(x)) {
    x <- seq(length=nbrOfLoci);
  } else {
    disallow <- c("Inf");
    x <- Arguments$getDoubles(x, length=length2, disallow=disallow);
  }

  # Argument 'alphaTCN':
  alphaTCN <- Arguments$getDouble(alphaTCN, range=c(0,1));

  # Argument 'alphaDH':
  alphaDH <- Arguments$getDouble(alphaDH, range=c(0,1));

  # Argument 'undoTCN':
  undoTCN <- Arguments$getDouble(undoTCN, range=c(0,Inf));

  # Argument 'undoDH':
  undoDH <- Arguments$getDouble(undoDH, range=c(0,Inf));


  # Argument 'flavor':
  flavor <- match.arg(flavor);
  knownFlavors <- c("tcn,dh", "tcn&dh", "sqrt(tcn),dh", "sqrt(tcn)&dh");
  if (!is.element(flavor, knownFlavors)) {
    throw("Segmentation flavor is not among the supported ones (", paste(sprintf("\"%s\"", knownFlavors), collapse=", "), "): ", flavor);
  }

  # Argument 'joinSegments':
  joinSegments <- Arguments$getLogical(joinSegments);

  # Argument 'knownSegments':
  if (is.null(knownSegments)) {
    knownSegments <- data.frame(chromosome=integer(0), start=integer(0), end=integer(0));
  } else {
    if (!joinSegments) {
##      warning("Argument 'knownSegments' should only be specified if argument 'joinSegments' is TRUE.");
    }
  }

  if (!is.data.frame(knownSegments)) {
    throw("Argument 'knownSegments' is not a data.frame: ", class(knownSegments)[1]);
  }

  if (!all(is.element(c("chromosome", "start", "end"), colnames(knownSegments)))) {
    throw("Argument 'knownSegments' does not have the required column names: ", hpaste(colnames(knownSegments)));
  }


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


  verbose && enter(verbose, "Segmenting paired tumor-normal signals using Paired PSCBS");

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


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Call genotypes
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # If muN is missing, call genotypes from betaN
  if (is.null(muN)) {
    verbose && enter(verbose, "Calling genotypes from normal allele B fractions");
    verbose && str(verbose, betaN);
    muN <- aroma.light::callNaiveGenotypes(betaN, censorAt=c(0,1));
    verbose && cat(verbose, "Called genotypes:");
    verbose && str(verbose, muN);
    verbose && print(verbose, table(muN));
    verbose && exit(verbose);
  }



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Normalize betaT using betaN (TumorBoost normalization)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (tbn) {
    verbose && enter(verbose, "Normalizing betaT using betaN (TumorBoost)");
    betaTN <- aroma.light::normalizeTumorBoost(betaT=betaT, betaN=betaN, muN=muN, preserveScale=TRUE);
    verbose && cat(verbose, "Normalized BAFs:");
    verbose && str(verbose, betaTN);

    # Assert that no missing values where introduced
    keep <- (is.finite(betaT) & is.finite(betaN) & is.finite(muN));
    if (any(is.na(betaTN[keep]))) {
      throw("Internal error: normalizeTumorBoost() introduced missing values.");
    }
    rm(keep);
    verbose && exit(verbose);
  } else {
    betaTN <- betaT;
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Setup up data");
  data <- data.frame(chromosome=chromosome, x=x, CT=CT, betaT=betaT, betaTN=betaTN, betaN=betaN, muN=muN);
  verbose && str(verbose, data);
  rm(chromosome, x, CT, betaT, betaTN, betaN, muN); # Not needed anymore
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Drop data points without known genomic positions, because that
  # is what DNAcopy::CNA() will do otherwise.  At the end, we will 
  # undo this such that the returned 'data' object is complete.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ok <- (!is.na(data$chromosome) & !is.na(data$x));
  if (any(!ok)) {
    verbose && enter(verbose, "Dropping loci with unknown locations");
    verbose && cat(verbose, "Number of loci dropped: ", sum(!ok));
    data <- data[ok,,drop=FALSE];
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
  o <- order(data$chromosome, data$x, decreasing=FALSE, na.last=TRUE);
  # Any change?
  if (any(o != seq(along=o))) {
    data <- data[o,,drop=FALSE];
  }
  rm(o); # Not needed anymore
  verbose && str(verbose, data);
  verbose && exit(verbose);

  # Attach 'index' (guaranteed to be ordered)
  data$index <- seq(length=nrow(data));


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Multiple chromosomes?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify all chromosomes, excluding missing values
  chromosomes <- sort(unique(data$chromosome), na.last=NA);
  nbrOfChromosomes <- length(chromosomes);
  if (nbrOfChromosomes > 1) {
    verbose && enter(verbose, "Segmenting multiple chromosomes");
    verbose && cat(verbose, "Number of chromosomes: ", nbrOfChromosomes);

    fitList <- list();
    for (kk in seq(length=nbrOfChromosomes)) {
      chromosomeKK <- chromosomes[kk];
      chrTag <- sprintf("Chr%02d", chromosomeKK);
      verbose && enter(verbose, sprintf("Chromosome #%d ('%s') of %d", kk, chrTag, nbrOfChromosomes));

      # Extract subset of data and parameters for this chromosome
      dataKK <- subset(data, chromosome == chromosomeKK);
      verbose && str(verbose, dataKK);
      fields <- attachLocally(dataKK, fields=c("CT", "betaT", "betaTN", "betaN", "muN", "chromosome", "x"));
      rm(dataKK); # Not needed anymore

      knownSegmentsKK <- NULL;
      if (!is.null(knownSegments)) {
        knownSegmentsKK <- subset(knownSegments, chromosome == chromosomeKK);
        verbose && cat(verbose, "Known segments:");
        verbose && print(verbose, knownSegmentsKK);
      }

      fit <- segmentByPairedPSCBS(CT=CT, betaT=betaTN, 
                betaN=betaN, muN=muN, 
                chromosome=chromosome, x=x,
                tbn=FALSE, joinSegments=joinSegments,
                knownSegments=knownSegmentsKK,
                alphaTCN=alphaTCN, alphaDH=alphaDH,
                undoTCN=undoTCN, undoDH=undoDH,
                flavor=flavor,
                ..., verbose=verbose);

      # Sanity checks
      if (nrow(knownSegmentsKK) == 0) {
        stopifnot(nrow(fit$data) == length(CT));
        stopifnot(all.equal(fit$data$CT, CT));
        stopifnot(all.equal(fit$data$muN, muN));
      }

      # Updata betaT (which is otherwise equals betaTN)
      fit$data$betaT <- betaT;

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
  # Subset 'knownSegments'
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Keeping only current chromosome for 'knownSegments'");

  currChromosome <- chromosome[1];
  verbose && cat(verbose, "Chromosome: ", currChromosome);

  knownSegments <- subset(knownSegments, chromosome == currChromosome);
  nbrOfSegments <- nrow(knownSegments);

  verbose && cat(verbose, "Known segments for this chromosome:");
  verbose && print(verbose, knownSegments);

  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Sanity checks
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Here 'knownSegments' should specify at most a single chromosome
  uChromosomes <- sort(unique(knownSegments$chromosome));
  if (length(uChromosomes) > 1) {
    throw("INTERNAL ERROR: Argument 'knownSegments' specifies more than one chromosome: ", hpaste(uChromosomes));
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup input data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && cat(verbose, "alphaTCN: ", alphaTCN);
  verbose && cat(verbose, "alphaDH: ", alphaDH);
  verbose && cat(verbose, "Number of loci: ", nbrOfLoci);

  # SNPs are identifies as those loci that have non-missing 'betaTN' & 'muN'
  isSnp <- (!is.na(data$betaTN) & !is.na(data$muN));
  nbrOfSnps <- sum(isSnp);
  verbose && cat(verbose, "Number of SNPs: ", nbrOfSnps);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Calculate decrease-of-heterozygosity signals (DHs)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Calculating DHs");
  # DH is by definition only defined for heterozygous SNPs.  
  # For simplicity, we set it to be NA for non-heterozygous loci.
  isHet <- isSnp & (data$muN == 1/2);
  verbose && printf(verbose, "Number of heterozygous SNPs: %d (%.2f%%)\n",
                                     sum(isHet), 100*sum(isHet)/nbrOfSnps);
  naValue <- as.double(NA);
  rho <- rep(naValue, length=nbrOfLoci);
  rho[isHet] <- 2*abs(data$betaTN[isHet]-1/2);
  verbose && cat(verbose, "Normalized DHs:");
  verbose && str(verbose, rho);
  data$rho <- rho;
  rm(rho); # Not needed anymore
  verbose && exit(verbose);



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # 1a. Identification of change points in total copy numbers
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Identification of change points by total copy numbers");

  fields <- attachLocally(data, fields=c("CT", "chromosome", "x", "index"));

  # Physical positions of loci
  fit <- segmentByCBS(CT, 
                      chromosome=chromosome, x=x, index=index,
                      joinSegments=joinSegments,
                      knownSegments=knownSegments,
                      alpha=alphaTCN, undo=undoTCN, ...,
                      verbose=verbose);
  verbose && str(verbose, fit);

  rm(list=fields); # Not needed anymore

  # Sanity check
  if (nrow(knownSegments) == 0) {
    stopifnot(nrow(fit$data) == nrow(data));
    stopifnot(all(fit$data$chromosome == data$chromosome));
    stopifnot(all(fit$data$x == data$x));
    stopifnot(all(fit$data$index == data$index));
    stopifnot(all.equal(fit$data$y, data$CT));
  }

  tcnSegments <- fit$output;
  tcnSegRows <- fit$segRows;
  rm(fit);

  # Sanity checks
  stopifnot(all(tcnSegRows[,1] <= tcnSegRows[,2], na.rm=TRUE));
  stopifnot(all(tcnSegRows[-nrow(tcnSegRows),2] < tcnSegRows[-1,1], na.rm=TRUE));

  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # 1b. Restructure TCN segmentation results
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Restructure TCN segmentation results");
  # Drop dummy columns
  keep <- setdiff(colnames(tcnSegments), c("sampleName"));
  tcnSegments <- tcnSegments[,keep,drop=FALSE];

  # Tag fields by TCN
  names <- names(tcnSegments);
  # Adding 'tcn' prefix to column names
  names <- sprintf("tcn%s", capitalize(names));
  names <- gsub("tcnChromosome", "chromosome", names, fixed=TRUE);
  names(tcnSegments) <- names;
  rm(names);
  verbose && print(verbose, tcnSegments);

  nbrOfSegs <- nrow(tcnSegments);
  verbose && cat(verbose, "Number of TCN segments: ", nbrOfSegs);

  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # 2a. Identification of additional change points using DH
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # For each segment independently, segment decrease of heterozygousity (DH)
  # using CBS. By definition, only heterozygous SNPs are used.

  dhSegRows <- NULL;
  tcnSegsExpanded <- NULL;

  # For each TCN segment...
  segs <- vector("list", length=nbrOfSegs);
  for (kk in seq(length=nbrOfSegs)) {
    tcnId <- kk;

    xStart <- tcnSegments[kk,"tcnStart"];
    xEnd <- tcnSegments[kk,"tcnEnd"];
    regionTag <- sprintf("[%g,%g]", xStart, xEnd);
    verbose && enter(verbose, sprintf("Total CN segment #%d (%s) of %d", kk, regionTag, nbrOfSegs));

    nbrOfTCNLociKK <- tcnSegments[kk,"tcnNbrOrLoci"];
    verbose && cat(verbose, "Number of TCN loci in segment: ", nbrOfTCNLociKK);

    rowStart <- tcnSegRows[kk,1];
    rowEnd <- tcnSegRows[kk,2];

    # Empty segment or a segment separator?
    isSplitter <- (is.na(rowStart) && is.na(rowEnd));
    if (isSplitter) {
      verbose && cat(verbose, "No signals to segment. Just a \"splitter\" segment. Skipping.");

      # Sanity check
      stopifnot(kk > 1);

      # Add a splitter segment
      segT <- segs[[kk-1]];
      segT <- segT[as.integer(NA),];
      keys <- colnames(tcnSegments);
      segT[,keys] <- tcnSegments[kk,keys];
      segT[,"tcnId"] <- tcnId;
      segT[,"dhId"] <- 1L;
      segT[,c("tcnNbrOfSNPs", "tcnNbrOfHets", "dhNbrOfLoci")] <- 0L;
      segT[,"dhStart"] <- xStart;
      segT[,"dhEnd"] <- xEnd;
      segs[[kk]] <- segT;
      verbose && print(verbose, segT);

      # Add a splitter to TCN and DH segment row matrix
      segRowsT <- dhSegRows[as.integer(NA),];
      dhSegRows <- rbind(dhSegRows, segRowsT);

      segRowsT <- tcnSegsExpanded[as.integer(NA),];
      tcnSegsExpanded <- rbind(tcnSegsExpanded, segRowsT);

      verbose && exit(verbose);
      next;
    }

    # Extract locus data for TCN segment
    rows <- rowStart:rowEnd;
    dataKK <- data[rows,,drop=FALSE];
    nbrOfLociKK <- nrow(dataKK);

    # Sanity check
    stopifnot(sum(!is.na(dataKK$CT)) == nbrOfTCNLociKK);
    gammaT <- tcnSegments[kk,"tcnMean"];

##    if (nrow(knownSegments) == 0) {
##      verbose && print(verbose, all.equal(mean(dataKK$CT, na.rm=TRUE), gammaT, tolerance=tol));
##      stopifnot(all.equal(mean(dataKK$CT, na.rm=TRUE), gammaT, tolerance=tol));
##    }

    verbose && cat(verbose, "Locus data for TCN segment:");
    verbose && str(verbose, dataKK);

    verbose && cat(verbose, "Number of loci: ", nbrOfLociKK);
    nbrOfSnpsKK <- sum(!is.na(dataKK$muN));
    verbose && printf(verbose, "Number of SNPs: %d (%.2f%%)\n", 
                                  nbrOfSnpsKK, 100*nbrOfSnpsKK/nbrOfLociKK);
    nbrOfHetsKK <- sum(!is.na(dataKK$muN) & dataKK$muN == 1/2);
    verbose && printf(verbose, "Number of heterozygous SNPs: %d (%.2f%%)\n",
                                  nbrOfHetsKK, 100*nbrOfHetsKK/nbrOfSnpsKK);


    verbose && enter(verbose, "Segmenting DH signals");
    fields <- attachLocally(dataKK, fields=c("chromosome", "x", "rho", "index"));

    # Since segments in 'knownSegments' has already been used in the TCN
    # segmentation, they are not needed in the DH segmentation.
    currChromosome <- chromosome[1];
    verbose && cat(verbose, "Chromosome: ", currChromosome);
    knownSegmentsT <- data.frame(chromosome=currChromosome, start=xStart, end=xEnd);

    fit <- segmentByCBS(rho, 
                        chromosome=chromosome, x=x,
                        joinSegments=joinSegments,
                        knownSegments=knownSegmentsT,
                        alpha=alphaDH, undo=undoDH, ...,
                        verbose=verbose);
    verbose && str(verbose, fit);
    dhSegments <- fit$output;
    dhSegRowsKK <- fit$segRows;

    verbose && cat(verbose, "DH segmentation (locally-indexed) rows:");
    verbose && print(verbose, dhSegRowsKK);
    verbose && str(verbose, index);

    # Remap to genome-wide indices
    for (cc in 1:2) {
      dhSegRowsKK[,cc] <- index[dhSegRowsKK[,cc]];
    }

    verbose && cat(verbose, "DH segmentation rows:");
    verbose && print(verbose, dhSegRowsKK);

    # Not needed anymore
    rm(list=fields);
    rm(fit);
    verbose && exit(verbose);

    # Drop dummy columns
    keep <- setdiff(colnames(dhSegments), c("sampleName", "chromosome"));
    dhSegments <- dhSegments[,keep,drop=FALSE];

    # Tag fields by DH
    names <- names(dhSegments);
    # Adding 'dh' prefix to column names
    names <- sprintf("dh%s", capitalize(names));
    names(dhSegments) <- names;
    rm(names);

    # Special case: If there where not enough data to segment DH...
    if (nrow(dhSegments) == 0) {
      dhSegments <- dhSegments[as.integer(NA),,drop=FALSE];
      dhSegRowsKK <- dhSegRowsKK[as.integer(NA),,drop=FALSE];
    }

    verbose && cat(verbose, "DH segmentation table:");
    verbose && print(verbose, dhSegments);
    verbose && print(verbose, dhSegRowsKK);

    # Expand the TCN segmentation result data frame
    rows <- rep(kk, times=nrow(dhSegments));
    verbose && cat(verbose, "Rows:");
    verbose && print(verbose, rows);
    tcnSegmentsKK <- tcnSegments[rows,,drop=FALSE];
    tcnSegRowsKK <- tcnSegRows[rows,,drop=FALSE];
    # Sanity check
    stopifnot(nrow(tcnSegmentsKK) == nrow(dhSegments));
    stopifnot(nrow(tcnSegRowsKK) == nrow(dhSegments));
    stopifnot(is.na(tcnSegRowsKK[,1]) || is.na(dhSegRowsKK[,1]) || (tcnSegRowsKK[,1] <= dhSegRowsKK[,1]));
    stopifnot(is.na(tcnSegRowsKK[,2]) || is.na(dhSegRowsKK[,2]) || (dhSegRowsKK[,2] <= tcnSegRowsKK[,2]));
    verbose && cat(verbose, "TCN segmentation rows:");
    verbose && print(verbose, tcnSegRowsKK);
    stopifnot(all(tcnSegRowsKK[,1] == tcnSegRowsKK[1,1], na.rm=TRUE));
    stopifnot(all(tcnSegRowsKK[,2] == tcnSegRowsKK[1,2], na.rm=TRUE));

    verbose && cat(verbose, "TCN and DH segmentation rows:");
    verbose && print(verbose, tcnSegRowsKK);
    verbose && print(verbose, dhSegRowsKK);
    verbose && print(verbose, tcnSegsExpanded);

    # Append
    tcnSegsExpanded <- rbind(tcnSegsExpanded, tcnSegRowsKK);
    verbose && cat(verbose, "TCN segmentation (expanded) rows:");
    verbose && print(verbose, tcnSegsExpanded);
    rownames(tcnSegsExpanded) <- NULL;

    dhSegRows <- rbind(dhSegRows, dhSegRowsKK);
    rownames(dhSegRows) <- NULL;

    verbose && cat(verbose, "TCN and DH segmentation rows:");
    verbose && print(verbose, tcnSegRows);
    verbose && print(verbose, dhSegRows);
    verbose && print(verbose, tcnSegsExpanded);

    # Sanity checks
    stopifnot(all(tcnSegRows[,1] <= tcnSegRows[,2], na.rm=TRUE));
    stopifnot(all(tcnSegRows[-nrow(tcnSegRows),2] < tcnSegRows[-1,1], na.rm=TRUE));
    stopifnot(all(dhSegRows[,1] <= dhSegRows[,2], na.rm=TRUE));
    stopifnot(all(dhSegRows[-nrow(dhSegRows),2] < dhSegRows[-1,1], na.rm=TRUE));
    stopifnot(all(tcnSegsExpanded[,1] <= tcnSegsExpanded[,2], na.rm=TRUE));
    stopifnot(all(tcnSegsExpanded[,1] <= dhSegRows[,1], na.rm=TRUE));
    stopifnot(all(tcnSegsExpanded[,2] >= dhSegRows[,2], na.rm=TRUE));
##    if (!all(tcnSegsExpanded[-nrow(tcnSegsExpanded),2] < tcnSegsExpanded[-1,1], na.rm=TRUE)) {
##      stopifnot(all(tcnSegsExpanded[-nrow(tcnSegsExpanded),2] < tcnSegsExpanded[-1,1], na.rm=TRUE));
##    }


    # Sanity check
    stopifnot(nrow(dhSegRows) == nrow(tcnSegsExpanded));

    # Append information on number of SNPs and hets in CN region
    tcnSegmentsKK <- cbind(
      tcnSegmentsKK, 
      tcnNbrOfSNPs=nbrOfSnpsKK,
      tcnNbrOfHets=nbrOfHetsKK
    );
    verbose && cat(verbose, "Total CN segmentation table (expanded):");
    verbose && print(verbose, tcnSegmentsKK);

    # Sanity check
    stopifnot(nrow(tcnSegmentsKK) == nrow(dhSegments));

    # Combine TCN and DH segmentation results
    tcndhSegments <- cbind(
      tcnId=rep(kk, times=nrow(dhSegments)),
      dhId=seq(length=nrow(dhSegments)),
      tcnSegmentsKK,
      dhSegments
    );

    segs[[kk]] <- tcndhSegments;

    verbose && cat(verbose, "(TCN,DH) segmentation for one total CN segment:");
    verbose && print(verbose, segs[[kk]]);

    verbose && exit(verbose);    
  } # for (kk ...)

  segs <- Reduce(rbind, segs);
  rownames(segs) <- NULL;

  # Sanity check
  stopifnot(nrow(dhSegRows) == nrow(tcnSegsExpanded));
  rownames(tcnSegRows) <- rownames(dhSegRows) <- NULL;

  stopifnot(all(tcnSegRows[,1] <= tcnSegRows[,2], na.rm=TRUE));
  stopifnot(all(tcnSegRows[-nrow(tcnSegRows),2] < tcnSegRows[-1,1], na.rm=TRUE));
  stopifnot(all(dhSegRows[,1] <= dhSegRows[,2], na.rm=TRUE));
  stopifnot(all(dhSegRows[-nrow(dhSegRows),2] < dhSegRows[-1,1], na.rm=TRUE));
  stopifnot(all(tcnSegsExpanded[,1] <= tcnSegsExpanded[,2], na.rm=TRUE));
##  stopifnot(all(tcnSegsExpanded[-nrow(tcnSegsExpanded),2] < tcnSegsExpanded[-1,1], na.rm=TRUE));

  # Move 'chromosome' column to the first column
  idx <- match("chromosome", names(segs));
  idxs <- c(idx, seq(length=ncol(segs))[-idx]);
  segs <- segs[,idxs,drop=FALSE];
  verbose && print(verbose, segs);

  verbose && enter(verbose, "Calculating (C1,C2) per segment");
  # Append (C1,C2) estimates
  tcn <- segs$tcnMean;
  dh <- segs$dhMean;
  C1 <- 1/2*(1-dh)*tcn;
  C2 <- tcn - C1;
  segs <- cbind(segs, c1Mean=C1, c2Mean=C2);
  verbose && exit(verbose);

  nbrOfSegs <- nrow(segs);
  verbose && cat(verbose, "Number of segments: ", nbrOfSegs);

  verbose && exit(verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Create result object
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  params <- list(
    alphaTCN = alphaTCN,
    alphaDH = alphaDH,
    flavor = flavor,
    tbn = tbn,
    joinSegments = joinSegments,
    knownSegments = knownSegments,
    seed = seed
  );

  # Should we drop attributes? /HB 2010-09-24
  stopifnot(all(data$index == seq(length=nrow(data))));
  data$index <- NULL; # Drop, because it is guaranteed to be ordered
  class(data) <- c("PairedPSCNData", class(data));

  class(segs) <- c("PairedPSCNSegments", class(segs));

  fit <- list(
    data = data,
    output = segs,
    tcnSegRows = tcnSegsExpanded,
    dhSegRows = dhSegRows,
    params = params
  );

  class(fit) <- c("PairedPSCBS", "PSCBS", "AbstractCBS");

  # Update 
  if (is.element(flavor, c("tcn&dh", "sqrt(tcn)&dh"))) {
    fit$params$flavor <- gsub("&", ",", flavor, fixed=TRUE); # AD HOC.
    fit <- postsegmentTCN(fit, verbose=verbose);

    # Sanity check
    CT <- fit$data$CT;
    tcnSegRows <- fit$tcnSegRows;
    dhSegRows <- fit$dhSegRows;
    for (jj in 1:nrow(tcnSegRows)) {
      tcnSegRowJJ <- unlist(tcnSegRows[jj,,drop=TRUE], use.names=FALSE);
      dhSegRowJJ <- unlist(dhSegRows[jj,,drop=TRUE], use.names=FALSE);
      stopifnot(
        is.na(tcnSegRowJJ[1]) || is.na(dhSegRowJJ[1]) || 
        # A TCN segment must start at or before a DH segment...
        (tcnSegRowJJ[1] <= dhSegRowJJ[1]) ||
        # ...unless there was an outlier at the left edge.
        (is.na(CT[dhSegRowJJ[1]]) && (tcnSegRowJJ[1] - 1L <= dhSegRowJJ[1]))
      );
      stopifnot(
        is.na(tcnSegRowJJ[2]) || is.na(dhSegRowJJ[2]) ||
        # A TCN segment must end at or after a DH segment...
        (dhSegRowJJ[2] <= tcnSegRowJJ[2]) ||
        # ...unless there was an outlier at the right edge.
        (is.na(CT[dhSegRowJJ[2]]) && (dhSegRowJJ[2] <= tcnSegRowJJ[2] + 1L))
      );
    } # for (jj ...)
    rm(CT, tcnSegRows, dhSegRows);  # Clean up
  }

  verbose && print(verbose, head(as.data.frame(fit)));
  verbose && print(verbose, tail(as.data.frame(fit)));

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Return results
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  fit;
}) # segmentByPairedPSCBS()



setMethodS3("segmentByPairedPSCBS", "data.frame", function(CT, ...) {
  # To please R CMD check
  data <- CT;


  segmentByPairedPSCBS(CT=data$CT, betaT=data$betaT, betaN=data$betaN, 
                   muN=data$muN, chromosome=data$chromosome, x=data$x, ...);
})


############################################################################
# HISTORY:
# 2011-10-21
# o Now segmentByPairedCBS() handles forced separators in 'knownSegments'.
# 2011-10-20
# o CLEANUP: Dropped a stray debug output message in segmentByPairedPSCBS().
# o Replaced argument 'knownCPs' with 'knownSegments' for  segmentByCBS().
# 2011-10-02
# o Added segmentByPairedPSCBS() for data.frame such that the locus-level
#   data arguments can also be passed via a data.frame. 
# 2011-09-04
# o ROBUSTNESS: Added drop=FALSE to matrix subsettings.
# o CLEANUP: Removed all references to/usage of DNAcopy fields, which
#   are no longer part of the CBS class.
# 2011-09-03
# o Updated code to not use deprecated argument 'columnNamesFlavor'
#   of segmentByCBS().
# 2011-08-08
# o BUG FIX: If dropSegmentationOutliers() would drop an outlier next to
#   a change point, such that total copy-number signal would become NA,
#   then the sanity checks that TCN segments always overlaps DH segments
#   would fail.  Now the sanity checks are aware of this special case.
# o Moved the sanity checks that tests the TCN and DH "segRows" from the
#   bootstrapTCNandDHByRegion() to segmentByPairedPSCBS().  This is the
#   first step to fix a failure in the sanity checks that could happend
#   iff one first run dropSegmentationOutliers().
# 2011-07-15
# o DOCUMENTATION: Added a section to help("segmentByPairedPSCBS") on
#   the importance of doing a whole-genome PSCBS segmentations if 
#   calling AB and LOH states afterward.
# 2011-07-14
# o DOCUMENTATION: Added to the help that arguments betaT, betaN and muN
#   may contain NAs for non-polymorphic loci.
# o BUG FIX/ROBUSTNESS: In some cases, the segmentation table would 
#   contain column names with incorrect capitalization, e.g. "tcnnbrOfLoci"
#   instead of "tcnNbrOfLoci".  This would cause several downstream 
#   methods to give an error.  The reason for this is that the Hmisc
#   package, if loaded after R.utils, overrides capitalize() in R.utils
#   with another (buggy?) capitalize() function.  To avoid this, we
#   now everywhere specify explicitly that we want to the one in R.utils.
# 2011-07-06
# o DOCUMENTATION: The description of argument 'chromosome' for 
#   segmentByPairedPSCBS() did not describe how to segment multiple
#   chromosomes in one call.
# 2011-07-05
# o BUG FIX: Output fields 'tcnNbrOfSNPs'and 'tcnNbrOfHets' were mistakenly
#   labelled as 'tcnNbrOr...'.  Thanks Christine Ho at UC Berkeley for
#   reporting on this.
# 2011-06-28
# o DOCUMENTATION: Clarified that argument 'CT' should be tumor copy
#   number ratios relative to the normal.
# 2011-06-14
# o CONVENTION: Changed the column names of returned data frames. 
#   They now follow the camelCase naming convention and are shorter.
# 2011-05-29
# o Renamed options to reflect new package name.
# 2011-04-07
# o ROBUSTNESS: Added validation of the different 'tcnSegRows' and
#   'dhSegRows' calculations in segmentByPairedPSCBS().  This helps
#   us narrow down a bug in postsegmentTCN().
# 2010-12-09
# o BUG FIX: When there were multiple chromsomes processed by
#   segmentByPairedPSCBS(), then the returned data object would
#   contain 'betaT' identical to 'betaTN'.
# 2010-12-02
# o Now segmentByPairedPSCBS() uses option "psCBS/sanityChecks/tolerance".
# 2010-11-30
# o Now segmentByPairedPSCBS() returns data frames 'tcnLociToExclude'
#   and 'dhLociToExclude'.
# o BUG FIX: Argument 'flavor' of segmentByPairedPSCBS() would be ignored
#   if multiple chromsomomes were segmented.
# 2010-11-28
# o BUG FIX: Iff argument 'chromosome' to segmentByPairedPSCBS() was of
#   length greater than one and specified exactly one unique chromosome,
#   then exception "Number of elements in argument 'chromosome' should
#   be exactly 8712 not 86209 value(s)" would be thrown.
# 2010-11-27
# o BUG FIX: segmentByPairedPSCBS() would not accept missing values in
#   argument 'chromosome'.
# o Now arguments '...' of segmentByPairedPSCBS() are passed to
#   the two segmentByCBS() calls.
# 2010-11-22
# o BUG FIX: segmentByPairedPSCBS() would not subset the correct set of
#   DH signals if there were some missing values in TCN.
# 2010-11-21
# o Changed the default to flavor="tch&dh".
# o Added support for flavors "tcn&dh", which, contrary to "tcn,dh",
#   enforces TCN and DH to have the same change points.
# o Now segmentByPairedPSCBS() also returns minor and major copy numbers
#   for each segment.
# o Forgot to return arguments 'joinSegments' & 'knownCPs' in 'params'.
# 2010-11-20
# o Now it is possible to specify the boundaries of the regions to be
#   segmented as known change points via argument 'knownCPs'.
# o Added argument 'joinSegments' to segmentByPairedPSCBS() in order to 
#   specify if neighboring segments should be joined or not.
# o Now segmentByCBS() allows for unknown genomic positions.
# o Now segmentByCBS() allows also for missing total CN signals.
# 2010-11-16
# o BUG FIX: In the rare cases where two loci at the same positions are
#   split up into two neighboring segments, then segmentByPairedPSCBS()
#   would fail to infer which they were if and only if the loci were not
#   ordered along the genome.  This could happen with for instance
#   Affymetrix GenomeWideSNP_6 data.
# o DOCUMENTATION: Clarified the form of argument 'muN', and added
#   references to papers and cross links to more internal methods.
# 2010-11-04
# o BUG FIX: There was a stray/debug stop() statement left in  
#   segmentByPairedPSCBS() causing an "error" in the rare case 
#   when loci that have the same physical locations are split
#   into two different segments.
# 2010-11-02
# o Added arguments 'undoTCN' and 'undoDH' to segmentByPairedPSCBS().
# o BUG FIX: Arguments 'alphaTCN' and 'alphaDH' of segmentByPairedPSCBS() 
#   were not used when more than one chromosome were segmented.
# 2010-10-25
# o BUG FIX: Now the correct set of loci are extracted from each TCN
#   segment, in the rare case that two neighboring TCN segments have
#   the same end points.
# 2010-10-18
# o Added arguments 'alphaTCN' and 'alphaDH' to segmentByPairedPSCBS().
# o Now segmentByPairedPSCBS() can segment multiple chromosomes.
# 2010-10-17
# o Added argument 'tbn' to segmentByPairedPSCBS() specifying whether
#   TumorBoostNormalization should be applied or not.
# 2010-10-10
# o The default is now to segment TCN on the original scale, not the sqrt().
# o Added flavor "sqrt(tcn),dh", which is segments sqrt(TCN) and then DH,
#   as original proposed by ABO.
# 2010-10-03
# o CLEAN UP: Now segmentByPairedPSCBS() is making use of argument
#   'chromosome' of segmentByCBS().
# 2010-10-02
# o Argument 'chromosome' default to 0 and have to be a finite integer.
# 2010-09-24
# o Now the 'data' field returned is a data.frame (no longer a list).
# o Now the 'chromosome' field of the data field is expanded to have the
#   same number of elements as the other locus fields.
# 2010-09-18
# o Added argument 'chromosome' to segmentByPairedPSCBS(), which, if given,
#   adds a chromosome column to the data and segmentation results.
# 2010-09-08
# o Now segmentByPairedPSCBS() also returns the TumorBoost normalized data.
#   This also means that plot() for PairedPSCBS no longer has to 
#   recalculate them.
# 2010-09-04
# o Added drawLevels() for PairedPSCBS.
# o Added as.data.frame() and print() for PairedPSCBS.
# 2010-09-03
# o Added plot().
# 2010-07-09
# o The segmentByPairedPSCBS() method was written completely from scratch.
# o Created.
############################################################################
