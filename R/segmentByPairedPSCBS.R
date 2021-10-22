###########################################################################/**
# @RdocDefault segmentByPairedPSCBS
# @alias segmentByPairedPSCBS.data.frame
# @alias segmentByPairedPSCBS.PairedPSCBS
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
#   \item{CT}{A @numeric @vector of J tumor total copy number (TCN)
#        ratios in [0,+@Inf) (due to noise, small negative values are
#        also allowed).  The TCN ratios are typically scaled such that
#        copy-neutral diploid loci have a mean of two.}
#   \item{thetaT, thetaN}{(alternative) As an alternative to specifying
#        tumor TCN \emph{ratios} relative to the match normal by
#        argument \code{CT}, on may specify total tumor and normal
#        signals seperately, in which case the TCN ratios \code{CT} are
#        calculated as \eqn{CT = 2*thetaT/thetaN}.}
#   \item{betaT}{A @numeric @vector of J tumor allele B fractions (BAFs)
#        in [0,1] (due to noise, values may be slightly outside as well)
#        or @NA for non-polymorphic loci.}
#   \item{betaN}{A @numeric @vector of J matched normal BAFs in [0,1]
#        (due to noise, values may be slightly outside as well) or @NA
#        for non-polymorphic loci.}
#   \item{muN}{An optional @numeric @vector of J genotype calls in
#        \{0,1/2,1\} for AA, AB, and BB, respectively,
#        and @NA for non-polymorphic loci.
#        If not given, they are estimated from the normal BAFs using
#        @see "aroma.light::callNaiveGenotypes" as described in [2].}
#   \item{rho}{(alternative to \code{betaT} and \code{betaN}/\code{muN})
#        A @numeric @vector of J decrease-of-heterozygosity signals (DHs)
#        in [0,1] (due to noise, values may be slightly larger than one
#        as well).  By definition, DH should be @NA for homozygous loci
#        and for non-polymorphic loci.}
#   \item{chromosome}{(Optional) An @integer scalar (or a @vector of length J),
#        which can be used to specify which chromosome each locus belongs to
#        in case multiple chromosomes are segments.
#        This argument is also used for annotation purposes.}
#   \item{x}{Optional @numeric @vector of J genomic locations.
#            If @NULL, index locations \code{1:J} are used.}
#   \item{alphaTCN, alphaDH}{The significance levels for segmenting total
#        copy numbers (TCNs) and decrease-in-heterozygosity signals (DHs),
#        respectively.}
#   \item{undoTCN, undoDH}{Non-negative @numerics.  If greater than 0,
#        then a cleanup of segmentions post segmentation is done.
#        See argument \code{undo} of @see "segmentByCBS" for more
#        details.}
#   \item{avgTCN, avgDH}{A @character string specifying how to calculating
#         segment mean levels \emph{after} change points have been
#         identified.}
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
#     not share loci.  See @see "findLargeGaps" and @see "gapsToSegments".}
#   \item{dropMissingCT}{If @TRUE, loci for which 'CT' is missing
#     are dropped, otherwise not.}
#   \item{seed}{An (optional) @integer specifying the random seed to be
#     set before calling the segmentation method.  The random seed is
#     set to its original state when exiting.  If @NULL, it is not set.}
#   \item{verbose}{See @see "R.utils::Verbose".}
#   \item{preserveScale}{\emph{Defunct - gives an error is specified.}}
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
#   fail to provide valid calls if done chromosome by chromosome.
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
# \section{Paired PSCBS with only genotypes}{
#   If allele B fractions for the matched normal (\code{betaN}) are
#   not available, but genotypes (\code{muN}) are, then it is possible
#   to run a version of Paired PSCBS where TumorBoost normalization
#   of the tumor allele B fractions is skipped.  In order for this
#   to work, argument \code{tbn} must be set to @FALSE.
# }
#
# @examples "../incl/segmentByPairedPSCBS.Rex"
#
# @author "HB"
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
#   To segment tumor total copy numbers and allele B fractions
#   \emph{without} a matched normal, see @see "segmentByNonPairedPSCBS".
#
#   To segment total copy-numbers, or any other unimodal signals,
#   see @see "segmentByCBS".
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("segmentByPairedPSCBS", "default", function(CT, thetaT=NULL, thetaN=NULL, betaT=NULL, betaN=NULL, muN=NULL, rho=NULL, chromosome=0, x=NULL, alphaTCN=0.009, alphaDH=0.001, undoTCN=0, undoDH=0, ..., avgTCN=c("mean", "median"), avgDH=c("mean", "median"), flavor=c("tcn&dh", "tcn,dh", "sqrt(tcn),dh", "sqrt(tcn)&dh", "tcn"), tbn=is.null(rho), joinSegments=TRUE, knownSegments=NULL, dropMissingCT=TRUE, seed=NULL, verbose=FALSE, preserveScale=FALSE) {
  # WORKAROUND: If Hmisc is loaded after R.utils, it provides a buggy
  # capitalize() that overrides the one we want to use. Until PSCBS
  # gets a namespace, we do the following workaround. /HB 2011-07-14
  capitalize <- R.utils::capitalize

  # To please R CMD check
  index <- NULL; rm(list="index")

  # Settings for sanity checks
  tol <- getOption("PSCBS/sanityChecks/tolerance", 0.0005)

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'thetaT' & 'thetaN':
  if (!is.null(thetaT) && !is.null(thetaN)) {
    thetaT <- Arguments$getDoubles(thetaT, disallow=disallow)
    nbrOfLoci <- length(thetaT)
    length2 <- rep(nbrOfLoci, times=2L)
    thetaN <- Arguments$getDoubles(thetaN, length=length2, disallow=disallow)
    CT <- 2 * thetaT / thetaN
  } else if (!is.null(thetaT) || !is.null(thetaN)) {
    stop("Either argument 'CT' needs to be specified or *both* of arguments 'thetaT' and 'thetaN'")
  }

  # Argument 'CT':
  disallow <- c("Inf")
  CT <- Arguments$getDoubles(CT, disallow=disallow)
  nbrOfLoci <- length(CT)
  length2 <- rep(nbrOfLoci, times=2L)


  # Argument 'betaT':
  if (!is.null(betaT)) {
    betaT <- Arguments$getDoubles(betaT, length=length2, disallow="Inf")
  }

  # Argument 'betaN':
  if (!is.null(betaN)) {
    betaN <- Arguments$getDoubles(betaN, length=length2, disallow="Inf")
  }

  # Argument 'muN':
  if (!is.null(muN)) {
    muN <- Arguments$getDoubles(muN, length=length2, range=c(0,1), disallow="Inf")
    if (all(is.na(muN)) == nbrOfLoci) {
      stop(sprintf("All genotypes ('muN') are NAs: %d (100%%) out of %d", nbrOfLoci, nbrOfLoci))
    }
  }

  # Argument 'rho':
  if (!is.null(rho)) {
    rho <- Arguments$getDoubles(rho, range=c(0,Inf), length=length2, disallow="Inf")
  }

  if (is.null(muN)) {
    if (is.null(betaN) && is.null(rho)) {
      stop("If argument 'muN' is not given, then either 'betaN' or 'rho' must be.")
    }
  }

  # Argument 'tbn':
  tbn <- Arguments$getLogical(tbn)
  if (!is.null(tbn)) {
    if (tbn) {
      if (is.null(betaT)) {
        stop("Cannot do TumorBoost normalization (tbn=TRUE) without tumor BAFs ('betaT').")
      }
      if (is.null(betaN)) {
        stop("Cannot do TumorBoost normalization (tbn=TRUE) with normal BAFs ('betaN').")
      }
    }
  }

  # Argument 'chromosome':
  if (is.null(chromosome)) {
    chromosome <- 0L
  } else {
    disallow <- c("Inf")
    chromosome <- Arguments$getIntegers(chromosome, range=c(0,Inf), disallow=disallow)
    if (length(chromosome) > 1) {
      chromosome <- Arguments$getIntegers(chromosome, length=length2, disallow=disallow)
    }
  }

  # Argument 'x':
  if (is.null(x)) {
    x <- seq_len(nbrOfLoci)
  } else {
    disallow <- c("Inf")
    x <- Arguments$getDoubles(x, length=length2, disallow=disallow)
  }

  # Argument 'alphaTCN':
  alphaTCN <- Arguments$getDouble(alphaTCN, range=c(0,1))

  # Argument 'alphaDH':
  alphaDH <- Arguments$getDouble(alphaDH, range=c(0,1))

  # Argument 'undoTCN':
  undoTCN <- Arguments$getDouble(undoTCN, range=c(0,Inf))

  # Argument 'undoDH':
  undoDH <- Arguments$getDouble(undoDH, range=c(0,Inf))

  # Argument 'avgTCN' & 'avgDH':
  avgTCN <- match.arg(avgTCN)
  avgDH <- match.arg(avgDH)

  # Argument 'flavor':
  flavor <- match.arg(flavor)
  knownFlavors <- eval(formals(segmentByPairedPSCBS.default)$flavor, enclos = baseenv())
  if (!is.element(flavor, knownFlavors)) {
    stop("Segmentation flavor is not among the supported ones (", paste(sprintf("\"%s\"", knownFlavors), collapse=", "), "): ", flavor)
  }

  # Argument 'joinSegments':
  joinSegments <- Arguments$getLogical(joinSegments)

  # Argument 'knownSegments':
  if (is.null(knownSegments)) {
    knownSegments <- data.frame(chromosome=integer(0), start=integer(0), end=integer(0))
  } else {
    if (!joinSegments) {
##      warning("Argument 'knownSegments' should only be specified if argument 'joinSegments' is TRUE.")
    }
  }

  if (!is.data.frame(knownSegments)) {
    stop("Argument 'knownSegments' is not a data.frame: ", class(knownSegments)[1])
  }

  if (!all(is.element(c("chromosome", "start", "end"), colnames(knownSegments)))) {
    stop("Argument 'knownSegments' does not have the required column names: ", hpaste(colnames(knownSegments)))
  }

  # Argument 'dropMissingCT':
  dropMissingCT <- Arguments$getLogical(dropMissingCT)
  if (!dropMissingCT) {
    if (is.element(flavor, c("tcn&dh", "sqrt(tcn)&dh"))) {
      stop("Missing values in 'CT' are (currently) not supported by the chosen 'flavor': ", flavor)
    }
  }


  # Argument 'seed':
  if (!is.null(seed)) {
    seed <- Arguments$getIntegers(seed)
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
    pushState(verbose)
    on.exit(popState(verbose))
  }

  # Argument 'preserveScale' is deprecated
  if (!missing(preserveScale)) {
    .Defunct(msg = "Argument 'preserveScale' for segmentByPairedPSCBS() is defunct; as of PSCBS 0.64.0 (Mar 2018) it is effectively fixed to FALSE, which has been the default since PSCBS 0.50.0 (Oct 2015). To avoid this error, do not specify 'preserveScale' when calling segmentByPairedPSCBS().")
  }
  
  verbose && enter(verbose, "Segmenting paired tumor-normal signals using Paired PSCBS")

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Call genotypes?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Are genotype calls muN missing and can they be called?
  if (is.null(muN) && !is.null(betaN)) {
    verbose && enter(verbose, "Calling genotypes from normal allele B fractions")
    verbose && str(verbose, betaN)
    muN <- callNaiveGenotypes(betaN, censorAt=c(0,1))
    verbose && cat(verbose, "Called genotypes:")
    verbose && str(verbose, muN)
    verbose && print(verbose, table(muN))
    # Assert proper calls
    muN <- Arguments$getDoubles(muN, length=length2, range=c(0,1), disallow="Inf")
    # Sanity check
    if (all(is.na(muN))) {
      stop(sprintf("All genotypes ('muN') called from the normal allele B fractions ('betaN') are NAs: %d (100%%) out of %d", nbrOfLoci, nbrOfLoci))
    }
    verbose && exit(verbose)
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Normalize betaT using betaN (TumorBoost normalization)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (tbn) {
    verbose && enter(verbose, "Normalizing betaT using betaN (TumorBoost)")
    betaTN <- normalizeTumorBoost(betaT=betaT, betaN=betaN, muN=muN, preserveScale=FALSE)
    verbose && cat(verbose, "Normalized BAFs:")
    verbose && str(verbose, betaTN)

    # Assert that no missing values where introduced
    keep <- (is.finite(betaT) & is.finite(betaN) & is.finite(muN))
    if (anyNA(betaTN[keep])) {
      stop("Internal error: normalizeTumorBoost() introduced missing values.")
    }
    # Not needed anymore
    keep <- NULL
    verbose && exit(verbose)
  } else {
    betaTN <- betaT
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Setup up data")
  data <- data.frame(chromosome=chromosome, x=x, CT=CT)
  if (!is.null(thetaT)) {
    data$thetaT <- thetaT
    data$thetaN <- thetaN
  }
  if (!is.null(betaT)) data$betaT <- betaT
  if (!is.null(betaTN)) data$betaTN <- betaTN
  if (!is.null(betaN)) data$betaN <- betaN
  if (!is.null(muN)) data$muN <- muN
  if (!is.null(rho)) data$rho <- rho
  verbose && str(verbose, data)
  # Not needed anymore
  chromosome <- x <- CT <- thetaT <- thetaN <- betaT <- betaTN <- betaN <- muN <- rho <- NULL

  # Sanity check
  .stop_if_not(nrow(data) == nbrOfLoci)
  verbose && exit(verbose)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Drop data points without known genomic positions, because that
  # is what DNAcopy::CNA() will do otherwise.  At the end, we will
  # undo this such that the returned 'data' object is complete.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ok <- (!is.na(data$chromosome) & !is.na(data$x))
  if (any(!ok)) {
    verbose && enter(verbose, "Dropping loci with unknown locations")
    verbose && cat(verbose, "Number of loci dropped: ", sum(!ok))
    data <- data[ok,,drop=FALSE]
    nbrOfLoci <- nrow(data)
    verbose && exit(verbose)
  }
  ok <- NULL # Not needed anymore

  # Sanity check
  .stop_if_not(nrow(data) == nbrOfLoci)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Drop loci for which CT is missing (regardless of betaT)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (dropMissingCT) {
    ok <- (!is.na(data$CT))
    if (any(!ok)) {
      verbose && enter(verbose, "Dropping loci for which TCNs are missing")
      verbose && cat(verbose, "Number of loci dropped: ", sum(!ok))
      data <- data[ok,,drop=FALSE]
      nbrOfLoci <- nrow(data)
      verbose && exit(verbose)
    }
    ok <- NULL # Not needed anymore

    # Sanity check
    .stop_if_not(nrow(data) == nbrOfLoci)
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Reorder data points along the genome, because that is what
  # DNAcopy::segment() will return.  At the end, we will undo
  # the sort such that the returned 'data' object is always in
  # the same order and number of loci as the input data.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Ordering data along genome")
  o <- order(data$chromosome, data$x, decreasing=FALSE, na.last=TRUE)
  # Any change?
  if (any(o != seq_along(o))) {
    data <- data[o,,drop=FALSE]
  }
  o <- NULL # Not needed anymore
  verbose && str(verbose, data)
  verbose && exit(verbose)

  # Attach 'index' (guaranteed to be ordered)
  data$index <- seq_len(nrow(data))

  # Sanity check
  .stop_if_not(nrow(data) == nbrOfLoci)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Assert no missing values in (chromosome, x, CT)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Sanity check
  ok <- (!is.na(data$chromosome) & !is.na(data$x))
  if (!all(ok)) {
    stop("INTERNAL ERROR: Detected (chromosome, x) with missing values also after filtering.")
  }

  # Sanity check
  if (dropMissingCT) {
    ok <- (!is.na(data$CT))
    if (!all(ok)) {
      stop("INTERNAL ERROR: Detected TCN with missing values also after filtering.")
    }
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Multiple chromosomes?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify all chromosomes, excluding missing values
  chromosomes <- sort(unique(data$chromosome), na.last=NA)
  nbrOfChromosomes <- length(chromosomes)
  if (nbrOfChromosomes > 1) {
    verbose && enter(verbose, "Segmenting multiple chromosomes")
    verbose && cat(verbose, "Number of chromosomes: ", nbrOfChromosomes)

    # Generate random seeds?
    seeds <- NULL
    if (!is.null(seed)) {
      randomSeed("set", seed=seed, kind="L'Ecuyer-CMRG")
      on.exit(randomSeed("reset"), add=TRUE)
      verbose && printf(verbose, "Random seed temporarily set (seed=c(%s), kind=\"L'Ecuyer-CMRG\")\n", paste(seed, collapse=", "))
      seeds <- randomSeed("advance", n=nbrOfChromosomes)
      verbose && printf(verbose, "Produced %d seeds from this stream for future usage\n", length(seeds))
    }

    fitList <- listenv()
    for (kk in seq_len(nbrOfChromosomes)) {
      chromosomeKK <- chromosomes[kk]
      chrTag <- sprintf("Chr%02d", chromosomeKK)
      verbose && enter(verbose, sprintf("Chromosome #%d ('%s') of %d", kk, chrTag, nbrOfChromosomes))

      seedKK <- seeds[[kk]]

      # Extract subset of data and parameters for this chromosome
      dataKK <- subset(data, chromosome == chromosomeKK)
      verbose && str(verbose, dataKK)
      fields <- attachLocally(dataKK, fields=c("CT", "thetaT", "thetaN", "betaT", "betaTN", "betaN", "muN", "rho", "chromosome", "x"))
      dataKK <- NULL # Not needed anymore

      knownSegmentsKK <- NULL
      if (!is.null(knownSegments)) {
        knownSegmentsKK <- subset(knownSegments, chromosome == chromosomeKK)
        verbose && cat(verbose, "Known segments:")
        verbose && print(verbose, knownSegmentsKK)
      }

      fitList[[chrTag]] %<-% {
        fit <- segmentByPairedPSCBS(CT=CT, thetaT=thetaT, thetaN=thetaN,
                  betaT=betaTN, betaN=betaN, muN=muN, rho=rho,
                  chromosome=chromosome, x=x,
                  tbn=FALSE, joinSegments=joinSegments,
                  knownSegments=knownSegmentsKK,
                  alphaTCN=alphaTCN, alphaDH=alphaDH,
                  undoTCN=undoTCN, undoDH=undoDH,
                  avgTCN=avgTCN, avgDH=avgDH,
                  flavor=flavor,
                  ...,
                  seed=seedKK,
                  verbose=verbose)

        # Sanity checks
        if (nrow(knownSegmentsKK) == 0) {
          .stop_if_not(nrow(fit$data) == length(CT))
          .stop_if_not(all.equal(fit$data$CT, CT))
          .stop_if_not(all.equal(fit$data$muN, muN))
        }

        # Update betaT (which is otherwise equals betaTN)
        fit$data$betaT <- betaT

        verbose && print(verbose, head(as.data.frame(fit)))
        verbose && print(verbose, tail(as.data.frame(fit)))

        fit
      } %seed% TRUE %label% sprintf("segmentByPairedPSCBS-%s", chrTag)  ## fitList[[chrTag]] <- ...

      rm(list=fields) # Not needed anymore
      verbose && exit(verbose)
    } # for (kk ...)

    verbose && enter(verbose, "Merging (independently) segmented chromosome")
    fitList <- as.list(fitList)
    ## former Reduce() w/ append(..., addSplit = TRUE)
    fit <- do.call(c, args = c(fitList, addSplit = TRUE))
    fitList <- NULL # Not needed anymore
    verbose && str(verbose, fit)
    verbose && exit(verbose)

    # Update parameters that otherwise may be incorrect
    fit$params$tbn <- tbn
    fit$params$seed <- seed

    segs <- as.data.frame(fit)
    if (nrow(segs) < 6) {
      verbose && print(verbose, segs)
    } else {
      verbose && print(verbose, head(segs))
      verbose && print(verbose, tail(segs))
    }

    verbose && exit(verbose)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Return results
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    return(fit)
  } # if (nbrOfChromosomes > 1)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Subset 'knownSegments'
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Keeping only current chromosome for 'knownSegments'")

  currChromosome <- data$chromosome[1]
  verbose && cat(verbose, "Chromosome: ", currChromosome)

  knownSegments <- subset(knownSegments, chromosome == currChromosome)
  nbrOfSegments <- nrow(knownSegments)

  verbose && cat(verbose, "Known segments for this chromosome:")
  verbose && print(verbose, knownSegments)

  verbose && exit(verbose)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Sanity checks
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Here 'knownSegments' should specify at most a single chromosome
  uChromosomes <- sort(unique(knownSegments$chromosome))
  if (length(uChromosomes) > 1) {
    stop("INTERNAL ERROR: Argument 'knownSegments' specifies more than one chromosome: ", hpaste(uChromosomes))
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Assert no missing values in (chromosome, x, CT)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Sanity check
  ok <- (!is.na(data$chromosome) & !is.na(data$x))
  if (!all(ok)) {
    stop("INTERNAL ERROR: Detected (chromosome, x) with missing values also after filtering.")
  }

  # Sanity check
  if (dropMissingCT) {
    ok <- (!is.na(data$CT))
    if (!all(ok)) {
      stop("INTERNAL ERROR: Detected TCN with missing values also after filtering.")
    }
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup input data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && cat(verbose, "alphaTCN: ", alphaTCN)
  verbose && cat(verbose, "alphaDH: ", alphaDH)
  verbose && cat(verbose, "Number of loci: ", nbrOfLoci)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Calculate decrease-of-heterozygosity signals (DHs)?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (is.null(data$rho)) {
    verbose && enter(verbose, "Calculating DHs")
    # SNPs are identifies as those loci that have non-missing 'betaTN' & 'muN'
    isSnp <- (!is.na(data$betaTN) & !is.na(data$muN))
    nbrOfSnps <- sum(isSnp)
    verbose && cat(verbose, "Number of SNPs: ", nbrOfSnps)

    # DH is by definition only defined for heterozygous SNPs.
    # For simplicity, we set it to be NA for non-heterozygous loci.
    isHet <- isSnp & (data$muN == 1/2)
    verbose && printf(verbose, "Number of heterozygous SNPs: %d (%.2f%%)\n",
                                       sum(isHet), 100*sum(isHet)/nbrOfSnps)
    rho <- rep(NA_real_, times=nbrOfLoci)
    rho[isHet] <- 2*abs(data$betaTN[isHet]-1/2)
    verbose && cat(verbose, "Normalized DHs:")
    verbose && str(verbose, rho)
    data$rho <- rho
    isSnp <- isHet <- rho <- NULL # Not needed anymore
    verbose && exit(verbose)
  }
  ## Sanity check
  .stop_if_not(!is.null(data$rho))



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Generate random seeds?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  seeds <- NULL
  if (!is.null(seed)) {
    randomSeed("set", seed=seed, kind="L'Ecuyer-CMRG")
    on.exit(randomSeed("reset"), add=TRUE)
    verbose && printf(verbose, "Random seed temporarily set (seed=c(%s), kind=\"L'Ecuyer-CMRG\")\n", paste(seed, collapse=", "))
    seeds <- randomSeed("advance", n=2L) ## For TCN and DH
    names(seeds) <- c("TCN", "DH")
    verbose && printf(verbose, "Produced %d seeds from this stream for future usage\n", length(seeds))
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # 1a. Identification of change points in total copy numbers
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Identification of change points by total copy numbers")

  fields <- attachLocally(data, fields=c("CT", "thetaT", "thetaN", "chromosome", "x", "index"))

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Assert no missing values in (chromosome, x, CT)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Sanity check
  ok <- (!is.na(data$chromosome) & !is.na(data$x))
  if (!all(ok)) {
    stop("INTERNAL ERROR: Detected (chromosome, x) with missing values also after filtering.")
  }

  # Sanity check
  if (dropMissingCT) {
    ok <- (!is.na(data$CT))
    if (!all(ok)) {
      stop("INTERNAL ERROR: Detected CT with missing values also after filtering.")
    }
  }


  # Segment TCN ratios
  # Calculate tumor-normal TCN ratios?
  fit <- segmentByCBS(CT,
                      chromosome=chromosome, x=x, index=index,
                      joinSegments=joinSegments,
                      knownSegments=knownSegments,
                      alpha=alphaTCN, undo=undoTCN, ...,
                      seed=seeds[["TCN"]],
                      verbose=verbose)
  verbose && str(verbose, fit)

  rm(list=fields) # Not needed anymore

  # Sanity check
  if (nrow(knownSegments) == 0) {
    .stop_if_not(nrow(fit$data) == nrow(data))
    .stop_if_not(all(fit$data$chromosome == data$chromosome))
    .stop_if_not(all(fit$data$x == data$x))
    .stop_if_not(all(fit$data$index == data$index))
    .stop_if_not(all.equal(fit$data$y, data$CT))
  }

  tcnSegments <- fit$output
  tcnSegRows <- fit$segRows
  fit <- NULL # Not needed anymore

  # Sanity checks
  .stop_if_not(all(tcnSegRows[,1] <= tcnSegRows[,2], na.rm=TRUE))
  .stop_if_not(all(tcnSegRows[-nrow(tcnSegRows),2] < tcnSegRows[-1,1], na.rm=TRUE))

  verbose && exit(verbose)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # 1b. Restructure TCN segmentation results
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Restructure TCN segmentation results")
  # Drop dummy columns
  keep <- setdiff(colnames(tcnSegments), c("sampleName"))
  tcnSegments <- tcnSegments[,keep,drop=FALSE]

  # Tag fields by TCN
  names <- names(tcnSegments)
  # Adding 'tcn' prefix to column names
  names <- sprintf("tcn%s", capitalize(names))
  names <- gsub("tcnChromosome", "chromosome", names, fixed=TRUE)
  names(tcnSegments) <- names
  verbose && print(verbose, tcnSegments)

  nbrOfSegs <- nrow(tcnSegments)
  verbose && cat(verbose, "Number of TCN segments: ", nbrOfSegs)

  verbose && exit(verbose)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # 2a. Identification of additional change points using DH
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # For each segment independently, segment decrease of heterozygousity (DH)
  # using CBS. By definition, only heterozygous SNPs are used.

  if (flavor == "tcn") {
    verbose && enter(verbose, "TCN-only segmentation")

    tcnSegsExpanded <- tcnSegRows
    dhSegRows <- tcnSegRows

    # Segments
    segs <- tcnSegments
    segs[,"tcnId"] <- seq_len(nbrOfSegs)
    segs[,"dhId"] <- rep(1L, times=nbrOfSegs)
    segs[,c("tcnNbrOfSNPs", "tcnNbrOfHets", "dhNbrOfLoci")] <- 0L
    segs[,"dhStart"] <- segs[,"tcnStart"]
    segs[,"dhEnd"] <- segs[,"tcnEnd"]

    # For each TCN segment...
    for (kk in seq_len(nbrOfSegs)) {
      tcnId <- kk

      xStart <- tcnSegments[kk,"tcnStart"]
      xEnd <- tcnSegments[kk,"tcnEnd"]
      regionTag <- sprintf("[%10g,%10g]", xStart, xEnd)
      verbose && enter(verbose, sprintf("Total CN segment #%d (%s) of %d", kk, regionTag, nbrOfSegs))

      # Empty segment?
      rowStart <- tcnSegRows[kk,1]
      rowEnd <- tcnSegRows[kk,2]

      # Empty segment or a segment separator?
      isEmptySegment <- (is.na(rowStart) && is.na(rowEnd))

      # Nothing to do?
      if (isEmptySegment) {
        verbose && exit(verbose)
        next
      }

      nbrOfTCNLociKK <- tcnSegments[kk,"tcnNbrOfLoci"]
      verbose && cat(verbose, "Number of TCN loci in segment: ", nbrOfTCNLociKK)
      rows <- seq(from=rowStart, length.out=nbrOfTCNLociKK)
      dataKK <- data[rows,,drop=FALSE]
      nbrOfLociKK <- nrow(dataKK)

      verbose && cat(verbose, "Locus data for TCN segment:")
      verbose && str(verbose, dataKK)

      verbose && cat(verbose, "Number of loci: ", nbrOfLociKK)
      hasDH <- !is.null(dataKK$rho)
      if (hasDH) {
        isSnpKK <- !is.na(dataKK$rho)
        isHetsKK <- (isSnpKK & (dataKK$rho > 0))
      } else {
        isSnpKK <- !is.na(dataKK$muN)
        isHetsKK <- (isSnpKK & (dataKK$muN == 1/2))
      }
      nbrOfSnpsKK <- sum(isSnpKK)
      nbrOfHetsKK <- sum(isHetsKK)
      verbose && printf(verbose, "Number of SNPs: %d (%.2f%%)\n",
                                    nbrOfSnpsKK, 100*nbrOfSnpsKK/nbrOfLociKK)

      verbose && printf(verbose, "Number of heterozygous SNPs: %d (%.2f%%)\n",
                                    nbrOfHetsKK, 100*nbrOfHetsKK/nbrOfSnpsKK)

      segs[kk,"tcnNbrOfSNPs"] <- nbrOfSnpsKK
      segs[kk,"tcnNbrOfHets"] <- nbrOfHetsKK
      segs[kk,"dhNbrOfLoci"] <- nbrOfHetsKK

      # Adjust 'dhRows[kk,]'
      rows <- rows[isHetsKK]
      rows <- range(rows, na.rm=TRUE)
      dhSegRows[kk,] <- rows

      # Sanity check
      if (nbrOfHetsKK > 0) {
        .stop_if_not(all(dhSegRows[kk,1] <= dhSegRows[kk,2], na.rm=TRUE))
      }

      # Calculate dhMean
      rhoKK <- dataKK[["rho"]][isHetsKK]
      segs[kk,"dhMean"] <- mean(rhoKK, na.rm=TRUE)

      verbose && exit(verbose)
    } # for (kk ...)

    # Reorder segmentation columns
    keys <- c("tcnId", "dhId", colnames(tcnSegments))
    keys <- c(keys, setdiff(colnames(segs), keys))
    segs <- segs[,keys]

    verbose && exit(verbose)
  } else {
    dhSegRows <- NULL
    tcnSegsExpanded <- NULL

    # For each TCN segment...
    segs <- vector("list", length=nbrOfSegs)
    for (kk in seq_len(nbrOfSegs)) {
      tcnId <- kk

      xStart <- tcnSegments[kk,"tcnStart"]
      xEnd <- tcnSegments[kk,"tcnEnd"]
      regionTag <- sprintf("[%10g,%10g]", xStart, xEnd)
      verbose && enter(verbose, sprintf("Total CN segment #%d (%s) of %d", kk, regionTag, nbrOfSegs))

      # Empty segment?
      rowStart <- tcnSegRows[kk,1]
      rowEnd <- tcnSegRows[kk,2]

      # Empty segment or a segment separator?
      isEmptySegment <- (is.na(rowStart) && is.na(rowEnd))
      isSplitter <- (isEmptySegment && is.na(xStart) && is.na(xEnd))
      isEmptySegment <- (isEmptySegment & !isSplitter)

      if (isSplitter) {
        verbose && cat(verbose, "No signals to segment. Just a \"splitter\" segment. Skipping.")

        # Sanity check
        .stop_if_not(kk >= 1)

        # Add a splitter segment
        segT <- segs[[kk-1]]
        segT <- segT[NA_integer_,]
        keys <- colnames(tcnSegments)
        segT[,keys] <- tcnSegments[kk,keys]
        segT[,"tcnId"] <- tcnId
        segT[,"dhId"] <- 1L
        segT[,c("tcnNbrOfSNPs", "tcnNbrOfHets", "dhNbrOfLoci")] <- 0L
        segT[,"dhStart"] <- xStart
        segT[,"dhEnd"] <- xEnd
        segs[[kk]] <- segT
        verbose && print(verbose, segT)

        # Add a splitter to TCN and DH segment row matrix
        segRowsT <- dhSegRows[NA_integer_,]
        dhSegRows <- rbind(dhSegRows, segRowsT)

        segRowsT <- tcnSegsExpanded[NA_integer_,]
        tcnSegsExpanded <- rbind(tcnSegsExpanded, segRowsT)

        verbose && exit(verbose)
        next
      } # if (isSplitter)


      nbrOfTCNLociKK <- tcnSegments[kk,"tcnNbrOfLoci"]
      verbose && cat(verbose, "Number of TCN loci in segment: ", nbrOfTCNLociKK)

      # Sanity check
      .stop_if_not(!isEmptySegment || (isEmptySegment && (nbrOfTCNLociKK == 0)))

      if (nbrOfTCNLociKK > 0) {
        # Extract locus data for TCN segment
        rows <- rowStart:rowEnd
    ##    if (nrow(knownSegments) == 0) {
    ##      gammaT <- tcnSegments[kk,"tcnMean"]
    ##      verbose && print(verbose, all.equal(mean(dataKK$CT, na.rm=TRUE), gammaT, tolerance=tol))
    ##      .stop_if_not(all.equal(mean(dataKK$CT, na.rm=TRUE), gammaT, tolerance=tol))
    ##    }
      } else {
        rows <- integer(0)
      } # if (nbrOfTCNLociKK > 0)

      dataKK <- data[rows,,drop=FALSE]
      nbrOfLociKK <- nrow(dataKK)

      # Sanity check
      .stop_if_not(length(dataKK$CT) == nbrOfTCNLociKK)
    ##  .stop_if_not(sum(!is.na(dataKK$CT)) == nbrOfTCNLociKK)

      verbose && cat(verbose, "Locus data for TCN segment:")
      verbose && str(verbose, dataKK)

      verbose && cat(verbose, "Number of loci: ", nbrOfLociKK)

      hasDH <- !is.null(dataKK$rho)
      if (hasDH) {
        isSnpKK <- !is.na(dataKK$rho)
        isHetsKK <- (isSnpKK & (dataKK$rho > 0))
      } else {
        isSnpKK <- !is.na(dataKK$muN)
        isHetsKK <- (isSnpKK & (dataKK$muN == 1/2))
      }
      nbrOfSnpsKK <- sum(isSnpKK)
      nbrOfHetsKK <- sum(isHetsKK)
      verbose && printf(verbose, "Number of SNPs: %d (%.2f%%)\n",
                                    nbrOfSnpsKK, 100*nbrOfSnpsKK/nbrOfLociKK)
      verbose && printf(verbose, "Number of heterozygous SNPs: %d (%.2f%%)\n",
                                    nbrOfHetsKK, 100*nbrOfHetsKK/nbrOfSnpsKK)

      # Since segments in 'knownSegments' has already been used in the TCN
      # segmentation, they are not needed in the DH segmentation.
      currChromosome <- data$chromosome[1]
      verbose && cat(verbose, "Chromosome: ", currChromosome)
      knownSegmentsT <- data.frame(chromosome=currChromosome, start=xStart, end=xEnd)

      verbose && enter(verbose, "Segmenting DH signals")
      fields <- attachLocally(dataKK, fields=c("chromosome", "x", "rho", "index"))

      fit <- segmentByCBS(rho,
                          chromosome=chromosome, x=x,
                          joinSegments=joinSegments,
                          knownSegments=knownSegmentsT,
                          alpha=alphaDH, undo=undoDH, ...,
                          seed=seeds[["DH"]],
                          verbose=verbose)
      verbose && str(verbose, fit)
      dhSegments <- fit$output
      dhSegRowsKK <- fit$segRows

      verbose && cat(verbose, "DH segmentation (locally-indexed) rows:")
      verbose && print(verbose, dhSegRowsKK)
      verbose && str(verbose, index)

      # Remap to genome-wide indices
      for (cc in 1:2) {
        dhSegRowsKK[,cc] <- index[dhSegRowsKK[,cc]]
      }

      verbose && cat(verbose, "DH segmentation rows:")
      verbose && print(verbose, dhSegRowsKK)

      # Not needed anymore
      rm(list=fields)
      fit <- NULL
      verbose && exit(verbose)

      # Drop dummy columns
      keep <- setdiff(colnames(dhSegments), c("sampleName", "chromosome"))
      dhSegments <- dhSegments[,keep,drop=FALSE]

      # Tag fields by DH
      names <- names(dhSegments)
      # Adding 'dh' prefix to column names
      names <- sprintf("dh%s", capitalize(names))
      names(dhSegments) <- names

      # Special case: If there where not enough data to segment DH...
      if (nrow(dhSegments) == 0) {
        dhSegments <- dhSegments[NA_integer_,,drop=FALSE]
        dhSegRowsKK <- dhSegRowsKK[NA_integer_,,drop=FALSE]
      }

      verbose && cat(verbose, "DH segmentation table:")
      verbose && print(verbose, dhSegments)
      verbose && print(verbose, dhSegRowsKK)

      # Expand the TCN segmentation result data frame
      rows <- rep(kk, times=nrow(dhSegments))
      verbose && cat(verbose, "Rows:")
      verbose && print(verbose, rows)
      tcnSegmentsKK <- tcnSegments[rows,,drop=FALSE]
      tcnSegRowsKK <- tcnSegRows[rows,,drop=FALSE]
      # Sanity check
      .stop_if_not(nrow(tcnSegmentsKK) == nrow(dhSegments))
      .stop_if_not(nrow(tcnSegRowsKK) == nrow(dhSegments))
      .stop_if_not(all(is.na(tcnSegRowsKK[,1]) | is.na(dhSegRowsKK[,1]) | (tcnSegRowsKK[,1] <= dhSegRowsKK[,1])))
      .stop_if_not(all(is.na(tcnSegRowsKK[,2]) | is.na(dhSegRowsKK[,2]) | (dhSegRowsKK[,2] <= tcnSegRowsKK[,2])))
      verbose && cat(verbose, "TCN segmentation rows:")
      verbose && print(verbose, tcnSegRowsKK)
      .stop_if_not(all(tcnSegRowsKK[,1] == tcnSegRowsKK[1,1], na.rm=TRUE))
      .stop_if_not(all(tcnSegRowsKK[,2] == tcnSegRowsKK[1,2], na.rm=TRUE))

      verbose && cat(verbose, "TCN and DH segmentation rows:")
      verbose && print(verbose, tcnSegRowsKK)
      verbose && print(verbose, dhSegRowsKK)
      verbose && print(verbose, tcnSegsExpanded)

      # Append
      tcnSegsExpanded <- rbind(tcnSegsExpanded, tcnSegRowsKK)
      verbose && cat(verbose, "TCN segmentation (expanded) rows:")
      verbose && print(verbose, tcnSegsExpanded)
      rownames(tcnSegsExpanded) <- NULL

      dhSegRows <- rbind(dhSegRows, dhSegRowsKK)
      rownames(dhSegRows) <- NULL

      verbose && cat(verbose, "TCN and DH segmentation rows:")
      verbose && print(verbose, tcnSegRows)
      verbose && print(verbose, dhSegRows)
      verbose && print(verbose, tcnSegsExpanded)

      # Sanity checks
      .stop_if_not(all(tcnSegRows[,1] <= tcnSegRows[,2], na.rm=TRUE))
      .stop_if_not(all(tcnSegRows[-nrow(tcnSegRows),2] < tcnSegRows[-1,1], na.rm=TRUE))
      .stop_if_not(all(dhSegRows[,1] <= dhSegRows[,2], na.rm=TRUE))
      .stop_if_not(all(dhSegRows[-nrow(dhSegRows),2] < dhSegRows[-1,1], na.rm=TRUE))
      .stop_if_not(all(tcnSegsExpanded[,1] <= tcnSegsExpanded[,2], na.rm=TRUE))
      .stop_if_not(all(tcnSegsExpanded[,1] <= dhSegRows[,1], na.rm=TRUE))
      .stop_if_not(all(tcnSegsExpanded[,2] >= dhSegRows[,2], na.rm=TRUE))
  ##    if (!all(tcnSegsExpanded[-nrow(tcnSegsExpanded),2] < tcnSegsExpanded[-1,1], na.rm=TRUE)) {
  ##      .stop_if_not(all(tcnSegsExpanded[-nrow(tcnSegsExpanded),2] < tcnSegsExpanded[-1,1], na.rm=TRUE))
  ##    }


      # Sanity check
      .stop_if_not(nrow(dhSegRows) == nrow(tcnSegsExpanded))

      # Append information on number of SNPs and hets in CN region
      tcnSegmentsKK <- cbind(
        tcnSegmentsKK,
        tcnNbrOfSNPs=nbrOfSnpsKK,
        tcnNbrOfHets=nbrOfHetsKK
      )
      verbose && cat(verbose, "Total CN segmentation table (expanded):")
      verbose && print(verbose, tcnSegmentsKK)

      # Sanity check
      .stop_if_not(nrow(tcnSegmentsKK) == nrow(dhSegments))

      # Combine TCN and DH segmentation results
      tcndhSegments <- cbind(
        tcnId=rep(kk, times=nrow(dhSegments)),
        dhId=seq_len(nrow(dhSegments)),
        tcnSegmentsKK,
        dhSegments
      )

      segs[[kk]] <- tcndhSegments

      verbose && cat(verbose, "(TCN,DH) segmentation for one total CN segment:")
      verbose && print(verbose, segs[[kk]])

      verbose && exit(verbose)
    } # for (kk ...)

    segs <- Reduce(rbind, segs)
    rownames(segs) <- NULL
  } # if (flavor == "tcn")

  # Sanity check
  .stop_if_not(nrow(dhSegRows) == nrow(tcnSegsExpanded))
  rownames(tcnSegRows) <- rownames(dhSegRows) <- NULL

  .stop_if_not(all(tcnSegRows[,1] <= tcnSegRows[,2], na.rm=TRUE))
  .stop_if_not(all(tcnSegRows[-nrow(tcnSegRows),2] < tcnSegRows[-1,1], na.rm=TRUE))
  if (flavor != "tcn") {
    .stop_if_not(all(dhSegRows[,1] <= dhSegRows[,2], na.rm=TRUE))
  }
  .stop_if_not(all(dhSegRows[-nrow(dhSegRows),2] < dhSegRows[-1,1], na.rm=TRUE))
  .stop_if_not(all(tcnSegsExpanded[,1] <= tcnSegsExpanded[,2], na.rm=TRUE))
##  .stop_if_not(all(tcnSegsExpanded[-nrow(tcnSegsExpanded),2] < tcnSegsExpanded[-1,1], na.rm=TRUE))

  # Move 'chromosome' column to the first column
  idx <- match("chromosome", names(segs))
  idxs <- c(idx, seq_len(ncol(segs))[-idx])
  segs <- segs[,idxs,drop=FALSE]
  verbose && print(verbose, segs)

  verbose && enter(verbose, "Calculating (C1,C2) per segment")
  # Append (C1,C2) estimates
  tcn <- segs$tcnMean
  dh <- segs$dhMean
  C1 <- 1/2*(1-dh)*tcn
  C2 <- tcn - C1
  segs <- cbind(segs, c1Mean=C1, c2Mean=C2)
  verbose && exit(verbose)

  nbrOfSegs <- nrow(segs)
  verbose && cat(verbose, "Number of segments: ", nbrOfSegs)

  verbose && exit(verbose)

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
  )

  # Should we drop attributes? /HB 2010-09-24
  .stop_if_not(all(data$index == seq_len(nrow(data))))
  data$index <- NULL # Drop, because it is guaranteed to be ordered
  class(data) <- c("PairedPSCNData", class(data))

  class(segs) <- c("PairedPSCNSegments", class(segs))

  fit <- list(
    data = data,
    output = segs,
    tcnSegRows = tcnSegsExpanded,
    dhSegRows = dhSegRows,
    params = params
  )

  class(fit) <- c("PairedPSCBS", "PSCBS", "AbstractCBS")


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Update?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (avgTCN != "mean" || avgDH != "mean") {
    verbose && enter(verbose, "Updating mean level using different estimator")
    verbose && cat(verbose, "TCN estimator: ", avgTCN)
    verbose && cat(verbose, "DH estimator: ", avgDH)
    fit <- updateMeans(fit, avgTCN=avgTCN, avgDH=avgDH, verbose=less(verbose, 20))
    verbose && exit(verbose)
  }

  if (is.element(flavor, c("tcn&dh", "sqrt(tcn)&dh"))) {
    fit$params$flavor <- gsub("&", ",", flavor, fixed=TRUE) # AD HOC.
    fit <- postsegmentTCN(fit, verbose=verbose)

    # Sanity check
    CT <- fit$data$CT
    tcnSegRows <- fit$tcnSegRows
    dhSegRows <- fit$dhSegRows
    for (jj in 1:nrow(tcnSegRows)) {
      tcnSegRowJJ <- unlist(tcnSegRows[jj,,drop=TRUE], use.names=FALSE)
      dhSegRowJJ <- unlist(dhSegRows[jj,,drop=TRUE], use.names=FALSE)
      .stop_if_not(
        is.na(tcnSegRowJJ[1]) || is.na(dhSegRowJJ[1]) ||
        # A TCN segment must start at or before a DH segment...
        (tcnSegRowJJ[1] <= dhSegRowJJ[1]) ||
        # ...unless there was an outlier at the left edge.
        (is.na(CT[dhSegRowJJ[1]]) && (tcnSegRowJJ[1] - 1L <= dhSegRowJJ[1]))
      )
      .stop_if_not(
        is.na(tcnSegRowJJ[2]) || is.na(dhSegRowJJ[2]) ||
        # A TCN segment must end at or after a DH segment...
        (dhSegRowJJ[2] <= tcnSegRowJJ[2]) ||
        # ...unless there was an outlier at the right edge.
        (is.na(CT[dhSegRowJJ[2]]) && (dhSegRowJJ[2] <= tcnSegRowJJ[2] + 1L))
      )
    } # for (jj ...)
    # Not needed anymore
    CT <- tcnSegRows <- dhSegRows <- NULL
  }

  verbose && print(verbose, head(as.data.frame(fit)))
  verbose && print(verbose, tail(as.data.frame(fit)))

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Return results
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  fit
}) # segmentByPairedPSCBS()



setMethodS3("segmentByPairedPSCBS", "data.frame", function(CT, ...) {
  # To please R CMD check
  data <- CT

  segmentByPairedPSCBS(CT=data$CT, thetaT=data$thetaT, thetaN=data$thetaN,
                       betaT=data$betaT, betaN=data$betaN,
                       muN=data$muN, rho=data$rho,
                       chromosome=data$chromosome, x=data$x, ...)
})



setMethodS3("segmentByPairedPSCBS", "PairedPSCBS", function(...) {
  resegment(...)
}) # segmentByPairedPSCBS()
