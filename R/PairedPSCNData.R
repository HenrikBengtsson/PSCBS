###########################################################################/**
# @RdocClass PairedPSCNData
#
# @title "The PairedPSCNData class"
#
# \description{
#  @classhierarchy
#
#  A PairedPSCNData object holds paired tumor-normal parent-specific 
#  copy number data.
# }
# 
# @synopsis
#
# \arguments{
#   \item{CT}{A @numeric @vector of J tumor total copy number (TCN)
#        ratios in [0,+@Inf) (due to noise, small negative values are
#        also allowed).  The TCN ratios are typically scaled such that
#        copy-neutral diploid loci have a mean of two.}
#   \item{CN}{An optional @numeric @vector of J normal TCN ratios.}
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
#   \item{isSNP}{An optional @logical @vector of length J specifying
#        whether each locus is a SNP or not (non-polymorphic loci).}
#   \item{chromosome}{(Optional) An @integer scalar (or a @vector of length J),
#        which can be used to specify which chromosome each locus belongs to
#        in case multiple chromosomes are segments.
#        This argument is also used for annotation purposes.}
#   \item{x}{Optional @numeric @vector of J genomic locations.
#            If @NULL, index locations \code{1:J} are used.}
#   \item{...}{Optional named locus-specific signal @vectors of length J.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
#*/###########################################################################
setConstructorS3("PairedPSCNData", function(chromosome=NULL, x=NULL, isSNP=NULL, muN=NULL, CT=NULL, betaT=NULL, CN=NULL, betaN=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!is.null(chromosome)) {
    # Is first argument special?
    if (is.data.frame(chromosome)) {
      data <- chromosome;

      chromosome <- data$chromosome;
      x <- data$x;
      isSNP <- data$isSNP;
      muN <- data$muN;
      CT <- data$CT;
      betaT <- data$betaT;
      CN <- data$CN;
      betaN <- data$betaN;
    }

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Required arguments
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Argument 'CT':
    disallow <- c("Inf");
    CT <- Arguments$getDoubles(CT, disallow=disallow);
  
    nbrOfLoci <- length(CT);
    length2 <- rep(nbrOfLoci, times=2);
  
    # Argument 'betaT':
    betaT <- Arguments$getDoubles(betaT, length=length2, disallow="Inf");
  

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Mutually optional arguments (that are not validated in superclass)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Argument 'betaN':
    if (!is.null(betaN)) {
      betaN <- Arguments$getDoubles(betaN, length=length2, disallow="Inf");
    }
  
    if (is.null(betaN) && is.null(muN)) {
      throw("If argument 'betaN' is not given, then 'muN' must be.");
    }
  

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Optional arguments
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Argument 'CN':
    if (!is.null(CN)) {
      disallow <- c("Inf");
      CN <- Arguments$getDoubles(CN, length=length2, disallow=disallow);
    }

    # Argument 'chromosome':
    if (is.null(chromosome)) {
      chromosome <- 0L; 
    } else {
      disallow <- c("Inf");
      chromosome <- Arguments$getIntegers(chromosome, range=c(0,Inf), disallow=disallow);
    }
  
    if (length(chromosome) == 1) {
      chromosome <- rep(chromosome, times=nbrOfLoci);
    } else {
      chromosome <- Arguments$getIntegers(chromosome, length=length2, disallow=disallow);
    }
  }

  extend(AbstractPSCNData(chromosome=chromosome, x=x, isSNP=isSNP, mu=muN, CT=CT, CN=CN, betaT=betaT, betaN=betaN, ...), "PairedPSCNData");
})


setMethodS3("as.PairedPSCNData", "PairedPSCNData", function(this, ...) {
  this;
})

setMethodS3("as.PairedPSCNData", "data.frame", function(this, ...) {
  data <- this;
  PairedPSCNData(data, ...);
})

setMethodS3("as.PairedPSCNData", "NonPairedPSCNData", function(T, N, ...) {
  dataT <- T;
  dataN <- N;

  for (key in c("chromosome", "x", "isSNP")) {
    vT <- dataT[[key]];
    vN <- dataN[[key]];
    if (is.null(vT) || is.null(vN)) {
      next;
    }
    stopifnot(all.equal(vT, vN));
  } # for (key ...)

  names <- colnames(dataT);
  idxs <- which(is.element(names, c("C", "beta", "mu")));
  names[idxs] <- sprintf("%sT", names[idxs])
  colnames(dataT) <- names;

  names <- colnames(dataN);
  idxs <- which(is.element(names, c("C", "beta", "mu")));
  names[idxs] <- sprintf("%sN", names[idxs]);
  colnames(dataN) <- names;
  
  data <- cbind(dataT, dataN);
  names <- colnames(data);
  keep <- !duplicated(names);
  data <- data[,keep,drop=FALSE];

  PairedPSCNData(data, ...);
})


setMethodS3("extractNonPairedPSCNData", "PairedPSCNData", function(this, which=c("T", "N"), ...) {
  # Argument 'which':
  which <- match.arg(which);

  # Identify fields to keep
  fields0 <- c("chromosome", "x", "isSNP");
  fieldsW <- sprintf("%s%s", c("mu", "C", "beta"), which);
  fields <- c(fields0, fieldsW);

  # Extract subset of fields
  data <- getLocusData(this);
  names <- colnames(data);
  keep <- names[is.element(names, fields)];
  data <- data[,keep, drop=FALSE];

  # Rename fields
  names <- colnames(data);
  idxs <- match(fieldsW, names);
  idxs <- na.omit(idxs);
  names[idxs] <- gsub("N$", "", names[idxs]);
  colnames(data) <- names;

  # Return
  NonPairedPSCNData(data);
}, protected=TRUE)



setMethodS3("callNaiveGenotypes", "PairedPSCNData", function(this, force=FALSE, ..., verbose=FALSE) {
  data <- this;

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Calling genotypes from normal allele B fractions");

  muN <- data$muN;
  if (!force && !is.null(muN)) {
    verbose && cat(verbose, "Already done. Skipping.");
    verbose && exit(verbose);
    return(data);
  }

  dataN <- extractNonPairedPSCNData(this, which="N");
  dataN <- callNaiveGenotypes(dataN, ..., verbose=less(verbose, 5));
  muN <- dataN$mu;

  # Not needed anymore
  rm(dataN);  

  data$muN <- muN;
  verbose && exit(verbose);

  data;
})


setMethodS3("normalizeTumorBoost", "PairedPSCNData", function(this, preserveScale=TRUE, force=FALSE, ..., verbose=FALSE) {
  data <- this;

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Normalizing tumor BAFs using normal BAFs (TumorBoost)");

  betaTN <- data$betaTN;
  if (!force && !is.null(betaTN)) {
    verbose && cat(verbose, "Already done. Skipping.");
    verbose && exit(verbose);
    return(data);
  }

  betaT <- data$betaT;
  betaN <- data$betaN;
  muN <- data$muN;

  betaTN <- aroma.light::normalizeTumorBoost(betaT=betaT, betaN=betaN, muN=muN, preserveScale=preserveScale);
  verbose && cat(verbose, "Normalized tumor BAFs:");
  verbose && str(verbose, betaTN);

  # Assert that no missing values where introduced
  keep <- (is.finite(betaT) & is.finite(betaN) & is.finite(muN));
  if (any(is.na(betaTN[keep]))) {
    throw("Internal error: normalizeTumorBoost() introduced missing values.");
  }
  rm(keep);

  data$betaTN <- betaTN;
  verbose && exit(verbose);

  data;
}) # normalizeTumorBoost()


setMethodS3("getTotalCopyNumbers", "PairedPSCNData", function(this, what=c("CT", "2*CT/CN"), ...) {
  data <- this;

  # Argument 'what':
  what <- match.arg(what);

  CT <- data$CT;

  if (what == "2*CT/CN") {
    CN <- data$CN;
    CT <- 2*CT/CN;
  }

  CT;
})

setMethodS3("getTCNs", "PairedPSCNData", function(...) {
  getTotalCopyNumbers(...);
})


setMethodS3("callSegmentationOutliers", "PairedPSCNData", function(y, ...) {
  # To please R CMD check
  data <- y;

  callSegmentationOutliers(y=data$CT, chromosome=data$chromosome, x=data$x, ...);
}) # callSegmentationOutliers()


setMethodS3("dropSegmentationOutliers", "PairedPSCNData", function(CT, ...) {
  # To please R CMD check
  data <- CT;

  isOutlier <- callSegmentationOutliers(data, ...);
  naValue <- as.double(NA);
  data$CT[isOutlier] <- naValue;
  data;
}) # dropSegmentationOutliers()


setMethodS3("segmentByPairedPSCBS", "PairedPSCNData", function(CT, ...) {
  # To please R CMD check
  data <- CT;

  segmentByPairedPSCBS(CT=data$CT, betaT=data$betaT, betaN=data$betaN, 
                   muN=data$muN, chromosome=data$chromosome, x=data$x, ...); 
}) # segmentByPairedPSCBS()



############################################################################
# HISTORY:
# 2012-02-29
# o Added extractNonPairedPSCNData() for PairedPSCNData.
# o Added callSegmentationOutliers() and dropSegmentationOutliers() 
#   for PairedPSCNData.
# o Added segmentByPairedPSCBS() for PairedPSCNData.
# o Added as.PairedPSCNData().
# o Added callNaiveGenotypes() for PairedPSCNData.
# o Added normalizeTumorBoost() for PairedPSCNData.
# o Created.
############################################################################
