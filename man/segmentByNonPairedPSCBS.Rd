%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  segmentByNonPairedPSCBS.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{segmentByNonPairedPSCBS}
\alias{segmentByNonPairedPSCBS.default}
\alias{segmentByNonPairedPSCBS}
 \alias{segmentByNonPairedPSCBS.data.frame}
 \alias{segmentByNonPairedPSCBS.PairedPSCBS}
 \alias{segmentByNonPairedPSCBS}

 \title{Segment total copy numbers and allele B fractions using the Non-paired PSCBS method}

 \description{
  Segment total copy numbers and allele B fractions using the Non-paired PSCBS method [1].
  This method does not requires matched normals.
  This is a low-level segmentation method.
  It is intended to be applied to one tumor sample at the time.
 }

 \usage{\method{segmentByNonPairedPSCBS}{default}(CT, betaT, ..., flavor=c("tcn", "tcn&dh", "tcn,dh", "sqrt(tcn),dh", "sqrt(tcn)&dh"), tauA=NA, tauB=1 - tauA, verbose=FALSE)}

 \arguments{
   \item{CT}{A \code{\link[base]{numeric}} \code{\link[base]{vector}} of J tumor total copy number (TCN)
        ratios in [0,+\code{\link[base:is.finite]{Inf}}) (due to noise, small negative values are
        also allowed).  The TCN ratios are typically scaled such that
        copy-neutral diploid loci have a mean of two.}
   \item{betaT}{A \code{\link[base]{numeric}} \code{\link[base]{vector}} of J tumor allele B fractions (BAFs)
        in [0,1] (due to noise, values may be slightly outside as well)
        or \code{\link[base]{NA}} for non-polymorphic loci.}
   \item{...}{Additional arguments passed to \code{\link{segmentByPairedPSCBS}}().}
   \item{flavor}{A \code{\link[base]{character}} specifying what type of segmentation and
     calling algorithm to be used.}
   \item{tauA, tauB}{...}
   \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
 }

 \value{
   Returns the segmentation results as a \code{\link{NonPairedPSCBS}} object.
 }

 \details{
   Internally \code{\link{segmentByPairedPSCBS}}() is used for segmentation.
   This segmentation method does \emph{not} support weights.
 }

 \section{Reproducibility}{
   The "DNAcopy::segment" implementation of CBS uses approximation
   through random sampling for some estimates.  Because of this,
   repeated calls using the same signals may result in slightly
   different results, unless the random seed is set/fixed.
 }

 \section{Whole-genome segmentation is preferred}{
   Although it is possible to segment each chromosome independently
   using Paired PSCBS, we strongly recommend to segment whole-genome
   (TCN,BAF) data at once.  The reason for this is that downstream
   CN-state calling methods, such as the AB and the LOH callers,
   performs much better on whole-genome data.  In fact, they may
   fail to provide valid calls if done chromsome by chromosome.
 }

 \section{Missing and non-finite values}{
   The total copy number signals as well as any optional positions
   must not contain missing values, i.e. \code{\link[base]{NA}}s or \code{\link[base:is.finite]{NaN}}s.
   If there are any, an informative error is thrown.
   Allele B fractions may contain missing values, because such are
   interpreted as representing non-polymorphic loci.

   None of the input signals may have infinite values, i.e. -\code{\link[base:is.finite]{Inf}} or +\code{\link[base:is.finite]{Inf}}.
   If so, an informative error is thrown.
 }

 \section{Paired PSCBS with only genotypes}{
   If allele B fractions for the matched normal (\code{betaN}) are
   not available, but genotypes (\code{muN}) are, then it is possible
   to run a version of Paired PSCBS where TumorBoost normalization
   of the tumor allele B fractions is skipped.  In order for this
   to work, argument \code{tbn} must be set to \code{\link[base:logical]{FALSE}}.
 }

 \examples{

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Load SNP microarray data
# (note to package developers: this example data set may
#  be replaced in a future release of the package)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
pathname <- system.file("data-ex/PairedPSCBS,exData,chr01.Rbin", package="PSCBS")
data <- R.utils::loadObject(pathname)
str(data)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Paired PSCBS segmentation
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Drop single-locus outliers
dataS <- dropSegmentationOutliers(data)

# Speed up example by segmenting fewer loci
dataS <- dataS[seq(from=1, to=nrow(data), by=10),]

str(dataS)

R.oo::attachLocally(dataS)

# Paired PSCBS segmentation
fit <- segmentByPairedPSCBS(CT, betaT=betaT, betaN=betaN,
                            chromosome=chromosome, x=x,
                            seed=0xBEEF, verbose=-10)
print(fit)

# Plot results
plotTracks(fit)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Bootstrap segment level estimates
# (used by the AB caller, which, if skipped here,
#  will do it automatically)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fit <- bootstrapTCNandDHByRegion(fit, B=100, verbose=-10)
print(fit)
plotTracks(fit)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Calling segments in allelic balance (AB)
# NOTE: Ideally, this should be done on whole-genome data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Explicitly estimate the threshold in DH for calling AB
# (which be done by default by the caller, if skipped here)
deltaAB <- estimateDeltaAB(fit, flavor="qq(DH)", verbose=-10)
print(deltaAB)
## [1] 0.1657131

fit <- callAB(fit, delta=deltaAB, verbose=-10)
print(fit)
plotTracks(fit)

# Even if not explicitly specified, the estimated
# threshold parameter is returned by the caller
stopifnot(fit$params$deltaAB == deltaAB)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Calling segments in loss-of-heterozygosity (LOH)
# NOTE: Ideally, this should be done on whole-genome data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Explicitly estimate the threshold in C1 for calling LOH
# (which be done by default by the caller, if skipped here)
deltaLOH <- estimateDeltaLOH(fit, flavor="minC1|nonAB", verbose=-10)
print(deltaLOH)
## [1] 0.625175

fit <- callLOH(fit, delta=deltaLOH, verbose=-10)
print(fit)
plotTracks(fit)

# Even if not explicitly specified, the estimated
# threshold parameter is returned by the caller
stopifnot(fit$params$deltaLOH == deltaLOH)
}

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \references{
  [1] A.B. Olshen, H. Bengtsson, P. Neuvial, P.T. Spellman, R.A. Olshen, V.E. Seshan, \emph{Parent-specific copy number in paired tumor-normal studies using circular binary segmentation}, Bioinformatics, 2011
 \cr
  [2] H. Bengtsson, P. Neuvial and T.P. Speed, \emph{TumorBoost: Normalization of allele-specific tumor copy numbers from a single pair of tumor-normal genotyping microarrays}, BMC Bioinformatics, 2010
 \cr
 }

 \seealso{
   To segment paired tumor-normal total copy numbers and allele B fractions,
   see \code{\link{segmentByPairedPSCBS}}().

   To segment total copy numbers, or any other unimodal signals,
   see \code{\link{segmentByCBS}}().
 }


\keyword{IO}
