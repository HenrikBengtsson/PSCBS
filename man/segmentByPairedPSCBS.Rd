%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  segmentByPairedPSCBS.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{segmentByPairedPSCBS}
\alias{segmentByPairedPSCBS.default}
\alias{segmentByPairedPSCBS}
 \alias{segmentByPairedPSCBS.data.frame}
 \alias{segmentByPairedPSCBS.PairedPSCBS}
 \alias{segmentByPairedPSCBS}

 \title{Segment total copy numbers and allele B fractions using the Paired PSCBS method}

 \description{
  Segment total copy numbers and allele B fractions using the Paired PSCBS method [1].
  This method requires matched normals.
  This is a low-level segmentation method.
  It is intended to be applied to one tumor-normal sample at the time.
 }

 \usage{\method{segmentByPairedPSCBS}{default}(CT, betaT, betaN=NULL, muN=NULL, chromosome=0, x=NULL, alphaTCN=0.009, alphaDH=0.001, undoTCN=Inf, undoDH=Inf, ..., flavor=c("tcn&dh", "tcn,dh", "sqrt(tcn),dh", "sqrt(tcn)&dh"), tbn=TRUE, joinSegments=TRUE, knownSegments=NULL, seed=NULL, verbose=FALSE)}

 \arguments{
   \item{CT}{A \code{\link[base]{numeric}} \code{\link[base]{vector}} of J tumor total tumor copy number (TCN)
        ratios in [0,+\code{\link[base:is.finite]{Inf}}) (due to noise, small negative values are
        also allowed).  The TCN ratios are typically scaled such that
        copy-neutral diploid loci have a mean of two.}
   \item{betaT}{A \code{\link[base]{numeric}} \code{\link[base]{vector}} of J tumor allele B fractions (BAFs)
        in [0,1] (due to noise, values may be slightly outside as well)
        or \code{\link[base]{NA}} for non-polymorphic loci.}
   \item{betaN}{A \code{\link[base]{numeric}} \code{\link[base]{vector}} of J matched normal BAFs in [0,1]
        (due to noise, values may be slightly outside as well) or \code{\link[base]{NA}}
        for non-polymorphic loci.}
   \item{muN}{An optional \code{\link[base]{numeric}} \code{\link[base]{vector}} of J genotype calls in
        \{0,1/2,1\} for AA, AB, and BB, respectively,
        and \code{\link[base]{NA}} for non-polymorphic loci.
        If not given, they are estimated from the normal BAFs using
        \code{\link[aroma.light]{callNaiveGenotypes}} as described in [2].}
   \item{chromosome}{(Optional) An \code{\link[base]{integer}} scalar (or a \code{\link[base]{vector}} of length J),
        which can be used to specify which chromosome each locus belongs to
        in case multiple chromosomes are segments.
        This argument is also used for annotation purposes.}
   \item{x}{Optional \code{\link[base]{numeric}} \code{\link[base]{vector}} of J genomic locations.
            If \code{\link[base]{NULL}}, index locations \code{1:J} are used.}
   \item{alphaTCN, alphaDH}{The significance levels for segmenting total
        copy numbers (TCNs) and decrease-in-heterozygosity signals (DHs),
        respectively.}
   \item{undoTCN, undoDH}{Non-negative \code{\link[base]{numeric}}s.  If less than +\code{\link[base:is.finite]{Inf}},
        then a cleanup of segmentions post segmentation is done.
        See argument \code{undo} of \code{\link{segmentByCBS}}() for more
        details.}
   \item{...}{Additional arguments passed to \code{\link{segmentByCBS}}().}
   \item{flavor}{A \code{\link[base]{character}} specifying what type of segmentation and
     calling algorithm to be used.}
   \item{tbn}{If \code{\link[base:logical]{TRUE}}, \code{betaT} is normalized before segmentation
     using the TumorBoost method [2], otherwise not.}
   \item{joinSegments}{If \code{\link[base:logical]{TRUE}}, there are no gaps between neighboring
     segments.
     If \code{\link[base:logical]{FALSE}}, the boundaries of a segment are defined by the support
     that the loci in the segments provides, i.e. there exist a locus
     at each end point of each segment.  This also means that there
     is a gap between any neighboring segments, unless the change point
     is in the middle of multiple loci with the same position.
     The latter is what \code{DNAcopy::segment()} returns.
   }
   \item{knownSegments}{Optional \code{\link[base]{data.frame}} specifying
     \emph{non-overlapping} known segments.  These segments must
     not share loci.}
   \item{seed}{An (optional) \code{\link[base]{integer}} specifying the random seed to be
     set before calling the segmentation method.  The random seed is
     set to its original state when exiting.  If \code{\link[base]{NULL}}, it is not set.}
   \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
 }

 \value{
   Returns the segmentation results as a \code{\link{PairedPSCBS}} object.
 }

 \details{
   Internally \code{\link{segmentByCBS}}() is used for segmentation.
   The Paired PSCBS segmentation method does \emph{not} support weights.
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

# Order by chromosome and position
o <- order(data$chromosome, data$position)
data <- data[o,]
str(data)
R.oo::attachLocally(data)
x <- position
J <- length(x)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Paired PSCBS segmentation
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Drop single-locus outliers
CTs <- dropSegmentationOutliers(CT, chromosome=1, x=x, verbose=-10)

# Paired PSCBS segmentation
fit <- segmentByPairedPSCBS(CTs, betaT=betaT, betaN=betaN,
                            chromosome=1, x=x,
                            seed=0xBEEF, verbose=-10)
print(fit)

# Plot results
plotTracks(fit)

# Sanity check
stopifnot(nbrOfSegments(fit) == 12)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Bootstrap segment level estimates
# (used by the AB caller, which, if skipped here,
#  will do it automatically)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fit <- bootstrapTCNandDHByRegion(fit, verbose=-10)
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
   Internally, \code{\link[aroma.light]{callNaiveGenotypes}} is used to
   call naive genotypes, \code{\link[aroma.light]{normalizeTumorBoost}} is
   used for TumorBoost normalization, and \code{\link{segmentByCBS}}() is used
   to segment TCN and DH separately.

   To segment total copy-numbers, or any other unimodal signals,
   see \code{\link{segmentByCBS}}().
 }


\keyword{IO}
