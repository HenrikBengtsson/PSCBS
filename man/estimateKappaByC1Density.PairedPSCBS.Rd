%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  PairedPSCBS.estimateKappa.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 \name{estimateKappaByC1Density.PairedPSCBS}
\alias{estimateKappaByC1Density.PairedPSCBS}
\alias{PairedPSCBS.estimateKappaByC1Density}
\alias{estimateKappaByC1Density,PairedPSCBS-method}

 \title{Estimate global background in segmented copy numbers}

 \description{
  Estimate global background in segmented copy numbers based on the location of peaks in a weighted
  density estimator of the minor copy number mean levels.
  The global background, here called \eqn{\kappa},
  may have multiple origins where normal contamination is one,
  but not necessarily the only one.
 }

 \usage{
\method{estimateKappaByC1Density}{PairedPSCBS}(this, typeOfWeights=c("dhNbrOfLoci", "sqrt(dhNbrOfLoci)"), adjust=1, minDensity=0.2, ..., verbose=FALSE)
}

 \arguments{
   \item{typeOfWeights}{A \code{\link[base]{character}} string specifying how weights
    are calculated.}
   \item{adjust}{A \code{\link[base]{numeric}} scale factor specifying the size of
    the bandwidth parameter used by the density estimator.}
   \item{minDensity}{A non-negative \code{\link[base]{numeric}} threshold specifying
    the minimum density a peak should have in order to consider
    it a peak.}
   \item{...}{Not used.}
   \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
 }

 \value{
   Returns the background estimate as a \code{\link[base]{numeric}} scalar.
 }

 \section{Algorithm}{
  \itemize{
  \item Retrieve segment-level minor copy numbers and corresponding weights:
   \enumerate{
    \item Grabs the segment-level C1 estimates.
    \item Calculate segment weights.
          The default (\code{typeOfWeights="dhNbrOfLoci"}) is to use
          weights proportional to the number of heterozygous SNPs.
          An alternative (\code{typeOfWeights="sqrt(dhNbrOfLoci)"}) is
          to use the square root of those counts.
   }

  \item Identify subset of regions with C1=0:
   \enumerate{
    \item Estimates the weighted empirical density function
          (truncated at zero below).  Tuning parameter 'adjust'.
    \item Find the first two peaks
          (with a density greater than tuning parameter 'minDensity').
    \item Assumes that the two peaks corresponds to C1=0 and C1=1.
    \item Defines threshold Delta0.5 as the center location between
          these two peaks.
   }

  \item Estimate the normal contamination:
   \enumerate{
    \item For all segments with C1 < Delta0.5, calculate the weighted
          median of their C1:s.
    \item Let kappa be the above weighted median.
          This is the estimated background.
   }
  }
 }

 \author{Henrik Bengtsson}

 \seealso{
   Instead of calling this method explicitly, it is recommended
   to use the \code{\link[PSCBS:estimateKappa.PairedPSCBS]{*estimateKappa}()} method.
 }


\keyword{internal}
\keyword{methods}
