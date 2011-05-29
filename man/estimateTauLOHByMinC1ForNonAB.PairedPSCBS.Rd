%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  PairedPSCBS.estimateTauLOH.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{estimateTauLOHByMinC1ForNonAB.PairedPSCBS}
\alias{estimateTauLOHByMinC1ForNonAB.PairedPSCBS}
\alias{PairedPSCBS.estimateTauLOHByMinC1ForNonAB}
\alias{estimateTauLOHByMinC1ForNonAB.PairedPSCBS}
\alias{estimateTauLOHByMinC1ForNonAB,PairedPSCBS-method}


\title{Estimate a threshold for calling LOH from DH}

\description{
 Estimate a threshold for calling LOH from DH based on the location of guessed C1=0 and C1=1 peaks.
}

\usage{\method{estimateTauLOHByMinC1ForNonAB}{PairedPSCBS}(this, midpoint=1/2, maxC=3, ..., verbose=FALSE)}

\arguments{
  \item{midpoint}{A \code{\link[base]{numeric}} scalar in [0,1] specifying the relative
   position of the midpoint between the estimated locations of
   C1=0 and C1=1 mean parameters.}
  \item{maxC}{Maximum total copy number of a segment in order to
   be included in the initial set of segments.}
  \item{...}{Not used.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
  Returns the estimated LOH treshold as a \code{\link[base]{numeric}} scalar.
}

\details{
  This method requires that calls for allelic balances already have
  been me made, cf. \code{\link[PSCBS:callAllelicBalance.PairedPSCBS]{*callAllelicBalance}()}.
}

\section{Algorithm}{
 \itemize{
  \item Grabs the segment-level C1 estimates.
  \item Calculate segment weights proportional to the number of heterozygous SNPs.
  \item Estimate the C1=1 location as the weighted median C1 for segments that have been called to be in allelic balance.
  \item Estimate the C1=0 location as the smallest C1 among segments that are not in allelic balance.
  \item Let the LOH threshold be the midpoint of the estimates C1=0 and C1=1 locations.
 }
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  Instead of calling this method explicitly, it is recommended
  to use the \code{\link[PSCBS:estimateTauLOH.PairedPSCBS]{*estimateTauLOH}()} method.
}


\keyword{internal}
\keyword{methods}
