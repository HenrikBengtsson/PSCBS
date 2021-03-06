%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  callSegmentationOutliers.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{callSegmentationOutliers}
\alias{callSegmentationOutliers}

\alias{callSegmentationOutliers.default}
\alias{callSegmentationOutliers.data.frame}
\alias{dropSegmentationOutliers}
\alias{dropSegmentationOutliers.default}
\alias{dropSegmentationOutliers.data.frame}

\title{Calls/drops single-locus outliers along the genome}

\description{
 Calls/drops single-locus outliers along the genome that have a signal that differ significantly from the
 neighboring loci.
}

\usage{
 \method{callSegmentationOutliers}{default}(y, chromosome=0, x=NULL, method="DNAcopy::smooth.CNA", ...,
  verbose=FALSE)
 \method{callSegmentationOutliers}{data.frame}(y, ...)
 \method{dropSegmentationOutliers}{default}(y, ...)
 \method{dropSegmentationOutliers}{data.frame}(y, ...)
}

\arguments{
  \item{y}{A \code{\link[base]{numeric}} \code{\link[base]{vector}} of J genomic signals to be segmented.}
  \item{chromosome}{(Optional) An \code{\link[base]{integer}} scalar
      (or a \code{\link[base]{vector}} of length J contain a unique value).
      Only used for annotation purposes.}
  \item{x}{Optional \code{\link[base]{numeric}} \code{\link[base]{vector}} of J genomic locations.
           If \code{\link[base]{NULL}}, index locations \code{1:J} are used.}
  \item{method}{A \code{\link[base]{character}} string specifying the method
       used for calling outliers.}
  \item{...}{Additional arguments passed to internal outlier
       detection method.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
  \code{callSegmentationOutliers()} returns a \code{\link[base]{logical}} \code{\link[base]{vector}} of length J.
  \code{dropSegmentationOutliers()} returns an object of the same type
  as argument \code{y}, where the signals for which outliers were called
  have been set to \code{\link[base]{NA}}.
}

\section{Missing and non-finite values}{
  Signals as well as genomic positions may contain missing
  values, i.e. \code{\link[base]{NA}}s or \code{\link[base:is.finite]{NaN}}s.  By definition, these cannot
  be outliers.
}

\author{Henrik Bengtsson}

\seealso{
  Internally \code{\link[DNAcopy]{smooth.CNA}} is utilized to identify
  the outliers.
}


\keyword{methods}
\keyword{IO}
