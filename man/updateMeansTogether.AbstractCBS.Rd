%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AbstractCBS.HCLUST.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{updateMeansTogether.AbstractCBS}
\alias{updateMeansTogether.AbstractCBS}
\alias{AbstractCBS.updateMeansTogether}
\alias{updateMeansTogether,AbstractCBS-method}
\alias{updateMeansTogether.CBS}
\alias{updateMeansTogether.PairedPSCBS}

\title{Updates the CN mean levels jointly in sets of segments}

\description{
 Updates the CN mean levels jointly in sets of segments as if they were one large segment.
 The locus-level data is not updated/modified.
}

\usage{
\method{updateMeansTogether}{AbstractCBS}(...)
}

\arguments{
 \item{idxList}{A \code{\link[base]{list}}, where each element is an \code{\link[base]{integer}} \code{\link[base]{vector}}
   specifying segment indices of segments for which the mean levels
   should be calculated jointly.}
 \item{...}{Not used.}
 \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
  Returns an object of the same class.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  This method is utilized by \code{\link[PSCBS:pruneByHClust.AbstractCBS]{*pruneByHClust}()}.
}


\keyword{internal}
\keyword{methods}
