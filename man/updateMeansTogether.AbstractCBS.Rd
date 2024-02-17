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
 \item{...}{Not used.}
}

\value{
  Returns an object of the same class.
}

\author{Henrik Bengtsson}

\seealso{
  This method is utilized by \code{\link[PSCBS:pruneByHClust.AbstractCBS]{*pruneByHClust}()}.
}


\keyword{internal}
\keyword{methods}
