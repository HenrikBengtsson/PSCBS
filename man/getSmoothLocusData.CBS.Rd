%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  CBS.SMOOTH.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{getSmoothLocusData.CBS}
\alias{getSmoothLocusData.CBS}
\alias{CBS.getSmoothLocusData}
\alias{getSmoothLocusData,CBS-method}

\title{Gets smoothed locus-level data}

\description{
 Gets smoothed locus-level data.
}

\usage{
\method{getSmoothLocusData}{CBS}(fit, by, ...)
}

\arguments{
  \item{fit}{An \code{\link{CBS}} object.}
  \item{by}{A \code{\link[base]{numeric}} scalar specifying the bin size.}
  \item{...}{Not used.}
}

\value{
  Returns a \code{\link[base]{data.frame}} where the
  first two columns are 'chromosome' and 'x' (position)
  and the remaining ones are the smoothed locus-level data.
}

\author{Henrik Bengtsson}

\seealso{
  For more information see \code{\link{CBS}}.
}


\keyword{internal}
\keyword{methods}
