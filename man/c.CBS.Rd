%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  CBS.RESTRUCT.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{c.CBS}
\alias{c.CBS}
\alias{CBS.c}
\alias{c,CBS-method}
\alias{c.PSCBS}

\title{Concatenates segmentation results}

\description{
  Concatenates segmentation results.
}

\usage{
\method{c}{CBS}(..., addSplit=TRUE)
}

\arguments{
 \item{\dots}{One or more \code{\link{AbstractCBS}} objects to be combined.}
 \item{addSplit}{If \code{\link[base:logical]{TRUE}}, a "divider" is added between chromosomes.}
}

\value{
  Returns an \code{\link{AbstractCBS}} object of the same class in \dots.
}

\author{Henrik Bengtsson}

\seealso{
  For more information see \code{\link{CBS}}.
}
\keyword{internal}
\keyword{methods}
