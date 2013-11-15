%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  CBS.joinSegments.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{joinSegments.CBS}
\alias{joinSegments.CBS}
\alias{CBS.joinSegments}
\alias{joinSegments,CBS-method}

\title{Joins neighboring segments such that there is no gap in between them}

\description{
 Joins neighboring segments such that there is no gap in between them.
 For instance, consider two neighboring segments [x1,x2] and [x3,x4]
 with x1 < x2 < x3 < x4.  After join the segments, they are
 [x1,x23] and [x23,x4] where x23 = (x2 + x3)/2.
}

\usage{
\method{joinSegments}{CBS}(fit, range=NULL, verbose=FALSE, ...)
}

\arguments{
  \item{range}{(optional) A \code{\link[base]{numeric}} \code{\link[base]{vector}} of length two.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
  \item{...}{Not used.}
}

\value{
  Returns an updated \code{\link{CBS}} object.
}

\details{
  This function assumes only chromosome exists.
  If more, an error will be thrown.
}

\author{Henrik Bengtsson}



\keyword{internal}
\keyword{methods}
\keyword{IO}
