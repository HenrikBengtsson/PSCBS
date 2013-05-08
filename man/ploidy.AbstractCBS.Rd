%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AbstractCBS.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{ploidy.AbstractCBS}
\alias{ploidy.AbstractCBS}
\alias{AbstractCBS.ploidy}
\alias{ploidy,AbstractCBS-method}
\alias{AbstractCBS.ploidy<-}
\alias{ploidy<-.AbstractCBS}
\alias{ploidy<-,AbstractCBS-method}

\alias{AbstractCBS.setPloidy}
\alias{setPloidy.AbstractCBS}
\alias{setPloidy,AbstractCBS-method}

\alias{AbstractCBS.adjustPloidyScale}
\alias{adjustPloidyScale.AbstractCBS}
\alias{adjustPloidyScale,AbstractCBS-method}

\alias{adjustPloidyScale.PairedPSCBS}
\alias{adjustPloidyScale}
\alias{ploidy}
\alias{ploidy<-}
\alias{setPloidy}

\title{Gets and sets ploidy}

\description{
 Gets and sets ploidy.
}

\usage{
  \method{ploidy}{AbstractCBS}(fit, ...)
  \method{ploidy}{AbstractCBS}(fit) <- value
}

\arguments{
  \item{fit}{An \code{\link{AbstractCBS}} object.}
  \item{value}{An \code{\link[base]{integer}} (in \eqn{1,2,\ldots}) specifying the genome ploidy .}
  \item{...}{Not used.}
}

\value{
  Returns (invisibly) an updated object.
}

\author{Henrik Bengtsson}

\seealso{
  For more information see \code{\link{AbstractCBS}}..
}


\keyword{internal}
\keyword{methods}
