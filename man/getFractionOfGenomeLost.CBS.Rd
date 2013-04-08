%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  CBS.CALL.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getFractionOfGenomeLost.CBS}
\alias{getFractionOfGenomeLost.CBS}
\alias{CBS.getFractionOfGenomeLost}
\alias{getFractionOfGenomeLost,CBS-method}
\alias{CBS.getFractionOfGenomeGained}
\alias{getFractionOfGenomeGained.CBS}
\alias{getFractionOfGenomeGained,CBS-method}

\alias{CBS.getFractionOfGenomeAltered}
\alias{getFractionOfGenomeAltered.CBS}
\alias{getFractionOfGenomeAltered,CBS-method}

\alias{CBS.getFGL}
\alias{getFGL.CBS}
\alias{getFGL,CBS-method}

\alias{CBS.getFGG}
\alias{getFGG.CBS}
\alias{getFGG,CBS-method}

\alias{CBS.getFGA}
\alias{getFGA.CBS}
\alias{getFGA,CBS-method}


\title{Calculates the fraction of the genome lost, gained, or aberrant either way}

\description{
 Calculates the fraction of the genome lost, gained, or aberrant either way (in sense of total copy numbers),
 using definitions closely related to those presented in [1].
}

\usage{
\method{getFractionOfGenomeLost}{CBS}(fit, ...)
}

\arguments{
 \item{...}{Not used.}
}

\value{
 Returns a \code{\link[base]{double}} in [0,1].
}

\author{Henrik Bengtsson}

\references{
  [1] Fridlyand et al. \emph{Breast tumor copy number aberration
      phenotypes and genomic instability}, BMC Cancer, 2006. \cr
}

\seealso{
  Internally, \code{\link[PSCBS:getCallStatistics.CBS]{*getCallStatistics}()} is used.
  For more information see \code{\link{CBS}}.
}


\keyword{internal}
\keyword{methods}
