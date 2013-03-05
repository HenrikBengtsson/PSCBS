%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  PairedPSCBS.estimateKappa.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{estimateKappa.PairedPSCBS}
\alias{estimateKappa.PairedPSCBS}
\alias{PairedPSCBS.estimateKappa}
\alias{estimateKappa,PairedPSCBS-method}

\title{Estimate global background in segmented copy numbers}

\description{
 Estimate global background in segmented copy numbers.
 The global background, here called \eqn{\kappa},
 may have multiple origins where normal contamination is one,
 but not necessarily the only one.
}

\usage{
\method{estimateKappa}{PairedPSCBS}(this, flavor=c("density(C1)"), ...)
}

\arguments{
  \item{flavor}{A \code{\link[base]{character}} string specifying which type of
   estimator to use.}
  \item{...}{Additional arguments passed to the estimator.}
}

\value{
  Returns the background estimate as a \code{\link[base]{numeric}} scalar.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  Internally, one of the following methods are used:
  \code{\link[PSCBS:estimateKappaByC1Density.PairedPSCBS]{*estimateKappaByC1Density}()}.
}
\keyword{internal}
\keyword{methods}
