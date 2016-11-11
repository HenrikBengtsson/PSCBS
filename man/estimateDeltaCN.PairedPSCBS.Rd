%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  PairedPSCBS.callCopyNeutral.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{estimateDeltaCN.PairedPSCBS}
\alias{estimateDeltaCN.PairedPSCBS}
\alias{PairedPSCBS.estimateDeltaCN}
\alias{estimateDeltaCN,PairedPSCBS-method}
\alias{estimateDeltaCN}
\alias{estimateDeltaCN.CBS}

\title{Estimates the length of one total copy-number (TCN) unit}

\description{
 Estimates the length of one total copy-number (TCN) unit
}

\usage{
\method{estimateDeltaCN}{PairedPSCBS}(fit, scale=1, kappa=estimateKappa(fit), ...)
}

\arguments{
 \item{scale}{A \code{\link[base]{numeric}} scale factor in (0,Inf) used for rescaling
  (multiplying) the final estimate with.}
 \item{...}{Not used.}
}

\value{
 Returns a positive scalar \code{\link[base]{numeric}}.
}

\details{
  For parent-specific copy-number (PSCN) data, the TCN unit length is
  estimated as \eqn{(1-kappa)/2}, where \eqn{kappa} is estimated from
  data (by \code{\link[PSCBS]{estimateKappa}}).

  For total copy-number (TCN) data (only),
}

\author{Henrik Bengtsson}

\seealso{
  For more information see \code{\link{PairedPSCBS}}.
}


\keyword{internal}
\keyword{methods}
