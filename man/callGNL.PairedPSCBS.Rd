%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  PairedPSCBS.callGNL.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{callGNL.PairedPSCBS}
\alias{callGNL.PairedPSCBS}
\alias{PairedPSCBS.callGNL}
\alias{callGNL,PairedPSCBS-method}
\alias{PairedPSCBS.callGainNeutralLoss}
\alias{callGainNeutralLoss.PairedPSCBS}
\alias{callGainNeutralLoss,PairedPSCBS-method}


\title{Calls segments that are gained, copy neutral, or lost}

\description{
 Calls segments that are gained, copy neutral, or lost, where copy neutral means having a total copy number
 that corresponds to the ploidy of the genome.
}

\usage{
\method{callGNL}{PairedPSCBS}(fit, flavor=c("TCN|AB"), ..., minSize=1, force=FALSE)
}

\arguments{
  \item{flavor}{A \code{\link[base]{character}} string specifying which type of
   call to use.}
  \item{...}{Additional arguments passed to the caller.}
  \item{minSize}{An optional \code{\link[base]{integer}} specifying the minimum number
   of data points in order to call a segments.  If fewer data points,
   then the call is set to \code{\link[base]{NA}} regardless.}
  \item{force}{If \code{\link[base:logical]{FALSE}}, and copy-neutral calls already exits,
   then nothing is done, otherwise the calls are done.}
}

\value{
  Returns a \code{\link{PairedPSCBS}} object with added calls.
}

\author{Henrik Bengtsson}

\seealso{
  Internally, one of the following methods are used:
  \code{callGNLByTCNofAB()}.
}
\keyword{internal}
\keyword{methods}
