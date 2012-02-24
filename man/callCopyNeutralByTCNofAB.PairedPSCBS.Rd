%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  PairedPSCBS.callCN.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{callCopyNeutralByTCNofAB.PairedPSCBS}
\alias{callCopyNeutralByTCNofAB.PairedPSCBS}
\alias{PairedPSCBS.callCopyNeutralByTCNofAB}
\alias{callCopyNeutralByTCNofAB.PairedPSCBS}
\alias{callCopyNeutralByTCNofAB,PairedPSCBS-method}


\title{Calls regions that are copy neutral}

\description{
 Calls regions that are copy neutral from the total copy numbers (TCNs) of segments
 in allelic balance (AB).
}

\usage{\method{callCopyNeutralByTCNofAB}{PairedPSCBS}(fit, ..., force=FALSE, verbose=FALSE)}

\arguments{
  \item{fit}{A PairedPSCBS fit object as returned by
    \code{\link[PSCBS]{segmentByPairedPSCBS}}.}
  \item{...}{Additional arguments passed to
    \code{\link{findNeutralCopyNumberState}}().}
  \item{force}{If \code{\link[base:logical]{TRUE}}, an already called object is skipped, otherwise not.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
  Returns a \code{\link{PairedPSCBS}} fit object where a column
  with the copy-neutral call.
}

%% examples "../incl/callCopyNeutralByTCNofAB.PairedPSCBS.Rex"

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}


\keyword{internal}
\keyword{methods}
