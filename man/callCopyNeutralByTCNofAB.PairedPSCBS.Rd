%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  PairedPSCBS.callCopyNeutral.R
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

\usage{\method{callCopyNeutralByTCNofAB}{PairedPSCBS}(fit, delta=estimateDeltaCN(fit), alpha=0.05, ..., force=FALSE, verbose=FALSE)}

\arguments{
  \item{fit}{A PairedPSCBS fit object as returned by
    \code{\link[PSCBS]{segmentByPairedPSCBS}}.}
  \item{delta}{A non-negative \code{\link[base]{double}} specifying the width of the
    "acceptance" region.
    Defaults to half of the distance between two integer TCN states,
    i.e. 1/2.  This argument should be shrunken as a function of
    the amount of the normal contaminator.}
  \item{alpha}{A \code{\link[base]{double}} in [0,0.5] specifying the significance level
    of the confidence intervals used.}
  \item{...}{Additional arguments passed to
             \code{\link[PSCBS:calcStatsForCopyNeutralABs.PairedPSCBS]{*calcStatsForCopyNeutralABs}()}.}
  \item{force}{If \code{\link[base:logical]{TRUE}}, an already called object is skipped, otherwise not.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
  Returns a \code{\link{PairedPSCBS}} fit object where a column
  with the copy-neutral call.
}

\details{
  ...
}

%% examples "../incl/callCopyNeutralByTCNofAB.PairedPSCBS.Rex"

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}


\keyword{internal}
\keyword{methods}
