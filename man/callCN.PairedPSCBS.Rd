%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  PairedPSCBS.callCN.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{callCN.PairedPSCBS}
\alias{callCN.PairedPSCBS}
\alias{PairedPSCBS.callCN}
\alias{callCN.PairedPSCBS}
\alias{callCN,PairedPSCBS-method}

\alias{PairedPSCBS.callCopyNeutral}
\alias{callCopyNeutral.PairedPSCBS}
\alias{callCopyNeutral,PairedPSCBS-method}


\title{Calls segments that are copy neutral}

\description{
 Calls segments that are copy neutral, i.e. that have a total copy number that corresponds
 to the ploidy of the genome.
}

\usage{\method{callCN}{PairedPSCBS}(fit, flavor=c("TCN|AB"), ..., minSize=1, force=FALSE)}

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
  Returns a \code{\link{PairedPSCBS}} object with copy-neutral calls.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  Internally, one of the following methods are used:
  \code{\link[PSCBS:callCopyNeutralByTCNofAB.PairedPSCBS]{*callCopyNeutralByTCNofAB}()}.
}
\keyword{internal}
\keyword{methods}
