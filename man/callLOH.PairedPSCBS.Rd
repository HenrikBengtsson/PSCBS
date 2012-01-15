%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  PairedPSCBS.callLOH.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{callLOH.PairedPSCBS}
\alias{callLOH.PairedPSCBS}
\alias{PairedPSCBS.callLOH}
\alias{callLOH.PairedPSCBS}
\alias{callLOH,PairedPSCBS-method}


\title{Calls segments that are in LOH}

\description{
 Calls segments that are in LOH, i.e. that have "zero" minor copy number.
}

\usage{\method{callLOH}{PairedPSCBS}(fit, flavor=c("SmallC1", "LargeDH"), ..., minSize=1, xorCalls=TRUE, force=FALSE)}

\arguments{
  \item{flavor}{A \code{\link[base]{character}} string specifying which type of
   call to use.}
  \item{...}{Additional arguments passed to the caller.}
  \item{minSize}{An optional \code{\link[base]{integer}} specifying the minimum number
   of data points in order to call a segments.  If fewer data points,
   then the call is set to \code{\link[base]{NA}} regardless.}
  \item{xorCalls}{If \code{\link[base:logical]{TRUE}}, a region already called AB, will
   for consistency never be called LOH, resulting in either an LOH
   call set to \code{\link[base:logical]{FALSE}} or \code{\link[base]{NA}} (as explained below).}
  \item{force}{If \code{\link[base:logical]{FALSE}}, and allelic-balance calls already exits,
   then nothing is done, otherwise the calls are done.}
}

\value{
  Returns a \code{\link{PairedPSCBS}} object with LOH calls.
}

\section{AB and LOH consistency}{
  Biologically, a segment can not be both in allelic balance (AB) and
  in loss-of-heterozygosity (LOH) at the same time.
  To avoid reporting such inconsistencies, the LOH caller will,
  if argument \code{xorCalls=TRUE}, never report a segment to be in
  LOH if it is already called to be in AB.
  However, regardless of of the AB call, a segment is still always
  tested for LOH, to check weather the LOH caller is consistent with the
  AB caller or not.  Thus, in order to distinguish the case where
  the AB caller and LOH caller agree from when they disagree,
  we report either (AB,LOH)=(TRUE,FALSE) or (TRUE,NA).  The former is
  reported when they are consistent, and the latter when they are not,
  or when the LOH caller could not call it.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  Internally, one of the following methods are used:
  \code{\link[PSCBS:callLowC1ByC1.PairedPSCBS]{*callLowC1ByC1}()},
  \code{\link[PSCBS:callExtremeAllelicImbalanceByDH.PairedPSCBS]{*callExtremeAllelicImbalanceByDH}()}.
}
\keyword{internal}
\keyword{methods}
