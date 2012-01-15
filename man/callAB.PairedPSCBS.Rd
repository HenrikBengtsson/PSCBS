%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  PairedPSCBS.callAB.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{callAB.PairedPSCBS}
\alias{callAB.PairedPSCBS}
\alias{PairedPSCBS.callAB}
\alias{callAB.PairedPSCBS}
\alias{callAB,PairedPSCBS-method}

\alias{PairedPSCBS.callAllelicBalance}
\alias{callAllelicBalance.PairedPSCBS}
\alias{callAllelicBalance,PairedPSCBS-method}


\title{Calls segments that are in allelic balance}

\description{
 Calls segments that are in allelic balance, i.e. that have equal minor and major copy numbers.
}

\usage{\method{callAB}{PairedPSCBS}(fit, flavor=c("DeltaAB*"), ..., minSize=1, xorCalls=TRUE, force=FALSE)}

\arguments{
  \item{flavor}{A \code{\link[base]{character}} string specifying which type of
   call to use.}
  \item{...}{Additional arguments passed to the caller.}
  \item{minSize}{An optional \code{\link[base]{integer}} specifying the minimum number
   of data points in order to call a segments.  If fewer data points,
   then the call is set to \code{\link[base]{NA}} regardless.}
  \item{xorCalls}{If \code{\link[base:logical]{TRUE}}, a region already called LOH, will
   for consistency never be called AB, resulting in either an AB
   call set to \code{\link[base:logical]{FALSE}} or \code{\link[base]{NA}} (as explained below).}
  \item{force}{If \code{\link[base:logical]{FALSE}}, and allelic-balance calls already exits,
   then nothing is done, otherwise the calls are done.}
}

\value{
  Returns a \code{\link{PairedPSCBS}} object with allelic-balance calls.
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
  or when the AB caller could not call it.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  Internally, one of the following methods are used:
  \code{\link[PSCBS:callAllelicBalanceByDH.PairedPSCBS]{*callAllelicBalanceByDH}()}.
}
\keyword{internal}
\keyword{methods}
