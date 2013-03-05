%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  PairedPSCBS.estimateDeltaLOH.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{estimateDeltaLOH.PairedPSCBS}
\alias{estimateDeltaLOH.PairedPSCBS}
\alias{PairedPSCBS.estimateDeltaLOH}
\alias{estimateDeltaLOH,PairedPSCBS-method}

\title{Estimate a threshold for calling LOH from DH}

\description{
 Estimate a threshold for calling LOH from DH to be used by the \code{\link[PSCBS:callLOH.PairedPSCBS]{*callLOH}()} method.
}

\usage{
\method{estimateDeltaLOH}{PairedPSCBS}(this, flavor=c("minC1|nonAB"), ..., max=Inf, verbose=FALSE)
}

\arguments{
  \item{scale}{An optional \code{\link[base]{numeric}} scale factor.}
  \item{flavor}{A \code{\link[base]{character}} string specifying which type of
   estimator to use.}
  \item{...}{Additional arguments passed to the estimator.}
  \item{max}{(Optional) The maxium estimate allowed. If greater than
   this value, the estimate will be truncated.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
  Returns the threshold estimate as a \code{\link[base]{numeric}} scalar or -\code{\link[base:is.finite]{Inf}}.
  In case it is not possible to estimate the LOH threshold, then
  -\code{\link[base:is.finite]{Inf}} is returned.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  Internally, one of the following methods are used:
  \code{\link[PSCBS:estimateDeltaLOHByMinC1ForNonAB.PairedPSCBS]{*estimateDeltaLOHByMinC1ForNonAB}()}.
}
\keyword{internal}
\keyword{methods}
