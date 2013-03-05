%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  PairedPSCBS.estimateDeltaAB.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{estimateDeltaAB.PairedPSCBS}
\alias{estimateDeltaAB.PairedPSCBS}
\alias{PairedPSCBS.estimateDeltaAB}
\alias{estimateDeltaAB,PairedPSCBS-method}

\title{Estimate a threshold for calling allelic balance from DH}

\description{
 Estimate a threshold for calling allelic balance from DH to be used by the \code{\link[PSCBS:callAB.PairedPSCBS]{*callAB}()} method.
}

\usage{
\method{estimateDeltaAB}{PairedPSCBS}(this, scale=NULL, flavor=c("qq(DH)", "q(DH)", "mad(hBAF)", "median(DH)"), ..., max=Inf, verbose=FALSE)
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
  Returns the threshold estimate as a \code{\link[base]{numeric}} scalar.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  Internally, one of the following methods are used:
  \code{\link[PSCBS:estimateDeltaABBySmallDH.PairedPSCBS]{*estimateDeltaABBySmallDH}()},
  \code{\link[PSCBS:estimateStdDevForHeterozygousBAF.PairedPSCBS]{*estimateStdDevForHeterozygousBAF}()},
  \code{\link[PSCBS:estimateMeanForDH.PairedPSCBS]{*estimateMeanForDH}()}, and
  \code{\link[PSCBS:estimateHighDHQuantileAtAB.PairedPSCBS]{*estimateHighDHQuantileAtAB}()}.
}
\keyword{internal}
\keyword{methods}
