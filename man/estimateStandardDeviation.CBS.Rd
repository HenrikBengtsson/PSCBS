%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  CBS.EXTS.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{estimateStandardDeviation.CBS}
\alias{estimateStandardDeviation.CBS}
\alias{CBS.estimateStandardDeviation}
\alias{estimateStandardDeviation.CBS}
\alias{estimateStandardDeviation,CBS-method}


\title{Estimates the whole-genome standard deviation of the signals}

\description{
 Estimates the whole-genome standard deviation of the signals.
}

\usage{\method{estimateStandardDeviation}{CBS}(fit, chromosomes=NULL, method=c("diff", "res", "abs"), estimator=c("mad", "sd"), na.rm=TRUE, weights=NULL, ...)}

\arguments{
 \item{chromosomes}{An optional \code{\link[base]{vector}} specifying the subset of
   chromosomes used for the estimate.  If \code{\link[base]{NULL}}, all chromosomes are used.}
 \item{method}{A \code{\link[base]{character}} string specifying the method used.}
 \item{estimator}{A \code{\link[base]{character}} string or a \code{\link[base]{function}} specifying the
   internal estimator.}
 \item{na.rm}{If \code{\link[base:logical]{TRUE}}, missing values are dropped, otherwise not.}
 \item{weights}{An optional \code{\link[base]{double}} \code{\link[base]{vector}} of \code{\link[PSCBS:nbrOfLoci.CBS]{*nbrOfLoci}()}
   non-negative weights.}
 \item{...}{Not used.}
}

\value{
 Returns a non-negative \code{\link[base]{numeric}} scale.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{CBS}}.
}


\keyword{internal}
\keyword{methods}
