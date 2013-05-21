%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  weightedQuantile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{weightedQuantile}
\alias{weightedQuantile.default}
\alias{weightedQuantile}

\title{Weighted Quantile Value}

\usage{
\method{weightedQuantile}{default}(x, w, probs=c(0, 0.25, 0.5, 0.75, 1), na.rm=TRUE,
  method=c("wtd.quantile"), ...)
}

\description{
  Computes a weighted quantile of a numeric vector.
}

\arguments{
  \item{x}{a \code{\link[base]{numeric}} \code{\link[base]{vector}} containing the values whose weighted
           quantile is to be computed.}
  \item{w}{a numeric \code{\link[base]{vector}} of weights the same length as
           \code{x} giving the weights to use for each element of \code{x}.
           Negative weights are treated as zero weights.
           Default value is equal weight to all values.}
  \item{probs}{a \code{\link[base]{numeric}} \code{\link[base]{vector}} of quantiles in [0,1] to be retrieved.}
  \item{na.rm}{a \code{\link[base]{logical}} value indicating whether \code{\link[base]{NA}} values in
           \code{x} should be stripped before the computation proceeds,
           or not.}
  \item{method}{If \code{"wtd.quantile"}, then \code{\link[Hmisc]{wtd.quantile}}
           of the \pkg{Hmisc} package is used.
           No other methods are currently supported.}
  \item{...}{Additional arguments passed to the estimator.}
}

\value{
  Returns the weighted quantile.
}

\author{Henrik Bengtsson}

\seealso{
  Internally the following functions may be used:
  \code{\link[stats]{quantile}} (if no weights are specified), or
  \code{\link[Hmisc]{wtd.quantile}}.
  For a weighted median estimator, \code{\link[matrixStats]{weightedMedian}}
  of the \pkg{matrixStats} package.
}




\keyword{univar}
\keyword{robust}
\keyword{internal}
