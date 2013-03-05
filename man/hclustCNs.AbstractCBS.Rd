%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AbstractCBS.HCLUST.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{hclustCNs.AbstractCBS}
\alias{hclustCNs.AbstractCBS}
\alias{AbstractCBS.hclustCNs}
\alias{hclustCNs,AbstractCBS-method}

\title{Performs a hierarchical clustering of the CN mean levels}

\description{
 Performs a hierarchical clustering of the CN mean levels.
}

\usage{
\method{hclustCNs}{AbstractCBS}(fit, size=NULL, distMethod="euclidean", hclustMethod="ward", ..., verbose=FALSE)
}

\arguments{
 \item{size}{Argument passed to \code{\link[PSCBS:sampleCNs.AbstractCBS]{*sampleCNs}()}.}
 \item{distMethod, hclustMethod}{Argument \code{method} for
   \code{\link[stats]{dist}} and "stats::hclust", respectively.}
 \item{...}{Not used.}
 \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
  Returns a \code{hclust} object as returned by \code{\link[stats]{hclust}}.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  This method is utilized by \code{\link[PSCBS:pruneByHClust.AbstractCBS]{*pruneByHClust}()}.
}


\keyword{internal}
\keyword{methods}
