%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AbstractCBS.HCLUST.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{pruneByHClust.AbstractCBS}
\alias{pruneByHClust.AbstractCBS}
\alias{AbstractCBS.pruneByHClust}
\alias{pruneByHClust,AbstractCBS-method}

\title{Prunes the CN profile by pruning and merging through hierarchical clustering}

\description{
 Prunes the CN profile by pruning and merging through hierarchical clustering.
}

\usage{
\method{pruneByHClust}{AbstractCBS}(fit, ..., size=NULL, distMethod="euclidean", hclustMethod="ward", merge=TRUE, update=TRUE, verbose=FALSE)
}

\arguments{
 \item{...}{Arguments passed to \code{\link[stats]{cutree}},
   particularly either of thresholds \code{h} or \code{k}.}
 \item{size, distMethod, hclustMethod}{Arguments (as well as
   some of \code{...}) passed to \code{\link[PSCBS:hclustCNs.AbstractCBS]{*hclustCNs}()}.}
 \item{merge}{If \code{\link[base:logical]{TRUE}}, consecutive segments that belong to the
   same PSCN cluster will be merged into one large segment.}
 \item{update}{If \code{\link[base:logical]{TRUE}}, segment means are updated afterwards, otherwise not.}
 \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
  Returns a pruned object of the same class.
}

\examples{\dontrun{
 fitP <- pruneByHClust(fit, h=0.25);
}}

\author{Henrik Bengtsson}


\keyword{internal}
\keyword{methods}
