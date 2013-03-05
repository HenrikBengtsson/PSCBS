%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AbstractCBS.RESTRUCT.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{dropChangePoints.AbstractCBS}
\alias{dropChangePoints.AbstractCBS}
\alias{AbstractCBS.dropChangePoints}
\alias{dropChangePoints,AbstractCBS-method}

 \title{Drops zero or more change points}

 \description{
   Drops zero or more change points, which is done by dropping one change point at the
   time using \code{\link[PSCBS:dropChangePoint.AbstractCBS]{*dropChangePoint}()}
   and recalculating the segment statistics at the end.

   \emph{NOTE: This method only works if there is only one chromosome.}
 }

 \usage{
\method{dropChangePoints}{AbstractCBS}(fit, idxs, update=TRUE, ...)
}

 \arguments{
  \item{idxs}{An \code{\link[base]{integer}} \code{\link[base]{vector}} specifying the change points to be dropped.}
  \item{update}{If \code{\link[base:logical]{TRUE}}, segment statistics are updated.}
  \item{...}{Other arguments passed to \code{\link[PSCBS:dropChangePoint.AbstractCBS]{*dropChangePoint}()}
             and \code{\link[PSCBS:updateMeans.AbstractCBS]{*updateMeans}()}.}
 }

 \value{
   Returns an \code{\link{AbstractCBS}} of the same class with
   \code{length(idxs)} segments.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
   For more information see \code{\link{AbstractCBS}}.
 }
\keyword{internal}
\keyword{methods}
