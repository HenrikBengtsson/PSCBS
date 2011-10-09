%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AbstractCBS.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{nbrOfChangePoints.AbstractCBS}
\alias{nbrOfChangePoints.AbstractCBS}
\alias{AbstractCBS.nbrOfChangePoints}
\alias{nbrOfChangePoints.AbstractCBS}
\alias{nbrOfChangePoints,AbstractCBS-method}


 \title{Gets the number of change points}

 \description{
   Gets the number of change points, which is defined as the number of segments minus
   the number of chromosomes.
 }

 \usage{\method{nbrOfChangePoints}{AbstractCBS}(fit, ...)}

 \arguments{
  \item{splitters, ...}{Arguments passed to \code{\link[PSCBS:getSegments.AbstractCBS]{*getSegments}()}.}
 }

 \value{
   Returns an \code{\link[base]{integer}}.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
   \code{\link[PSCBS:nbrOfSegments.AbstractCBS]{*nbrOfSegments}()}
   \code{\link[PSCBS:nbrOfChromosomes.AbstractCBS]{*nbrOfChromosomes}()}
   For more information see \code{\link{AbstractCBS}}.
 }
\keyword{internal}
\keyword{methods}
