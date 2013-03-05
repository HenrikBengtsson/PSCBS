%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AbstractCBS.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{as.data.frame.AbstractCBS}
\alias{as.data.frame.AbstractCBS}
\alias{AbstractCBS.as.data.frame}
\alias{as.data.frame,AbstractCBS-method}

 \title{Gets the table of segments}

 \description{
  Gets the table of segments.
 }

 \usage{
\method{as.data.frame}{AbstractCBS}(x, ...)
}

 \arguments{
   \item{...}{Not used.}
 }

 \value{
   Returns a \code{\link[base]{data.frame}}, where each row corresponds to
   a unique segment.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
   Utilizes \code{\link[PSCBS:getSegments.AbstractCBS]{*getSegments}()}.
   For more information see \code{\link{AbstractCBS}}..
 }
\keyword{internal}
\keyword{methods}
