%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  PSCBS.IO.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 \name{writeSegments.PSCBS}
\alias{writeSegments.PSCBS}
\alias{PSCBS.writeSegments}
\alias{writeSegments,PSCBS-method}

 \title{Writes the table of segments to file}

 \description{
  Writes the table of segments to file.
 }

 \usage{
\method{writeSegments}{PSCBS}(fit, name=getSampleName(fit), tags=NULL, ext="tsv", path=NULL, addHeader=TRUE, createdBy=NULL, sep="	", nbrOfDecimals=4, splitters=FALSE, overwrite=FALSE, skip=FALSE, ...)
}

 \arguments{
   \item{name, tags}{Name and optional tags part of the filename}.
   \item{path}{The directory where the file will be written.}
   \item{addHeader}{If \code{\link[base:logical]{TRUE}}, header comments are written.}
   \item{createdBy}{A header comment of whom created the file.}
   \item{splitters}{If \code{\link[base:logical]{TRUE}}, each chromosome is separated by a row
     of missing values.}
   \item{overwrite, skip}{If an output file already exists, these
     arguments specifies what should happen.}
   \item{...}{Additional arguments pass to \code{getSegments()}.}
 }

 \value{
   Returns the pathname of the the file written.
 }

 \author{Henrik Bengtsson}

 \seealso{
   Utilizes \code{\link[PSCBS:getSegments.PSCBS]{*getSegments}()}.
   For more information see \code{\link{PSCBS}}..
 }


\keyword{internal}
\keyword{methods}
