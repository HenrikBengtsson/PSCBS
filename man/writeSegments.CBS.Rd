%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  CBS.IO.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 \name{writeSegments.CBS}
\alias{writeSegments.CBS}
\alias{CBS.writeSegments}
\alias{writeSegments.CBS}
\alias{writeSegments,CBS-method}


 \title{Writes the table of segments to file}

 \description{
  Writes the table of segments to file.
 }

 \usage{\method{writeSegments}{CBS}(fit, filename=sprintf("\%s.tsv", getSampleName(fit)), path=NULL, addHeader=TRUE, createdBy=NULL, sep="	", nbrOfDecimals=4, splitters=FALSE, overwrite=FALSE, skip=FALSE, ...)}

 \arguments{
   \item{filename, path}{The filename and the path of the file to be written.}
   \item{addHeader}{If \code{\link[base:logical]{TRUE}}, header comments are written.}
   \item{createdBy}{A header comment of whom created the file.}
   \item{splitters}{If \code{\link[base:logical]{TRUE}}, each chromosome is separated by a row
     of missing values.}
   \item{overwrite, skip}{If an output file already exists, these
     arguments specifies what should happen.}
   \item{...}{Not used.}
 }

 \value{
   Returns the pathname of the the file written.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
   Utilizes \code{\link[PSCBS:getSegments.CBS]{*getSegments}()}.
   For more information see \code{\link{CBS}}..
 }


\keyword{internal}
\keyword{methods}
