%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AbstractCBS.RESTRUCT.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 \name{Restructuring AbstractCBS objects}
\alias{Restructuring AbstractCBS objects}
\title{Restructuring AbstractCBS objects}

 \alias{RestructPSCBS}

 \description{
   This page describes available methods for restructuring an
   \code{\link{AbstractCBS}} object.

   \itemize{
     \item \code{\link[PSCBS:append.AbstractCBS]{*append}()}
           - Appends one \code{\link{AbstractCBS}} to another.

     \item \code{\link[PSCBS:extractChromosomes.AbstractCBS]{*extractChromosomes}()} /
           \code{\link[PSCBS:extractChromosome.AbstractCBS]{*extractChromosome}()}
           - Extracts an \code{\link{AbstractCBS}} with the specified chromosomes.

     \item \code{\link[PSCBS:extractRegions.AbstractCBS]{*extractRegions}()} /
           \code{\link[PSCBS:extractRegion.AbstractCBS]{*extractRegion}()}
           - Extracts an \code{\link{AbstractCBS}} with the specified regions
             each of a certain size, where a region is defined as a
             connected set of segments.

     \item \code{\link[PSCBS:dropRegions.AbstractCBS]{*dropRegions}()} /
           \code{\link[PSCBS:dropRegion.AbstractCBS]{*dropRegion}()}
           - Drops specified regions and returns an \code{\link{AbstractCBS}}
             without them.

     \item \code{\link[PSCBS:dropChangePoint.AbstractCBS]{*dropChangePoint}()} /
           \code{\link[PSCBS:mergeTwoSegments.AbstractCBS]{*mergeTwoSegments}()}
           - Drops a change point by merging two neighboring segments
             and recalculating the statistics for the merged segment
             before returning an \code{\link{AbstractCBS}}.
   }

   All of the above methods are implemented for \code{\link{CBS}} and
   \code{\link{PairedPSCBS}} objects.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
   For more information see \code{\link{AbstractCBS}}.
 }


\keyword{documentation}
\keyword{internal}
