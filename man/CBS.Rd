%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  CBS.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{CBS}
\docType{class}
\alias{CBS}


 \title{The CBS class}

 \description{
   A CBS object holds results from the
   Circular Binary Segmentation (CBS) method
   for \emph{one} sample for one or more chromosomes.

  Package:  PSCBS \cr
\bold{Class CBS}\cr

\code{list}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[PSCBS]{AbstractCBS}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\emph{\code{CBS}}\cr

\bold{Directly known subclasses:}\cr
\cr

public abstract static class \bold{CBS}\cr
extends \emph{\link[PSCBS]{AbstractCBS}}\cr


 }

 \usage{CBS(...)}

 \arguments{
  \item{...}{Arguments passed to the constructor of \code{\link{AbstractCBS}}.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{\link[PSCBS:append.CBS]{append}} \tab Appends one segmentation result to another.\cr
 \tab \code{as} \tab  -\cr
 \tab \code{\link[PSCBS:estimateStandardDeviation.CBS]{estimateStandardDeviation}} \tab Estimates the whole-genome standard deviation of the signals.\cr
 \tab \code{getSegmentSizes} \tab  -\cr
 \tab \code{\link[PSCBS:plotTracks.CBS]{plotTracks}} \tab Plots copy numbers along the genome.\cr
 \tab \code{\link[PSCBS:pruneBySdUndo.CBS]{pruneBySdUndo}} \tab Prune the CBS profile by dropping change points that are too small.\cr
 \tab \code{segmentByCBS} \tab  -\cr
 \tab \code{\link[PSCBS:writeSegments.CBS]{writeSegments}} \tab Writes the table of segments to file.\cr
}


 \bold{Methods inherited from AbstractCBS}:\cr
all.equal, append, as.data.frame, dropChangePoint, dropRegion, dropRegions, extractChromosome, extractChromosomes, extractCNs, extractRegions, extractSegments, getChromosomes, getLocusData, getSampleName, getSegments, getSegmentSizes, load, mergeThreeSegments, mergeTwoSegments, nbrOfChangePoints, nbrOfChromosomes, nbrOfLoci, nbrOfSegments, plotTracks, print, pruneByHClust, resegment, sampleCNs, sampleName, sampleName<-, save, setLocusData, setSampleName, setSegments, tileChromosomes, updateMeans

 \bold{Methods inherited from list}:\cr
all.equal, as.data.frame, attachLocally, callHooks, relist, within


 }

 \section{Difference to DNAcopy object}{
   A CBS object is similar to DNAcopy objects with the major
   difference that a CBS object holds only one sample, whereas
   a DNAcopy object can hold more than one sample.
 }

 \section{See also}{
  The \code{\link{segmentByCBS}}() method returns an object of this class.
 }

 \author{Henrik Bengtsson}
\keyword{classes}
