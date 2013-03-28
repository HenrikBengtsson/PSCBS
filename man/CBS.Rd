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

 \usage{
CBS(...)
}

 \arguments{
  \item{...}{Arguments passed to the constructor of \code{\link{AbstractCBS}}.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{all.equal} \tab  -\cr
 \tab \code{\link[PSCBS:append.CBS]{append}} \tab Appends one segmentation result to another.\cr
 \tab \code{as} \tab  -\cr
 \tab \code{as.character} \tab  -\cr
 \tab \code{\link[PSCBS:as.data.frame.CBS]{as.data.frame}} \tab Gets the table of segments.\cr
 \tab \code{\link[PSCBS:as.DNAcopy.CBS]{as.DNAcopy}} \tab Coerces a CBS object to a DNAcopy object.\cr
 \tab \code{\link[PSCBS:callAmplifications.CBS]{callAmplifications}} \tab Calls (focal) amplifications.\cr
 \tab \code{callArms} \tab  -\cr
 \tab \code{\link[PSCBS:callGainsAndLosses.CBS]{callGainsAndLosses}} \tab Calls gains and losses.\cr
 \tab \code{\link[PSCBS:callOutliers.CBS]{callOutliers}} \tab Calls outliers.\cr
 \tab \code{drawCentromeres} \tab  -\cr
 \tab \code{drawChangePoints} \tab  -\cr
 \tab \code{drawChromosomes} \tab  -\cr
 \tab \code{drawLevels} \tab  -\cr
 \tab \code{\link[PSCBS:estimateStandardDeviation.CBS]{estimateStandardDeviation}} \tab Estimates the whole-genome standard deviation of the signals.\cr
 \tab \code{extractCallsByLocus} \tab  -\cr
 \tab \code{extractChromosomes} \tab  -\cr
 \tab \code{extractCNs} \tab  -\cr
 \tab \code{\link[PSCBS:extractSegmentMeansByLocus.CBS]{extractSegmentMeansByLocus}} \tab Extracts segments means at each locus.\cr
 \tab \code{extractSegments} \tab  -\cr
 \tab \code{extractTotalCNs} \tab  -\cr
 \tab \code{\link[PSCBS:getCallStatistics.CBS]{getCallStatistics}} \tab Calculates various call statistics per chromosome.\cr
 \tab \code{getCallStatisticsByArms} \tab  -\cr
 \tab \code{getChromosomeRanges} \tab  -\cr
 \tab \code{getFGA} \tab  -\cr
 \tab \code{getFGG} \tab  -\cr
 \tab \code{getFGL} \tab  -\cr
 \tab \code{getFractionOfGenomeAltered} \tab  -\cr
 \tab \code{getFractionOfGenomeGained} \tab  -\cr
 \tab \code{\link[PSCBS:getFractionOfGenomeLost.CBS]{getFractionOfGenomeLost}} \tab Calculates the fraction of the genome lost, gained, or aberrant either way.\cr
 \tab \code{getLocusData} \tab  -\cr
 \tab \code{getLocusSignalNames} \tab  -\cr
 \tab \code{getSegments} \tab  -\cr
 \tab \code{getSegmentTrackPrefixes} \tab  -\cr
 \tab \code{getSignalType} \tab  -\cr
 \tab \code{highlightArmCalls} \tab  -\cr
 \tab \code{highlightCalls} \tab  -\cr
 \tab \code{highlightLocusCalls} \tab  -\cr
 \tab \code{isSegmentSplitter} \tab  -\cr
 \tab \code{isWholeChromosomeGained} \tab  -\cr
 \tab \code{isWholeChromosomeLost} \tab  -\cr
 \tab \code{joinSegments} \tab  -\cr
 \tab \code{\link[PSCBS:mergeNonCalledSegments.CBS]{mergeNonCalledSegments}} \tab Merge neighboring segments that are not called.\cr
 \tab \code{mergeTwoSegments} \tab  -\cr
 \tab \code{nbrOfAmplifications} \tab  -\cr
 \tab \code{nbrOfGains} \tab  -\cr
 \tab \code{nbrOfLosses} \tab  -\cr
 \tab \code{plot} \tab  -\cr
 \tab \code{\link[PSCBS:plotTracks.CBS]{plotTracks}} \tab Plots copy numbers along the genome.\cr
 \tab \code{plotTracksManyChromosomes} \tab  -\cr
 \tab \code{\link[PSCBS:pruneBySdUndo.CBS]{pruneBySdUndo}} \tab Prune the CBS profile by dropping change points that are too small.\cr
 \tab \code{resegment} \tab  -\cr
 \tab \code{segmentByCBS} \tab  -\cr
 \tab \code{seqOfSegmentsByDP} \tab  -\cr
 \tab \code{shiftTCN} \tab  -\cr
 \tab \code{signalType} \tab  -\cr
 \tab \code{signalType<-} \tab  -\cr
 \tab \code{subset} \tab  -\cr
 \tab \code{tileChromosomes} \tab  -\cr
 \tab \code{updateBoundaries} \tab  -\cr
 \tab \code{updateMeans} \tab  -\cr
 \tab \code{updateMeansTogether} \tab  -\cr
 \tab \code{writeLocusData} \tab  -\cr
 \tab \code{\link[PSCBS:writeSegments.CBS]{writeSegments}} \tab Writes the table of segments to file.\cr
}


 \bold{Methods inherited from AbstractCBS}:\cr
all.equal, append, as.data.frame, drawChangePoints, drawKnownSegments, dropChangePoint, dropChangePoints, dropRegion, dropRegions, extractChromosome, extractChromosomes, extractCNs, extractRegions, extractSegments, getChromosomeOffsets, getChromosomeRanges, getChromosomes, getLocusData, getLocusSignalNames, getMeanEstimators, getSampleName, getSegments, getSegmentSizes, getSegmentTrackPrefixes, load, mergeThreeSegments, mergeTwoSegments, nbrOfChangePoints, nbrOfChromosomes, nbrOfLoci, nbrOfSegments, plotTracks, print, pruneByDP, pruneByHClust, renameChromosomes, report, resegment, resetSegments, sampleCNs, sampleName, sampleName<-, save, seqOfSegmentsByDP, setLocusData, setMeanEstimators, setSampleName, setSegments, shiftTCN, tileChromosomes, updateMeans

 \bold{Methods inherited from list}:\cr
all.equal, as.data.frame, attachLocally, callHooks, relist, within

 "public"
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
