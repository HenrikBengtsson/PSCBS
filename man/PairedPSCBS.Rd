%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  PairedPSCBS.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{PairedPSCBS}
\docType{class}
\alias{PairedPSCBS}


 \title{The PairedPSCBS class}

 \description{
  Package:  PSCBS \cr
\bold{Class PairedPSCBS}\cr

\code{list}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[PSCBS]{AbstractCBS}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[PSCBS]{PSCBS}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\emph{\code{PairedPSCBS}}\cr

\bold{Directly known subclasses:}\cr
\cr

public abstract static class \bold{PairedPSCBS}\cr
extends \emph{\link[PSCBS]{PSCBS}}\cr



  A PairedPSCBS is an object containing the results from the
  Paired PSCBS method.
 }

 \usage{PairedPSCBS(fit=list(), ...)}

 \arguments{
   \item{fit}{A \code{\link[base]{list}} structure containing the Paired PSCBS results.}
   \item{...}{Not used.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{applyByRegion} \tab  -\cr
 \tab \code{arrowsC1C2} \tab  -\cr
 \tab \code{arrowsDeltaC1C2} \tab  -\cr
 \tab \code{bootstrapCIs} \tab  -\cr
 \tab \code{bootstrapTCNandDHByRegion} \tab  -\cr
 \tab \code{calcStatsForCopyNeutralABs} \tab  -\cr
 \tab \code{\link[PSCBS:callAB.PairedPSCBS]{callAB}} \tab Calls segments that are in allelic balance.\cr
 \tab \code{callABandHighAI} \tab  -\cr
 \tab \code{callABandLowC1} \tab  -\cr
 \tab \code{\link[PSCBS:callAllelicBalanceByDH.PairedPSCBS]{callAllelicBalanceByDH}} \tab Calls segments that are in allelic balance.\cr
 \tab \code{\link[PSCBS:callCopyNeutral.PairedPSCBS]{callCopyNeutral}} \tab Calls segments that have a neutral total copy number.\cr
 \tab \code{\link[PSCBS:callCopyNeutralByTCNofAB.PairedPSCBS]{callCopyNeutralByTCNofAB}} \tab Calls regions that are copy neutral.\cr
 \tab \code{callExtremeAllelicImbalanceByDH} \tab  -\cr
 \tab \code{callGainNeutralLoss} \tab  -\cr
 \tab \code{\link[PSCBS:callGNL.PairedPSCBS]{callGNL}} \tab Calls segments that are gained, copy neutral, or lost.\cr
 \tab \code{callGNLByTCNofAB} \tab  -\cr
 \tab \code{callGNLByTCNofABv1} \tab  -\cr
 \tab \code{\link[PSCBS:callLOH.PairedPSCBS]{callLOH}} \tab Calls segments that are in LOH.\cr
 \tab \code{callLowC1ByC1} \tab  -\cr
 \tab \code{callNTCN} \tab  -\cr
 \tab \code{\link[PSCBS:callROH.PairedPSCBS]{callROH}} \tab Calls segments that are in ROH.\cr
 \tab \code{callROHOneSegment} \tab  -\cr
 \tab \code{drawChangePointsC1C2} \tab  -\cr
 \tab \code{drawConfidenceBands} \tab  -\cr
 \tab \code{drawLevels} \tab  -\cr
 \tab \code{\link[PSCBS:estimateDeltaAB.PairedPSCBS]{estimateDeltaAB}} \tab Estimate a threshold for calling allelic balance from DH.\cr
 \tab \code{\link[PSCBS:estimateDeltaABBySmallDH.PairedPSCBS]{estimateDeltaABBySmallDH}} \tab Estimate a threshold for calling allelic balance from DH.\cr
 \tab \code{estimateDeltaCN} \tab  -\cr
 \tab \code{\link[PSCBS:estimateDeltaLOH.PairedPSCBS]{estimateDeltaLOH}} \tab Estimate a threshold for calling LOH from DH.\cr
 \tab \code{\link[PSCBS:estimateDeltaLOHByMinC1ForNonAB.PairedPSCBS]{estimateDeltaLOHByMinC1ForNonAB}} \tab Estimate a threshold for calling LOH from DH.\cr
 \tab \code{estimateHighDHQuantileAtAB} \tab  -\cr
 \tab \code{\link[PSCBS:estimateKappa.PairedPSCBS]{estimateKappa}} \tab Estimate global background in segmented copy numbers.\cr
 \tab \code{\link[PSCBS:estimateKappaByC1Density.PairedPSCBS]{estimateKappaByC1Density}} \tab Estimate global background in segmented copy numbers.\cr
 \tab \code{estimateMeanForDH} \tab  -\cr
 \tab \code{estimateStdDevForHeterozygousBAF} \tab  -\cr
 \tab \code{extractC1C2} \tab  -\cr
 \tab \code{extractCallsByLocus} \tab  -\cr
 \tab \code{extractCNs} \tab  -\cr
 \tab \code{extractDeltaC1C2} \tab  -\cr
 \tab \code{extractDhSegment} \tab  -\cr
 \tab \code{extractLocusLevelC1C2} \tab  -\cr
 \tab \code{extractLocusLevelTCN} \tab  -\cr
 \tab \code{\link[PSCBS:extractMinorMajorCNs.PairedPSCBS]{extractMinorMajorCNs}} \tab Extract minor and major copy-number mean levels per segment.\cr
 \tab \code{extractSegments} \tab  -\cr
 \tab \code{\link[PSCBS:extractTCNAndDHs.PairedPSCBS]{extractTCNAndDHs}} \tab Extract TCN and DH mean levels per segment.\cr
 \tab \code{getChromosomeOffsets} \tab  -\cr
 \tab \code{getChromosomeRanges} \tab  -\cr
 \tab \code{getLocusData} \tab  -\cr
 \tab \code{linesC1C2} \tab  -\cr
 \tab \code{linesDeltaC1C2} \tab  -\cr
 \tab \code{\link[PSCBS:mergeTwoSegments.PairedPSCBS]{mergeTwoSegments}} \tab Merge two neighboring segments.\cr
 \tab \code{plot} \tab  -\cr
 \tab \code{plotC1C2} \tab  -\cr
 \tab \code{plotDeltaC1C2} \tab  -\cr
 \tab \code{\link[PSCBS:plotTracks.PairedPSCBS]{plotTracks}} \tab Plots parental specific copy numbers along the genome.\cr
 \tab \code{plotTracks1} \tab  -\cr
 \tab \code{plotTracks2} \tab  -\cr
 \tab \code{plotTracksManyChromosomes} \tab  -\cr
 \tab \code{pointsC1C2} \tab  -\cr
 \tab \code{pointsDeltaC1C2} \tab  -\cr
 \tab \code{postsegmentTCN} \tab  -\cr
 \tab \code{resegment} \tab  -\cr
 \tab \code{segmentByNonPairedPSCBS} \tab  -\cr
 \tab \code{segmentByPairedPSCBS} \tab  -\cr
 \tab \code{seqOfSegmentsByDP} \tab  -\cr
 \tab \code{shiftTCN} \tab  -\cr
 \tab \code{tileChromosomes} \tab  -\cr
 \tab \code{updateMeans} \tab  -\cr
 \tab \code{updateMeansTogether} \tab  -\cr
}


 \bold{Methods inherited from PSCBS}:\cr
append, as.data.frame, drawChangePoints, extractChromosomes, getLocusData, getLocusSignalNames, getSegmentTrackPrefixes, isSegmentSplitter, writeSegments

 \bold{Methods inherited from AbstractCBS}:\cr
all.equal, append, as.data.frame, drawChangePoints, drawKnownSegments, dropChangePoint, dropChangePoints, dropRegion, dropRegions, extractChromosome, extractChromosomes, extractCNs, extractRegions, extractSegments, getChromosomeOffsets, getChromosomeRanges, getChromosomes, getLocusData, getLocusSignalNames, getMeanEstimators, getSampleName, getSegments, getSegmentSizes, getSegmentTrackPrefixes, load, mergeThreeSegments, mergeTwoSegments, nbrOfChangePoints, nbrOfChromosomes, nbrOfLoci, nbrOfSegments, plotTracks, print, pruneByDP, pruneByHClust, renameChromosomes, report, resegment, resetSegments, sampleCNs, sampleName, sampleName<-, save, seqOfSegmentsByDP, setLocusData, setMeanEstimators, setSampleName, setSegments, shiftTCN, tileChromosomes, updateMeans

 \bold{Methods inherited from list}:\cr
all.equal, as.data.frame, attachLocally, callHooks, relist, within

 "public"
 }

 \author{Henrik Bengtsson}

 \seealso{
   The \code{\link{segmentByPairedPSCBS}}() method returns an object of this class.
 }
\keyword{classes}
