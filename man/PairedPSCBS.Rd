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
 \tab \code{\link[PSCBS:callAB.PairedPSCBS]{callAB}} \tab Calls segments that are in allelic balance.\cr
 \tab \code{\link[PSCBS:callCopyNeutral.PairedPSCBS]{callCopyNeutral}} \tab Calls segments that have a neutral total copy number.\cr
 \tab \code{\link[PSCBS:callGNL.PairedPSCBS]{callGNL}} \tab Calls segments that are gained, copy neutral, or lost.\cr
 \tab \code{callGNLByTCNofAB} \tab  -\cr
 \tab \code{callGainNeutralLoss} \tab  -\cr
 \tab \code{\link[PSCBS:callLOH.PairedPSCBS]{callLOH}} \tab Calls segments that are in LOH.\cr
 \tab \code{callNTCN} \tab  -\cr
 \tab \code{\link[PSCBS:callROH.PairedPSCBS]{callROH}} \tab Calls segments that are in ROH.\cr
 \tab \code{\link[PSCBS:estimateDeltaAB.PairedPSCBS]{estimateDeltaAB}} \tab Estimate a threshold for calling allelic balance from DH.\cr
 \tab \code{\link[PSCBS:estimateDeltaLOH.PairedPSCBS]{estimateDeltaLOH}} \tab Estimate a threshold for calling LOH from DH.\cr
 \tab \code{\link[PSCBS:estimateKappa.PairedPSCBS]{estimateKappa}} \tab Estimate global background in segmented copy numbers.\cr
 \tab \code{extractCNs} \tab  -\cr
 \tab \code{hasBootstrapSummaries} \tab  -\cr
 \tab \code{\link[PSCBS:plotTracks.PairedPSCBS]{plotTracks}} \tab Plots parental specific copy numbers along the genome.\cr
 \tab \code{segmentByNonPairedPSCBS} \tab  -\cr
 \tab \code{segmentByPairedPSCBS} \tab  -\cr
 \tab \code{seqOfSegmentsByDP} \tab  -\cr
}


 \bold{Methods inherited from PSCBS}:\cr
as.data.frame, c, drawChangePoints, extractChromosomes, extractWIG, getLocusData, getLocusSignalNames, getSegmentTrackPrefixes, isLocallyPhased, isSegmentSplitter, normalizeTotalCNs, writeSegments

 \bold{Methods inherited from AbstractCBS}:\cr
adjustPloidyScale, all.equal, append, as.data.frame, clearCalls, drawChangePoints, drawKnownSegments, dropChangePoint, dropChangePoints, dropRegion, dropRegions, extractCNs, extractChromosome, extractChromosomes, extractRegions, extractSegments, extractWIG, getChangePoints, getChromosomeOffsets, getChromosomeRanges, getChromosomes, getLocusData, getLocusSignalNames, getMeanEstimators, getSampleName, getSegmentSizes, getSegmentTrackPrefixes, getSegments, load, mergeThreeSegments, mergeTwoSegments, nbrOfChangePoints, nbrOfChromosomes, nbrOfLoci, nbrOfSegments, normalizeTotalCNs, ploidy, ploidy<-, plotTracks, print, pruneByDP, pruneByHClust, renameChromosomes, report, resegment, resetSegments, sampleCNs, sampleName, sampleName<-, save, seqOfSegmentsByDP, setLocusData, setMeanEstimators, setPloidy, setSampleName, setSegments, shiftTCN, tileChromosomes, updateMeans, writeWIG

 \bold{Methods inherited from list}:\cr
Ops,nonStructure,vector-method, Ops,structure,vector-method, Ops,vector,nonStructure-method, Ops,vector,structure-method, all.equal, as.data.frame, attachLocally, averageQuantile, callHooks, coerce,ANY,list-method, coerce,Matrix,vector-method, coerce,abIndex,vector-method, coerce,dgCMatrix,vector-method, coerce,dsCMatrix,vector-method, coerce,list,indMatrix-method, coerce,sparseVector,vector-method, contents, ggplot_add, latex, normalizeAverage, normalizeDifferencesToAverage, normalizeQuantileRank, normalizeQuantileSpline, plotDensity, relist, scale_type, type.convert, within, wrap


}

\author{Henrik Bengtsson}

\seealso{
  The \code{\link{segmentByPairedPSCBS}}() method returns an object of this class.
}
\keyword{classes}
