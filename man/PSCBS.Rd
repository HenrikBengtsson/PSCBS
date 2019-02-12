%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  PSCBS.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{PSCBS}
\docType{class}
\alias{PSCBS}


\title{The PSCBS class}

\description{
 Package:  PSCBS \cr
\bold{Class PSCBS}\cr

\code{list}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[PSCBS]{AbstractCBS}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\emph{\code{PSCBS}}\cr

\bold{Directly known subclasses:}\cr
\emph{\link[PSCBS]{NonPairedPSCBS}}, \emph{\link[PSCBS]{PairedPSCBS}}\cr

public abstract static class \bold{PSCBS}\cr
extends \emph{\link[PSCBS]{AbstractCBS}}\cr



 A PSCBS is an object containing results from parent-specific copy-number
 (PSCN) segmentation.
}

\usage{PSCBS(fit=list(), ...)}

\arguments{
  \item{fit}{A \code{\link[base]{list}} structure containing the PSCN segmentation results.}
  \item{...}{Not used.}
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{c} \tab  -\cr
 \tab \code{isLocallyPhased} \tab  -\cr
 \tab \code{normalizeTotalCNs} \tab  -\cr
 \tab \code{\link[PSCBS:writeSegments.PSCBS]{writeSegments}} \tab Writes the table of segments to file.\cr
}


 \bold{Methods inherited from AbstractCBS}:\cr
adjustPloidyScale, all.equal, append, as.data.frame, clearCalls, drawChangePoints, drawKnownSegments, dropChangePoint, dropChangePoints, dropRegion, dropRegions, extractCNs, extractChromosome, extractChromosomes, extractRegions, extractSegments, extractWIG, getChangePoints, getChromosomeOffsets, getChromosomeRanges, getChromosomes, getLocusData, getLocusSignalNames, getMeanEstimators, getSampleName, getSegmentSizes, getSegmentTrackPrefixes, getSegments, load, mergeThreeSegments, mergeTwoSegments, nbrOfChangePoints, nbrOfChromosomes, nbrOfLoci, nbrOfSegments, normalizeTotalCNs, ploidy, ploidy<-, plotTracks, print, pruneByDP, pruneByHClust, renameChromosomes, report, resegment, resetSegments, sampleCNs, sampleName, sampleName<-, save, seqOfSegmentsByDP, setLocusData, setMeanEstimators, setPloidy, setSampleName, setSegments, shiftTCN, tileChromosomes, updateMeans, writeWIG

 \bold{Methods inherited from list}:\cr
Ops,nonStructure,vector-method, Ops,structure,vector-method, Ops,vector,nonStructure-method, Ops,vector,structure-method, all.equal, as.data.frame, attachLocally, averageQuantile, callHooks, coerce,ANY,list-method, coerce,Matrix,vector-method, coerce,abIndex,vector-method, coerce,dgCMatrix,vector-method, coerce,dsCMatrix,vector-method, coerce,list,indMatrix-method, coerce,sparseVector,vector-method, contents, ggplot_add, latex, normalizeAverage, normalizeDifferencesToAverage, normalizeQuantileRank, normalizeQuantileSpline, plotDensity, relist, scale_type, type.convert, within, wrap


}

\author{Henrik Bengtsson}

\seealso{
  \code{\link{PairedPSCBS}}.
}
\keyword{classes}
