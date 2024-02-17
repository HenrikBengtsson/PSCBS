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

public abstract class \bold{CBS}\cr
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
 \tab \code{as} \tab  -\cr
 \tab \code{c} \tab  -\cr
 \tab \code{estimateStandardDeviation} \tab  -\cr
 \tab \code{plotTracks} \tab  -\cr
 \tab \code{pruneBySdUndo} \tab  -\cr
 \tab \code{segmentByCBS} \tab  -\cr
 \tab \code{seqOfSegmentsByDP} \tab  -\cr
 \tab \code{writeSegments} \tab  -\cr
}


 \bold{Methods inherited from AbstractCBS}:\cr
adjustPloidyScale, all.equal, as.data.frame, clearCalls, drawChangePoints, drawKnownSegments, dropChangePoint, dropChangePoints, dropRegion, dropRegions, extractCNs, extractChromosome, extractChromosomes, extractRegions, extractSegments, extractWIG, getChangePoints, getChromosomeOffsets, getChromosomeRanges, getChromosomes, getLocusData, getLocusSignalNames, getMeanEstimators, getSampleName, getSegmentSizes, getSegmentTrackPrefixes, getSegments, mergeThreeSegments, mergeTwoSegments, nbrOfChangePoints, nbrOfChromosomes, nbrOfLoci, nbrOfSegments, normalizeTotalCNs, ploidy, ploidy<-, plotTracks, print, pruneByDP, pruneByHClust, renameChromosomes, report, resegment, resetSegments, sampleCNs, sampleName, sampleName<-, seqOfSegmentsByDP, setLocusData, setMeanEstimators, setPloidy, setSampleName, setSegments, shiftTCN, tileChromosomes, updateMeans, writeWIG

 \bold{Methods inherited from list}:\cr
Ops,nonStructure,vector-method, Ops,structure,vector-method, Ops,vector,nonStructure-method, Ops,vector,structure-method, all.equal, as.data.frame, attachLocally, callHooks, coerce,ANY,list-method, relist, type.convert, within


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
