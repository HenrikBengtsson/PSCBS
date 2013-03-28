%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AbstractCBS.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{AbstractCBS}
\docType{class}
\alias{AbstractCBS}


 \title{The AbstractCBS class}

 \description{
  Package:  PSCBS \cr
\bold{Class AbstractCBS}\cr

\code{list}\cr
\code{~~|}\cr
\code{~~+--}\emph{\code{AbstractCBS}}\cr

\bold{Directly known subclasses:}\cr
\emph{\link[PSCBS]{CBS}}, \emph{\link[PSCBS]{NonPairedPSCBS}}, \emph{\link[PSCBS]{PairedPSCBS}}, \emph{\link[PSCBS]{PSCBS}}\cr

public abstract static class \bold{AbstractCBS}\cr
extends list\cr



  All CBS-style segmentation results extend this class, e.g.
  \code{\link{CBS}} and \code{\link{PairedPSCBS}}.
 }

 \usage{
AbstractCBS(fit=list(), sampleName=fit$sampleName, ...)
}

 \arguments{
   \item{fit}{A \code{\link[base]{list}} structure containing the segmentation results.}
   \item{sampleName}{A \code{\link[base]{character}} string.}
   \item{...}{Not used.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{all.equal} \tab  -\cr
 \tab \code{\link[PSCBS:append.AbstractCBS]{append}} \tab Appends one segmentation result to another.\cr
 \tab \code{\link[PSCBS:as.data.frame.AbstractCBS]{as.data.frame}} \tab Gets the table of segments.\cr
 \tab \code{drawChangePoints} \tab  -\cr
 \tab \code{drawKnownSegments} \tab  -\cr
 \tab \code{dropChangePoint} \tab  -\cr
 \tab \code{\link[PSCBS:dropChangePoints.AbstractCBS]{dropChangePoints}} \tab Drops zero or more change points.\cr
 \tab \code{dropRegion} \tab  -\cr
 \tab \code{\link[PSCBS:dropRegions.AbstractCBS]{dropRegions}} \tab Drops chromosomal regions (a connected set of segments).\cr
 \tab \code{extractChromosome} \tab  -\cr
 \tab \code{extractChromosomes} \tab  -\cr
 \tab \code{extractCNs} \tab  -\cr
 \tab \code{extractRegion} \tab  -\cr
 \tab \code{extractRegions} \tab  -\cr
 \tab \code{extractSegment} \tab  -\cr
 \tab \code{extractSegments} \tab  -\cr
 \tab \code{getChromosomeOffsets} \tab  -\cr
 \tab \code{getChromosomeRanges} \tab  -\cr
 \tab \code{\link[PSCBS:getChromosomes.AbstractCBS]{getChromosomes}} \tab Gets the set of chromosomes.\cr
 \tab \code{\link[PSCBS:getLocusData.AbstractCBS]{getLocusData}} \tab Gets the locus-level data.\cr
 \tab \code{getLocusSignalNames} \tab  -\cr
 \tab \code{getMeanEstimators} \tab  -\cr
 \tab \code{\link[PSCBS:getSampleName.AbstractCBS]{getSampleName}} \tab Gets the name of the sample segmented.\cr
 \tab \code{\link[PSCBS:getSegments.AbstractCBS]{getSegments}} \tab Gets the segments.\cr
 \tab \code{getSegmentSizes} \tab  -\cr
 \tab \code{getSegmentTrackPrefixes} \tab  -\cr
 \tab \code{\link[PSCBS:hclustCNs.AbstractCBS]{hclustCNs}} \tab Performs a hierarchical clustering of the CN mean levels.\cr
 \tab \code{\link[PSCBS:load.AbstractCBS]{load}} \tab Loads an AbstractCBS object from file.\cr
 \tab \code{\link[PSCBS:mergeThreeSegments.AbstractCBS]{mergeThreeSegments}} \tab Merge a segment and its two flanking segments.\cr
 \tab \code{\link[PSCBS:mergeTwoSegments.AbstractCBS]{mergeTwoSegments}} \tab Merge two neighboring segments.\cr
 \tab \code{\link[PSCBS:nbrOfChangePoints.AbstractCBS]{nbrOfChangePoints}} \tab Gets the number of change points.\cr
 \tab \code{\link[PSCBS:nbrOfChromosomes.AbstractCBS]{nbrOfChromosomes}} \tab Gets the number of chromosomes.\cr
 \tab \code{\link[PSCBS:nbrOfLoci.AbstractCBS]{nbrOfLoci}} \tab Gets the number of loci.\cr
 \tab \code{\link[PSCBS:nbrOfSegments.AbstractCBS]{nbrOfSegments}} \tab Gets the number of segments.\cr
 \tab \code{\link[PSCBS:plotTracks.AbstractCBS]{plotTracks}} \tab Plots the segmentation result along the genome.\cr
 \tab \code{print} \tab  -\cr
 \tab \code{\link[PSCBS:pruneByDP.AbstractCBS]{pruneByDP}} \tab Prunes the CN profile using dynamical programming.\cr
 \tab \code{\link[PSCBS:pruneByHClust.AbstractCBS]{pruneByHClust}} \tab Prunes the CN profile by pruning and merging through hierarchical clustering.\cr
 \tab \code{renameChromosomes} \tab  -\cr
 \tab \code{\link[PSCBS:report.AbstractCBS]{report}} \tab Generates a report of the segmentation results.\cr
 \tab \code{resegment} \tab  -\cr
 \tab \code{\link[PSCBS:resetSegments.AbstractCBS]{resetSegments}} \tab Reset the segments.\cr
 \tab \code{sampleCNs} \tab  -\cr
 \tab \code{sampleName} \tab  -\cr
 \tab \code{sampleName<-} \tab  -\cr
 \tab \code{\link[PSCBS:save.AbstractCBS]{save}} \tab Saves an AbstractCBS object to file.\cr
 \tab \code{seqOfSegmentsByDP} \tab  -\cr
 \tab \code{setLocusData} \tab  -\cr
 \tab \code{setMeanEstimators} \tab  -\cr
 \tab \code{\link[PSCBS:setSampleName.AbstractCBS]{setSampleName}} \tab Sets the name of the sample segmented.\cr
 \tab \code{setSegments} \tab  -\cr
 \tab \code{shiftTCN} \tab  -\cr
 \tab \code{tileChromosomes} \tab  -\cr
 \tab \code{\link[PSCBS:updateMeans.AbstractCBS]{updateMeans}} \tab Updates the CN mean levels for each segment independently.\cr
 \tab \code{\link[PSCBS:updateMeansTogether.AbstractCBS]{updateMeansTogether}} \tab Updates the CN mean levels jointly in sets of segments.\cr
}


 \bold{Methods inherited from list}:\cr
all.equal, as.data.frame, attachLocally, callHooks, relist, within

 "public"
 }

 \author{Henrik Bengtsson}
\keyword{classes}
