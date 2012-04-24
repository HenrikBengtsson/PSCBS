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

 \usage{AbstractCBS(fit=list(), sampleName=fit$sampleName, ...)}

 \arguments{
   \item{fit}{A \code{\link[base]{list}} structure containing the segmentation results.}
   \item{sampleName}{A \code{\link[base]{character}} string.}
   \item{...}{Not used.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{\link[PSCBS:append.AbstractCBS]{append}} \tab Appends one segmentation result to another.\cr
 \tab \code{extractCNs} \tab  -\cr
 \tab \code{\link[PSCBS:getChromosomes.AbstractCBS]{getChromosomes}} \tab Gets the set of chromosomes.\cr
 \tab \code{\link[PSCBS:getLocusData.AbstractCBS]{getLocusData}} \tab Gets the locus-level data.\cr
 \tab \code{\link[PSCBS:getSegments.AbstractCBS]{getSegments}} \tab Gets the segments.\cr
 \tab \code{getSegmentSizes} \tab  -\cr
 \tab \code{\link[PSCBS:load.AbstractCBS]{load}} \tab Loads an AbstractCBS object from file.\cr
 \tab \code{\link[PSCBS:mergeThreeSegments.AbstractCBS]{mergeThreeSegments}} \tab Merge a segment and its two flanking segments.\cr
 \tab \code{\link[PSCBS:nbrOfChangePoints.AbstractCBS]{nbrOfChangePoints}} \tab Gets the number of change points.\cr
 \tab \code{\link[PSCBS:nbrOfChromosomes.AbstractCBS]{nbrOfChromosomes}} \tab Gets the number of chromosomes.\cr
 \tab \code{\link[PSCBS:nbrOfLoci.AbstractCBS]{nbrOfLoci}} \tab Gets the number of loci.\cr
 \tab \code{\link[PSCBS:nbrOfSegments.AbstractCBS]{nbrOfSegments}} \tab Gets the number of segments.\cr
 \tab \code{\link[PSCBS:plotTracks.AbstractCBS]{plotTracks}} \tab Plots the segmentation result along the genome.\cr
 \tab \code{sampleCNs} \tab  -\cr
 \tab \code{\link[PSCBS:save.AbstractCBS]{save}} \tab Saves an AbstractCBS object to file.\cr
}


 \bold{Methods inherited from list}:\cr
all.equal, as.data.frame, attachLocally, callHooks, relist, within


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
