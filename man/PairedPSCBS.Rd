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

public abstract class \bold{PairedPSCBS}\cr
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
 \tab \code{\link[PSCBS:callLOH.PairedPSCBS]{callLOH}} \tab Calls segments that are in LOH.\cr
 \tab \code{\link[PSCBS:estimateDeltaAB.PairedPSCBS]{estimateDeltaAB}} \tab Estimate a threshold for calling allelic balance from DH.\cr
 \tab \code{\link[PSCBS:estimateDeltaLOH.PairedPSCBS]{estimateDeltaLOH}} \tab Estimate a threshold for calling LOH from DH.\cr
 \tab \code{\link[PSCBS:estimateKappa.PairedPSCBS]{estimateKappa}} \tab Estimate global background in segmented copy numbers.\cr
 \tab \code{\link[PSCBS:plotTracks.PairedPSCBS]{plotTracks}} \tab Plots parental specific copy numbers along the genome.\cr
}


 \bold{Methods inherited from PSCBS}:\cr
append, as.data.frame, extractChromosomes, getLocusData

 \bold{Methods inherited from AbstractCBS}:\cr
all.equal, append, as.data.frame, dropChangePoint, dropRegion, dropRegions, extractChromosome, extractChromosomes, extractRegions, getChromosomes, getLocusData, getSampleName, getSegments, mergeTwoSegments, nbrOfChangePoints, nbrOfChromosomes, nbrOfLoci, nbrOfSegments, plotTracks, print, sampleName, sampleName<-, setLocusData, setSampleName, setSegments, tileChromosomes, updateMeans

 \bold{Methods inherited from list}:\cr
all.equal, as.data.frame, attachLocally, callHooks, relist, within


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
   The \code{\link{segmentByPairedPSCBS}}() method returns an object of this class.
 }
\keyword{classes}
