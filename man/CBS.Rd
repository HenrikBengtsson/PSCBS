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
\code{~~+--}\code{CBS}\cr

\bold{Directly known subclasses:}\cr
\cr

public class \bold{CBS}\cr
extends list\cr


 }

 \usage{CBS(...)}

 \arguments{
  \item{...}{Not used.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{append} \tab  -\cr
 \tab \code{as} \tab  -\cr
 \tab \code{as.character} \tab  -\cr
 \tab \code{as.data.frame} \tab  -\cr
 \tab \code{\link[PSCBS:as.DNAcopy.CBS]{as.DNAcopy}} \tab Coerces a CBS object to a DNAcopy object.\cr
 \tab \code{drawCentromeres} \tab  -\cr
 \tab \code{drawLevels} \tab  -\cr
 \tab \code{\link[PSCBS:estimateStandardDeviation.CBS]{estimateStandardDeviation}} \tab Estimates the whole-genome standard deviation of the signals.\cr
 \tab \code{extractByChromosomes} \tab  -\cr
 \tab \code{getCallStatistics} \tab  -\cr
 \tab \code{getChromosomeRanges} \tab  -\cr
 \tab \code{getChromosomes} \tab  -\cr
 \tab \code{getFGA} \tab  -\cr
 \tab \code{getFGG} \tab  -\cr
 \tab \code{getFGL} \tab  -\cr
 \tab \code{getFractionOfGenomeAltered} \tab  -\cr
 \tab \code{getFractionOfGenomeGained} \tab  -\cr
 \tab \code{getFractionOfGenomeLost} \tab  -\cr
 \tab \code{getSignalType} \tab  -\cr
 \tab \code{highlightCalls} \tab  -\cr
 \tab \code{isWholeChromosomeGained} \tab  -\cr
 \tab \code{isWholeChromosomeLost} \tab  -\cr
 \tab \code{nbrOfAmplifications} \tab  -\cr
 \tab \code{nbrOfChromosomes} \tab  -\cr
 \tab \code{nbrOfGains} \tab  -\cr
 \tab \code{nbrOfLoci} \tab  -\cr
 \tab \code{nbrOfLosses} \tab  -\cr
 \tab \code{nbrOfSegments} \tab  -\cr
 \tab \code{plot} \tab  -\cr
 \tab \code{\link[PSCBS:plotTracks.CBS]{plotTracks}} \tab Plots copy numbers along the genome.\cr
 \tab \code{print} \tab  -\cr
 \tab \code{signalType} \tab  -\cr
 \tab \code{writeLocusData} \tab  -\cr
 \tab \code{writeSegments} \tab  -\cr
}


 \bold{Methods inherited from list}:\cr
all.equal, as.data.frame, attachLocally, callHooks, relist, within


 }

 \section{Difference to DNAcopy object}{
   A CBS object is similar to DNAcopy objects with the major
   difference that a CBS object holds only one sample, whereas
   a DNAcopy object can hold more than one sample.
 }

 \section{See also}{
   \code{\link{segmentByCBS}}()
 }

 \author{Henrik Bengtsson}
\keyword{classes}
