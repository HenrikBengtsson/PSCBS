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
\link[PSCBS]{CBS}, \link[PSCBS]{PairedPSCBS}, \emph{\link[PSCBS]{PSCBS}}\cr

public abstract class \bold{AbstractCBS}\cr
extends list\cr



  All CBS-style segmentation results extend this class, e.g.
  \code{\link{CBS}} and \code{\link{PairedPSCBS}}.
 }

 \usage{AbstractCBS(fit=list(), ...)}

 \arguments{
   \item{fit}{A \code{\link[base]{list}} structure containing the segmentation results.}
   \item{...}{Not used.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{append} \tab  -\cr
 \tab \code{\link[PSCBS:as.data.frame.AbstractCBS]{as.data.frame}} \tab Gets the table of segments.\cr
 \tab \code{\link[PSCBS:getChromosomes.AbstractCBS]{getChromosomes}} \tab Gets the set of chromosomes.\cr
 \tab \code{\link[PSCBS:getSegments.AbstractCBS]{getSegments}} \tab Gets the segments.\cr
 \tab \code{\link[PSCBS:nbrOfChromosomes.AbstractCBS]{nbrOfChromosomes}} \tab Gets the number of chromosomes.\cr
 \tab \code{\link[PSCBS:nbrOfSegments.AbstractCBS]{nbrOfSegments}} \tab Gets the number of segments.\cr
 \tab \code{plotTracks} \tab  -\cr
}


 \bold{Methods inherited from list}:\cr
all.equal, as.data.frame, attachLocally, callHooks, relist, within


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
