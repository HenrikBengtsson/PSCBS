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
\emph{\link[PSCBS]{CBS}}, \emph{\link[PSCBS]{NonPairedPSCBS}}, \emph{\link[PSCBS]{PSCBS}}, \emph{\link[PSCBS]{PairedPSCBS}}\cr

public abstract class \bold{AbstractCBS}\cr
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
 \tab \code{adjustPloidyScale} \tab  -\cr
 \tab \code{extractCNs} \tab  -\cr
 \tab \code{getChangePoints} \tab  -\cr
 \tab \code{getChromosomes} \tab  -\cr
 \tab \code{getLocusData} \tab  -\cr
 \tab \code{getSegmentSizes} \tab  -\cr
 \tab \code{getSegments} \tab  -\cr
 \tab \code{mergeThreeSegments} \tab  -\cr
 \tab \code{nbrOfChangePoints} \tab  -\cr
 \tab \code{nbrOfChromosomes} \tab  -\cr
 \tab \code{nbrOfLoci} \tab  -\cr
 \tab \code{nbrOfSegments} \tab  -\cr
 \tab \code{normalizeTotalCNs} \tab  -\cr
 \tab \code{ploidy} \tab  -\cr
 \tab \code{ploidy<-} \tab  -\cr
 \tab \code{plotTracks} \tab  -\cr
 \tab \code{sampleCNs} \tab  -\cr
 \tab \code{writeWIG} \tab  -\cr
}


 \bold{Methods inherited from list}:\cr
Ops,nonStructure,vector-method, Ops,structure,vector-method, Ops,vector,nonStructure-method, Ops,vector,structure-method, all.equal, as.data.frame, attachLocally, callHooks, coerce,ANY,list-method, relist, type.convert, within


}

\author{Henrik Bengtsson}


\keyword{classes}
\keyword{internal}
