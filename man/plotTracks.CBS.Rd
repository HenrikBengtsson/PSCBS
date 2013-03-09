%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  CBS.PLOT.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 \name{plotTracks.CBS}
\alias{plotTracks.CBS}
\alias{CBS.plotTracks}
\alias{plotTracks,CBS-method}

 \title{Plots copy numbers along the genome}

 \description{
  Plots copy numbers along the genome for one or more chromosomes.
  Each type of track is plotted in its own panel.
 }

 \usage{
\method{plotTracks}{CBS}(x, scatter=TRUE, pch=20, col="gray", meanCol="purple", cex=1, grid=FALSE, Clim="auto", xScale=1e-06, Clab="auto", ..., byIndex=FALSE, mar=NULL, add=FALSE)
}

 \arguments{
   \item{x}{A result object returned by \code{\link{segmentByCBS}}().}
   \item{pch}{The type of points to use.}
   \item{Clim}{The range of copy numbers.}
   \item{xScale}{The scale factor used for genomic positions.}
   \item{...}{Not used.}
   \item{add}{If \code{\link[base:logical]{TRUE}}, the panels plotted are added to the existing plot,
     otherwise a new plot is created.}
 }

 \value{
   Returns nothing.
 }

 \author{Henrik Bengtsson}



\keyword{internal}
\keyword{methods}
\keyword{IO}
