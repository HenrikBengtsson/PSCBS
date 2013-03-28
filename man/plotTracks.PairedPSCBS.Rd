%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  PairedPSCBS.PLOT.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{plotTracks.PairedPSCBS}
\alias{plotTracks.PairedPSCBS}
\alias{PairedPSCBS.plotTracks}
\alias{plotTracks,PairedPSCBS-method}
\alias{PairedPSCBS.plotTracks1}
\alias{plotTracks1.PairedPSCBS}
\alias{plotTracks1,PairedPSCBS-method}

\alias{PairedPSCBS.plotTracks2}
\alias{plotTracks2.PairedPSCBS}
\alias{plotTracks2,PairedPSCBS-method}

\alias{PairedPSCBS.plotTracksManyChromosomes}
\alias{plotTracksManyChromosomes.PairedPSCBS}
\alias{plotTracksManyChromosomes,PairedPSCBS-method}

\alias{plot}

\title{Plots parental specific copy numbers along the genome}

\description{
 Plots parental specific copy numbers along the genome for one or more chromosomes.
 It is possible to specify what type of tracks to plot.
 Each type of track is plotted in its own panel.
}

\usage{
\method{plotTracks}{PairedPSCBS}(fit, ...)
}

\arguments{
  \item{x}{A result object returned by \code{\link{segmentByPairedPSCBS}}().}
  \item{tracks}{A \code{\link[base]{character}} \code{\link[base]{vector}} specifying what types of tracks to plot.}
  \item{scatter}{A \code{\link[base]{character}} \code{\link[base]{vector}} specifying which of the tracks should
    have scatter plot.}
  \item{calls}{A \code{\link[base]{character}} \code{\link[base]{vector}} of regular expression identifying
    call labels to be highlighted in the panels.}
  \item{pch}{The type of the scatter points, if any.}
  \item{col}{The color of the scatter points, if any.}
  \item{cex}{The size of the scatter points, if any.}
  \item{changepoints}{If \code{\link[base:logical]{TRUE}}, changepoints are drawn as vertical lines.}
  \item{grid}{If \code{\link[base:logical]{TRUE}}, horizontal lines are displayed.}
  \item{quantiles}{A \code{\link[base]{numeric}} \code{\link[base]{vector}} in [0,1] specifying the quantiles
     of the confidence bands to be drawn, if any.}
  \item{xlim}{(Optional) The genomic range to plot.}
  \item{Clim}{The range of copy numbers.}
  \item{Blim}{The range of allele B fractions (BAFs) and
    decrease of heterozygosity (DHs).}
  \item{xScale}{The scale factor used for genomic positions.}
  \item{...}{Not used.}
  \item{add}{If \code{\link[base:logical]{TRUE}}, the panels plotted are added to the existing plot,
    otherwise a new plot is created.}
  \item{subplots}{If \code{\link[base:logical]{TRUE}}, then subplots are automatically setup.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
  Returns nothing.
}

\author{Henrik Bengtsson}



\keyword{internal}
\keyword{methods}
\keyword{IO}
