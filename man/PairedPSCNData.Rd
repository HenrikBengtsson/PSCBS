%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  PairedPSCNData.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{PairedPSCNData}
\docType{class}
\alias{PairedPSCNData}


 \title{The PairedPSCNData class}

 \description{
  Package:  PSCBS \cr
\bold{Class PairedPSCNData}\cr

\code{data.frame}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[PSCBS]{AbstractCNData}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[PSCBS]{AbstractPSCNData}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{PairedPSCNData}\cr

\bold{Directly known subclasses:}\cr
\cr

public class \bold{PairedPSCNData}\cr
extends \link[PSCBS]{AbstractPSCNData}\cr



  A PairedPSCNData object holds paired tumor-normal parent-specific
  copy number data.
 }

 \usage{PairedPSCNData(chromosome=NULL, x=NULL, isSNP=NULL, muN=NULL, CT=NULL, betaT=NULL, CN=NULL, betaN=NULL, ...)}

 \arguments{
   \item{CT}{A \code{\link[base]{numeric}} \code{\link[base]{vector}} of J tumor total copy number (TCN)
        ratios in [0,+\code{\link[base:is.finite]{Inf}}) (due to noise, small negative values are
        also allowed).  The TCN ratios are typically scaled such that
        copy-neutral diploid loci have a mean of two.}
   \item{CN}{An optional \code{\link[base]{numeric}} \code{\link[base]{vector}} of J normal TCN ratios.}
   \item{betaT}{A \code{\link[base]{numeric}} \code{\link[base]{vector}} of J tumor allele B fractions (BAFs)
        in [0,1] (due to noise, values may be slightly outside as well)
        or \code{\link[base]{NA}} for non-polymorphic loci.}
   \item{betaN}{A \code{\link[base]{numeric}} \code{\link[base]{vector}} of J matched normal BAFs in [0,1]
        (due to noise, values may be slightly outside as well) or \code{\link[base]{NA}}
        for non-polymorphic loci.}
   \item{muN}{An optional \code{\link[base]{numeric}} \code{\link[base]{vector}} of J genotype calls in
        \{0,1/2,1\} for AA, AB, and BB, respectively,
        and \code{\link[base]{NA}} for non-polymorphic loci.
        If not given, they are estimated from the normal BAFs using
        \code{\link[aroma.light]{callNaiveGenotypes}} as described in [2].}
   \item{isSNP}{An optional \code{\link[base]{logical}} \code{\link[base]{vector}} of length J specifying
        whether each locus is a SNP or not (non-polymorphic loci).}
   \item{chromosome}{(Optional) An \code{\link[base]{integer}} scalar (or a \code{\link[base]{vector}} of length J),
        which can be used to specify which chromosome each locus belongs to
        in case multiple chromosomes are segments.
        This argument is also used for annotation purposes.}
   \item{x}{Optional \code{\link[base]{numeric}} \code{\link[base]{vector}} of J genomic locations.
            If \code{\link[base]{NULL}}, index locations \code{1:J} are used.}
   \item{...}{Optional named locus-specific signal \code{\link[base]{vector}}s of length J.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{as} \tab  -\cr
 \tab \code{as.PairedPSCNData} \tab  -\cr
 \tab \code{callNaiveGenotypes} \tab  -\cr
 \tab \code{callSegmentationOutliers} \tab  -\cr
 \tab \code{dropSegmentationOutliers} \tab  -\cr
 \tab \code{getTCNs} \tab  -\cr
 \tab \code{getTotalCopyNumbers} \tab  -\cr
 \tab \code{normalizeTumorBoost} \tab  -\cr
 \tab \code{segmentByPairedPSCBS} \tab  -\cr
}


 \bold{Methods inherited from AbstractPSCNData}:\cr
getSNPFields, isSNP

 \bold{Methods inherited from AbstractCNData}:\cr
findLargeGaps, getChromosomes, getLocusData, hasKnownPositions, nbrOfChromosomes, nbrOfLoci, orderAlongGenome

 \bold{Methods inherited from data.frame}:\cr
$<-, [, [[, [[<-, [<-, aggregate, anyDuplicated, as.data.frame, as.list, as.matrix, as.NonPairedPSCNData, as.PairedPSCNData, attachLocally, by, callSegmentationOutliers, cbind, dim, dimnames, dimnames<-, droplevels, dropSegmentationOutliers, duplicated, edit, findLargeGaps, format, formula, gapsToSegments, head, is.na, Math, mean, merge, na.exclude, na.omit, Ops, plot, plotDensity, print, prompt, rbind, row.names, row.names<-, rowsum, segmentByCBS, segmentByPairedPSCBS, split, split<-, stack, str, subset, summary, Summary, t, tail, transform, unique, unstack, unwrap, within, wrap, writeDataFrame


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
