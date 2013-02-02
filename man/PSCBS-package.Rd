%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  999.package.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{PSCBS-package}
\alias{PSCBS-package}

\docType{package}

\title{Package PSCBS}



 \description{
   Segmentation of allele-specific DNA copy number data and detection of regions with abnormal copy number within each parental chromosome.  Both tumor-normal paired and tumor-only analyses are supported..

   This package should be considered to be in an alpha or beta phase.
   You should expect the API to be changing over time.
 }

 \section{Requirements}{
   This package requires external packages
   R.utils (>= 1.19.5), utils, R.methodsS3 (>= 1.4.2), R.oo (>= 1.11.7), R.utils (>=
1.19.5), DNAcopy (>= 1.27.2), R.rsp (>= 0.8.2),
   and also suggests DNAcopy (>= 1.27.2), aroma.light (>= 1.22.0), matrixStats (>=
0.6.2), R.cache (>= 0.6.5), Hmisc (>= 3.10-1), ggplot2 (>=
0.9.3).
 }

 \section{Installation and updates}{
   To install this package, use \code{install.packages("PSCBS")}.
 }

 \section{To get started}{
   To get started, see:
   \enumerate{
     \item \code{\link{segmentByCBS}}() - segments total copy-numbers, or any
           other unimodal genomic signals, using the CBS method [3,4].
     \item \code{\link{segmentByPairedPSCBS}}() - segments allele-specific
           tumor signal from a tumor with a matched normal
           using the Paired PSCBS method [1,2].
   }
 }

 \section{How to cite}{
   Please use [1] and [2] to cite when using Paired PSCBS,
   and [3] and [4] when using CBS.
 }

 \author{
  A.B. Olshen, H. Bengtsson, P. Neuvial, P.T. Spellman, R.A. Olshen, V.E. Seshan.
 }

 \section{License}{
  GPL (>= 2).
 }

 \references{
  [1] A.B. Olshen, H. Bengtsson, P. Neuvial, P.T. Spellman, R.A. Olshen, V.E. Seshan, \emph{Parent-specific copy number in paired tumor-normal studies using circular binary segmentation}, Bioinformatics, 2011
 \cr
  [2] H. Bengtsson, P. Neuvial and T.P. Speed, \emph{TumorBoost: Normalization of allele-specific tumor copy numbers from a single pair of tumor-normal genotyping microarrays}, BMC Bioinformatics, 2010
 \cr
  [3] A.B. Olshen, E.S. Venkatraman (aka Venkatraman E. Seshan), R. Lucito and M. Wigler, \emph{Circular binary segmentation for the analysis of array-based DNA copy number data}, Biostatistics, 2004
 \cr
  [4] E.S. Venkatraman and A.B. Olshen, \emph{A faster circular binary segmentation algorithm for the analysis of array CGH data}, Bioinformatics, 2007
 \cr
 }
\keyword{package}
