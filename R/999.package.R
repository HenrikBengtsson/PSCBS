#########################################################################/**
# @RdocPackage "PSCBS"
#
# \description{
#   @eval "packageDescription('PSCBS')$Description".
#
#   This package should be considered to be in an alpha or beta phase.
#   You should expect the API to be changing over time.
# }
#
# \section{Requirements}{
#   This package requires external packages
#   @eval "hpaste(unlist(packageDescription('PSCBS')[c('Depends', 'Imports')]))",	
#   and also suggests @eval "packageDescription('PSCBS')$Suggests".
# }
#
# \section{Installation and updates}{
#   To install this package, use \code{install.packages("PSCBS")}.
# }
#
# \section{To get started}{
#   To get started, see:
#   \enumerate{
#     \item @see "segmentByCBS" - segments total copy-numbers, or any
#           other unimodal genomic signals, using the CBS method [3,4].
#     \item @see "segmentByPairedPSCBS" - segments allele-specific 
#           tumor signal from a tumor with a matched normal
#           using the Paired PSCBS method [1,2].
#   }
# }
# 
# \section{How to cite}{
#   Please use [1] and [2] to cite when using Paired PSCBS,
#   and [3] and [4] when using CBS.
# }
#
# \author{
#  @eval "packageDescription('PSCBS')$Author".
# }
#
# \section{License}{
#  @eval "packageDescription('PSCBS')$License".
# }
# 
# \references{
#  [1] @include "../incl/OlshenA_etal_2011.Rd" \cr
#  [2] @include "../incl/BengtssonH_etal_2010.Rd" \cr 
#  [3] @include "../incl/OlshenVenkatraman_2004.Rd" \cr
#  [4] @include "../incl/VenkatramanOlshen_2007.Rd" \cr
# }
#*/#########################################################################

