###########################################################################/**
# @set "class=AbstractCBS"
# @RdocMethod plotTracks
#
# @title "Plots the segmentation result along the genome"
#
# \description{
#   @get "title".
# }
# 
# @synopsis
#
# \arguments{
#  \item{...}{...}
# }
#
# \value{
#   Returns nothing.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################  
setMethodS3("plotTracks", "AbstractCBS", abstract=TRUE);


setMethodS3("tileChromosomes", "AbstractCBS", abstract=TRUE, protected=TRUE);


setMethodS3("drawChangePoints", "AbstractCBS", abstract=TRUE, protected=TRUE);


############################################################################
# HISTORY:
# 2011-12-03
# o Added drawChangePoints().
# 2011-10-02
# o Created.
############################################################################
