#############################################################
# PROBLEMATIC CASE
# Title: All chromosome is in LOH
# Created by: Henrik Bengtsson
# Created on: 2011-07-08
#
# Description:
# The following example illustrates a case when the AB
# caller fails, because Delta_AB is estimated too large,
# resulting in all segments being called AB and none LOH.
#
# Use case:
# This could happend when PSCBS runs on a single chromosome
# that is in LOH.  It is unlikely to happen when processing
# whole-genome data.
#############################################################

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Load SNP microarray data
# (note to package developers: this example data set may
#  be replaced in a future release of the package)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
pathname <- system.file("data-ex/PairedPSCBS,exData,chr01.Rbin", package="PSCBS")
data <- R.utils::loadObject(pathname)

# Subset data (such that all regions will be in LOH)
data <- subset(data, position <= 125e6)

str(data)
R.oo::attachLocally(data)
x <- position


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Paired PSCBS segmentation
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Drop single-locus outliers
CTs <- dropSegmentationOutliers(CT, chromosome=1, x=x, verbose=-10)

# Paired PSCBS segmentation
fit <- segmentByPairedPSCBS(CTs, betaT=betaT, betaN=betaN,
                            chromosome=1, x=x, 
                            seed=0xBEEF, verbose=-10)

fit <- callAB(fit, verbose=-10)
fit <- callLOH(fit, verbose=-10)
print(fit)
str(fit$params) # deltaAB=1.17, deltaLowC1=-Inf
plotTracks(fit)
