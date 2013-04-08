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
verbose <- Arguments$getVerbose(-10*interactive(), timestamp=TRUE)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Load SNP microarray data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
data("paired.chr01", package="PSCBS.data")
data <- paired.chr01

# Subset data (such that all regions will be in LOH)
data <- subset(data, x <= 125e6)

str(data)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Paired PSCBS segmentation
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Drop single-locus outliers
dataS <- dropSegmentationOutliers(data)

# Paired PSCBS segmentation
fit <- segmentByPairedPSCBS(dataS, seed=0xBEEF, verbose=verbose)

fit <- callAB(fit, verbose=verbose)
fit <- callLOH(fit, verbose=verbose)
print(fit)
str(fit$params) # deltaAB=1.17, deltaLowC1=-Inf
plotTracks(fit)
