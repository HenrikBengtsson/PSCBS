library("PSCBS")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Load SNP microarray data
# (note to package developers: this example data set may
#  be replaced in a future release of the package)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
pathname <- system.file("data-ex/PairedPSCBS,exData,chr01.Rbin", package="PSCBS")
data <- R.utils::loadObject(pathname)
str(data)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Paired PSCBS segmentation
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Drop single-locus outliers
dataS <- dropSegmentationOutliers(data)

# Run light-weight tests by default
if (Sys.getenv("_R_CHECK_FULL_") == "") {
  # Use only every 5th data point
  dataS <- dataS[seq(from=1, to=nrow(data), by=5),]
  # Number of segments (for assertion)
  nSegs <- 3L
  # Number of bootstrap samples (see below)
  B <- 100L
} else {
  # Full tests
  nSegs <- 12L
  B <- 1000L
}

str(dataS)

R.oo::attachLocally(dataS)

# Paired PSCBS segmentation
fit <- segmentByPairedPSCBS(CT, betaT=betaT, betaN=betaN,
                            chromosome=chromosome, x=x, 
                            seed=0xBEEF, verbose=-10)
print(fit)

# Plot results
plotTracks(fit)

# Sanity check
stopifnot(nbrOfSegments(fit) == nSegs)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Bootstrap segment level estimates
# (used by the AB caller, which, if skipped here,
#  will do it automatically)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fit <- bootstrapTCNandDHByRegion(fit, B=B, verbose=-10)
print(fit)
plotTracks(fit)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Calling segments in allelic balance (AB)
# NOTE: Ideally, this should be done on whole-genome data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Explicitly estimate the threshold in DH for calling AB
# (which be done by default by the caller, if skipped here)
deltaAB <- estimateDeltaAB(fit, flavor="qq(DH)", verbose=-10)
print(deltaAB)
## [1] 0.1657131

fit <- callAB(fit, delta=deltaAB, verbose=-10)
print(fit)
plotTracks(fit)

# Even if not explicitly specified, the estimated 
# threshold parameter is returned by the caller
stopifnot(fit$params$deltaAB == deltaAB)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Calling segments in loss-of-heterozygosity (LOH)
# NOTE: Ideally, this should be done on whole-genome data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Explicitly estimate the threshold in C1 for calling LOH
# (which be done by default by the caller, if skipped here)
deltaLOH <- estimateDeltaLOH(fit, flavor="minC1|nonAB", verbose=-10)
print(deltaLOH)
## [1] 0.625175

fit <- callLOH(fit, delta=deltaLOH, verbose=-10)
print(fit)
plotTracks(fit)

# Even if not explicitly specified, the estimated 
# threshold parameter is returned by the caller
stopifnot(fit$params$deltaLOH == deltaLOH)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Calling segments with run of homozygosity (ROH)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fit <- callROH(fit, verbose=-10)
print(fit)
plotTracks(fit)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Calling segments that are gained, copy neutral, and lost
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fit <- callGNL(fit, verbose=-10)
print(fit)
plotTracks(fit)

