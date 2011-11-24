library("PSCBS")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Load SNP microarray data
# (note to package developers: this example data set may
#  be replaced in a future release of the package)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
pathname <- system.file("data-ex/PairedPSCBS,exData,chr01.Rbin", package="PSCBS")
data <- R.utils::loadObject(pathname)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Paired PSCBS segmentation
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Drop single-locus outliers
dataS <- dropSegmentationOutliers(data)
str(dataS)

fig <- 1;


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# (a) Don't segment the centromere (and force a separator)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
knownSegments <- data.frame(
  chromosome = c(        1,  1,         1),
  start      = c(     -Inf, NA, 141510003),
  end        = c(120992603, NA,      +Inf)
)


# Paired PSCBS segmentation
fit <- segmentByPairedPSCBS(dataS, knownSegments=knownSegments, 
                            seed=0xBEEF, verbose=-10)
print(fit)

# Plot results
devSet(list(fit, "tracks"));
plotTracks(fit)

# Sanity check
stopifnot(nbrOfSegments(fit) == 12)

fit1 <- fit


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# (b) Segment also the centromere (which will become NAs)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
knownSegments <- data.frame(
  chromosome = c(        1,         1,         1),
  start      = c(     -Inf, 120992604, 141510003),
  end        = c(120992603, 141510002,      +Inf)
)


# Paired PSCBS segmentation
fit <- segmentByPairedPSCBS(dataS, knownSegments=knownSegments, 
                            seed=0xBEEF, verbose=-10)
print(fit)

# Plot results
devSet(list(fit, "tracks"));
plotTracks(fit)

# Sanity check [TO FIX: See above]
stopifnot(nbrOfSegments(fit) == 12)

fit2 <- fit


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# (c) Do not segment the centromere (without a separator)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
knownSegments <- data.frame(
  chromosome = c(        1,         1),
  start      = c(     -Inf, 141510003),
  end        = c(120992603,      +Inf)
)

# Paired PSCBS segmentation
fit <- segmentByPairedPSCBS(dataS, knownSegments=knownSegments, 
                            seed=0xBEEF, verbose=-10)
print(fit)

# Plot results
devSet(list(fit, "tracks"));
plotTracks(fit)

# Sanity check
stopifnot(nbrOfSegments(fit) == 11)

fit3 <- fit
