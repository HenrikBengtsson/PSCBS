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

# Run light-weight tests by default
if (Sys.getenv("_R_CHECK_FULL_") == "") {
  # Use only every 5th data point
  dataS <- dataS[seq(from=1, to=nrow(data), by=5),]
  # Number of segments (for assertion)
  nSegs <- 4L
} else {
  # Full tests
  nSegs <- 13L
}

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
dev.set(2L)
plotTracks(fit)
abline(v=c(knownSegments$start, knownSegments$end)/1e6, lty=3)

# Sanity check
stopifnot(nbrOfSegments(fit) == nSegs)

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
dev.set(3L)
plotTracks(fit)
abline(v=c(knownSegments$start, knownSegments$end)/1e6, lty=3)

# Sanity check [TO FIX: See above]
stopifnot(nbrOfSegments(fit) == nSegs)

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
dev.set(4L)
plotTracks(fit)
abline(v=c(knownSegments$start, knownSegments$end)/1e6, lty=3)

# Sanity check
stopifnot(nbrOfSegments(fit) == nSegs-1L)

fit3 <- fit


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# (d) Skip the identification of new change points
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
knownSegments <- data.frame(
  chromosome = c(        1,         1),
  start      = c(     -Inf, 141510003),
  end        = c(120992603,      +Inf)
)

# Paired PSCBS segmentation
fit <- segmentByPairedPSCBS(dataS, knownSegments=knownSegments, 
                            undoTCN=Inf, undoDH=Inf,
                            seed=0xBEEF, verbose=-10)
print(fit)

# Plot results
dev.set(5L)
plotTracks(fit)
abline(v=c(knownSegments$start, knownSegments$end)/1e6, lty=3)

# Sanity check
stopifnot(nbrOfSegments(fit) == nrow(knownSegments))

fit4 <- fit
