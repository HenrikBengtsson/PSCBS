
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Load SNP microarray data
# (note to package developers: this example data set may
#  be replaced in a future release of the package)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
pathname <- system.file("data-ex/PairedPSCBS,exData,chr01.Rbin", package="PSCBS")
data <- R.utils::loadObject(pathname)

# Order by chromosome and position
o <- order(data$chromosome, data$position)
data <- data[o,]
str(data)
R.oo::attachLocally(data)
x <- position
J <- length(x)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Paired PSCBS segmentation
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Drop single-locus outliers
CTs <- dropSegmentationOutliers(CT, chromosome=1, x=x, verbose=-10)


fig <- 1;


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# (a) Don't segment the centromere (and force a separator)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
knownSegments <- data.frame(
  chromosome = c(        1,  1,         1),
  start      = c(     x[1], NA, 141510003),
  end        = c(120992603, NA,      x[J])
)


# Paired PSCBS segmentation
fit <- segmentByPairedPSCBS(CTs, betaT=betaT, betaN=betaN,
                            chromosome=1, x=x, 
                            knownSegments=knownSegments, 
                            seed=0xBEEF, verbose=-10)
print(fit)

# Plot results
devSet(list(fit, "tracks"));
plotTracks(fit)

devSet(list(fit, "(C1,C2)"), width=4.5);
plotC1C2(fit)
linesC1C2(fit)

# Sanity check
stopifnot(nbrOfSegments(fit) == 12)

fit1 <- fit


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# (b) Segment also the centromere (which will become NAs)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
knownSegments <- data.frame(
  chromosome = c(        1,         1,         1),
  start      = c(     x[1], 120992604, 141510003),
  end        = c(120992603, 141510002,      x[J])
)


# Paired PSCBS segmentation
fit <- segmentByPairedPSCBS(CTs, betaT=betaT, betaN=betaN,
                            chromosome=1, x=x, 
                            knownSegments=knownSegments, 
                            seed=0xBEEF, verbose=-10)
print(fit)

# Plot results
devSet(list(fit, "tracks"));
plotTracks(fit)

devSet(list(fit, "(C1,C2)"), width=4.5);
plotC1C2(fit)
linesC1C2(fit)

# Sanity check [TO FIX: See above]
stopifnot(nbrOfSegments(fit) == 12)

fit2 <- fit


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# (c) Do not segment the centromere (without a separator)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
knownSegments <- data.frame(
  chromosome = c(        1,         1),
  start      = c(     x[1], 141510003),
  end        = c(120992603,       x[J])
)

# Paired PSCBS segmentation
fit <- segmentByPairedPSCBS(CTs, betaT=betaT, betaN=betaN,
                            chromosome=1, x=x, 
                            knownSegments=knownSegments, 
                            seed=0xBEEF, verbose=-10)
print(fit)

# Plot results
devSet(list(fit, "tracks"));
plotTracks(fit)

devSet(list(fit, "(C1,C2)"), width=4.5);
plotC1C2(fit)
linesC1C2(fit)

# Sanity check
stopifnot(nbrOfSegments(fit) == 11)

fit3 <- fit

