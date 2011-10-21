
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



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# (a) Do not segment the centromere
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
knownSegments <- data.frame(
  chromosome = c(1L, 1L),
  start      = c(x[1], 141510003),
  end        = c(120992603, x[J])
)

# Paired PSCBS segmentation
fit <- segmentByPairedPSCBS(CTs, betaT=betaT, betaN=betaN,
                            chromosome=1, x=x, 
                            knownSegments=knownSegments, 
                            seed=0xBEEF, verbose=-10)
print(fit)

# Plot results
plotTracks(fit)

# Sanity check
stopifnot(nbrOfSegments(fit) == 11)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# (b) Segment the centromere but specify it
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
knownSegments <- data.frame(
  chromosome = c(1L, 1L, 1L),
  start      = c(x[1], 120992604, 141510003),
  end        = c(120992603, 141510002, x[J])
)


# Paired PSCBS segmentation
fit <- segmentByPairedPSCBS(CTs, betaT=betaT, betaN=betaN,
                            chromosome=1, x=x, 
                            knownSegments=knownSegments, 
                            seed=0xBEEF, verbose=-10)
print(fit)

# Plot results
plotTracks(fit)

# Sanity check
stopifnot(nbrOfSegments(fit) == 12)

