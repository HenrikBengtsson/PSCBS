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
str(dataS)

# Paired PSCBS segmentation
fit <- segmentByPairedPSCBS(dataS, seed=0xBEEF, verbose=-10)
print(fit)

tB0 <- system.time({
  fitB <- bootstrapTCNandDHByRegion(fit, B=100L);
});
print(tB0);

# Verbose output adds a substantial overhead
tB1 <- system.time({
  fitB <- bootstrapTCNandDHByRegion(fit, B=100L, verbose=-10);
});
print(tB1);
print(tB1/tB0);
