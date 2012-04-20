# This test script calls a report generator which requires
# the 'ggplot2' package, which in turn will require packages
# 'colorspace', 'dichromat', 'munsell', 'reshape2' and 'scales'.

# Only run this test in full testing mode
if (Sys.getenv("_R_CHECK_FULL_") != "") {
library("PSCBS")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Load SNP microarray data
# (note to package developers: this example data set may
#  be replaced in a future release of the package)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
pathname <- system.file("data-ex/PairedPSCBS,exData,chr01.Rbin", package="PSCBS")
data <- R.utils::loadObject(pathname)
str(data)

data <- data.frame(chromosome=data$chromosome, x=data$x, y=data$CT)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# CBS segmentation
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Drop single-locus outliers
dataS <- dropSegmentationOutliers(data)

# Speed up example by segmenting fewer loci
dataS <- dataS[seq(from=1, to=nrow(data), by=5),]

str(dataS)

gaps <- findLargeGaps(dataS, minLength=2e6)
knownSegments <- gapsToSegments(gaps)

# Paired PSCBS segmentation
fit <- segmentByCBS(dataS, knownSegments=knownSegments,
                            seed=0xBEEF, verbose=-10)
signalType(fit) <- "ratio"

# Fake a multi-chromosome segmentation
fit1 <- fit
fit2 <- renameChromosomes(fit, from=1, to=2)
fit <- append(fit1, fit2)

report(fit, sampleName="CBS", studyName="CBS-Ex", verbose=-10)

} # if (Sys.getenv("_R_CHECK_FULL_"))