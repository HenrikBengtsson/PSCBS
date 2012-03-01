library("PSCBS")


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Load SNP microarray data
# (note to package developers: this example data set may
#  be replaced in a future release of the package)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
pathname <- system.file("data-ex/PairedPSCBS,exData,chr01.Rbin", package="PSCBS")
data <- R.utils::loadObject(pathname)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Paired PSCN data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
data <- PairedPSCNData(data);
print(head(data));

data <- callNaiveGenotypes(data);
print(head(data));

data <- normalizeTumorBoost(data);
print(head(data));

data <- isSNP(data);
print(head(data));

data <- dropSegmentationOutliers(data);
print(head(data));

idxs <- hasKnownPositions(data);
str(idxs);
