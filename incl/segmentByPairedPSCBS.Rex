library("R.utils");
library("aroma.light");


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Load data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
pathname <- system.file("data-ex/PairedPSCBS,exData,chr01.Rbin", package="PSCBS");
data <- R.utils::loadObject(pathname);

# Subset to speed up example
#data0 <- data;
#keep <- seq(from=1, to=nrow(data), by=5);
#data <- data[keep,,drop=FALSE];
# Order
o <- order(data$position);
data <- data[o,];
str(data);
R.oo::attachLocally(data);
x <- position;



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Paired PSCBS segmentation
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Drop single-locus outliers
CTs <- dropSegmentationOutliers(CT, chromosome=1, x=x, verbose=-10);

# Paired PSCBS segmentation
fit <- segmentByPairedPSCBS(CTs, betaT=betaT, betaN=betaN,
                            chromosome=1, x=x, verbose=-10);
print(fit);

# Plot results
plotTracks(fit);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Bootstrap segment level estimates
# (used by the AB caller)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fit <- bootstrapTCNandDHByRegion(fit, verbose=-10);
print(fit);
plotTracks(fit);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Calling segments in allelic balance (AB)
# NOTE: Ideally, this should be done on whole-genome data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Explicitly estimate the threshold in DH for calling AB
deltaAB <- estimateDeltaAB(fit, flavor="qq(DH)", verbose=-10);
print(deltaAB);
## [1] 

fit <- callAB(fit, delta=deltaAB, flavor="DeltaAB*", verbose=-10);
print(fit);
plotTracks(fit);



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Calling segments in loss-of-heterozygosity (LOH)
# NOTE: Ideally, this should be done on whole-genome data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Explicitly estimate the threshold in C1 for calling LOH
deltaLOH <- estimateDeltaLOH(fit, flavor="minC1|nonAB", verbose=-10);
print(deltaLOH);
## [1] 

fit <- callLOH(fit, delta=deltaLOH, flavor="SmallC1", verbose=-10);
print(fit);
plotTracks(fit);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Estimate background signals
# NOTE: Ideally, this should be done on whole-genome data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
kappa <- estimateKappa(fit);
print(kappa);
## [1] ...

