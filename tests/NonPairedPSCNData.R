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


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Split up in tumor and normal non-paired PSCN data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dataT <- extractNonPairedPSCNData(data, "T");
print(head(dataT));

dataN <- extractNonPairedPSCNData(data, "N");
print(head(dataN));


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Merge tumor and normal into paired PSCN data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
data2 <- as.PairedPSCNData(T=dataT, N=dataN);
print(head(data2))

# Validate
stopifnot(all.equal(data2, data));
