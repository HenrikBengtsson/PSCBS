# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Simulating copy-number data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
set.seed(0xbeef);

# Number of loci
J <- 1000;

mu <- double(J);
mu[200:300] <- mu[200:300] + 1;
mu[650:800] <- mu[650:800] - 1;
eps <- rnorm(J, sd=1/2);
y <- mu + eps;
x <- sort(runif(length(y), max=length(y))) * 1e5;
w <- runif(J);
w[650:800] <- 0.001;


xlab <- "Position (Mb)";
ylim <- c(-3,3);
xMb <- x/1e6;
plot(xMb,y, pch=20, col="#aaaaaa", xlab=xlab, ylim=ylim);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Segment
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fit <- segmentByCBS(y, x=x);
print(fit);
drawLevels(fit, col="red", lwd=2, xScale=1e-6);



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# TESTS
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fit <- segmentByCBS(y, x=x, seed=0xbeef);

# Reverse the ordering and segment 
fitR <- segmentByCBS(rev(y), x=rev(x), seed=0xbeef);
stopifnot(all.equal(fitR$output, fit$output));
##stopifnot(all.equal(fitR$data[,2], fit$data[,2]));

# Reverse, but preserve ordering of 'data' object
fitRP <- segmentByCBS(rev(y), x=rev(x), preserveOrder=TRUE);
stopifnot(all.equal(fitRP$output, fit$output));
##stopifnot(all.equal(fitRP$data[,2], rev(x)));


# Change points inbetween data points at the same locus
x[650:654] <- x[649];
fitC <- segmentByCBS(rev(y), x=rev(x), preserveOrder=TRUE, seed=0xbeef);
## stopifnot(!is.null(fitC$lociNotPartOfSegment));
##stopifnot(fitC$data$maploc[unlist(fitC$lociNotPartOfSegment)] == x[649]);


# Allow for some missing values in signals
y[450] <- NA;
fitD <- segmentByCBS(y, x=x, seed=0xbeef);


# Allow for some missing genomic annotations
x[495] <- NA;
fitD <- segmentByCBS(y, x=x, seed=0xbeef, verbose=-10);


# Emulate a centromere
x[650:699] <- NA;
fit <- segmentByCBS(y, x=x, seed=0xbeef);
xMb <- x/1e6;
plot(xMb,y, pch=20, col="#aaaaaa", xlab=xlab, ylim=ylim);
drawLevels(fit, col="red", lwd=2, xScale=1e-6);

fitC <- segmentByCBS(y, x=x, joinSegments=TRUE, seed=0xbeef, verbose=-10);
drawLevels(fitC, col="blue", lwd=2, xScale=1e-6);
drawLevels(fit, col="red", lwd=2, xScale=1e-6);
