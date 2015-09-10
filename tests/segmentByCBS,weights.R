library("PSCBS")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Simulating copy-number data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
set.seed(0xBEEF)

# Number of loci
J <- 1000

x <- sort(runif(J, max=J)) * 1e5

mu <- double(J)
mu[200:300] <- mu[200:300] + 1
mu[350:400] <- NA # centromere
mu[650:800] <- mu[650:800] - 1
eps <- rnorm(J, sd=1/2)
y <- mu + eps

outliers <- seq(from=1L, to=J, length.out=0.2*J)
y[outliers] <- y[outliers] + 1.5

w <- rep(1.0, times=J)
w[outliers] <- 0.01

data <- data.frame(chromosome=1L, x=x, y=y)
dataW <- cbind(data, w=w)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Single-chromosome segmentation
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
par(mar=c(2,3,0.2,1)+0.1)
# Segment without weights
fit <- segmentByCBS(data)
sampleName(fit) <- "CBS_Example"
print(fit)
plotTracks(fit)
## Highlight outliers (they pull up the mean levels)
points(x[outliers]/1e6, y[outliers], col="purple")

# Segment with weights
fitW <- segmentByCBS(dataW)
sampleName(fitW) <- "CBS_Example (weighted)"
print(fitW)
drawLevels(fitW, col="red")

legend("topright", bg="white", legend=c("outliers", "non-weighted CBS", "weighted CBS"), col=c("purple", "purple", "red"), lwd=c(NA,3,3), pch=c(1,NA,NA))

## Assert that weighted segment means are less biased
dmean <- getSegments(fit)$mean - getSegments(fitW)$mean
cat("Segment mean differences:\n")
print(dmean)
stopifnot(all(dmean > 0, na.rm=TRUE))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Multi-chromosome segmentation
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
data2 <- data
data2$chromosome <- 2L
data <- rbind(data, data2)
dataW <- cbind(data, w=w)

par(mar=c(2,3,0.2,1)+0.1)
# Segment without weights
fit <- segmentByCBS(data)
sampleName(fit) <- "CBS_Example"
print(fit)
plotTracks(fit, Clim=c(-3,3))

# Segment with weights
fitW <- segmentByCBS(dataW)
sampleName(fitW) <- "CBS_Example (weighted)"
print(fitW)
drawLevels(fitW, col="red")

legend("topright", bg="white", legend=c("outliers", "non-weighted CBS", "weighted CBS"), col=c("purple", "purple", "red"), lwd=c(NA,3,3), pch=c(1,NA,NA))

## Assert that weighted segment means are less biased
dmean <- getSegments(fit)$mean - getSegments(fitW)$mean
cat("Segment mean differences:\n")
print(dmean)
stopifnot(all(dmean > 0, na.rm=TRUE))
