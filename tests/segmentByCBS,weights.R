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

par(mar=c(1.7,1,0.2,1)+0.1)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Segmentation without weights
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fit <- segmentByCBS(y, x=x)
sampleName(fit) <- "CBS_Example"
print(fit)
plotTracks(fit)
## Highlight outliers (they pull up the mean levels)
points(x[outliers]/1e6, y[outliers], col="purple")


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Segmentation with weights
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fitW <- segmentByCBS(y, x=x, w=w)
sampleName(fitW) <- "CBS_Example (weighted)"
print(fitW)
drawLevels(fitW, col="red")

legend("topright", legend=c("outliers", "non-weighted CBS", "weighted CBS"), col=c("purple", "purple", "red"), lwd=c(NA,3,3), pch=c(1,NA,NA))
