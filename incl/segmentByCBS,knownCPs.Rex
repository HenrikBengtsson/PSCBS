
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Simulating copy-number data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
set.seed(0xBEEF)

# Number of loci
J <- 1000

mu <- double(J)
mu[200:300] <- mu[200:300] + 1
mu[350:400] <- NA # centromere
mu[650:800] <- mu[650:800] - 1
eps <- rnorm(J, sd=1/2)
y <- mu + eps
x <- sort(runif(length(y), max=length(y))) * 1e5
w <- runif(J)
w[650:800] <- 0.001


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Segmentation
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fit <- segmentByCBS(y, x=x)
print(fit)
plotTracks(fit)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Segmentation with some known change points
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
knownSegments <- data.frame(
  chromosome=c(    0,   0),
  start     =x[c(  1, 401)],
  end       =x[c(349,   J)]+c(0,10e6)
)
fit2 <- segmentByCBS(y, x=x, knownSegments=knownSegments, verbose=TRUE)
print(fit2)
plotTracks(fit2)
abline(v=c(knownSegments$start, knownSegments$end)/1e6, lty=3)


knownSegments <- data.frame(
  chromosome=c(    0,  0,   0),
  start     =x[c(  1, NA, 401)],
  end       =x[c(349, NA,   J)]+c(0,10e6)
)
fit2b <- segmentByCBS(y, x=x, knownSegments=knownSegments, verbose=TRUE)
print(fit2b)
plotTracks(fit2b)
abline(v=c(knownSegments$start, knownSegments$end)/1e6, lty=3)


# Segment just the (empty) centromere. All statistics will be NAs.
knownSegments <- data.frame(
  chromosome=c(    0),
  start     =x[c(350)],
  end       =x[c(400)]
)
fit3 <- segmentByCBS(y, x=x, knownSegments=knownSegments, verbose=TRUE)
print(fit3)



# If one specify the (empty) centromere as a segment, then its
# estimated statistics will be NAs.  [TODO /HB 2011-10-19]
knownSegments <- data.frame(
  chromosome=c(    0,   0,   0),
  start     =x[c(  1, 350, 401)],
  end       =x[c(349, 400,   J)]
)
fit4 <- segmentByCBS(y, x=x, knownSegments=knownSegments, verbose=TRUE)
print(fit4)
plotTracks(fit4)
abline(v=c(knownSegments$start, knownSegments$end)/1e6, lty=3)

