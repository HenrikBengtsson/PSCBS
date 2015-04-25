library("PSCBS")

# BUG FIX: PSCBS GitHub Issue #6
gaps <- findLargeGaps(chromosome=rep(0,10), x=1:10, minLength=2)
print(gaps)
stopifnot(is.data.frame(gaps))
stopifnot(nrow(gaps) == 0L)

# Simulating copy-number data
set.seed(0xBEEF)

# Simulate CN data
J <- 1000
mu <- double(J)
mu[200:300] <- mu[200:300] + 1
mu[350:400] <- NA # centromere
mu[650:800] <- mu[650:800] - 1
eps <- rnorm(J, sd=1/2)
y <- mu + eps
x <- seq(from=1, to=100e6, length.out=J)

data <- data.frame(chromosome=0L, x=x)

gaps <- findLargeGaps(x=x, minLength=1e6)
print(gaps)
stopifnot(is.data.frame(gaps))
stopifnot(nrow(gaps) == 0L)

gaps <- findLargeGaps(data, minLength=1e6)
print(gaps)
stopifnot(is.data.frame(gaps))
stopifnot(nrow(gaps) == 0L)

gaps <- findLargeGaps(x=x2, minLength=1e6)
print(gaps)
stopifnot(is.data.frame(gaps))
stopifnot(nrow(gaps) == 1L)

## Add missing values
x2 <- x
x2[30e6 < x & x < 50e6] <- NA
gaps <- findLargeGaps(x=x2, minLength=1e6)
print(gaps)
stopifnot(is.data.frame(gaps))
stopifnot(nrow(gaps) == 1L)
