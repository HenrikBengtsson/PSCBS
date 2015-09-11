library("PSCBS")
set.seed(0xBEEF)
J <- 1000
mu <- double(J)
mu[200:300] <- mu[200:300] + 1
mu[350:400] <- NA
mu[650:800] <- mu[650:800] - 1
eps <- rnorm(J, sd=1/2)
y <- mu + eps
x <- sort(runif(length(y), max=length(y))) * 1e5

data1 <- data.frame(chromosome=1L, x=x, y=y)
data2 <- data.frame(chromosome=2L, x=x, y=y)
dataM <- rbind(data1, data2)

fit1 <- segmentByCBS(data1)
print(fit1)
fit2 <- segmentByCBS(data2)
print(fit2)
fitM <- segmentByCBS(dataM)
print(fitM)

fit1p <- pruneBySdUndo(fit1)
print(fit1p)
fit2p <- pruneBySdUndo(fit2)
print(fit2p)
## FIXME: Issue #20
## fitMp <- pruneBySdUndo(fitM)
## print(fitMp)
