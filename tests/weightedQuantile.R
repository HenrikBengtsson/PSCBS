library("PSCBS")
library("stats")

message("weightedQuantile() ...")

if (requireNamespace("Hmisc")) {
  message(" - assert identical results to Hmisc::wtd.quantile()")
  wtd.quantile <- Hmisc::wtd.quantile
  for (kk in 1:100) {
    n <- 5L + sample.int(995, size = 1L)
    x <- rnorm(n, mean = 0.0, sd = 1.0)
    w <- runif(n, min = 0.5, max = 2.0) ## Non-normalized weights
    probs <- c(0.0, 0.25, 0.50, 0.75, 1.0)
    q0 <- wtd.quantile(x, weights = w, probs = probs, normwt = TRUE)
    q <- weightedQuantile(x, w = w, probs = probs)
    if (!isTRUE(all.equal(q, q0))) {
      print(q0)
      print(q)
      stopifnot(all.equal(q, q0))
    }
  }
}

message("weightedQuantile() ... DONE")
