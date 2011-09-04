
# Coerce an CBS object to a DNAcopy object
fitD <- as.DNAcopy(fit)

# Coerce an DNAcopy object to a CBS object
fitC <- as.CBS(fitD)

# Sanity check
fitD2 <- as.DNAcopy(fit)
stopifnot(all.equal(fitD2, fitD))

fitC2 <- as.CBS(fitD2)
stopifnot(all.equal(fitC2, fitC))
