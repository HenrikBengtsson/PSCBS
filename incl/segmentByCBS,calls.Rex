
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# CALLS
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Call gains and losses by segments
fitC <- callGainsAndLosses(fit)

# Call amplifications by segments
fitC <- callAmplifications(fitC)

# Call outliers by loci
fitC <- callOutliers(fitC)
