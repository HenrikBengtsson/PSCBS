library("revdepcheck")
if (isTRUE(as.logical(toupper(Sys.getenv("_R_CHECK_REVDEP_RESET_", "FALSE")))))
  revdep_reset()

options(warn = 1)

availableCores <- function() {
  getenv <- function(name) as.integer(Sys.getenv(name, NA_character_))
  getopt <- function(name) as.integer(getOption(name, NA_integer_))
  if (is.finite(n <- getopt("mc.cores") + 1L)) return(n)
  if (is.finite(n <- getopt("Ncpus") + 1L)) return(n)
  if (is.finite(n <- getenv("PBS_NUM_PPN"))) return(n)
  if (is.finite(n <- getenv("SLURM_CPUS_PER_TASK"))) return(n)
  if (is.finite(n <- getenv("NSLOTS"))) return(n)
  1L
}

revdep_check(bioc = TRUE, num_workers = availableCores(),
             timeout = as.difftime(10, units = "mins"), quiet = FALSE)
