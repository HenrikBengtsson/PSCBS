library("revdepcheck")
options(warn = 1L)

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

check <- function() {
  if (file_test("-f", p <- Sys.getenv("R_CHECK_ENVIRON", "~/.R/check.Renviron"))) {
    cat(sprintf("R CMD check will use env vars from %s\n", sQuote(p)))
    cat(sprintf("To disable, set 'R_CHECK_ENVIRON=false' (a fake pathname)\n"))
  }
  
  envs <- grep("^_R_CHECK_", names(Sys.getenv()), value = TRUE)
  if (length(envs) > 0L) {
    cat(sprintf("Detected _R_CHECK_* env vars that will affect R CMD check: %s\n",
                paste(sQuote(envs), collapse = ", ")))
  }
  
  revdep_check(bioc = TRUE, num_workers = availableCores(),
               timeout = as.difftime(20, units = "mins"), quiet = FALSE)
}


args <- base::commandArgs()
if ("--reset" %in% args) {
  revdep_reset()
} else if ("--todo" %in% args) {
  print(revdep_todo())
} else if ("--add" %in% args) {
  pos <- which("--add" == args)
  pkgs <- args[seq(from = pos + 1L, to = length(args))]
  revdep_add(packages = pkgs)
  cat(sprintf("* %s\n", revdep_todo()))
} else if ("--broken" %in% args) {
  revdep_add_broken()
  cat(sprintf("* %s\n", revdep_todo()))
} else {
  check()
}
