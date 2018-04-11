append <- function(...) UseMethod("append")
setMethodS3("append", "default", function(...) {
  base::append(...)
})
