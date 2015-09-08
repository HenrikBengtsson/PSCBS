setMethodS3("extractWIG", "CBS", function(fit, graphType=c("bar", "points", "line"), nbrOfDecimals=4L, Clim=NULL, colors=c(negative="231,41,138", positive="117,112,179"), ...) {
  # Argument 'graphType':
  graphType <- match.arg(graphType)

  # Argument 'nbrOfDecimals':
  nbrOfDecimals <- Arguments$getInteger(nbrOfDecimals);

  data <- getSegments(fit, splitter=FALSE)
  data <- data[,c("chromosome", "start", "end", "mean")]
  data$chromosome <- sprintf("chr%d", data$chromosome)
  colnames(data) <- c("chrom", "chromStart", "chromEnd", "signal")

  ## Round / truncate
  for (ff in c("chromStart", "chromEnd")) {
    data[[ff]] <- as.integer(round(data[[ff]], digits=0L))
  }

  # Round mean levels
  if (!is.null(nbrOfDecimals)) {
    data[["signal"]] <- round(data[["signal"]], digits=nbrOfDecimals);
  }


  ## Track information
  track <- list(
    type="wiggle_0",
    name=sampleName(fit),
    description="Circular Binary Segmentation from PSCBS::segmentByCBS()",
    graphType=graphType,
    visibility="full",
    maxHeightPixels="128:96:64",
    yLineOnOff="on",
    autoScale="true"
  )
  if (is.na(track$name)) track$name <- "Unknown sample"

  if (!is.null(Clim)) {
    track$viewLimits <- sprintf("%g:%g", Clim[1], Clim[2])
  }

  if (!is.null(colors)) {
    if (!is.null(names(colors))) colors <- colors[c("negative", "positive")]
    track$color <- colors[["negative"]]
    track$altColor <- colors[["positive"]]
  }

  attr(data, "track") <- track

  data
}, protected=TRUE)


# \references{
#  [1] Wiggle Track Format (WIG), UCSC Genome Browser
#      \url{http://genome.ucsc.edu/goldenPath/help/wiggle.html}
# }
setMethodS3("writeWIG", "AbstractCBS", function(fit, name=getSampleName(fit), tags=NULL, ext="wig", path=NULL, overwrite=FALSE, skip=FALSE, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'name' and 'tags':
  name <- Arguments$getCharacter(name);
  tags <- Arguments$getCharacters(tags);

  # Argument 'ext':
  ext <- Arguments$getCharacter(ext);

  # Arguments 'path':
  path <- Arguments$getWritablePath(path);


  fullname <- paste(c(name, tags), collapse=",");
  filename <- sprintf("%s.%s", fullname, ext);
  pathname <- Arguments$getWritablePathname(filename, path=path, mustNotExist=(!overwrite && !skip));

  # File already exists?
  if (isFile(pathname)) {
    # Skip?
    if (skip) {
      return(pathname);
    }

    # Overwrite!
    file.remove(pathname);
  }

  ## Write file (atomically)
  pathnameT <- pushTemporaryFile(pathname)

  bed <- extractWIG(fit, ...)

  ## Generate 'track' definition line
  track <- attr(bed, "track")
  track <- lapply(track, FUN=function(value) {
    if (is.character(value)) value <- dQuote(value)
    value
  })
  track <- unlist(track, use.names=TRUE)
  track <- sprintf("%s=%s", names(track), track)
  track <- paste(track, collapse=" ")
  track <- sprintf("track %s", track)


  cat(track, "\n", sep="", file=pathnameT)
  write.table(bed, file=pathnameT,
              col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE,
              append=TRUE)

  pathname <- popTemporaryFile(pathnameT)

  pathname
})


############################################################################
# HISTORY:
# 2015-09-08
# o Added extractWIG() and writeWIG() for CBS objects.
############################################################################
