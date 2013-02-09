setMethodS3("bootstrapTCNandDHByRegion", "PairedPSCBS", function(fit, B=1000L, statsFcn=function(x) quantile(x, probs=c(0.025, 0.050, 0.95, 0.975), na.rm=TRUE), by=c("betaTN", "betaT"), ..., seed=NULL, .debug=FALSE, force=FALSE, verbose=FALSE) {
  # Settings for sanity checks
  tol <- getOption("PSCBS/sanityChecks/tolerance", 0.0005);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'B':
  B <- Arguments$getInteger(B, range=c(1,Inf));

  # Argument 'by':
  by <- match.arg(by);

  # Argument 'seed':
  if (!is.null(seed)) {
    seed <- Arguments$getInteger(seed);
  }

  # Argument '.debug':
  .debug <- Arguments$getLogical(.debug);

  # Argument 'force':
  force <- Arguments$getLogical(force);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }



  verbose && enter(verbose, "Resample (TCN,DH) signals and re-estimate mean levels");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Set the random seed
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!is.null(seed)) {
    verbose && enter(verbose, "Setting (temporary) random seed");
    oldRandomSeed <- NULL;
    if (exists(".Random.seed", mode="integer")) {
      oldRandomSeed <- get(".Random.seed", mode="integer");
    }
    on.exit({
      if (!is.null(oldRandomSeed)) {
        .Random.seed <<- oldRandomSeed;
      }
    }, add=TRUE);
    verbose && cat(verbose, "The random seed will be reset to its original state afterward.");
    verbose && cat(verbose, "Seed: ", seed);
    set.seed(seed);
    verbose && exit(verbose);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extract data and estimates
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  data <- getLocusData(fit);
  tcnSegRows <- fit$tcnSegRows;
  dhSegRows <- fit$dhSegRows;
  segs <- getSegments(fit);
  params <- fit$params;

  # Sanity checks
  stopifnot(all(!is.na(data$chromosome) & !is.na(data$x)));

  # Sanity checks
  if (!params$joinSegments) {
    throw("Cannot bootstrap TCN and DH by segments unless PSCNs are segmented using joinSegments=TRUE.");
  }
  if (regexpr(",", params$flavor, fixed=TRUE) != -1) {
    throw(sprintf("Cannot bootstrap TCN and DH by segments if PSCNs are segmented using flavor=\"%s\".", params$flavor));
  }
  # Sanity check (same as above, but just in case)
  stopifnot(all(segs$tcnStart == segs$dhStart, na.rm=TRUE));
  stopifnot(all(segs$tcnEnd == segs$dhEnd, na.rm=TRUE));



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Find estimates to be done
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  stats <- statsFcn(1);
  stopifnot(!is.null(names(stats)));
  nbrOfStats <- length(stats);
  statsNames <- names(stats);

  # Already done?
  tcnStatsNames <- sprintf("tcn_%s", names(stats));
  dhStatsNames <- sprintf("dh_%s", names(stats));
  c1StatsNames <- sprintf("c1_%s", names(stats));
  c2StatsNames <- sprintf("c2_%s", names(stats));
  allStatsNames <- c(tcnStatsNames, dhStatsNames, c1StatsNames, c2StatsNames);
  isDone <- is.element(allStatsNames, names(segs));
  names(isDone) <- allStatsNames;
  verbose && cat(verbose, "Already done?");
  verbose && print(verbose, isDone);

  if (!force && all(isDone)) {
    verbose && cat(verbose, "Already done. Skipping.");
    verbose && exit(verbose);
    return(fit);
  }


  # Get mean estimators
  estList <- getMeanEstimators(fit, c("tcn", "dh"));
  avgTCN <- estList$tcn;
  avgDH <- estList$dh;
 

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Get signals
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Get (x,TCN,BAF) data
  chromosome <- data$chromosome;
  x <- data$x;
  CT <- data$CT;
  betaT <- data[[by]];
  muN <- data$muN;



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Classify each locus as (i) heterozygous SNP, (ii) homozygous SNP,
  # or (iii) non-polymorphic loci
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  verbose && enter(verbose, "Identifying heterozygous & homozygous SNPs and non-polymorphic loci");
  nbrOfLoci <- length(muN);
  verbose && cat(verbose, "Number of loci: ", nbrOfLoci);

  # SNPs are identifies as those loci that have non-missing 'muN' (& betaTN')
  isSNP <- (!is.na(muN) & !is.na(betaT));
  snps <- which(isSNP);
  nonSNPs <- which(!isSNP);
  nbrOfSNPs <- sum(isSNP);
  nbrOfNonSNPs <- sum(!isSNP);
  verbose && cat(verbose, "Number of SNPs: ", nbrOfSNPs);
  verbose && cat(verbose, "Number of non-SNPs: ", nbrOfNonSNPs);

  # Sanity checks
  stopifnot(length(intersect(snps, nonSNPs)) == 0);

  # Heterozygous SNPs
  isHet <- isSNP & (muN == 1/2);
  hets <- which(isSNP &  isHet);
  homs <- which(isSNP & !isHet);
  nbrOfHets <- length(hets);
  nbrOfHoms <- length(homs);
  verbose && printf(verbose, "Number of heterozygous SNPs: %d (%.2f%%)\n",
                                      nbrOfHets, 100*nbrOfHets/nbrOfSNPs);
  verbose && printf(verbose, "Number of homozygous SNPs: %d (%.2f%%)\n",
                                      nbrOfHoms, 100*nbrOfHoms/nbrOfSNPs);

  # Sanity checks
  stopifnot(length(intersect(hets, homs)) == 0);
  stopifnot(nbrOfHets + nbrOfHoms == nbrOfSNPs);

  # Sanity checks
  stopifnot(length(isSNP) == nbrOfLoci);
  stopifnot(length(isHet) == nbrOfLoci);

  # Not needed anymore
  rm(muN);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Precalculate DH signals
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Calculate DHs for heterozygous SNPs
  rho <- 2*abs(betaT - 1/2);

  # DH is by definition only defined for heterozygous SNPs.  
  # For simplicity, we set it to be NA for non-heterozygous loci.
  rho[!isHet] <- NA;

  data$rho <- rho;

  rm(betaT); # Not needed anymore



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Resample (TCN,DH) within each segments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  nbrOfSegments <- nrow(segs);

  # Allocate JxBx4 matrix M of bootstrap means
  naValue <- as.double(NA);
  dim <- c(nbrOfSegments, B, 4);
  dimnames <- list(NULL, NULL, c("tcn", "dh", "c1", "c2"));
  M <- array(naValue, dim=dim, dimnames=dimnames);
  verbose && str(verbose, M);

  # Identify all loci with non-missing signals
  idxsCT <- which(!is.na(CT));
  idxsRho <- which(!is.na(rho));

  # Vectorized pre-adjustments
  for (key in c("tcnNbrOfLoci", "dhNbrOfLoci")) {
    counts <- segs[[key]];
    counts[is.na(counts)] <- 0L;
    segs[[key]] <- counts;
  }

  hasTcnLoci <- (is.finite(tcnSegRows[,1]) & is.finite(tcnSegRows[,2]));
  hasDhLoci <- (is.finite(dhSegRows[,1]) & is.finite(dhSegRows[,2]));

  # Identify "splitter" segments which have no data
  chrs <- segs[["chromosome"]];
  tcnIds <- segs[["tcnId"]];
  dhIds <- segs[["dhId"]];
  tcnMeans <- segs[["tcnMean"]];
  dhMeans <- segs[["dhMean"]];
  isSplitter <- (is.na(chrs) && is.na(tcnIds) && is.na(dhIds));
 
  # Get all segment indices except for "splitters"
  jjs <- seq(length=nbrOfSegments);
  jjs <- jjs[!isSplitter];

  for (jj in jjs) {
    chr <- chrs[jj];
    tcnId <- tcnIds[jj];
    dhId <- dhIds[jj];

    verbose && enter(verbose, sprintf("Segment #%d (chr %d, tcnId=%d, dhId=%d) of %d", jj, chr, tcnId, dhId, nbrOfSegments));

    # Sanity check
    if (.debug) stopifnot(!is.na(chr) && !is.na(tcnId) && !is.na(dhId));

    # Get the segment data
    segJJ <- segs[jj,,drop=FALSE];
    verbose && print(verbose, segJJ);

    nbrOfTCNs <- segJJ[,"tcnNbrOfLoci"];
    nbrOfDHs <- segJJ[,"dhNbrOfLoci"];
    verbose && cat(verbose, "Number of TCNs: ", nbrOfTCNs);
    verbose && cat(verbose, "Number of DHs: ", nbrOfDHs);
    if (.debug) {
      stopifnot(!is.na(nbrOfTCNs));
      stopifnot(!is.na(nbrOfDHs));
    }

    tcnSegRowJJ <- unlist(tcnSegRows[jj,]);
    dhSegRowJJ <- unlist(dhSegRows[jj,]);

    # Indices of all loci
    if (hasTcnLoci[jj]) {
      idxsAll <- tcnSegRowJJ[1]:tcnSegRowJJ[2];
    } else {
      idxsAll <- integer(0);
    }

    verbose && str(verbose, idxsAll);
    verbose && print(verbose, hpaste(idxsAll));
    verbose && str(verbose, idxsCT);
    verbose && print(verbose, hpaste(idxsCT));

    # Keep only loci with finite TCNs
    idxsAll <- intersect(idxsAll, idxsCT);
    verbose && str(verbose, idxsAll);
    verbose && print(verbose, hpaste(idxsAll));

    # Sanity check
    if (length(idxsAll) != nbrOfTCNs) {
      verbose && str(verbose, setdiff(idxsCT, idxsAll));
      throw("INTERNAL ERROR: length(idxsAll) != nbrOfTCNs: ", length(idxsAll), " != ", nbrOfTCNs);
    }


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Identify loci used to calculate DH means
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    verbose && enter(verbose, "Identify loci used to bootstrap DH means");

    if (hasDhLoci[jj]) {
      idxsDH <- dhSegRowJJ[1]:dhSegRowJJ[2];
      idxsDH <- intersect(idxsDH, hets);
      # Drop missing values
      idxsDH <- intersect(idxsDH, idxsRho);
    } else {
      idxsDH <- integer(0);
    }

    verbose && cat(verbose, "Heterozygous SNPs to resample for DH:");
    verbose && str(verbose, idxsDH);

    # Sanity check
    if (.debug) stopifnot(length(idxsDH) == nbrOfDHs);

    verbose && exit(verbose);



    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Identify loci used to calculate TCN means
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    verbose && enter(verbose, "Identify loci used to bootstrap TCN means");

    # Identify SNPs and non-SNPs
    idxsSNP <- intersect(snps, idxsAll);
    idxsNonSNP <- setdiff(idxsAll, idxsSNP);
    verbose && cat(verbose, "SNPs:");
    verbose && str(verbose, idxsSNP);
    verbose && cat(verbose, "Non-polymorphic loci:");
    verbose && str(verbose, idxsNonSNP);
    # Sanity check
    if (.debug) stopifnot(length(idxsSNP) + length(idxsNonSNP) == length(idxsAll));

    # Identify heterozygous and homozygous SNPs
    idxsHet <- intersect(idxsSNP, hets);
    idxsHom <- intersect(idxsSNP, homs);

    # Drop missing values
    idxsNonSNP <- intersect(idxsNonSNP, idxsCT);
    idxsHet <- intersect(idxsHet, idxsCT);
    idxsHom <- intersect(idxsHom, idxsCT);
    idxsHetNonDH <- setdiff(idxsHet, idxsDH);

    verbose && cat(verbose, "Heterozygous SNPs to resample for TCN:");
    verbose && str(verbose, idxsHet);
    verbose && cat(verbose, "Homozygous SNPs to resample for TCN:");
    verbose && str(verbose, idxsHom);
        verbose && cat(verbose, "Non-polymorphic loci to resample for TCN:");
    verbose && str(verbose, idxsNonSNP);
    verbose && cat(verbose, "Heterozygous SNPs with non-DH to resample for TCN:");
    verbose && str(verbose, idxsHetNonDH);
    # Note that length(idxsHetNonDH) may differ from zero in case CT is non-missing
    # but rho is missing, e.g. CT = sum(c(thetaA,thetaB), na.rm=TRUE) and
    # thetaB is missing. /HB 2010-12-01

    idxsTCN <- sort(unique(c(idxsHet, idxsHom, idxsNonSNP)));
    verbose && cat(verbose, "Loci to resample for TCN:");
    verbose && str(verbose, idxsTCN);

    # Sanity check
    if (.debug) {
      stopifnot(length(idxsHet) + length(idxsHom) + length(idxsNonSNP) == nbrOfTCNs);
      stopifnot(length(intersect(idxsDH, idxsHetNonDH)) == 0);
      stopifnot(length(idxsTCN) == nbrOfTCNs);
    }

    verbose && exit(verbose);


    # These numbers should be preserved when the resampling
    verbose && printf(verbose, "Number of (#hets, #homs, #nonSNPs): (%d,%d,%d)\n",
                      length(idxsHet), length(idxsHom), length(idxsNonSNP));
                                              

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Sanity checks
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (nbrOfTCNs > 0L) {
      # Sanity check
      ys <- CT[idxsTCN];
      mu <- avgTCN(ys, na.rm=TRUE);
      dMu <- (mu - tcnMeans[jj]);
      if (abs(dMu) > tol) {
        str(list(nbrOfTCNs=nbrOfTCNs, tcnNbrOfLoci=segJJ$tcnNbrOfLoci, mu=mu, tcnMean=tcnMeans[jj], dMu=dMu, "abs(dMu)"=abs(dMu), "range(x[units])"=range(x[idxsTCN])));
        stop("INTERNAL ERROR: Incorrect TCN mean!");
      }
    }

    shouldHaveDHs <- (nbrOfDHs > 0L && !is.na(dhMeans[jj]));
    if (shouldHaveDHs) {
      # Sanity check
      ys <- rho[idxsDH];
      mu <- avgDH(ys, na.rm=TRUE);
      dMu <- (mu - dhMeans[jj]);
      if (abs(dMu) > tol) {
        str(list(nbrOfDHs=nbrOfDHs, dhNbrOfLoci=segJJ$dhNbrOfLoci, mu=mu, dhMean=dhMeans[jj], dMu=dMu, "abs(dMu)"=abs(dMu), "range(x[units])"=range(x[idxsDH])));
        stop("INTERNAL ERROR: Incorrect DH mean!");
      }
    }


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Bootstrap while preserving (#hets, #homs, #nonSNPs)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    verbose && enter(verbose, "Bootstrapping while preserving (#hets, #homs, #nonSNPs)");
    verbose && cat(verbose, "Number of bootstrap samples: ", B);

    if (!shouldHaveDHs) {
      idxsHetNonDH <- idxsDH;
    }

    nHoms <- length(idxsHom);
    nNonSNPs <- length(idxsNonSNP);
    nHetNonDHs <- length(idxsHetNonDH);

    # Defaults
    idxsDHBB <- NULL;

    # Bootstrap B times
    for (bb in seq(length=B)) {
      # (1) Bootstrap DHs
      if (shouldHaveDHs) {
        # (a) Resample heterozygous SNPs (=> resampled DH units)
        idxsDHBB <- resample(idxsDH, size=nbrOfDHs, replace=TRUE);

        # Extract signals
        rhoBB <- rho[idxsDHBB];

        # Calculate bootstrap mean
        M[jj,bb,"dh"] <- avgDH(rhoBB, na.rm=TRUE);
      } # if (shouldHaveDHs)
  
      # (2) Bootstrap TCNs
      if (nbrOfTCNs > 0L) {
        # (a) Resample non-DH hets SNPs
        idxsHetNonDHBB <- resample(idxsHetNonDH, size=nHetNonDHs, replace=TRUE);
        idxsHetBB <- c(idxsDHBB, idxsHetNonDHBB);

        # (a) Resample homozygous SNPs
        idxsHomBB <- resample(idxsHom, size=nHoms, replace=TRUE);
  
        # (b) Resample non-SNPs
        idxsNonSNPBB <- resample(idxsNonSNP, size=nNonSNPs, replace=TRUE);
      
        # (c) Resampled TCN units
        idxsTCNBB <- c(idxsHetBB, idxsHomBB, idxsNonSNPBB);

        # Sanity check
        if (.debug) {
          stopifnot(length(intersect(idxsDHBB, idxsHetNonDHBB)) == 0L);
          stopifnot(length(intersect(idxsHetBB, idxsHomBB)) == 0L);
          stopifnot(length(intersect(idxsHetBB, idxsNonSNPBB)) == 0L);
          stopifnot(length(intersect(idxsHomBB, idxsNonSNPBB)) == 0L);
        }

        # Extract signals
        CTBB <- CT[idxsTCNBB];

        # Calculate bootstrap mean
        M[jj,bb,"tcn"] <- avgTCN(CTBB, na.rm=TRUE);
      } # if (nbrOfTCNs > 0L)
    } # (for bb ...)
    verbose && exit(verbose);

    verbose && exit(verbose);
  } # for (jj ...)

  verbose && cat(verbose, "Bootstrap means");
  verbose && str(verbose, M);

  # Sanity check
  stopifnot(all(!is.nan(M)));


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Calculate (C1,C2) bootstrap statistics
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Calculating (C1,C2) from (TCN,DH) bootstraps");
  C1 <- (1-M[,,"dh"]) * M[,,"tcn"] / 2;
  C2 <- M[,,"tcn"] - C1;
  M[,,"c1"] <- C1;
  M[,,"c2"] <- C2;
  verbose && str(verbose, M);

  # Sanity check
  stopifnot(all(!is.nan(M)));
  verbose && exit(verbose);

  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Calculate bootstrap statistics
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  verbose && enter(verbose, "Calculating (TCN,DH) bootstrap statistics");
  # Allocate JxQx4 matrix S
  naValue <- as.double(NA);
  dim <- dim(M);
  dimnames <- dimnames(M);
  dim[2] <- nbrOfStats;
  dimnames[[2]] <- statsNames;
  S <- array(naValue, dim=dim, dimnames=dimnames);
  verbose && str(verbose, S);

  fields <- dimnames(M)[[3]];
  for (kk in seq(along=fields)) {
    field <- fields[kk];
    verbose && enter(verbose, sprintf("Field #%d ('%s') of %d", kk, field, length(fields)));

    Mkk <- M[,,kk,drop=FALSE];  # An JxB matrix
    dim(Mkk) <- dim(Mkk)[-3];
    # Sanity check
    stopifnot(is.matrix(Mkk));
    stopifnot(nrow(Mkk) == nbrOfSegments);
    stopifnot(ncol(Mkk) == B);

    for (jj in seq(length=nbrOfSegments)) {
      verbose && enter(verbose, sprintf("Segment #%d of %d", jj, nbrOfSegments));

      Mkkjj <- Mkk[jj,,drop=TRUE]; # A vector of length B

      S[jj,,kk] <- statsFcn(Mkkjj);

      verbose && exit(verbose);
    } # for (jj ...)

    verbose && exit(verbose);
  } # for (jj ...)
  verbose && exit(verbose);
  verbose && cat(verbose, "Bootstrap statistics");
  verbose && str(verbose, S);



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Store
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Reshape JxQx4 array to Jx(4*Q) matrix
  T <- wrap(S, map=list(1,NA), sep="_");
  colnames(T) <- gsub("(.*)_(.*)", "\\2_\\1", colnames(T));

  # Append
  segs <- cbind(segs, T);

  # Drop previously estimated values
  dups <- duplicated(colnames(segs), fromLast=TRUE);
  if (any(dups)) {
    segs <- segs[,!dups, drop=FALSE];
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Statistical sanity checks
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  if (B >= 100L) {
    verbose && enter(verbose, "Statistical sanity checks (iff B >= 100)");

    # Find extreme quantiles
    probs <- dimnames(S)[[2]];
    probs <- gsub("%", "", probs, fixed=TRUE);
    probs <- as.double(probs) / 100;

    # Is it possible to check?
    if (any(probs < 0.10) && any(probs > 0.90)) {
      tryCatch({
        fields <- dimnames(S)[[3]];
        for (kk in seq(along=fields)) {
          field <- fields[kk];
          verbose && enter(verbose, sprintf("Field #%d ('%s') of %d", kk, field, length(fields)));
          
          # Bootstrap statistics
          Skk <- S[,,kk, drop=FALSE];
          dim(Skk) <- dim(Skk)[-3];
  
          # Sanity checks
          stopifnot(is.matrix(Skk));
  
          range <- Skk[,c(1,ncol(Skk)),drop=FALSE];
    
          # Segmentation means
          key <- sprintf("%sMean", field);
          segMean <- segs[[key]];

          # Segmentation counts
          cfield <- sprintf("%sNbrOfLoci", ifelse(field == "tcn", "tcn", "dh"));
          counts <- segs[,cfield,drop=TRUE];
  
          if (verbose) {
            for (rr in seq_len(length(segMean))) {
              verbose && printf(verbose, "Seg %3d. mean=%g, range=[%g,%g], n=%d\n", rr, segMean[rr], range[rr,1], range[rr,2], counts[rr]);
            } # for (rr ...)
          }
  
          # Compare only segments with enough data points
          keep <- (counts > 1L);
          range <- range[keep,,drop=FALSE];
          segMean <- segMean[keep];  

          # Sanity checks
          stopifnot(all(range[,2] + tol >= range[,1], na.rm=TRUE));
          stopifnot(all(segMean + tol >= range[,1], na.rm=TRUE));
          stopifnot(all(segMean - tol <= range[,2], na.rm=TRUE));
    
          verbose && exit(verbose);
        } # for (kk ...)
      }, error = function(ex) {
        # If an error, display the data, then throw the exception
        verbose && cat(verbose, "Tolerance (option 'PSCBS/sanityChecks/tolerance'): ", tol);
        verbose && print(verbose, segs);
        throw(ex);
      })
    } else {
      verbose && cat(verbose, "Skipping. Not enough quantiles: ",
                               paste(dimnames(S)[[2]], collapse=", "));
    }

    verbose && exit(verbose);
  } # if (B >= 100L)
  
  
  fitB <- fit;
  fitB$output <- segs;

  verbose && exit(verbose);

  fitB;
}, private=TRUE) # bootstrapTCNandDHByRegion()



##############################################################################
# HISTORY
# 2013-02-09
# o BUG FIX: bootstrapTCNandDHByRegion() for PairedPSCBS did not bootstrap
#   for all available loci when calculating total CNs statistics, iff
#   the segment had been called run-of-homozygosity (ROH).
#   Thanks to Oscar Rueda at the Cancer Research UK Cambridge Institute 
#   for reporting on this.
# 2013-02-07
# o Improved some verbose outputs of bootstrapTCNandDHByRegion().
# 2013-01-15
# o Now bootstrapTCNandDHByRegion() uses the params$avgTCN and params$avgDH
#   estimators, iff given.
# 2012-11-05
# o GENERALIZATION: Now bootstrapTCNandDHByRegion() works for more "flavors",
#   including the default ('tcn') used by segmentByNonPairedPSCBS().
# 2012-09-20
# o SPEEDUP: By default bootstrapTCNandDHByRegion() for PairedPSCBS no
#   longer do sanity checks within the bootstrap loop.  This significantly
#   speed up the method.  To run checks, use argument .debug=TRUE.
# 2012-02-26
# o BUG FIX: bootstrapTCNandDHByRegion() for PairedPSCBS would resample
#   from a subset of the intended TCNs, iff the DH mean was non-finite
#   while there were still heterozygous SNPs.  This introduced a bias in
#   the estimates, which was neglectable for large segments, but for very
#   small segments (a few loci) it could be relatively large.
# 2012-02-24
# o BUG FIX: bootstrapTCNandDHByRegion(..., force=TRUE) for PairedPSCBS
#   would give an error iff previous bootstrap estimates already existed.
# o Added argument 'force' to bootstrapTCNandDHByRegion().
# 2011-11-26
# o Now bootstrapTCNandDHByRegion() for PairedPSCBS preserves NAs for DH
#   and (C1,C2) quantiles, if the DH mean level is NA, which can happen
#   when a segment is called ROH.
# o An internal sanity check of bootstrapTCNandDHByRegion() for PairedPSCBS
#   would give an error if DH mean levels had been set to NA for segments
#   called ROH.
# 2011-11-24
# o BUG FIX: bootstrapTCNandDHByRegion() for PairedPSCBS would give
#   an error, if a segment did not have any TCN signals, which can
#   occur when known segments are specified for Paired PSCBS.
# 2011-08-08
# o Moved the sanity checks that tests the TCN and DH "segRows" from the
#   bootstrapTCNandDHByRegion() to segmentByPairedPSCBS().  This is the
#   first step to fix a failure in the sanity checks that could happend
#   iff one first run dropSegmentationOutliers().
# 2011-06-14
# o Updated code to recognize new column names.
# 2011-05-29
# o Renamed options to reflect new package name.
# 2010-12-03
# o BUG FIX: In rare cases the bootstrap sanity checks can indeed produce
#   an invalid 'range', more precisely where (range[,2] >= range[,1]) is
#   not true.  This can happen if there is no variation in the bootstrap
#   estimates.  Beause of this we allow for some tolerance.
# 2010-12-02
# o Now bootstrapTCNandDHByRegion() uses option
#   "psCBS/sanityChecks/tolerance".
# 2010-12-01
# o BUG FIX: bootstrapTCNandDHByRegion() did not always exclude the correct
#   loci.
# 2010-11-27
# o BUG FIX: bootstrapTCNandDHByRegion() would incorrectly include
#   non-polymorphic loci in the set of homozygous SNPs during resampling.
# 2010-11-26
# o BUG FIX: The statistical sanity checks of the bootstrap estimates would
#   give an error when only single-sided bootstrap confidence interval was
#   calculated.
# 2010-11-23
# o ROBUSTNESS: Added more sanity checks to bootstrapTCNandDHByRegion().
# o WORKAROUND: The precision of the mean levels of DNAcopy::segment()
#   is not great enough to always compare it to that of R's estimates.
# o BUG FIX: bootstrapTCNandDHByRegion() would give an error if there was
#   only one segment.
# 2010-11-22
# o BUG FIX: The DH segmentation and bootstrap incorrectly included
#   missing values, when subseting.
# o BUG FIX: Some sanity checks were incorrect.
# o BUG FIX: bootstrapTCNandDHByRegion() for PairedPSCBS would not correctly
#   detect if bootstrap results are already available.
# 2010-11-21
# o Added argument 'seed'.
# o Added bootstrapTCNandDHByRegion() for PairedPSCBS.
# o Created from PairedPSCBS.BOOT.R.
##############################################################################
