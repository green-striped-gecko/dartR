#' A utility function to check for potential sources of errors in genlight objects.
#'
#' This function checks for errors in verbosity, loc.metrics, ind.metrics, ploidy, metadata, populations, individuals and loci.
#'
#' @param x -- name of the genlight object containing the SNP data or tag presence/absence data (SilicoDArT) [required]
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default NULL]
#' @return The modified genlight object
#' @examples
#' test <- utils.check.gl(testset.gl)

utils.check.gl <- function(x, verbose = NULL) {
  #### SET VERBOSITY ####
  if (is.null(verbose)) {
    if (!is.null(x@other$verbose)) {
      verbose <- x@other$verbose
    } else {
      verbose <- 2
      x@other$verbose <- 2
    }
  }
  if (verbose < 0 | verbose > 5) {
    message(paste(
      warn(
        "  Warning: Parameter 'verbose' must be an integer between 0 [silent] and 5 [full report], set to 2"
      )
    ))
    verbose <- 2
    x@other$verbose <- 2
  }
  
  #### CHECK GENLIGHT OBJECT ####
  if (class(x) != "genlight") {
    stop(error("Fatal Error: genlight object required!"))
  }
  
  #### CHECK METRICS ####
  # Check if the x@other$loc.metrics slot exists, if not, create as a dataframe
  if (is.null(x@other$loc.metrics)) {
    x@other$loc.metrics <- as.data.frame(array(NA, nLoc(x)))
  }
  x@other$loc.metrics <- data.frame(x@other$loc.metrics)
  # Check if the x@other$loc.metrics.flags slot exists, if not, create a dataframe
  if (is.null(x@other$loc.metrics.flags)) {
    x@other$loc.metrics.flags <- as.data.frame(array(NA, 1))
  }
  # Check if the x@other$ind.metrics slot exists, if not, create a dataframe
  if (verbose >= 2) {
    message(report("  Checking for individual metrics"))
  }
  if (is.null(x@other$ind.metrics)) {
    if (verbose >= 1) {
      message(report("    Creating a slot for individual metrics"))
    }
    x@other$ind.metrics <- as.data.frame(array(NA, nInd(x)))
    x@other$ind.metrics$id <- indNames(x)
  } else{
    if (verbose >= 1) {
      message(report("    Individual metrics confirmed"))
    }
  }
  # Check for the locus metrics, and create if they do not exist.
  # Check for the locus metrics flags, and create if they do not exist.
  # Check for the verbosity flag, and create if it does not exist.
  
  # if (verbose >= 2){cat("  Checking locus metrics and flags\n")}
  # x <- utils.reset.flags(x,set=FALSE,verbose=0)
  # 
  # # Calculate locus metrics
  # 
  # if (verbose >= 2){cat("  Recalculating locus metrics\n")}
  # x <- gl.recalc.metrics(x,verbose=0)
  # 
  # # Check for monomorphic loci
  # if (verbose >= 2){cat("  Checking for monomorphic loci\n")}
  # x2 <- gl.filter.monomorphs(x,verbose=0)
  # if(nLoc(x2)==nLoc(x)){
  #   if (verbose >= 1){cat("    No monomorphic loci detected\n")}
  #   x@other$loc.metrics.flags$monomorphs <- TRUE
  # } else {
  #   if (verbose >= 1){cat("    Dataset containes monomorphic loci\n")}
  #   x@other$loc.metrics.flags$monomorphs <- FALSE
  # }
  
  
  
  #### CHECK FLAGS ####
  # Check AvgPIC
  # if (is.null(x@other$loc.metrics$AvgPIC)) {
  #   x@other$loc.metrics.flags$AvgPIC <- F
  #   if (verbose >= 3) {
  #     message(warn("  Locus metric AvgPIC does not exist"))
  #   }
  # } else{
  #   x@other$loc.metrics.flags$AvgPIC <- T
  # }
  # # Check OneRatioRef
  # if (is.null(x@other$loc.metrics$OneRatioRef)) {
  #   x@other$loc.metrics.flags$OneRatioRef <- F
  #   if (verbose >= 3) {
  #     message(warn("  Locus metric OneRatioRef does not exist"))
  #   }
  # } else{
  #   x@other$loc.metrics.flags$OneRatioRef <- T
  # }
  # # Check OneRatioSnp
  # if (is.null(x@other$loc.metrics$OneRatioSnp)) {
  #   x@other$loc.metrics.flags$OneRatioSnp <- F
  #   if (verbose >= 3) {
  #     message(warn("  Locus metric OneRatioSnp does not exist"))
  #   }
  # } else{
  #   x@other$loc.metrics.flags$OneRatioSnp <- T
  # }
  # # Check PICRef
  # if (is.null(x@other$loc.metrics$PICRef)) {
  #   x@other$loc.metrics.flags$PICRef <- F
  #   if (verbose >= 3) {
  #     message(warn("  Locus metric PICRef does not exist"))
  #   }
  # } else{
  #   x@other$loc.metrics.flags$PICRef <- T
  # }
  # # Check PICSnp
  # if (is.null(x@other$loc.metrics$PICSnp)) {
  #   x@other$loc.metrics.flags$PICSnp <- F
  #   if (verbose >= 3) {
  #     message(warn("  Locus metric PICSnp does not exist"))
  #   }
  # } else{
  #   x@other$loc.metrics.flags$PICSnp <- T
  # }
  # # Check CallRate
  # if (is.null(x@other$loc.metrics$CallRate)) {
  #   x@other$loc.metrics.flags$CallRate <- F
  #   if (verbose >= 3) {
  #     message(warn("  Locus metric CallRate does not exist"))
  #   }
  # } else{
  #   x@other$loc.metrics.flags$CallRate <- T
  # }
  # # Check FreqHomRef
  # if (is.null(x@other$loc.metrics$FreqHomRef)) {
  #   x@other$loc.metrics.flags$FreqHomRef <- F
  #   if (verbose >= 3) {
  #     message(warn("  Locus metric FreqHomRef does not exist"))
  #   }
  # } else{
  #   x@other$loc.metrics.flags$FreqHomRef <- T
  # }
  # # Check FreqHomSnp
  # if (is.null(x@other$loc.metrics$FreqHomSnp)) {
  #   x@other$loc.metrics.flags$FreqHomSnp <- F
  #   if (verbose >= 3) {
  #     message(warn("  Locus metric FreqHomSnp does not exist"))
  #   }
  # } else{
  #   x@other$loc.metrics.flags$FreqHomSnp <- T
  # }
  # # Check FreqHets
  # if (is.null(x@other$loc.metrics$FreqHets)) {
  #   x@other$loc.metrics.flags$FreqHets <- F
  #   if (verbose >= 3) {
  #     message(warn("  Locus metric FreqHets does not exist"))
  #   }
  # } else{
  #   x@other$loc.metrics.flags$FreqHets <- T
  # }
  # # Check monomorphs
  # if (is.null(x@other$loc.metrics$monomorphs)) {
  #   x@other$loc.metrics.flags$monomorphs <- F
  #   if (verbose >= 3) {
  #     message(warn("  Locus metric monomorphs does not exist"))
  #   }
  # } else{
  #   x@other$loc.metrics.flags$monomorphs <- T
  # }
  # # Check maf
  # if (is.null(x@other$loc.metrics$maf)) {
  #   x@other$loc.metrics.flags$maf <- F
  #   if (verbose >= 3) {
  #     message(warn("  Locus metric maf does not exist"))
  #   }
  # } else{
  #   x@other$loc.metrics.flags$maf <- T
  # }
  # # Check OneRatio
  # if (is.null(x@other$loc.metrics$OneRatio)) {
  #   x@other$loc.metrics.flags$OneRatio <- F
  #   if (verbose >= 3) {
  #     message(warn("  Locus metric OneRatio does not exist"))
  #   }
  # } else{
  #   x@other$loc.metrics.flags$OneRatio <- T
  # }
  # # Check PIC
  # if (is.null(x@other$loc.metrics$PIC)) {
  #   x@other$loc.metrics.flags$PIC <- F
  #   if (verbose >= 3) {
  #     message(warn("  Locus metric PIC does not exist"))
  #   }
  # } else{
  #   x@other$loc.metrics.flags$PIC <- T
  # }
  
  #### CHECK PLOIDY ####
  if (all(x@ploidy == 1)) {
    message(report("  Processing a Tag Presence/Absence (SilicoDArT) dataset"))
  } else if (all(x@ploidy == 2)) {
    message(report("  Processing a SNP dataset"))
  } else {
    stop(error(
      "Fatal Error: Ploidy must be universally 1 (fragment P/A data) or 2 (SNP data)!"
    ))
  }
  
  #### CHECK POPULATIONS ####
  # Set a population if none is specified (such as if the genlight object has been generated manually)
  if (is.null(pop(x)) |
      is.na(length(pop(x))) | length(pop(x)) <= 0) {
    if (verbose >= 2) {
      message(
        warn(
          "  Population assignments not detected, individuals assigned to a single population labelled 'pop1'"
        )
      )
    }
    pop(x) <- array("pop1", dim = nInd(x))
    pop(x) <- as.factor(pop(x))
  }
  
  #### CHECK INDIVIDUALS ####
  # Check duplicated individual names
  
  #### CHECK LOCI ####
  # Check that the number of values in the loc.metrics dataframe is the same as the number of loci
  if (nLoc(x) != nrow(x@other$loc.metrics)) {
    message(
      warn(
        "  The number of rows in the loc.metrics table does not match the number of loci! This is potentially a major problem if there is a mismatch of the loci with the metadata. Trace back to identify the cause."
      )
    )
  }
  # Check monomorphic loci
  if (verbose >= 2) {
    message(report("  Checking for monomorphic loci"))
  }
  x2 <- gl.filter.monomorphs(x, verbose = 0)
  if (nLoc(x2) == nLoc(x)) {
    if (verbose >= 1) {
      message(report("    No monomorphic loci detected"))
    }
    x@other$loc.metrics.flags$monomorphs <- TRUE
  } else {
    if (verbose >= 1) {
      message(
        warn(
          "  Warning: Dataset contains monomorphic loci which will be included in this function calculations"
        )
      )
    }
    x@other$loc.metrics.flags$monomorphs <- FALSE
  }
  
  #### CHECK METADATA ####
  #check if coordinates are in the right place and not misspelled
  # if (!is.null(x@other$latlon))
  #   x@other$latlong <- x@other$latlon
  # 
  # if (!is.null(x@other$latlong)) {
  #   if (!is.null(x@other$latlong$long))
  #     x@other$latlong$lon <- x@other$latlong$long
  # }
  # #remove misspelled columns if they exist...
  # x@other$latlon <- NULL
  # x@other$latlong$long <- NULL
  # if (verbose >= 2) {
  #   message(report("  Spelling of coordinates checked and changed if necessary"))
  # }
  
  message(report("Checking genlight object finalised\n"))
  
  return(x)
  
}
