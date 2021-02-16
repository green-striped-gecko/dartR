#' Filter loci on the basis of numeric information stored in $other$loc.metrics in a genlight \{adegenet\} object
#'
#' This script uses any field with numeric values stored in $other$loc.metrics to filter loci. The loci to keep can be within the upper and lower thresholds ("within") or outside of the upper and lower thresholds ("outside").
#' The fields that are included in dartR, and a short description, are found below. Optionally, the user can also set his/her own filter by adding a vector into $other$loc.metrics as shown in the example.
#' 
#' - SnpPosition - position (zero is position 1) in the sequence tag of the defined SNP variant base 
#' - CallRate - proportion of samples for which the genotype call is non‐ missing (that is, not “‐” ) 
#' - OneRatioRef - proportion of samples for which the genotype score is 0 
#' - OneRatioSnp - proportion of samples for which the genotype score is 2 
#' - FreqHomRef - proportion of samples homozygous for the Reference allele 
#' - FreqHomSnp - proportion of samples homozygous for the Alternate (SNP) allele 
#' - FreqHets - proportion of samples which score as heterozygous, that is, scored as 1 
#' - PICRef - polymorphism information content (PIC) for the Reference allele 
#' - PICSnp - polymorphism information content (PIC) for the SNP 
#' - AvgPIC - average of the polymorphism information content (PIC) of the Reference and SNP alleles 
#' - AvgCountRef - sum of the tag read counts for all samples, divided by the number of samples with non‐zero tag read counts, for the Reference allele row 
#' - AvgCountSnp - sum of the tag read counts for all samples, divided by the number of samples with non‐zero tag read counts, for the Alternate (SNP) allele row 
#' - RepAvg - proportion of technical replicate assay pairs for which the marker score is consistent 
#' 
#' @param x -- name of the genlight object containing the SNP data [required]
#' @param metric -- name of the metric to be used for filtering [required]
#' @param upper -- filter upper threshold [required]
#' @param lower -- filter lower threshold  [required]
#' @param keep -- whether keep loci within of upper and lower thresholds or keep loci outside of upper and lower thresholds [within]
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2, unless specified using gl.set.verbosity]
#' @return The reduced genlight dataset
#' @export
#' @author Luis Mijangos (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' # adding dummy data
#' test <- testset.gl
#' test$other$loc.metrics$test <- 1:nLoc(test)
#' result <- gl.filter.locmetric(x=test, metric= "test", upper=255, lower=200, keep= "within", verbose=3)

gl.filter.locmetric <- function(x, metric, upper, lower, keep="within", verbose=2) {
  
# TRAP COMMAND, SET VERSION
  
  funname <- match.call()[[1]]
  build <- "Jacob"

# SET VERBOSITY
  
  if (is.null(verbose)){ 
    if(!is.null(x@other$verbose)){ 
      verbose <- x@other$verbose
    } else { 
      verbose <- 2
    }
  } 
  
  if (verbose < 0 | verbose > 5){
    cat(paste("  Warning: Parameter 'verbose' must be an integer between 0 [silent] and 5 [full report], set to 2\n"))
    verbose <- 2
  }
  
# FLAG SCRIPT START
  
  if (verbose >= 1){
    if(verbose==5){
      cat("Starting",funname,"[ Build =",build,"]\n")
    } else {
      cat("Starting",funname,"\n")
    }
  }

# STANDARD ERROR CHECKING
  
  if(class(x)!="genlight") {
    stop("  Fatal Error: genlight object required!\n")
  }

  # Work around a bug in adegenet if genlight object is created by subsetting
      if (nLoc(x)!=nrow(x@other$loc.metrics)) { stop("The number of rows in the loc.metrics table does not match the number of loci in your genlight object!")  }

  # Set a population if none is specified (such as if the genlight object has been generated manually)
    if (is.null(pop(x)) | is.na(length(pop(x))) | length(pop(x)) <= 0) {
      if (verbose >= 2){ cat("  Population assignments not detected, individuals assigned to a single population labelled 'pop1'\n")}
      pop(x) <- array("pop1",dim = nInd(x))
      pop(x) <- as.factor(pop(x))
    }

  # Check for monomorphic loci
    tmp <- gl.filter.monomorphs(x, verbose=0)
    if ((nLoc(tmp) < nLoc(x)) & verbose >= 2) {cat("  Warning: genlight object contains monomorphic loci\n")}

# FUNCTION SPECIFIC ERROR CHECKING
# check whether the field exists in the genlight object
  if ( !(metric %in% colnames(x$other$loc.metrics)) ) {
    stop("  Fatal Error: name of the metric not found\n")
  }

# DO THE JOB
  
  field <- which(colnames(x@other$loc.metrics) == metric)
  
  if(keep=="within"){
  index <- which(x@other$loc.metrics[,field] >= lower & x@other$loc.metrics[,field] <= upper)
  x2 <- x[ ,index]
  x2@other$loc.metrics <- x@other$loc.metrics[index,]
}
  if(keep=="outside"){
  index <- which(x@other$loc.metrics[,field] <= lower & x@other$loc.metrics[,field] >= upper)
  x2 <- x[ ,index]
  x2@other$loc.metrics <- x@other$loc.metrics[index,]
  }
 
  if (verbose > 2) {
    cat("    Initial number of loci:", nLoc(x), "\n")
    cat("    Number of loci deleted:", nLoc(x) - nLoc(x2), "\n")
    cat("    Final number of loci:", nLoc(x2), "\n")
  }

# ADD TO HISTORY  
  nh <- length(x2@other$history)
  x2@other$history[[nh + 1]] <- match.call()
  
# FLAG SCRIPT END

  if (verbose > 0) {
    cat("Completed:",funname,"\n")
  }

  return(x2)
}  
