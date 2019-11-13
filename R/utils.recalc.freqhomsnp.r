#' A utility script to recalculate the the frequency of the homozygous alternate SNP by locus after some populations have been deleted
#'
#' The locus metadata supplied by DArT has FreqHomSnp included,
#' but the frequency of the homozygous alternate will change when some individuals are removed from the dataset. 
#' This script recalculates the FreqHomSnp and places these recalculated values in the appropriate place in the genlight object.
#' Note that the frequency of the homozygote alternate SNPS is calculated from the individuals that could be scored.
#'
#' @param x -- name of the genlight object containing the SNP data [required]
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @return The modified genlight object
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' @seealso \code{utils.recalc.metrics} for recalculating all metrics, \code{utils.recalc.callrate} for recalculating CallRate,
#' \code{utils.recalc.freqhomref} for recalculating frequency of homozygous reference, \code{utils.recalc.avgpic} for recalculating AvgPIC,
#' \code{utils.recalc.freqhet} for recalculating frequency of heterozygotes, \code{gl.recalc.maf} for recalculating minor allele frequency,
#' \code{gl.recalc.rdepth} for recalculating average read depth
#' @examples
#' #result <- utils.recalc.freqhomsnp(testset.gl)

utils.recalc.freqhomsnp <- function(x, verbose=NULL) {

# TIDY UP FILE SPECS
  
  build <- "Jacob"
  funname <- match.call()[[1]]
  hold <- x
  # Note draws upon and modifies the loc.metrics.flags for freqhomsnp.
  
# FLAG SCRIPT START
  # set verbosity
  if (is.null(verbose) & !is.null(x@other$verbose)) verbose=x@other$verbose
  if (is.null(verbose)) verbose=2
 
  
  if (verbose < 0 | verbose > 5){
    cat("  Warning: Parameter 'verbose' must be an integer between 0 [silent] and 5 [full report], set to 2\n")
    verbose <- 2
  }
  
  if (verbose >= 1){
    cat("Starting",funname,"\n")
  }
  
# STANDARD ERROR CHECKING
  
  if(class(x)!="genlight") {
    cat("  Fatal Error: genlight object required!\n"); stop("Execution terminated\n")
  }
  
  if (all(x@ploidy == 1)){
    stop("  Processing  Presence/Absence (SilicoDArT) data\n")
    data.type <- "SilicoDArT"
  } else if (all(x@ploidy == 2)){
    if (verbose >= 2){cat("  Processing a SNP dataset\n")}
    data.type <- "SNP"
  } else {
    stop("Fatal Error: Ploidy must be universally 1 (fragment P/A data) or 2 (SNP data)")
  }
  
  # Check monomorphs have been removed up to date
  if (x@other$loc.metrics.flags$monomorphs == FALSE){
    if (verbose >= 2){
      cat("  Warning: Dataset contains monomorphic loci which will be included in the ",funname," calculations\n")
    }  
  }

# FUNCTION SPECIFIC ERROR CHECKING

  if (is.null(x@other$loc.metrics$FreqHomSnp)) {
    x@other$loc.metrics$FreqHomSnp <- array(NA,nLoc(x))
    if (verbose >= 2){
      cat("  Locus metric FreqHomSnp does not exist, creating slot @other$loc.metrics$FreqHomSnp\n")
    }
  }

# DO THE JOB

     t <- as.matrix(x)
     if (verbose >= 2) {cat("  Recalculating locus metric freqHomSnp\n")}
     x@other$loc.metrics$FreqHomSnp <- colMeans(t==2, na.rm = T)
     x@other$loc.metrics.flags$FreqHomSnp <- TRUE
     
# FLAG SCRIPT END

  if (verbose > 0) {
    cat("Completed:",funname,"\n")
  }
   
   return(x)
}
