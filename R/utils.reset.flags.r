#' A utility script to reset to FALSE (or TRUE) the locus metric flags after some individuals or populations have been deleted.
#'
#' The locus metadata supplied by DArT has OneRatioRef, OneRatioSnp, PICRef, PICSnp, and AvgPIC included,
#' but the allelec composition will change when some individuals are removed from the dataset and so the initial statistics will
#' no longer apply. This script resets the locus metrics flags to FALSE to indicate that these statistics in the genlight object
#' are no longer current.
#' 
#' If the locus metrics do not exist then they are added to the genlight object but not populated. If the locus metrics flags 
#' do not exist, then they are added to the genlight object and set to FALSE (or TRUE).
#'
#' @param x -- name of the genlight object containing the SNP data or tag presence/absence data (SilicoDArT) [required]
#' @param set -- set the flags to TRUE or FALSE [FALSE]
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @return The modified genlight object
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' @seealso \code{utils.recalc.metrics} for recalculating all metrics, \code{utils.recalc.callrate} for recalculating CallRate,
#' \code{utils.recalc.freqhomref} for recalculating frequency of homozygous reference, \code{utils.recalc.freqhomsnp} for recalculating frequency of homozygous alternate,
#' \code{utils.recalc.freqhet} for recalculating frequency of heterozygotes, \code{gl.recalc.maf} for recalculating minor allele frequency,
#' \code{gl.recalc.rdepth} for recalculating average read depth
#' @examples
#' #result <- utils.reset.flags(testset.gl)

utils.reset.flags <- function(x, set=FALSE, verbose=2) {
  
# TIDY UP FILE SPECS
  
  build <- "Jacob"
  funname <- match.call()[[1]]
  
# FLAG SCRIPT START
  
  if (verbose < 0 | verbose > 5){
    cat("  Warning: Parameter 'verbose' must be an integer between 0 [silent] and 5 [full report], set to 2\n")
    verbose <- 2
  }
  
  if (verbose >= 1){
    cat("Starting",funname,"[ Build =",build,"]\n")
  }
  
# STANDARD ERROR CHECKING
  
  if(class(x)!="genlight") {
    stop("  Fatal Error: genlight object required!\n")
  }
  
  if (all(x@ploidy == 1)){
    if (verbose >= 2){cat("  Processing  Presence/Absence (SilicoDArT) data\n")}
    data.type <- "SilicoDArT"
  } else if (all(x@ploidy == 2)){
    if (verbose >= 2){cat("  Processing a SNP dataset\n")}
    data.type <- "SNP"
  } else {
    stop("Fatal Error: Ploidy must be universally 1 (fragment P/A data) or 2 (SNP data)")
  }
  
  # DO THE JOB
  
  if (all(x@ploidy == 2)){
    if(verbose >= 2){
      cat("  Resetting flags for AvgPIC, OneRatioRef, OneRatioSnp, PICRef, PICSnp, CallRate, maf, FreqHets, FreqHomRef, FreqHomSnp, monomorphs\n")
    }
    if (is.null(x@other$loc.metrics$AvgPIC)) {
      x@other$loc.metrics$AvgPIC <- array(NA,nLoc(x))
      if (verbose >= 3){
        cat("  Locus metric AvgPIC does not exist, creating slot @other$loc.metrics$AvgPIC\n")
      }
    }
    x@other$loc.metrics.flags$AvgPIC <- set
    
    if (is.null(x@other$loc.metrics$OneRatioRef)) {
      x@other$loc.metrics$OneRatioRef <- array(NA,nLoc(x))
      if (verbose >= 3){
        cat("  Locus metric OneRatioRef does not exist, creating slot @other$loc.metrics$OneRatioRef\n")
      }
    }
    x@other$loc.metrics.flags$OneRatioRef <- set
    
    if (is.null(x@other$loc.metrics$OneRatioSnp)) {
      x@other$loc.metrics$OneRatioSnp <- array(NA,nLoc(x))
      if (verbose >= 3){
        cat("  Locus metric OneRatioSnp does not exist, creating slot @other$loc.metrics$OneRatioSnp\n")
      }
    }
    x@other$loc.metrics.flags$OneRatioSnp <- set
    
    if (is.null(x@other$loc.metrics$PICRef)) {
      x@other$loc.metrics$PICRef <- array(NA,nLoc(x))
      if (verbose >= 3){
        cat("  Locus metric PICRef does not exist, creating slot @other$loc.metrics$PICRef\n")
      }
    }
    x@other$loc.metrics.flags$PICRef <- set
    
    if (is.null(x@other$loc.metrics$PICSnp)) {
      x@other$loc.metrics$PICSnp <- array(NA,nLoc(x))
      if (verbose >= 3){
        cat("  Locus metric PICSnp does not exist, creating slot @other$loc.metrics$PICSnp\n")
      }
    }
    x@other$loc.metrics.flags$PICSnp <- set
    
    if (is.null(x@other$loc.metrics$CallRate)) {
      x@other$loc.metrics$CallRate <- array(NA,nLoc(x))
      if (verbose >= 3){
        cat("  Locus metric CallRate does not exist, creating slot @other$loc.metrics$CallRate\n")
      }
    }
    x@other$loc.metrics.flags$CallRate <- set
    
    if (is.null(x@other$loc.metrics$FreqHomRef)) {
      x@other$loc.metrics$FreqHomRef <- array(NA,nLoc(x))
      if (verbose >= 3){
        cat("  Locus metric FreqHomRef does not exist, creating slot @other$loc.metrics$FreqHomRef\n")
      }
    }
    x@other$loc.metrics.flags$FreqHomRef <- set
    
    if (is.null(x@other$loc.metrics$FreqHomSnp)) {
      x@other$loc.metrics$FreqHomSnp <- array(NA,nLoc(x))
      if (verbose >= 3){
        cat("  Locus metric FreqHomSnp does not exist, creating slot @other$loc.metrics$FreqHomSnp\n")
      }
    }
    x@other$loc.metrics.flags$FreqHomSnp <- set
    
    if (is.null(x@other$loc.metrics$FreqHets)) {
      x@other$loc.metrics$FreqHets <- array(NA,nLoc(x))
      if (verbose >= 3){
        cat("  Locus metric FreqHets does not exist, creating slot @other$loc.metrics$FreqHets\n")
      }
    }
    x@other$loc.metrics.flags$FreqHets <- set
      
    x@other$loc.metrics.flags$monomorphs <- set
    x@other$loc.metrics.flags$maf <- set

  }
  
  if (all(x@ploidy == 1)){
    if(verbose >= 2){
      cat("  Resetting flags for CallRate, PIC, OneRatio, monomorphs\n")
    }
    
    if (is.null(x@other$loc.metrics$PIC)) {
      x@other$loc.metrics$PIC <- array(NA,nLoc(x))
      if (verbose >= 3){
        cat("  Locus metric PIC does not exist, creating slot @other$loc.metrics$PIC\n")
      }
    }
    x@other$loc.metrics.flags$PIC <- set
    
    if (is.null(x@other$loc.metrics$OneRatio)) {
      x@other$loc.metrics$OneRatio <- array(NA,nLoc(x))
      if (verbose >= 3){
        cat("  Locus metric OneRatio does not exist, creating slot @other$loc.metrics$OneRatioRef\n")
      }
    }
    x@other$loc.metrics.flags$OneRatio <- set
    
    if (is.null(x@other$loc.metrics$CallRate)) {
      x@other$loc.metrics$CallRate <- array(NA,nLoc(x))
      if (verbose >= 3){
        cat("  Locus metric CallRate does not exist, creating slot @other$loc.metrics$CallRate\n")
      }
    }
    x@other$loc.metrics.flags$CallRate <- set
    
    if (is.null(x@other$loc.metrics$Qpmr)) {
      x@other$loc.metrics$Qpmr <- array(NA,nLoc(x))
      if (verbose >= 3){
        cat("  Locus metric Qpmr does not exist, creating slot @other$loc.metrics$Qpmr\n")
      }
    }
    x@other$loc.metrics.flags$Qpmr <- set
    
    x@other$loc.metrics.flags$monomorphs <- set

  }
  
  # FLAG SCRIPT END
  
  #add to history
  nh <- length(x@other$history)
  x@other$history[[nh + 1]] <- match.call() 
  
  if (verbose > 0) {
    cat("Completed:", funname, "\n")
  }
  
  return(x)
  
}

# # Test script
# gl <- testset.gl
# tmp <- utils.reset.flags(gl)
# gl@other$loc.metrics.flags
# tmp@other$loc.metrics.flags
# 
# tmp <- utils.reset.flags(gl,set=TRUE)
# gl@other$loc.metrics.flags
# tmp@other$loc.metrics.flags
# 
# 
# gs <- testset.gs
# tmp <- utils.reset.flags(gs)
# gs@other$loc.metrics.flags
# tmp@other$loc.metrics.flags
#   
# tmp <- utils.reset.flags(gs, set=TRUE)
# gs@other$loc.metrics.flags
# tmp@other$loc.metrics.flags
  