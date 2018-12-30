#' A utility script to calculate the read depth by locus, if it does not already exist
#'
#' The locus metadata supplied by DArT does not have red depth included, so it is calculated and
#' added to the locus metrics by this script. It calculates the read depth by adding the counts of 
#' the sequence tags for the reference allele (AvgCountRef provided by DArT PL) to the counts of the sequence tags for the alternate
#' (snp) allele (AvgCountSnp).
#'
#' @param x -- name of the genlight object containing the SNP data [required]
#' @param v -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @return The modified genlight dataset
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' #f <- dartR:::utils.recalc.rdepth(testset.gl)


utils.recalc.rdepth <- function(x, v=2) {
  
  if(class(x)!="genlight") {
    cat("Fatal Error: genlight object required!\n"); stop("Execution terminated\n")
  }
  # Work around a bug in adegenet if genlight object is created by subsetting
  x@other$loc.metrics <- x@other$loc.metrics[1:nLoc(x),]
  
  if (v > 0) {
    cat("Starting gl.recalc.rdepth: Read depth\n")
  }
  
  if (v < 0 | v > 5){
    cat("Warning: Verbosity must take on an integer value between 0 and 5, set to 3\n")
    v <- 3
  }
  
  # Calculate Read Depth
  
  if (v >= 3) {cat("Calculating Read Depth\n")}
  if (is.null(x@other$loc.metrics$rdepth)) {
    if (v >= 3){
      cat("  Locus metric rdepth does not exist, creating slot @other$loc.metrics$rdepth\n")
      cat("  Calculating read depth as the sum of AvgCountRef and AvgCountSnp\n")
    }
    x@other$loc.metrics$rdepth <- array(NA,nLoc(x))
    for (i in 1:nLoc(x)){
      x@other$loc.metrics$rdepth[i] <- x@other$loc.metrics$AvgCountRef[i] + x@other$loc.metrics$AvgCountSnp[i]
    }
  } else {
    if (v >= 3){
      cat("  Warning: Read depth metric already exists, no action taken\n")
    }
  }

  if (v > 0) {
    cat("Completed: gl.recalc.rdepth\n\n")
  }
  
  return(x)
}  
