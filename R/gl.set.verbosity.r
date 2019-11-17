#' Set the default verbosity level
#'
#' dartR functions have a verbosity parameter that sets the level of reporting during the execution of the function. The verbosity level,
#' set by parameter 'verbose' can be one of verbose 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary;
#' 5, full report. The default value for verbosity is stored as a flag in the gl@loc.metrics.flags slot of each genelight object. This script 
#' sets the value of the flag.
#' 
#' @param  x name of the genlight object containing the SNP data, or the genind object containing the SilocoDArT data [required]
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @param set.verbosity -- set the default verbosity to be: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @return The genlight with the verbosity flag reset
#' @export
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' gl <- gl.set.verbosity(testset.gl, value=2)

gl.set.verbosity <- function(x, value=2, verbose=NULL) {

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
  
  if (verbose >= 2){
    if (all(x@ploidy == 1)){
      cat("  Processing Presence/Absence (SilicoDArT) data\n")
    } else if (all(x@ploidy == 2)){
      cat("  Processing a SNP dataset\n")
    } else {
      stop("Fatal Error: Ploidy must be universally 1 (fragment P/A data) or 2 (SNP data)")
    }
  }
  
# DO THE JOB
  
  if (verbose >= 2){
    cat(paste("  Setting default verbosity to",value,"\n"))
  }
  x@other$verbose <- value
  
# ADD TO HISTORY
  
  nh <- length(x@other$history)
  x@other$history[[nh + 1]] <- match.call()      
  
# FLAG SCRIPT END
  
  if (verbose >= 1) {
    cat("Completed:",funname,"\n")
  }
  
  return(x)
}
  
