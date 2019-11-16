#' Subsample n loci from a genlight object and return as a genlight object
#'
#' This is a support script, to subsample a genlight \{adegenet\} object based on loci. Two methods are used
#' to subsample, random and based on information content (avgPIC).
#'
#' @param x -- name of the genlight object containing the SNP data [required]
#' @param n -- number of loci to include in the subsample [required]
#' @param method -- "random", in which case the loci are sampled at random; or avgPIC, in which case the top n loci
#' ranked on information content (AvgPIC) are chosen [default 'random']
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @return A genlight object with n loci
#' @export
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' result <- gl.subsample.loci(testset.gl, n=200, method="pic")

gl.subsample.loci <- function(x, n, method="random", verbose=NULL) {

  # TRAP COMMAND, SET VERSION
  
  funname <- match.call()[[1]]
  build <- "Jacob"
  hold <- x
  
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
    stop("Fatal Error: genlight object required!\n")
  }
  
  if (verbose >= 2){
    if (all(x@ploidy == 1)){
      cat("  Processing  Presence/Absence (SilicoDArT) data\n")
      pic <- x@other$loc.metrics$PIC
    } else if (all(x@ploidy == 2)){
      cat("  Processing a SNP dataset\n")
      pic <- x@other$loc.metrics$AvgPIC
    } else {
      stop("Fatal Error: Ploidy must be universally 1 (fragment P/A data) or 2 (SNP data)")
    }
  }
  
  # Check monomorphs have been removed up to date
  if (x@other$loc.metrics.flags$monomorphs == FALSE){
    if (verbose >= 2){
      cat("  Warning: Dataset contains monomorphic loci which will be included in the",funname,"selections\n")
    }  
  }
  
# FUNCTION SPECIFIC ERROR CHECKING

  # To be added

# DO THE JOB
  
  if(method=="random") {
    if (verbose>=3){cat("  Subsampling at random",n,"loci from",class(x),"object","\n")}
    randsel <- sample(1:nLoc(x), n, replace = F)
    x.new <- x[, randsel]
    x.new@other$loc.metrics <- x@other$loc.metrics[randsel,]
    
    if (verbose>=3) {
      cat("     No. of loci retained =", ncol(x.new),"\n")
      }
  } else if (method=="PIC" | method=="pic"){
    x.new <- x[, order(-pic)]
    x.new <- x.new[,1:n]
    x.new@other$loc.metrics <- x.new@other$loc.metrics[1:n,]
    if (verbose>=3) {
      cat("  No. of loci retained =", ncol(x.new),"\n")
    }
  } else {
    stop("  Fatal Error: method must be 'random' or 'pic'\n");
  }

# FLAG SCRIPT END

  if (verbose > 0) {
    cat("Completed:",funname,"\n")
  }

  return(x.new)

}
