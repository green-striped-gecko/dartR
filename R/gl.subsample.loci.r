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

# TIDY UP FILE SPECS
  
  build <- "Jacob"
  funname <- match.call()[[1]]
  hold <- x
  # Note does not draw upon or modifies the loc.metrics.flags.
  
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

# test <- gl.subsample.loci(gl, 6, method="pic")
# as.matrix(test)[1:20,]
# as.matrix(x)[1:20,1:6]
# as.matrix(test)[1:20,1:6]
# 
# test <- gl.subsample.loci(gs, 6, method="pic")
# as.matrix(test)[1:20,]
# as.matrix(x)[1:20,1:6]
# as.matrix(test)[1:20,1:6]

