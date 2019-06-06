#' Remove monomorphic loci, including those with all NAs
#'
#' This script deletes monomorphic loci from a genlight \{adegenet\} object
#'
#' A DArT dataset will not have monomorphic loci, but they can arise when populations are deleted by assignment or by using
#' the delete option in gl.pop.recode(). Retaining monomorphic loci unnecessarily increases the size of the dataset.
#' @param x -- name of the input genlight object [required]
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @param pb -- display progress bar [FALSE]
#' @return A genlight object with monomorphic loci removed
#' @import utils
#' @importFrom plyr count
#' @export
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' gl <- gl.filter.monomorphs(testset.gl, verbose=3)

# Last amended 3-Feb-19

gl.filter.monomorphs <- function (x, verbose=2, pb=FALSE) {

# TIDY UP FILE SPECS

  funname <- match.call()[[1]]

# FLAG SCRIPT START

  if (verbose < 0 | verbose > 5){
    cat("  Warning: Parameter 'verbose' must be an integer between 0 [silent] and 5 [full report], set to 2\n")
    verbose <- 2
  }

  if (verbose > 0) {
    cat("Starting",funname,"\n")
  }

# STANDARD ERROR CHECKING
  
  if(class(x)!="genlight") {
    cat("  Fatal Error: genlight object required!\n"); stop("Execution terminated\n")
  }

  if (nLoc(x)!=nrow(x@other$loc.metrics)) {
    stop("The number of rows in the loc.metrics table does not match the number of loci in your genlight object!")
  }

  # Set a population if none is specified (such as if the genlight object has been generated manually)
    if (is.null(pop(x)) | is.na(length(pop(x))) | length(pop(x)) <= 0) {
      if (verbose >= 2){ cat("  Population assignments not detected, individuals assigned to a single population labelled 'pop1'\n")}
      pop(x) <- factor(rep("pop1", nInd(x)))
    }

# DO THE JOB
  mml <- !(colMeans(as.matrix(x), na.rm=T) %%2==0)

#####  
  
    if (verbose > 2) {
    if (pb) {cat("\n")}
    cat("  Polymorphic loci:", sum(mml, na.rm=T), "\n")
    cat("  Monomorphic loci:", sum(!mml, na.rm=T), "\n")
    if (is.na(sum(mml))) {
    cat("  Loci with no scores (all NA):" , sum(is.na(mml)) ,"\n")
  }
}
  #Treat all na loci as monomorphic
  mml[is.na(mml)]<- FALSE
# Write the polymorphic loci to a new genlight object
  if (verbose > 1) {cat("  Deleting monomorphic loci and loci with all NA scores\n")}
  
  x <- x[,mml]
  x@other$loc.metrics <- x@other$loc.metrics[mml,]
  
# FLAG SCRIPT END

  if (verbose > 0)  cat("Completed:",funname,"\n")
  #add to history
  nh <- length(x@other$history)
  x@other$history[[nh + 1]] <- match.call()
return (x)
}
