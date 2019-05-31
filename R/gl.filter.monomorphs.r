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

  # Work around a bug in adegenet if genlight object is created by subsetting
    x@other$loc.metrics <- x@other$loc.metrics[1:nLoc(x),]

  # Set a population if none is specified (such as if the genlight object has been generated manually)
    if (is.null(pop(x)) | is.na(length(pop(x))) | length(pop(x)) <= 0) {
      if (verbose >= 2){ cat("  Population assignments not detected, individuals assigned to a single population labelled 'pop1'\n")}
      pop(x) <- array("pop1",dim = nLoc(x))
      pop(x) <- as.factor(pop(x))
    }

# DO THE JOB
  
# Create a vector to hold test results
  a <- vector(mode="logical", length=nLoc(x))
  for (i in 1:nLoc(x)) {a[i] <- NA}
# Set up the progress counter
  if (verbose > 1 && pb == TRUE) {
    progress <- txtProgressBar(min=0, max=1, style=3, initial=0, label="Working ....")
    getTxtProgressBar(progress)
  }
# Identify polymorphic, monomorphic and 'all na' loci
  # Set a <- TRUE if monomorphic, or if all NAs
  xmat <-as.matrix(x)
  for (i in (1:nLoc(x))) {
    a[i] <- all(xmat[,i]==0,na.rm=TRUE) || all(xmat[,i]==2,na.rm=TRUE)
    if (all(is.na(xmat[,i]))) {a[i] <- NA}
    if (verbose > 1 && pb == TRUE) {setTxtProgressBar(progress, i/nLoc(x))}
  }
# Count the number of monomorphic loci (TRUE), polymorphic loci (FALSE) and loci with no scores (all.na)
  counts <- plyr::count(a)
  if (verbose > 2) {
    if (pb) {cat("\n")}
    cat("  Polymorphic loci:", counts[1,2], "\n")
    cat("  Monomorphic loci:", counts[2,2], "\n")
    if (is.na(counts[3,2])) {counts[3,2] <- 0}
    cat("  Loci with no scores (all NA):" , counts[3,2] ,"\n")
  }

  #Treat all na loci as monomorphic
  a[is.na(a)] <- TRUE
# Write the polymorphic loci to a new genlight object
  if (verbose > 1) {cat("  Deleting monomorphic loci and loci with all NA scores\n")}

  x <- x[,(a==FALSE)]
  x@other$loc.metrics <- x@other$loc.metrics[(a==FALSE),]
  
# FLAG SCRIPT END

  if (verbose > 0) {
    cat("Completed:",funname,"\n")
  }
  
return (x)
}
