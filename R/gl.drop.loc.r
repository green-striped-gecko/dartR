#' Remove specified loci from a genelight \{adegenet\} object
#'
#' This script deletes selected loci from the nominated dataset.
#' 
#' @param x -- name of the genlight object containing SNP genotypes or a genind object containing presence/absence data [required]
#' @param loc.list -- vector of loci names to be dropped.
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @return A genlight object with the reduced data
#' @export
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#'    gl <- gl.drop.loc(testset.gl, loc.list=c("100049687|12-A/G","100050106|50-G/A"))

# Last amended 3-Feb-19

gl.drop.loc <- function(x, loc.list, verbose=2){

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

  if (length(loc.list) == 0) {
    cat("Fatal Error: list of loci to drop required!\n"); stop("Execution terminated\n")
  }
  test <- loc.list%in%locNames(x)
  if (!all(test,na.rm=FALSE)) {
    cat("Fatal Error: some of the listed loci are not present in the dataset!\n"); stop("Execution terminated\n")
  }

# DO THE JOB

# REMOVE LOCI
  
  if (verbose >= 2) {
    cat("  Deleting selected loci", loc.list, "\n")
  }

  # Delete listed loci
  
  # Remove rows flagged for deletion
    index <- !locNames(x)%in%loc.list
    x <- x[,index]
    x@other$loc.metrics <- x@other$loc.metrics[index,]

# REPORT A SUMMARY
    
  if (verbose >= 3) {
    cat("Summary of recoded dataset\n")
    cat(paste("  No. of loci:",nLoc(x),"\n"))
    cat(paste("  No. of individuals:", nInd(x),"\n"))
    cat(paste("  No. of populations: ", length(levels(factor(pop(x)))),"\n"))
  }

# FLAG SCRIPT END

  if (verbose > 0) {
    cat("Completed:",funname,"\n")
  }
    #add to history
    nh <- length(x@other$history)
    x@other$history[[nh + 1]] <- match.call()  
  return(x)
}

