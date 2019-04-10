#' Assign a individual metric as pop in a genlight \{adegenet\} object
#'
#' Individuals are assigned to populations based on the individual/sample/specimen metrics file (csv)
#'  used with gl.read.dart(). 
#'
#' One might want to define the population structure in accordance with another classification, such as
#' using an individual metric (e.g. sex, male or female). This script discards the current population 
#' assignments and replaces them with new population assignments defined by a specified individual metric.
#' 
#' The script returns a genlight object with the new population assignments Note that the original population
#' assigments are lost.
#'
#' @param x -- name of the genlight object containing SNP genotypes [required]
#' @param as.pop -- specify the name of the individual metric to set as the pop variable. [required]
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @return A genlight object with the reassigned populations
#' @export
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#'    gl <- gl.reassign.pop(testset.gl, as.pop='sex')

# Last amended 2-Apr-19

gl.reassign.pop <- function (x, as.pop, verbose = 2) {
  
  # TIDY UP FILE SPECS
  
  funname <- match.call()[[1]]
  
  # FLAG SCRIPT START
  
  if (verbose < 0 | verbose > 5) {
    cat("  Warning: Parameter 'verbose' must be an integer between 0 [silent] and 5 [full report], set to 2\n")
    verbose <- 2
  }
  
  if (verbose > 0) {
    cat("Starting", funname, "\n")
  }
  
  # STANDARD ERROR CHECKING
  
  if (class(x) != "genlight") {
    cat("  Fatal Error: genlight object required!\n")
    stop("Execution terminated\n")
  }
  
  x@other$loc.metrics <- x@other$loc.metrics[1:nLoc(x), ]
  if (is.null(pop(x)) | is.na(length(pop(x))) | length(pop(x)) <= 
      0) {
    if (verbose >= 2) {
      cat("  Population assignments not detected, individuals assigned to a single population labelled 'pop1'\n")
    }
    pop(x) <- array("pop1", dim = nLoc(x))
    pop(x) <- as.factor(pop(x))
  }
  tmp <- gl.filter.monomorphs(x, verbose = 0)
  if ((nLoc(tmp) < nLoc(x)) & verbose >= 2) {
    cat("  Warning: genlight object contains monomorphic loci\n")
  }
  
  # SCRIPT SPECIFIC ERROR CHECKING
  
  if (!(as.pop %in% names(x@other$ind.metrics))) {
    cat("  Fatal Error: Specified individual metric", as.pop, "not present in the dataset\n")
    stop()
  }
  
  # DO THE JOB
  
  pop(x) <- as.matrix(x@other$ind.metrics[as.pop])
  if (verbose >= 3) {
    cat("  Setting population assignments to", as.pop,"\n")
  }
  
  if (verbose >= 3) {
      cat("  Summary of recoded dataset\n")
      cat(paste("    No. of loci:", nLoc(x), "\n"))
      cat(paste("    No. of individuals:", nInd(x), "\n"))
      cat(paste("    No. of populations: ", length(levels(factor(pop(x)))), "\n"))
      cat(paste("    No. of populations: ", length(levels(factor(pop.hold))), "\n"))
  }

  if (verbose > 0) {
    cat("Completed:",funname,"\n")
  }
  
  return <- x
}