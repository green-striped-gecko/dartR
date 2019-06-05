#' Define a new population in a genlight \{adegenet\} object on the basis of specified individuals 
#'
#' The script reassigns existing individuals to a new population and removes their existing population assignment
#' 
#' The script returns a genlight object with the new population assignment.
#'
#' @param x -- name of the genlight object containing SNP genotypes [required]
#' @param ind.list -- a list of individuals to be assigned to the new population [required]
#' @param new -- name of the new population
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @return A genlight object with the redefined population structure
#' @export
#' @author Arthur Georges (bugs? Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#'    gl <- gl.define.pop(testset.gl, ind.list=c("AA019073","AA004859"), new="newguys")

# Last amended 3-Feb-19

gl.define.pop <- function(x, ind.list, new, verbose=2){

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

  # Check for monomorphic loci
    tmp <- gl.filter.monomorphs(x, verbose=0)
    if ((nLoc(tmp) < nLoc(x)) & verbose >= 2) {cat("  Warning: genlight object contains monomorphic loci\n")}

# FUNCTION SPECIFIC ERROR CHECKING

  for (case in ind.list){
    if (!(case%in%indNames(x))){
      cat("Warning: Listed individual",case,"not present in the dataset -- ignored\n")
      ind.list <- ind.list[!(ind.list==case)]
    }
  }
  if (length(ind.list) == 0) {
    cat("Fatal Error: no individuals listed to assign to population",new,"\n"); stop("Execution terminated\n")
  }

# DO THE JOB

# ASSIGN INDIVIDUALS
  
  if (verbose >= 2) {
    cat("Processing",class(x),"object\n")
    cat("  Assigned listed individuals", ind.list,"to new population",new, "\n")
  }

  tmp <- as.character(pop(x))
  for (case in ind.list){
    tmp[indNames(x) == case] <- new
  }
  pop(x) <- as.factor(tmp)
  
# REPORT A SUMMARY
    
  if (verbose >= 3) {
    cat("Summary of recoded dataset\n")
    cat(paste("  No. of loci:",nLoc(x),"\n"))
    cat(paste("  No. of individuals:", nInd(x),"\n"))
    cat(paste("  No. of populations: ", length(levels(factor(pop(x)))),"\n"))
  }

# FLAG SCRIPT END

  if (verbose > 0) {
    cat("Completed:", funname, "\n")
  }
  #add to history
  nh <- length(x@other$history)
  x@other$history[[nh + 1]] <- match.call()
  
  return(x)
}

