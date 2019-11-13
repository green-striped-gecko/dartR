#' Remove specified populations from a genelight \{adegenet\} object
#'
#' Individuals are assigned to populations based on the specimen metadata data file (csv) used with gl.read.dart() or gs.read.dart(). 
#'
#' The script, having deleted populations, optionally identifies resultant monomorphic loci or loci
#' with all values missing and deletes them (using gl.filter.monomorphs.r). The script also optionally
#' recalculates statistics made redundant by the deletion of individuals from the dataset.
#' 
#' The script returns a genlight object with the new population assignments and the recalculated locus metadata.
#'
#' @param x -- name of the genlight object containing SNP genotypes or Tag P/A data (SilicoDArT) [required]
#' @param pop.list -- a list of populations to be removed [required]
#' @param as.pop -- assign another metric to represent population [default NULL]
#' @param recalc -- Recalculate the locus metadata statistics [default FALSE]
#' @param mono.rm -- Remove monomorphic loci [default FALSE]
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @return A genlight object with the reduced data
#' @export
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#'    gl <- gl.drop.pop(testset.gl, pop.list=c("EmsubRopeMata","EmvicVictJasp"))
#'    gl <- gl.drop.pop(testset.gl, pop.list=c("Male","Unknown"),as.pop="sex")
#' @seealso \code{\link{gl.filter.monomorphs}}
#' @seealso \code{\link{gl.recalc.metrics}}

gl.drop.pop <- function(x, pop.list, as.pop=NULL, recalc=FALSE, mono.rm=FALSE, verbose=NULL){

# TIDY UP FILE SPECS

  funname <- match.call()[[1]]
  build <- "Jacob"

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
  
  if (all(x@ploidy == 1)){
    if (verbose >= 2){cat("  Processing  Presence/Absence (SilicoDArT) data\n")}
    data.type <- "SilicoDArT"
  } else if (all(x@ploidy == 2)){
    if (verbose >= 2){cat("  Processing a SNP dataset\n")}
    data.type <- "SNP"
  } else {
    stop("Fatal Error: Ploidy must be universally 1 (fragment P/A data) or 2 (SNP data)")
  }
  
# FUNCTION SPECIFIC ERROR CHECKING
    
    # Set a population if none is specified (such as if the genlight object has been generated manually)
    if (is.null(pop(x)) | is.na(length(pop(x))) | length(pop(x)) <= 0) {
      if (verbose >= 2){ 
        cat("  Population assignments not detected, individuals assigned to a single population labelled 'pop1'\n")
      }
      pop(x) <- array("pop1",dim = nInd(x))
      pop(x) <- as.factor(pop(x))
    }
    
  # Assign the new population list if as.pop is specified
    pop.hold <- pop(x)
    if (!is.null(as.pop)){
      pop.hold <- pop(x)
      pop(x) <- as.matrix(x@other$ind.metrics[as.pop])
      if (verbose >= 3) {cat("  Temporarily setting population assignments to",as.pop,"as specified by the as.pop parameter\n")}
    }
    
  if (verbose >= 2) {
    cat("  Checking for presence of nominated populations\n")
  }
  for (case in pop.list){
    if (!(case%in%popNames(x))){
      cat("  Warning: Listed population",case,"not present in the dataset -- ignored\n")
      pop.list <- pop.list[!(pop.list==case)]
    }
  }
  if (length(pop.list) == 0) {
    stop("Fatal Error: no populations listed to drop!\n")
  }
    
# DO THE JOB

# REMOVE POPULATIONS
  
  if (verbose >= 2) {
    cat("  Deleting populations", pop.list, "\n")
  }

# Delete listed populations, recalculate relevant locus metadata and remove monomorphic loci
  
  # Remove rows flagged for deletion
    x2 <- x[!x$pop%in%pop.list]
    pop.hold <- pop.hold[!x$pop%in%pop.list]
    x <- x2
    
  # Reset the flags
    x <- utils.reset.flags(x, verbose=0)
    
  # Remove monomorphic loci
    if (mono.rm) {
      x <- gl.filter.monomorphs(x,verbose=verbose)
    }
  # Recalculate statistics
    if (recalc) {
      x <- gl.recalc.metrics(x,verbose=verbose)
    }

  # REPORT A SUMMARY
    
  if (verbose >= 3) {
    if (!is.null(as.pop)) {
      cat("  Summary of recoded dataset\n")
      cat(paste("    No. of loci:",nLoc(x),"\n"))
      cat(paste("    No. of individuals:", nInd(x),"\n"))
      cat(paste("    No. of levels of",as.pop,"remaining: ", length(levels(factor(pop(x)))),"\n"))
      cat(paste("    No. of populations: ", length(levels(factor(pop.hold))),"\n"))
    } else {
      cat("  Summary of recoded dataset\n")
      cat(paste("    No. of loci:",nLoc(x),"\n"))
      cat(paste("    No. of individuals:", nInd(x),"\n"))
      cat(paste("    No. of populations: ", length(levels(factor(pop(x)))),"\n"))
    }  
  }
  if (verbose >= 2) {
    if (!recalc) {
      cat("  Note: Locus metrics not recalculated\n")
    } else {
      cat("  Note: Locus metrics recalculated\n")
    }
    if (!mono.rm) {
      cat("  Note: Resultant monomorphic loci not deleted\n")
    } else{
      cat("  Note: Resultant monomorphic loci deleted\n")
    }
  }
  
  # Reassign the initial population list if as.pop is specified
    
    if (!is.null(as.pop)){
      pop(x) <- pop.hold
      if (verbose >= 3) {cat("  Resetting population assignments to initial state\n")}
    }

# ADD TO HISTORY
  nh <- length(x@other$history)
  x@other$history[[nh + 1]] <- match.call() 
    
# FLAG SCRIPT END

    if (verbose > 0) {
      cat("Completed:", funname, "\n")
    }
    
    return(x)
}