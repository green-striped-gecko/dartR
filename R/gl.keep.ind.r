#' Remove all but the specified individuals from a genelight \{adegenet\} object
#'
#' The script, having deleted individuals, optionally identifies resultant monomorphic loci or loci
#' with all values missing and deletes them (using gl.filter.monomorphs.r). The script also optionally
#' recalculates statistics made redundant by the deletion of individuals from the dataset.
#' 
#' The script returns a genlight object with the individuals deleted and, optionally, the recalculated locus metadata.
#'
#' @param x -- name of the genlight object containing SNP genotypes or a genind object containing presence/absence data [required]
#' @param ind.list -- a list of individuals to be removed [required]
#' @param recalc -- Recalculate the locus metadata statistics [default FALSE]
#' @param mono.rm -- Remove monomorphic loci [default FALSE]
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @return A genlight object with the reduced data
#' @export
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#'    gl <- gl.keep.ind(testset.gl, ind.list=c("AA019073","AA004859"))
#' @seealso \code{\link{gl.filter.monomorphs}} for when mono.rm=TRUE, \code{\link{gl.recalc.metrics}} for when recalc=TRUE
#' @seealso \code{\link{gl.drop.ind}} to drop rather than keep specified individuals

gl.keep.ind <- function(x, ind.list, recalc=FALSE, mono.rm=FALSE, verbose=2){

# TIDY UP FILE SPECS
  
  funname <- match.call()[[1]]
  build <- "Jacob"
  
# FLAG SCRIPT START
  
  if (verbose < 0 | verbose > 5){
    cat("  Warning: Parameter 'verbose' must be an integer between 0 [silent] and 5 [full report], set to 2\n")
    verbose <- 2
  }
  
  if (verbose >= 1){
    cat("Starting",funname,"[ Build =",build,"]\n")
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

  for (case in ind.list){
    if (!(case%in%indNames(x))){
      cat("  Warning: Listed individual",case,"not present in the dataset -- ignored\n")
      ind.list <- ind.list[!(ind.list==case)]
    }
  }
  if (length(ind.list) == 0) {
    cat("  Fatal Error: no individuals listed to keep!\n"); stop("Execution terminated\n")
  }

# DO THE JOB

  if (verbose >= 2) {
    cat("  Processing",class(x),"object\n")
    cat("    Deleteing all but the listed individuals", ind.list, "\n")
  }

  # Delete all but the listed individuals, recalculate relevant locus metadata and remove monomorphic loci
  
  # Remove rows flagged for deletion
    x <- x[x$ind.names%in%ind.list]
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
    cat("  Summary of recoded dataset\n")
    cat(paste("    No. of loci:",nLoc(x),"\n"))
    cat(paste("    No. of individuals:", nInd(x),"\n"))
    cat(paste("    No. of populations: ", length(levels(factor(pop(x)))),"\n"))
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
    
# ADD TO HISTORY
    nh <- length(x@other$history)
    x@other$history[[nh + 1]] <- match.call() 
    
# FLAG SCRIPT END

    if (verbose > 0) {
      cat("Completed: gl.keep.ind\n")
    }
  
    return(x)
}

