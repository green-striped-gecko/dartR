#' Remove specified individuals from a genelight \{adegenet\} object
#'
#' The script, having deleted individuals, optionally identifies resultant monomorphic loci or loci
#' with all values missing and deletes them (using gl.filter.monomorphs.r). The script also optionally
#' recalculates statistics made redundant by the deletion of individuals from the dataset.
#' 
#' The script returns a genlight object with the individuals deleted and, optionally, the recalculated locus metadata.
#'
#' @param x -- name of the genlight object containing SNP genotypes [required]
#' @param ind.list -- a list of individuals to be removed [required]
#' @param recalc -- Recalculate the locus metadata statistics [default FALSE]
#' @param mono.rm -- Remove monomorphic loci [default FALSE]
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @return A genlight object with the reduced data
#' @export
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#'    gl <- gl.drop.ind(testset.gl, ind.list=c("AA019073","AA004859"))
#' @seealso \code{\link{gl.filter.monomorphs}}
#' @seealso \code{\link{gl.recalc.metrics}}

gl.drop.ind <- function(x, ind.list, recalc=FALSE, mono.rm=FALSE, verbose=2){

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
      cat("Warning: Listed individual",case,"not present in the dataset -- ignored\n")
      ind.list <- ind.list[!(ind.list==case)]
    }
  }
  if (length(ind.list) == 0) {
    stop("Fatal Error: no individuals to drop!\n")
  }

# DO THE JOB

# REMOVE INDIVIDUALS
  
  if (verbose >= 2) {
    cat("  Deleting individuals", ind.list, "\n")
  }

# Delete listed individuals, recalculate relevant locus metadata and remove monomorphic loci
  
  # Remove rows flagged for deletion
    x <- x[!x$ind.names%in%ind.list]
    
  # Reset the flags
    x <- utils.reset.flags(x, verbose=verbose)
    
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
    cat("Summary of recoded dataset\n")
    cat(paste("  No. of loci:",nLoc(x),"\n"))
    cat(paste("  No. of individuals:", nInd(x),"\n"))
    cat(paste("  No. of populations: ", length(levels(factor(pop(x)))),"\n"))
  }
  if (verbose >= 2) {
    if (!recalc) {
      cat("Note: Locus metrics not recalculated\n")
    } else {
      cat("Note: Locus metrics recalculated\n")
    }
    if (!mono.rm) {
      cat("Note: Resultant monomorphic loci not deleted\n")
    } else{
      cat("Note: Resultant monomorphic loci deleted\n")
    }
  }
    
# FLAG SCRIPT END
    
  #add to history
    nh <- length(x@other$history)
    x@other$history[[nh + 1]] <- match.call()  

  if (verbose > 0) {
    cat("Completed:",funname,"\n")
  }
    
  return(x)
  
}

# # Test script
# 
# gl <- testset.gl
# tmp <- gl.drop.ind(testset.gl, ind.list=c("AA019073","AA004859"))
# nInd(gl)
# nInd(tmp)
# gl@other$loc.metrics.flags
# tmp@other$loc.metrics.flags
# 
# gl <- testset.gl
# tmp <- gl.drop.ind(testset.gl, ind.list=c("AA019073","AA004859"),mono.rm=TRUE)
# nInd(gl)
# nInd(tmp)
# gl@other$loc.metrics.flags
# tmp@other$loc.metrics.flags
# 
# gl <- testset.gl
# tmp <- gl.drop.ind(testset.gl, ind.list=c("AA019073","AA004859"),recalc=TRUE)
# nInd(gl)
# nInd(tmp)
# gl@other$loc.metrics.flags
# tmp@other$loc.metrics.flags
# 
# gs <- testset.gs
# tmp <- gl.drop.pop(gs, pop.list=c("EmsubRopeMata","EmvicVictJasp"))
# nPop(gs)
# nPop(tmp)
# gs@other$loc.metrics.flags
# tmp@other$loc.metrics.flags
# 
# gl <- testset.gl
# tmp <- gl.drop.pop(gs, pop.list=c("EmsubRopeMata","EmvicVictJasp"),mono.rm = TRUE)
# nPop(gs)
# nPop(tmp)
# gs@other$loc.metrics.flags
# tmp@other$loc.metrics.flags
# 
# gl <- testset.gl
# tmp <- gl.drop.pop(gs, pop.list=c("EmsubRopeMata","EmvicVictJasp"),mono.rm = FALSE, recalc = TRUE)
# nPop(gs)
# nPop(tmp)
# gs@other$loc.metrics.flags
# tmp@other$loc.metrics.flags

