#' Remove specified populations from a genelight \{adegenet\} object
#'
#' Individuals are assigned to populations based on the specimen metadata data file (csv) used with gl.read.dart(). 
#'
#' The script, having deleted populations, optionally identifies resultant monomorphic loci or loci
#' with all values missing and deletes them (using gl.filter.monomorphs.r). The script also optionally
#' recalculates statistics made redundant by the deletion of individuals from the dataset.
#' 
#' The script returns a genlight object with the new population assignments and the recalculated locus metadata.
#'
#' @param x -- name of the genlight object containing SNP genotypes or a genind object containing presence/absence data [required]
#' @param pop.list -- a list of populations to be removed [required]
#' @param recalc -- Recalculate the locus metadata statistics [default FALSE]
#' @param mono.rm -- Remove monomorphic loci [default TRUE]
#' @param v -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @return A genlight object with the reduced data
#' @export
#' @author Arthur Georges (glbugs@@aerg.canberra.edu.au)
#' @examples
#'    gl <- gl.drop.pop(testset.gl, pop.list=c("EmsubRopeMata","EmvicVictJasp"))
#' @seealso \code{\link{gl.filter.monomorphs}}
#' @seealso \code{\link{gl.recalc.metrics}}
#' 

gl.drop.pop <- function(x, pop.list, recalc=FALSE, mono.rm=TRUE, v=2){

  if(class(x)!="genlight") {
    cat("Fatal Error: genlight object required!\n"); stop("Execution terminated\n")
  }
  if (length(pop.list) == 0) {
    cat("Fatal Error: list of populations to drop required!\n"); stop("Execution terminated\n")
  }
  test <- pop.list%in%levels(pop(x))
  if (!all(test,na.rm=FALSE)) {
    cat("Fatal Error: some of the listed populations are not present in the dataset!\n"); stop("Execution terminated\n")
  }
  if (v >= 1) {
    cat("Starting gl.drop.pop: Deleting selected populations\n")
  }

# REMOVE POPULATIONS
  if (v >= 2) {
    cat("  Deleting populations", pop.list, "\n")
  }

# Delete listed populations, recalculate relevant locus metadata and remove monomorphic loci
  
  # Remove rows flagged for deletion
    x <- x[!x$pop%in%pop.list]
  # Remove monomorphic loci
    if (mono.rm) {x <- gl.filter.monomorphs(x,v=0)}
  # Recalculate statistics
    if (recalc) {gl.recalc.metrics(x,v=v)}

# REPORT A SUMMARY
  if (v >= 3) {
    cat("Summary of recoded dataset\n")
    cat(paste("  No. of loci:",nLoc(x),"\n"))
    cat(paste("  No. of individuals:", nInd(x),"\n"))
    cat(paste("  No. of populations: ", length(levels(factor(pop(x)))),"\n"))
  }
    if (v >= 2) {  
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
    if (v > 0) {
      cat("Completed gl.drop.pop\n\n")
    }
    
    return <- x
}

