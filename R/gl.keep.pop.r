#' Remove all but specified populations from a genelight \{adegenet\} object
#'
#' Individuals are assigned to populations based on the specimen metadata data file (csv) used with gl.read.dart(). 
#'
#' The script, having deleted the specified populations, optionally identifies resultant monomorphic loci or loci
#' with all values missing and deletes them (using gl.filter.monomorphs.r). The script also optionally
#' recalculates statistics made redundant by the deletion of individuals from the dataset.
#' 
#' The script returns a genlight object with the new population assignments and the recalculated locus metadata.
#'
#' @param gl -- name of the genlight object containing SNP genotypes or a genind object containing presence/absence data [required]
#' @param pop.list -- a list of populations to be kept [required]
#' @param recalc -- Recalculate the locus metadata statistics [default TRUE]
#' @param mono.rm -- Remove monomorphic loci [default TRUE]
#' @param v -- verbosity: 0, silent; 1, brief; 2, verbose [default 1]
#' @return A genlight object with the reduced data
#' @export
#' @author Arthur Georges (glbugs@@aerg.canberra.edu.au)
#' @examples
#' \dontrun{
#'    gl <- gl.keep.pop(testset.gl, pop.list=c("EmsubRopeMata","EmvicVictJasp"))
#' }
#' @seealso \code{\link{gl.filter.monomorphs}}
#' 

gl.keep.pop <- function(gl, pop.list, recalc=FALSE, mono.rm=TRUE, v=1){
x <- gl

  if(class(x)!="genlight") {
    cat("Fatal Error: genlight object required for gl.recode.pop.r!\n"); stop()
  }

# REMOVE POPULATIONS
  if (v==2) {
    cat("Processing",class(x),"object\n")
    cat("  Deleteing populations", pop.list, "\n")
  }

# Delete listed populations, recalculate relevant locus metadata and remove monomorphic loci
  
  # Remove rows flagged for deletion
    cat("Deleting populations flagged for deletion\n")
    x2 <- x[x$pop%in%pop.list]
  # Remove monomorphic loci
    if (mono.rm) {x2 <- gl.filter.monomorphs(x2,v=0)}
  # Recalculate statistics
    if (recalc) {
      gl <- gl.recalc.metrics(x2,v=v)
    }

# REPORT A SUMMARY
  if (v==2) {
    cat("Summary of recoded dataset\n")
    cat(paste("  No. of loci:",nLoc(x2),"\n"))
    cat(paste("  No. of individuals:", nInd(x2),"\n"))
    cat(paste("  No. of populations: ", length(levels(factor(pop(x2)))),"\n"))
  }
  if (v>=1) {  
    if (!recalc) {cat("Note: Locus metrics not recalculated\n")}
    if (!mono.rm) {cat("note: Resultant monomorphic loci not deleted\n")}
  }
    
    return <- x2
}

