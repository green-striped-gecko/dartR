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
#' @param x -- name of the genlight object containing SNP genotypes or a genind object containing presence/absence data [required]
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
#' @seealso \code{\link{gl.recalc.metrics}}
#' 

gl.keep.pop <- function(x, pop.list, recalc=FALSE, mono.rm=FALSE, v=1){

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
    x <- x[x$pop%in%pop.list]
  # Remove monomorphic loci
    if (mono.rm) {x <- gl.filter.monomorphs(x,v=0)}
  # Recalculate statistics
    if (recalc) {
      x <- utils.recalc.avgpic(x,v=v)
      x <- utils.recalc.callrate(x,v=v)
      x <- utils.recalc.freqhets(x,v=v)
      x <- utils.recalc.freqhomref(x,v=v)
      x <- utils.recalc.freqhomsnp(x,v=v)
    }

# REPORT A SUMMARY
  if (v==2) {
    cat("Summary of recoded dataset\n")
    cat(paste("  No. of loci:",nLoc(x),"\n"))
    cat(paste("  No. of individuals:", nInd(x),"\n"))
    cat(paste("  No. of populations: ", length(levels(factor(pop(x)))),"\n"))
  }
    if (v>0) {  
      if (!recalc) {
        cat("Note: Locus metrics not recalculated\n")
      } else {
        cat("Note: Locus metrics recalculated\n")
      }
      if (!mono.rm) {
        cat("note: Resultant monomorphic loci not deleted\n")
      } else {
        cat("note: Resultant monomorphic loci deleted\n")
      }
    }
    
    return <- x
}

