#' Remove specified individuals from a genelight \{adegenet\} object
#'
#' The script, having deleted individuals, optionally identifies resultant monomorphic loci or loci
#' with all values missing and deletes them (using gl.filter.monomorphs.r). The script also optionally
#' recalculates statistics made redundant by the deletion of individuals from the dataset.
#' 
#' The script returns a genlight object with the individuals deleted and, optionally, the recalculated locus metadata.
#'
#' @param x -- name of the genlight object containing SNP genotypes or a genind object containing presence/absence data [required]
#' @param ind.list -- a list of individuals to be removed [required]
#' @param recalc -- Recalculate the locus metadata statistics [default TRUE]
#' @param mono.rm -- Remove monomorphic loci [default TRUE]
#' @param v -- verbosity: 0, silent; 1, brief; 2, verbose [default 1]
#' @return A genlight object with the reduced data
#' @export
#' @author Arthur Georges (glbugs@@aerg.canberra.edu.au)
#' @examples
#' \dontrun{
#'    gl <- gl.drop.ind(testset.gl, ind.list=c("AA019073","AA004859"))
#' }
#' @seealso \code{\link{gl.filter.monomorphs}}
#' @seealso \code{\link{gl.recalc.metrics}}
#' 

gl.drop.ind <- function(x, ind.list, recalc=FALSE, mono.rm=FALSE, v=1){

  if(class(x)!="genlight") {
    cat("Fatal Error: genlight object required for gl.drop.ind\n"); stop()
  }

# REMOVE POPULATIONS
  if (v==2) {
    cat("Processing",class(x),"object\n")
    cat("  Deleteing individuals", pop.list, "\n")
  }

# Delete listed individualss, recalculate relevant locus metadata and remove monomorphic loci
  
  # Remove rows flagged for deletion
    cat("Deleting individuals flagged for deletion\n")
    x <- x[!x$ind.names%in%ind.list]
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

