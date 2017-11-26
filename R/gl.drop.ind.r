#' Remove specified individuals from a genelight \{adegenet\} object
#'
#' The script, having deleted individuals, optionally identifies resultant monomorphic loci or loci
#' with all values missing and deletes them (using gl.filter.monomorphs.r). The script also optionally
#' recalculates statistics made redundant by the deletion of individuals from the dataset.
#' 
#' The script returns a genlight object with the individuals deleted and, optionally, the recalculated locus metadata.
#'
#' @param gl -- name of the genlight object containing SNP genotypes or a genind object containing presence/absence data [required]
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
#' 

gl.drop.ind <- function(gl, ind.list, recalc=TRUE, mono.rm=TRUE, v=1){
x <- gl

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
    x2 <- x[!x$ind.names%in%ind.list]
  # Remove monomorphic loci
    if (mono.rm) {x2 <- gl.filter.monomorphs(x2,v=0)}
  # Recalculate statistics
        if (recalc) {
      x2 <- utils.recalc.avgpic(x2,v=v)
      x2 <- utils.recalc.callrate(x2,v=v)
      x2 <- utils.recalc.freqhets(x2,v=v)
      x2 <- utils.recalc.freqhomref(x2,v=v)
      x2 <- utils.recalc.freqhomsnp(x2,v=v)
    }

# REPORT A SUMMARY
  if (v==2) {
    cat("Summary of recoded dataset\n")
    cat(paste("  No. of loci:",nLoc(x2),"\n"))
    cat(paste("  No. of individuals:", nInd(x2),"\n"))
    cat(paste("  No. of populations: ", length(levels(factor(pop(x2)))),"\n"))
    if (!recalc) {cat("Note: Locus metrics not recalculated\n")}
    if (!mono.rm) {cat("note: Resultant monomorphic loci not deleted\n")}
  }
  if (v==1) {  
    if (!recalc) {cat("Note: Locus metrics not recalculated\n")}
    if (!mono.rm) {cat("note: Resultant monomorphic loci not deleted\n")}
  }
    return <- x2
}

