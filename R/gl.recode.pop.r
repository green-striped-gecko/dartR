#' Recode population assignments in a genelight object \{adegenet\}
#'
#' This script recodes population assignments and/or deletes populations from a DaRT genlight SNP file 
#' based on information provided in a csv population recode file.
#'
#' Individuals are assigned to populations based on the specimen metadata data file (csv) used with gl.read.dart(). 
#' Recoding can be used to amalgamate populations or to selectively delete or retain populations.
#'
#' The population recode file contains a list of populations in the genelight object as
#' the first column of the csv file, and the new population assignments in the second column of the csv file.
#' The keyword Delete used as a new population assignment will result in the associated specimen being dropped from the dataset.
#' 
#' The script, having deleted populations, optionally identifies resultant monomorphic loci or loci
#' with all values missing and deletes them (using gl.filter.monomorphs.r). The script also optionally
#' recalculates the locus metadata as appropriate.
#'
#' @param x -- name of the genlight object containing SNP genotypes or a genind object containing presence/absence data [required]
#' @param pop.recode -- name of the csv file containing the population reassignments [required]
#' @param recalc -- Recalculate the locus metadata statistics if any individuals are deleted in the filtering [default TRUE]
#' @param mono.rm -- Remove monomorphic loci [default TRUE]
#' @param v -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @return A genlight object with the recoded and reduced data
#' @export
#' @author Arthur Georges (bugs? Post to https://groups.google.com/d/forum/dartr)
#' @examples
#' \donttest{
#'    gl <- gl.recode.pop(gl, pop.recode="pop_recode_table_0.csv")
#' }
#' @seealso \code{\link{gl.filter.monomorphs}}
#' 
#'

gl.recode.pop <- function(x, pop.recode, recalc=TRUE, mono.rm=TRUE, v=1){

  if (v > 0) {
    cat("Starting gl.recode.pop: Recoding populations using lookup table",pop.recode,"\n")
  }
  
  if(class(x)!="genlight") {
    cat("Fatal Error: genlight object required for gl.recode.pop.r!\n"); stop("Execution terminated\n")
  }

# RECODE POPULATIONS
  if (v > 1) {
    cat("  Reassigning entities to populations as per ", pop.recode, "\n")
  }
  recode.table <- read.csv(pop.recode, stringsAsFactors=FALSE, header=FALSE);
# Error check
  if(length(unique(pop(x))) != length(unique(recode.table[,1]))) {
    cat("Fatal Error: Population names in data file are not the same as in the recode table\n"); stop("Execution terminated\n")
  }
# Apply the recode to the populations
  pop.list <- as.character(pop(x));
  ntr <- length(recode.table[,1])
  for (i in 1:nInd(x)) {
    for (j in 1:ntr) {
      if (pop.list[i]==recode.table[j,1]) {pop.list[i] <- recode.table[j,2]}
    }
  }
  pop(x) <- pop.list

# Remove rows flagged for deletion
  x2 <- x[!x$pop=="delete" & !x$pop=="Delete"]
  
  if (length(pop(x2))!=length(pop(x))) {
     if (v > 1) {
       cat("  Removing entities flagged for deletion in ", pop.recode, "\n")
     }  
     # Remove monomorphic loci
       if (mono.rm) {x2 <- gl.filter.monomorphs(x2,v=v)}
       if (recalc) {gl.recalc.metrics(x2,v=v)}
 }

  # REPORT A SUMMARY
  if (v > 2) {
    cat("Summary of recoded dataset\n")
    cat(paste("  No. of loci:",nLoc(x2),"\n"))
    cat(paste("  No. of individuals:", nInd(x2),"\n"))
    cat(paste("  No. of populations: ", length(levels(factor(pop(x2)))),"\n"))
  }
  if (v > 1) {  
    if (!recalc) {cat("  Note: Locus metrics not recalculated\n")}
    if (!mono.rm) {cat("  Note: Resultant monomorphic loci not deleted\n")}
  }
  if (v > 0) {
    cat("Completed gl.recode.pop\n\n")
  }
    return(x2)
}

