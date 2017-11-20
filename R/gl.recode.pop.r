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
#' The script, having deleted populations, identifies resultant monomorphic loci or loci
#' with all values missing and deletes them (using gl.filter.monomorphs.r). The script also recalculates the locus metadata as
#' appropriate.
#'
#' @param gl -- name of the genlight object containing SNP genotypes or a genind object containing presence/absence data [required]
#' @param pop.recode -- name of the csv file containing the population reassignments [required]
#' @param recalc -- Recalculate the locus metadata statistics if any individuals are deleted in the filtering [default TRUE]
#' @param v -- verbosity: 0, silent; 1, brief; 2, verbose [default 1]
#' @return A genlight object with the recoded and reduced data
#' @export
#' @author Arthur Georges (glbugs@@aerg.canberra.edu.au)
#' @examples
#' \dontrun{
#'    gl <- gl.recode.pop(gl, pop.recode="pop_recode_table_0.csv")
#' }
#' @seealso \code{\link{gl.filter.monomorphs}}
#' 
#'

gl.recode.pop <- function(gl, pop.recode, recalc=TRUE, v=1){
x <- gl

  if(class(x)!="genlight") {
    cat("Fatal Error: genlight object required for gl.recode.pop.r!\n"); stop()
  }

# RECODE POPULATIONS
  if (v==2) {
    cat("Processing",class(x),"object\n")
    cat("  Reassigning entities to populations as per ", pop.recode, "\n")
  }
  recode.table <- read.csv(pop.recode, stringsAsFactors=FALSE, header=FALSE);
# Error check
  if(length(unique(pop(x))) != length(unique(recode.table[,1]))) {
    cat("Fatal Error: Population names in data file are not the same as in the recode table\n"); stop()
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
     if (v==2) {
       cat("  Removing entities flagged for deletion in ", pop.recode, "\n")
     }  
     # Remove monomorphic loci
       x2 <- gl.filter.monomorphs(x2,v=0)
       if (recalc) {
         # Recalculate statistics
         x2 <- utils.recalc.avgpic(x2,v=v)
         x2 <- utils.recalc.callrate(x2,v=v)
         x2 <- utils.recalc.freqhets(x2,v=v)
         x2 <- utils.recalc.freqhomref(x2,v=v)
         x2 <- utils.recalc.freqhomsnp(x2,v=v)
       }
 }

# REPORT A SUMMARY
  if (v==2) {
    cat("Summary of recoded dataset\n")
    cat(paste("  No. of loci:",nLoc(x2),"\n"))
    cat(paste("  No. of individuals:", nInd(x2),"\n"))
    cat(paste("  No. of populations: ", length(levels(factor(pop(x2)))),"\n"))
  }

    return <- x2
}

