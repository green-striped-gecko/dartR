#' Create or edit a population re-assignment table
#' 
#' A script to edit population assignments in a genlight object, or to 
#' create a reassignment table taking the population assignments
#' from a genlight object, or to edit existing population assignments in
#' a pop.recode.table.
#' 
#' Genlight objects assign specimens to populations based on information in the
#' ind.metadata file provided when the genlight object is first generated.
#' Often one wishes to subset the data by deleting populations or to amalgamate
#' populations. This can be done with a pop.recode table with two columns. The
#' first column is the population assignment in the genlight object, the second
#' column provides the new assignment.
#' 
#' This script will input an existing reassignment table for editting and
#' optionally save it as a new table, or if the name of an input table is not
#' supplied, will generate a table using the population assignments in the 
#' parent genlight object.
#' 
#' The script, having deleted populations, optionally identifies resultant monomorphic loci or loci
#' with all values missing and deletes them (using gl.filter.monomorphs.r). The script also optionally
#' recalculates statistics made redundant by the deletion of individuals from the dataset.
#' 
#' The script returns a genlight object with the new population assignments and the recalculated locus metadata.
#' 
#' @param gl Name of the genlight object for which populations are to be reassigned.[required]
#' @param pop.recode Name of the file to output the new assignments [optional]
#' @param recalc -- Recalculate the locus metadata statistics if any individuals are deleted [default TRUE]
#' @param mono.rm -- Remove monomorphic loci [default TRUE]
#' @param v -- verbosity: 0, silent; 1, brief; 2, verbose [default 1]
#' @return An object of class ("genlight") with the revised population assignments
#' @import utils
#' @export
#' @author Arthur Georges (glbugs@aerg.canberra.edu.au)
#' @examples
#' \dontrun{
#' gl <- gl.edit.recode.pop(gl)
#' gl <- gl.edit.recode.pop(gl, pop.recode="pop.recode.table.csv")
#' gl <- gl.edit.recode.pop(gl, pop.recode="pop.recode.table.csv", 
#' pop.recode.new="recode.table.minus.hybrids.csv")
#' }
#
# Ammended Georges 29-Oct-16

gl.edit.recode.pop <- function(gl, pop.recode=NULL, recalc=TRUE, mono.rm=TRUE, v=1) {
  
# Take assignments from gl  

  cat("Extracting current pop assignments from the gl object\n")
  recode.table <- cbind(levels(pop(gl)),levels(pop(gl)))

# Create recode table for editting, and bring up the editor
    new <- as.matrix(edit(recode.table))
    new <- new[,1:2]

# Write out the recode table, if requested
  if (is.null(pop.recode)) {
      cat("No output table specified, recode table not written to disk\n")
  } else {
    cat(paste("Writing population recode table to: ",pop.recode,"\n"))
    write.table(new, file=pop.recode, sep=",", row.names=FALSE, col.names=FALSE)    
  }

# Apply the new assignments  
  pop.list <- as.character(pop(gl));
  ntr <- length(new[,1])
  for (i in 1:nInd(gl)) {
    for (j in 1:ntr) {
      if (pop.list[i]==new[j,1]) {pop.list[i] <- new[j,2]}
    }
  }
  # Assigning new populations to gl
  cat("Assigning new population names\n")
  pop(gl) <- pop.list
  
  # If there are populations to be deleted, then recalculate relevant locus metadata and remove monomorphic loci
  
  if ("delete" %in% gl$pop | "Delete" %in% gl$pop) {
    # Remove rows flagged for deletion
      cat("Deleting populations flagged for deletion\n")
      gl <- gl[!gl$pop=="delete" & !gl$pop=="Delete"]
    # Remove monomorphic loci
      if(mono.rm) {gl <- gl.filter.monomorphs(gl,v=v)}
    # Recalculate statistics
      if (recalc) {
        gl <- utils.recalc.avgpic(gl,v=v)
        gl <- utils.recalc.callrate(gl,v=v)
        gl <- utils.recalc.freqhets(gl,v=v)
        gl <- utils.recalc.freqhomref(gl,v=v)
        gl <- utils.recalc.freqhomsnp(gl,v=v)
    }  
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
  return(gl)
  
}
