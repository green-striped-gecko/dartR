#' Create or edit a individual (=specimen) names and create an recode_ind file
#' 
#' A script to edit individual names in a genlight object, or to 
#' create a reassignment table taking the individual labels
#' from a genlight object, or to edit existing individual labels in
#' an existing recode_ind file.
#' 
#' Renaming individuals may be required when there have been errors in labelling arising
#' in the process from sample to DArT files. There may be occasions where renaming
#' individuals is required for preparation of figures. Caution needs to be exercised
#' because of the potential for breaking the "chain of evidence" between the samples themselves
#' and the analyses. Recoding individuals can also be done with a recode table (csv).
#' 
#' This script will input an existing recode table for editting and
#' optionally save it as a new table, or if the name of an input table is not
#' supplied, will generate a table using the individual labels in the 
#' parent genlight object.
#' 
#' The script, having deleted individuals, optionally identifies resultant monomorphic loci or loci
#' with all values missing and deletes them (using gl.filter.monomorphs.r). The script also optionally
#' recalculates statistics made redundant by the deletion of individuals from the dataset.
#' 
#' The script returns a genlight object with the new individual labels and the recalculated locus metadata.
#' 
#' @param gl Name of the genlight object for which individuals are to be relabelled.[required]
#' @param ind.recode Name of the file to output the new assignments [optional]
#' @param recalc -- Recalculate the locus metadata statistics [default TRUE]
#' @param mono.rm -- Remove monomorphic loci [default TRUE]
#' @return An object of class ("genlight") with the revised individual labels
#' @param v -- v=0, silent; v=1, low verbosity; v=2, high verbosity [default 1]
#' @import utils
#' @export
#' @author Arthur Georges (glbugs@aerg.canberra.edu.au)
#' @examples
#' \donttest{
#' gl <- gl.edit.recode.ind(testset.gl)
#' gl <- gl.edit.recode.ind(testset.gl, ind.recode="ind.recode.table.csv")
#' gl <- gl.edit.recode.ind(testset.gl, ind.recode="ind.recode.table.csv")
#' }
#' #Ammended Georges 9-Mar-17

gl.edit.recode.ind <- function(gl, ind.recode=NULL, recalc=TRUE, mono.rm=TRUE, v=1) {
  
# Take assignments from gl  

  cat("Extracting current individual labels from the gl object\n")
  recode.table <- cbind(indNames(gl),indNames(gl))

# Create recode table for editting, and bring up the editor
    new <- as.matrix(edit(recode.table))
    new <- new[,1:2]

# Write out the recode table, if requested
  if (is.null(ind.recode)) {
      cat("No output table specified, recode table not written to disk\n")
  } else {
    cat(paste("Writing individual recode table to: ",ind.recode,"\n"))
    write.table(new, file=ind.recode, sep=",", row.names=FALSE, col.names=FALSE)    
  }

# Apply the new assignments  
  ind.list <- as.character(indNames(gl));
  ntr <- length(new[,1])
  for (i in 1:nInd(gl)) {
    for (j in 1:ntr) {
      if (ind.list[i]==new[j,1]) {ind.list[i] <- new[j,2]}
    }
  }
  # Assigning new populations to gl
  cat("Assigning new individual (=specimen) names\n")
  indNames(gl) <- ind.list
  
  # If there are individuals to be deleted, then recalculate relevant locus metadata and remove monomorphic loci
  
  if ("delete" %in% gl$ind.names | "Delete" %in% gl$ind.names) {
    # Remove rows flagged for deletion
    cat("Deleting individuals flagged for deletion\n")
    gl <- gl[!gl$ind.names=="delete" & !gl$ind.names=="Delete"]
    # Remove monomorphic loci
    if(mono.rm) {gl <- gl.filter.monomorphs(gl,v=v)}
    # Recalculate statistics
    if (recalc) {
      gl.recalc.metrics(gl,v=v)
    }
  }

  # REPORT A SUMMARY
  if (v==2) {
    cat("Summary of recoded dataset\n")
    cat(paste("  No. of loci:",nLoc(gl),"\n"))
    cat(paste("  No. of individuals:", nInd(gl),"\n"))
    cat(paste("  No. of populations: ", length(levels(factor(pop(gl)))),"\n"))
    if (!recalc) {cat("Note: Locus metrics not recalculated\n")}
    if (!mono.rm) {cat("note: Resultant monomorphic loci not deleted\n")}
  }
  if (v>=1) {  
    if (!recalc) {cat("Note: Locus metrics not recalculated\n")}
    if (!mono.rm) {cat("note: Resultant monomorphic loci not deleted\n")}
  }
  
  return(gl)
  
}
