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
#' @param x Name of the genlight object for which individuals are to be relabelled.[required]
#' @param ind.recode Name of the file to output the new assignments [optional]
#' @param recalc -- Recalculate the locus metadata statistics [default TRUE]
#' @param mono.rm -- Remove monomorphic loci [default TRUE]
#' @return An object of class ("genlight") with the revised individual labels
#' @param verbose -- verbose=0, silent; verbose=1, low verbosity; verbose=2, high verbosity [default 1]
#' @import utils
#' @export
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' \dontrun{
#' gl <- gl.edit.recode.ind(testset.gl)
#' gl <- gl.edit.recode.ind(testset.gl, ind.recode="ind.recode.table.csv")
#' gl <- gl.edit.recode.ind(testset.gl, ind.recode="ind.recode.table.csv")
#' }

# Last amended 3-Feb-19

gl.edit.recode.ind <- function(x, ind.recode=NULL, recalc=TRUE, mono.rm=TRUE, verbose=1) {

# TIDY UP FILE SPECS

  funname <- match.call()[[1]]

# FLAG SCRIPT START

  if (verbose < 0 | verbose > 5){
    cat("  Warning: Parameter 'verbose' must be an integer between 0 [silent] and 5 [full report], set to 2\n")
    verbose <- 2
  }

  if (verbose > 0) {
    cat("Starting",funname,"\n")
  }

# STANDARD ERROR CHECKING
  
  if(class(x)!="genlight") {
    cat("  Fatal Error: genlight object required!\n"); stop("Execution terminated\n")
  }

  # Work around a bug in adegenet if genlight object is created by subsetting
      if (nLoc(x)!=nrow(x@other$loc.metrics)) { stop("The number of rows in the loc.metrics table does not match the number of loci in your genlight object!")  }

  # Set a population if none is specified (such as if the genlight object has been generated manually)
    if (is.null(pop(x)) | is.na(length(pop(x))) | length(pop(x)) <= 0) {
      if (verbose >= 2){ cat("  Population assignments not detected, individuals assigned to a single population labelled 'pop1'\n")}
      pop(x) <- array("pop1",dim = nInd(x))
      pop(x) <- as.factor(pop(x))
    }

  # Check for monomorphic loci
    tmp <- gl.filter.monomorphs(x, verbose=0)
    if ((nLoc(tmp) < nLoc(x)) & verbose >= 2) {cat("  Warning: genlight object contains monomorphic loci\n")}

# FUNCTION SPECIFIC ERROR CHECKING

# DO THE JOB
  
# Take assignments from x

  cat("Extracting current individual labels from the x object\n")
  recode.table <- cbind(indNames(x),indNames(x))

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
  ind.list <- as.character(indNames(x));
  ntr <- length(new[,1])
  for (i in 1:nInd(x)) {
    for (j in 1:ntr) {
      if (ind.list[i]==new[j,1]) {ind.list[i] <- new[j,2]}
    }
  }
  # Assigning new populations to x
  cat("Assigning new individual (=specimen) names\n")
  indNames(x) <- ind.list
  
  # If there are individuals to be deleted, then recalculate relevant locus metadata and remove monomorphic loci
  
  if ("delete" %in% x$ind.names | "Delete" %in% x$ind.names) {
    # Remove rows flagged for deletion
    cat("Deleting individuals flagged for deletion\n")
    x <- x[!x$ind.names=="delete" & !x$ind.names=="Delete"]
    # Remove monomorphic loci
    if(mono.rm) {x <- gl.filter.monomorphs(x,verbose=verbose)}
    # Recalculate statistics
    if (recalc) {
      gl.recalc.metrics(x,verbose=verbose)
    }
  }

  # REPORT A SUMMARY
  if (verbose==2) {
    cat("Summary of recoded dataset\n")
    cat(paste("  No. of loci:",nLoc(x),"\n"))
    cat(paste("  No. of individuals:", nInd(x),"\n"))
    cat(paste("  No. of populations: ", length(levels(factor(pop(x)))),"\n"))
    if (!recalc) {cat("Note: Locus metrics not recalculated\n")}
    if (!mono.rm) {cat("note: Resultant monomorphic loci not deleted\n")}
  }
  if (verbose>=1) {
    if (!recalc) {cat("Note: Locus metrics not recalculated\n")}
    if (!mono.rm) {cat("note: Resultant monomorphic loci not deleted\n")}
  }
  
# FLAG SCRIPT END

  if (verbose > 0) {
    cat("Completed:",funname,"\n")
  }
  #add to history
  nh <- length(x@other$history)
  x@other$history[[nh + 1]] <- match.call()  
  return(x)
  
}
