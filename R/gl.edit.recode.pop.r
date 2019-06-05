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
#' @param x Name of the genlight object for which populations are to be reassigned.[required]
#' @param pop.recode Name of the file to output the new assignments [optional]
#' @param recalc -- Recalculate the locus metadata statistics if any individuals are deleted [default TRUE]
#' @param mono.rm -- Remove monomorphic loci [default TRUE]
#' @param verbose -- specify the level of verbosity: 0, silent, fatal errors only; 1, flag function begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @return An object of class ("genlight") with the revised population assignments
#' @import utils
#' @export
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' \dontrun{
#' gl <- gl.edit.recode.pop(testset.gl)
#' }

# Last amended 3-Feb-19

gl.edit.recode.pop <- function(x, pop.recode=NULL, recalc=FALSE, mono.rm=TRUE, verbose=2) {

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
    x@other$loc.metrics <- x@other$loc.metrics[1:nLoc(x),]

  # Set a population if none is specified (such as if the genlight object has been generated manually)
    if (is.null(pop(x)) | is.na(length(pop(x))) | length(pop(x)) <= 0) {
      if (verbose >= 2){ cat("  Population assignments not detected, individuals assigned to a single population labelled 'pop1'\n")}
      pop(x) <- array("pop1",dim = nLoc(x))
      pop(x) <- as.factor(pop(x))
    }

# FUNCTION SPECIFIC ERROR CHECKING

# DO THE JOB

# Take assignments from x

  if (verbose >= 2){cat("  Extracting current pop assignments from the x object\n")}
  recode.table <- cbind(levels(pop(x)),levels(pop(x)))

# Create recode table for editting, and bring up the editor
    new <- as.matrix(edit(recode.table))
    new <- new[,1:2]

# Write out the recode table, if requested
  if (is.null(pop.recode)) {
      cat("  Error: No output table specified, recode table not written to disk\n")
  } else {
    if (verbose >= 2){cat(paste("  Writing population recode table to: ",pop.recode,"\n"))}
    write.table(new, file=pop.recode, sep=",", row.names=FALSE, col.names=FALSE)    
  }

# Apply the new assignments  
  pop.list <- as.character(pop(x));
  ntr <- length(new[,1])
  for (i in 1:nInd(x)) {
    for (j in 1:ntr) {
      if (pop.list[i]==new[j,1]) {pop.list[i] <- new[j,2]}
    }
  }
  # Assigning new populations to x
  if (verbose >= 2){cat("  Assigning new population names\n")}
  pop(x) <- pop.list
  
  # If there are populations to be deleted, then recalculate relevant locus metadata and remove monomorphic loci
  
  if ("delete" %in% x$pop | "Delete" %in% x$pop) {
    # Remove rows flagged for deletion
      if (verbose >= 2){cat("Deleting populations flagged for deletion\n")}
      x <- x[!x$pop=="delete" & !x$pop=="Delete"]
    # Remove monomorphic loci
      if(mono.rm) {x <- gl.filter.monomorphs(x,verbose=verbose)}
    # Recalculate statistics
      if (recalc) {
        gl.recalc.metrics(x,verbose=verbose)
    }  
  }
  
  # REPORT A SUMMARY
  if (verbose >= 3) {
    cat("  Summary of recoded dataset\n")
    cat(paste("    No. of loci:",nLoc(x),"\n"))
    cat(paste("    No. of individuals:", nInd(x),"\n"))
    cat(paste("    No. of populations: ", length(levels(factor(pop(x)))),"\n"))
    if (!recalc) {cat("  Note: Locus metrics not recalculated\n")}
    if (!mono.rm) {cat("  Note: Resultant monomorphic loci not deleted\n")}
  }
  if (verbose>=1) {
    if (!recalc) {cat("  Note: Locus metrics not recalculated\n")}
    if (!mono.rm) {cat("  Note: Resultant monomorphic loci not deleted\n")}
  }

# FLAG SCRIPT END
  if (verbose >= 2) {cat("  Population assignments recoded\n")}
  if (verbose > 0) {
    cat("Completed:",funname,"\n")
  }
  #add to history
  nh <- length(x@other$history)
  x@other$history[[nh + 1]] <- match.call()  
  return(x)
  
}
