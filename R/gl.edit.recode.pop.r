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
#' Use outpath=getwd() or outpath="." when calling this function to direct output files to your working directory.
#' 
#' The script returns a genlight object with the new population assignments and the recalculated locus metadata.
#' 
#' @param x Name of the genlight object for which populations are to be reassigned.[required]
#' @param out.recode.file Name of the file to output the new individual labels [optional]
#' @param outpath -- path where to save the output file [default tempdir(), mandated by CRAN].
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

gl.edit.recode.pop <- function(x, pop.recode=NULL, out.recode.file=NULL, outpath=tempdir(), recalc=FALSE, mono.rm=FALSE, verbose=2) {

  # TIDY UP FILE SPECS
  
  funname <- match.call()[[1]]
  build <- "Jacob"
  if (!is.null(out.recode.file)){
    outfilespec <- file.path(outpath, out.recode.file)
  }
  
  # FLAG SCRIPT START
  
  if (verbose < 0 | verbose > 5){
    cat("  Warning: Parameter 'verbose' must be an integer between 0 [silent] and 5 [full report], set to 2\n")
    verbose <- 2
  }
  
  if (verbose >= 1){
    cat("Starting",funname,"[ Build =",build,"]\n")
  }
  
  # STANDARD ERROR CHECKING
  
  if(class(x)!="genlight") {
    stop("Fatal Error: genlight object required!\n")
  }
  
  if (all(x@ploidy == 1)){
    if (verbose >= 2){cat("  Processing  Presence/Absence (SilicoDArT) data\n")}
    data.type <- "SilicoDArT"
  } else if (all(x@ploidy == 2)){
    if (verbose >= 2){cat("  Processing a SNP dataset\n")}
    data.type <- "SNP"
  } else {
    stop("Fatal Error: Ploidy must be universally 1 (fragment P/A data) or 2 (SNP data)")
  }
  
  # FUNCTION SPECIFIC ERROR CHECKING
  
  if (is.null(pop(x)) | is.na(length(pop(x))) | length(pop(x)) <= 0) {
    stop("  Fatal Error: Population names not detected\n")
  }

# DO THE JOB

  # Store variables
  hold.nLoc <- nLoc(x)
  hold.nInd <- nInd(x)
  hold.nPop <- nPop(x)
  
  # Take assignments from x

  if (verbose >= 2){cat("  Extracting current pop assignments from the x object\n")}
  recode.table <- cbind(levels(pop(x)),levels(pop(x)))

# Create recode table for editting, and bring up the editor
    new <- as.matrix(edit(recode.table))
    new <- new[,1:2]

# Write out the recode table, if requested
  if (is.null(pop.recode)) {
      cat("  No output table specified, recode table not written to disk\n")
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
  if (verbose >= 2){cat("Assigning new population names\n")}
  pop(x) <- pop.list
  
  # If there are populations to be deleted, then recalculate relevant locus metadata and remove monomorphic loci
  
  if ("delete" %in% x$pop | "Delete" %in% x$pop) {
   # Remove populations flagged for deletion
    if (verbose >= 2){cat("Deleting individuals/samples flagged for deletion (Flagged 'Delete' or 'delete')\n")}
    x <- gl.drop.pop(x,pop.list=c('Delete','delete'),verbose=0)
  }
  
  # Recalculate statistics
  if (recalc) {
    x <- gl.recalc.metrics(x,verbose=0)
  } 
  
  #  Remove monomorphic loci
  if (mono.rm) {
    x <- gl.filter.monomorphs(x,verbose=0)
  }
  
# REPORT A SUMMARY
  
  if (verbose>=2) {
    cat("\n  Summary of recoded dataset\n")
    cat(paste("  Original No. of loci:",hold.nLoc,"\n"))
    cat(paste("    New No. of loci:",nLoc(x),"\n"))
    cat(paste("  Original No. of individuals:", hold.nInd,"\n"))
    cat(paste("    New No. of individuals:", nInd(x),"\n"))
    cat(paste("  Original No. of populations:", hold.nPop,"\n"))
    cat(paste("    New No. of populations:", nPop(x),"\n\n"))
    if (!recalc) {cat("Note: Locus metrics not recalculated\n")}
    if (!mono.rm) {cat("Note: Resultant monomorphic loci not deleted\n")}
  }
  

# FLAG SCRIPT END
  
  if (verbose >= 2) {cat("\n  Success: Populaton assignments recoded\n\n")}  
  
  if (verbose > 0) {
    cat("Completed:",funname,"\n")
  }
  #add to history
  nh <- length(x@other$history)
  x@other$history[[nh + 1]] <- match.call()  
  return(x)
  
}

# # Test script
# gl <- gl.edit.recode.pop(testset.gl,verbose=2)
# gl <- gl.edit.recode.ind(gs, out.recode.file="ind.recode.table.csv",verbose=2)

