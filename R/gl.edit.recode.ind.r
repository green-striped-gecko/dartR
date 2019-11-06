#' Create or edit a individual (=specimen) names, create an recode_ind file amd apply the changes to a genlight object.
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
#' Use outpath=getwd() or outpath="." when calling this function to direct output files to your working directory.
#' 
#' The script returns a genlight object with the new individual labels and the recalculated locus metadata.
#' 
#' @param x Name of the genlight object for which individuals are to be relabelled.[required]
#' @param out.recode.file Name of the file to output the new individual labels [optional]
#' @param outpath -- path where to save the output file [default tempdir(), mandated by CRAN].
#' @param recalc -- Recalculate the locus metadata statistics [default TRUE]
#' @param mono.rm -- Remove monomorphic loci [default TRUE]
#' @return An object of class ("genlight") with the revised individual labels
#' @param verbose -- verbose=0, silent; verbose=1, low verbosity; verbose=2, high verbosity [default 2]
#' @import utils
#' @export
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' \dontrun{
#' gl <- gl.edit.recode.ind(testset.gl)
#' gl <- gl.edit.recode.ind(testset.gl, out.recode.file="ind.recode.table.csv")
#' gl <- gl.edit.recode.ind(testset.gl, out.recode.file="ind.recode.table.csv")
#' }

gl.edit.recode.ind <- function(x, out.recode.file=NULL, outpath=tempdir(), recalc=FALSE, mono.rm=FALSE, verbose=2){

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
  
# DO THE JOB
  
# Store variables
  hold.nLoc <- nLoc(x)
  hold.nInd <- nInd(x)
  hold.nPop <- nPop(x)
  
# Take assignments from x

  if(verbose >= 2){cat("Extracting current individual labels from the x object\n")}
  recode.table <- cbind(indNames(x),indNames(x))

# Create recode table for editting, and bring up the editor
    new <- as.matrix(edit(recode.table))
    new <- new[,1:2]

# Write out the recode table, if requested
  if (is.null(out.recode.file)) {
    if(verbose >= 2){cat("No output table specified, recode table not written to disk\n")}
  } else {
    if(verbose >= 2){cat(paste("Writing individual recode table to: ",out.recode.file,"\n"))}
    write.table(new, file=out.recode.file, sep=",", row.names=FALSE, col.names=FALSE)    
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
  if(verbose >= 2){cat("Assigning new individual (=specimen) names\n")}
  indNames(x) <- ind.list
  
  # If there are populations to be deleted, then recalculate relevant locus metadata and remove monomorphic loci
  
  if ("delete" %in% indNames(x) | "Delete" %in% indNames(x)) {
    # Remove populations flagged for deletion
    if (verbose >= 2){cat("Deleting individuals/samples flagged for deletion (Flagged 'Delete' or 'delete')\n")}
    x <- gl.drop.ind(x,ind.list=c('Delete','delete'),verbose=0)
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
 
# ADD TO HISTORY 
    nh <- length(x@other$history)
    x@other$history[[nh + 1]] <- match.call() 
  
# FLAG SCRIPT END

  if (verbose > 0) {
    cat("Completed:",funname,"\n")
  }
    
return(x)
  
}
