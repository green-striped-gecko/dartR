#' Recode individual (=specimen = sample) labels in a genelight object \{adegenet\}
#'
#' This script recodes individual labels and/or deletes individuals from a DaRT genlight SNP file
#' based on a lookup table provided as a csv file.
#'
#' Renaming individuals may be required when there have been errors in labelling arising
#' in the process from sample to DArT files. There may be occasions where renaming
#' individuals is required for preparation of figures. When caution needs to be exercised
#' because of the potential for breaking the "chain of evidence" associated with the samples,
#' recoding individuals using a recode table (csv) can provide a clear record of the changes.
#' 
#' The script, having deleted individuals, optionally identifies resultant monomorphic loci or loci
#' with all values missing and deletes them (using gl.filter.monomorphs.r). The script also optionally
#' recalculates statistics made incorrect by the deletion of individuals from the dataset.
#' 
#' The script returns a genlight object with the new individual labels, the monomorphic loci optionally removed
#' and the optionally recalculated locus metadata.
#'
#' @param x -- name of the genlight object containing SNP genotypes [required]
#' @param ind.recode -- name of the csv file containing the individual relabelling [required]
#' @param recalc -- if TRUE, recalculate the locus metadata statistics if any individuals are deleted in the filtering [default FALSE]
#' @param mono.rm -- if TRUE, remove monomorphic loci [default FALSE]
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @return A genlight or genind object with the recoded and reduced data
#' @export
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' \dontrun{
#'    gl <- gl.recode.ind(testset.gl, ind.recode="testset_pop_recode.csv")
#' }
#' @seealso \code{\link{gl.filter.monomorphs}} for filtering monomorphs, \code{\link{gl.recalc.metrics}} for recalculating locus metrics,
#' \code{\link{gl.recode.pop}} for recoding populations

gl.recode.ind <- function(x, ind.recode, recalc=FALSE, mono.rm=FALSE, verbose=2){

# TIDY UP FILE SPECS

  #outfilespec <- file.path(outpath, outfile)
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
  # Set a population if none is specified (such as if the genlight object has been generated manually)
    if (is.na(length(pop(x))) | length(pop(x)) <= 0) {
      if (verbose >= 2){cat("  Population assignments not detected, assigning all individuals to one population, label 'pop1'\n")}
      pop(x) <- array("pop1",dim = nInd(x))
      pop(x) <- as.factor(pop(x))
    }
   # Check for monomorphic loci
    tmp <- gl.filter.monomorphs(x)
    if ((nLoc(tmp) < nLoc(x)) & verbose >= 2) {cat("  Warning: genlight object contains monomorphic loci\n")}

# DO THE JOB

  if (verbose >= 2){
    cat("  Relabelling individuals (=specimens) as per ", ind.recode, "\n")
    cat("    Reading lookup table\n")
  }
  recode.table <- read.csv(ind.recode, stringsAsFactors=FALSE, header=FALSE);
# Error check
  if(length(unique(indNames(x))) != length(unique(recode.table[,1]))) {
    cat("  Fatal Error: Individual names in data file are not the same as in the recode table\n"); stop()
  }
# Apply the recode to the individuals
  if (verbose >= 2){
    cat("    Applying the recoding\n")
  }
  ind.list <- as.character(indNames(x))
  ntr <- length(recode.table[,1])
  for (i in 1:nInd(x)) {
    for (j in 1:ntr) {
      if (ind.list[i]==recode.table[j,1]) {ind.list[i] <- recode.table[j,2]}
    }
  }
  indNames(x) <- ind.list

  # If there are individuals to be deleted, then recalculate relevant locus metadata and remove monomorphic loci
  
  if ("delete" %in% x$ind.names | "Delete" %in% x$ind.names) {
    # Remove rows flagged for deletion
      if (verbose >= 2){
        cat("    Deleting individuals or samples flagged for deletion\n")
      }
      x2 <- x[!x$ind.names=="delete" & !x$ind.names=="Delete"]
    #  Remove monomorphic loci
      if (mono.rm) {x2 <- gl.filter.monomorphs(x2,verbose=verbose)}
    # Recalculate statistics
      if (recalc) {x2 <- gl.recalc.metrics(x2,verbose=verbose)}
  } else {
    x2 <- x
  }

# REPORT A SUMMARY
  if (verbose >= 3) {
    cat("  Summary of recoded dataset\n")
    cat(paste("    No. of loci:",nLoc(x2),"\n"))
    cat(paste("    No. of individuals:", nInd(x2),"\n"))
    cat(paste("    No. of populations: ", length(levels(factor(pop(x2)))),"\n"))
  }
  if (verbose >= 2) {
    if (!recalc) {
      cat("  Note: Locus metrics not recalculated\n")
    } else {
      cat("  Note: Locus metrics recalculated\n")
    }
    if (!mono.rm) {
      cat("  Note: Resultant monomorphic loci not deleted\n")
    } else{
      cat("  Note: Resultant monomorphic loci deleted\n")
    }
  }

# FLAG SCRIPT END

  if (verbose > 0) {
    cat("Completed:",funname,"\n")
  }
  #add to history
  nh <- length(x2@other$history)
  x2@other$history[[nh + 1]] <- match.call()  
    return(x2)
}

