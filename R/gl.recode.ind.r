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
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2 or as specified using gl.set.verbosity]
#' @return A genlight or genind object with the recoded and reduced data
#' @export
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#'   file <- system.file("extdata","testset_pop_recode.csv", package="dartR")
#'   #gl <- gl.recode.ind(testset.gl, ind.recode=file, verbose=0)
#' @seealso \code{\link{gl.filter.monomorphs}} for filtering monomorphs, \code{\link{gl.recalc.metrics}} for recalculating locus metrics,
#' \code{\link{gl.recode.pop}} for recoding populations

gl.recode.ind <- function(x, ind.recode, recalc=FALSE, mono.rm=FALSE, verbose=NULL){

# TRAP COMMAND, SET VERSION
  
  funname <- match.call()[[1]]
  build <- "Jacob"
  
# SET VERBOSITY
  
  if (is.null(verbose)){ 
    if(!is.null(x@other$verbose)){ 
      verbose <- x@other$verbose
    } else { 
      verbose <- 2
    }
  } 
  
  if (verbose < 0 | verbose > 5){
    cat(paste("  Warning: Parameter 'verbose' must be an integer between 0 [silent] and 5 [full report], set to 2\n"))
    verbose <- 2
  }
  
# FLAG SCRIPT START
  
  if (verbose >= 1){
    if(verbose==5){
      cat("Starting",funname,"[ Build =",build,"]\n")
    } else {
      cat("Starting",funname,"\n")
    }
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
  
# SCRIPT SPECIFIC ERROR CHECKING
  
  recode.table <- read.csv(ind.recode, stringsAsFactors=FALSE, header=FALSE, stringsAsFactors = TRUE);
  if(length(unique(indNames(x))) != length(unique(recode.table[,1]))) {
    stop("  Fatal Error: Individual names in data file are not the same as in the recode table\n")
  }
  
# DO THE JOB

  if (verbose >= 2){
    cat("  Relabelling individuals (=specimens) as per ", ind.recode, "\n")
    cat("    Reading lookup table\n")
  }
  
  # Store variables
  hold.nLoc <- nLoc(x)
  hold.nInd <- nInd(x)
  hold.nPop <- nPop(x)

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
        cat("    Deleting individuals/samples flagged for deletion\n")
      }
      x <- gl.drop.ind(x,ind.list=c("Delete","delete"),verbose=0)
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
