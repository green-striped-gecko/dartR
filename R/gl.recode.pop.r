#' Recode population assignments in a genlight object \{adegenet\}
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
#' @param x -- name of the genlight object containing the SNP data [required]
#' @param pop.recode -- name of the csv file containing the population reassignments [required]
#' @param recalc -- Recalculate the locus metadata statistics if any individuals are deleted in the filtering [default FALSE]
#' @param mono.rm -- Remove monomorphic loci [default FALSE]
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2 or as specified using gl.set.verbosity]
#' @return A genlight object with the recoded and reduced data
#' @export
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#'   mfile <- system.file("extdata", "testset_pop_recode.csv", package="dartR")
#'   nPop(testset.gl)
#'   gl <- gl.recode.pop(gl, pop.recode=mfile, verbose=3)
#' @seealso \code{\link{gl.filter.monomorphs}}
#' 

gl.recode.pop <- function(x, pop.recode, recalc=TRUE, mono.rm=TRUE, verbose=NULL){

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
  
# FUNCTION SPECIFIC ERROR CHECKING

  if (is.null(pop(x)) | is.na(length(pop(x))) | length(pop(x)) <= 0) {
    stop("  Fatal Error: Population names not detected\n")
  }

  recode.table <- read.csv(pop.recode, stringsAsFactors=FALSE, header=FALSE)
  if(length(unique(pop(x))) != length(unique(recode.table[,1]))) {
    stop("Fatal Error: Population names in data file are not the same as in the recode table\n")
  }

# DO THE JOB

  if (verbose >= 2) {
    cat("  Reassigning entities to populations as per ", pop.recode, "\n")
  }
  
# Store variables
  hold.nLoc <- nLoc(x)
  hold.nInd <- nInd(x)
  hold.nPop <- nPop(x)

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
  
  if ("delete" %in% indNames(x) | "Delete" %in% indNames(x)) {
    if (verbose >= 2){cat("Deleting populations flagged for deletion (flagged 'Delete' or 'delete')\n")}
    x <- gl.drop.pop(x,pop.list=c("Delete","delete"),verbose=0)
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
  
  if (verbose>=3) {
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
