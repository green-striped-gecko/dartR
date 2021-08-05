#' @name gl.recode.ind
#' @title Recode individual (=specimen = sample) labels in a genelight object \{adegenet\}
#' @description 
#' This script recodes individual labels and/or deletes individuals from a DaRT genlight SNP file
#' based on a lookup table provided as a csv file.
#' @details
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
#'   file <- system.file("extdata","testset_ind_recode.csv", package="dartR")
#'   gl <- gl.recode.ind(testset.gl, ind.recode=file, verbose=3)
#' @seealso \code{\link{gl.filter.monomorphs}} for filtering monomorphs, \code{\link{gl.recalc.metrics}} for recalculating locus metrics,
#' \code{\link{gl.recode.pop}} for recoding populations

gl.recode.ind <- function(x, 
                          ind.recode, 
                          recalc = FALSE, 
                          mono.rm = FALSE, 
                          verbose = NULL){

# SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
# FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func=funname,build="Jackson",v=verbose)
  
# CHECK DATATYPE 
  datatype <- utils.check.datatype(x,verbose=verbose)
  
# SCRIPT SPECIFIC ERROR CHECKING
  # change stringsAsFactors=FALSE due to the new default in r
  recode.table <- read.csv(ind.recode,  header=FALSE, stringsAsFactors = FALSE)
  
  v1 <- unique(indNames(x))
  v2 <- unique(recode.table[,1])
  v1_v2 <- v1[!(v1 %in% v2)]
  l1 <- length(v1)
  l2 <- length(v2)
  
  if(l1 != l2) {
    stop(error("Fatal Error: Individuals do not agree in number with those listed in the recode table\n"))
  }
  if(!(length(v1_v2)==0)){
    stop(error("Fatal Error: Some individuals have no reassignment specified in the recode table:",v1_v2,"\n"))
  }
  
# DO THE JOB

  if (verbose >= 2){
    cat(report("  Relabelling individuals (=specimens) as per ", ind.recode, "\n"))
    cat(report("    Reading lookup table\n"))
  }
  
  # Store variables
  hold.nLoc <- nLoc(x)
  hold.nInd <- nInd(x)
  hold.nPop <- nPop(x)

# Apply the recode to the individuals
  if (verbose >= 2){
    cat(report("    Applying the recoding\n"))
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
        cat(report("    Deleting individuals/samples flagged for deletion\n"))
      }
    deletions <- indNames(x)[tolower(recode.table[,2])=="delete"]
    if (verbose==3){
      cat("  Dropping\n",paste(deletions,collapse=", "),"\n")
      cat("  A total of",length(deletions),"individuals dropped\n")
    }
    x <- gl.drop.ind(x,ind.list=c("Delete","delete"),verbose=0)
  } 
  
  # Remove monomorphic loci
  if(mono.rm){
    if(verbose >= 2){cat(report("  Deleting monomorphic loc\n"))}
    x <- gl.filter.monomorphs(x,verbose=0)
  } 
  # Check monomorphs have been removed
  if (x@other$loc.metrics.flags$monomorphs == FALSE){
    if (verbose >= 2){
      cat(warn("  Warning: Resultant dataset may contain monomorphic loci\n"))
    }  
  }
  
  # Recalculate statistics
  if (recalc) {
    x <- gl.recalc.metrics(x,verbose=0)
    if(verbose >= 2){cat(report("  Recalculating locus metrics\n"))}
  } else {
    if(verbose >= 2){
      cat(warn("  Locus metrics not recalculated\n"))
      x <- utils.reset.flags(x,verbose=0)
    }
  }
  
# REPORT A SUMMARY
  
  if (verbose>=2) {
    cat("  Summary of recoded dataset\n")
    cat(paste("  Original No. of loci:",hold.nLoc,"\n"))
    cat(paste("    New No. of loci:",nLoc(x),"\n"))
    cat(paste("  Original No. of individuals:", hold.nInd,"\n"))
    cat(paste("    New No. of individuals:", nInd(x),"\n"))
    cat(paste("  Original No. of populations:", hold.nPop,"\n"))
    cat(paste("    New No. of populations:", nPop(x),"\n"))
    if (!recalc) {cat(report("  Note: Locus metrics not recalculated\n"))}
    if (!mono.rm) {cat(report("  Note: Resultant monomorphic loci not deleted\n"))}
  }
  
# ADD TO HISTORY
    nh <- length(x@other$history)
    x@other$history[[nh + 1]] <- match.call()  

# FLAG SCRIPT END

  if (verbose > 0) {
    cat(report("Completed:",funname,"\n"))
  }
    
  return(x)
}
