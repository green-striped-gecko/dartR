#' Remove all but the specified individuals from a genelight \{adegenet\} object
#'
#' The script, having deleted individuals, optionally identifies resultant monomorphic loci or loci
#' with all values missing and deletes them (using gl.filter.monomorphs.r). The script also optionally
#' recalculates statistics made redundant by the deletion of individuals from the dataset.
#' 
#' The script returns a genlight object with the individuals deleted and, optionally, the recalculated locus metadata.
#'
#' @param x -- name of the genlight object containing SNP genotypes or a genind object containing presence/absence data [required]
#' @param ind.list -- a list of individuals to be removed [required]
#' @param recalc -- Recalculate the locus metadata statistics [default FALSE]
#' @param mono.rm -- Remove monomorphic loci [default FALSE]
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2 or as specified using gl.set.verbosity]
#' @return A genlight object with the reduced data
#' @export
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#'   # SNP data
#'     gl2 <- gl.keep.ind(testset.gl, ind.list=c("AA019073","AA004859"))
#'   # Tag P/A data
#'    gs2 <- gl.keep.ind(testset.gs, ind.list=c("AA020656","AA19077","AA004859"))
#'    
#' @seealso \code{\link{gl.filter.monomorphs}} for when mono.rm=TRUE, \code{\link{gl.recalc.metrics}} for when recalc=TRUE
#' @seealso \code{\link{gl.drop.ind}} to drop rather than keep specified individuals

gl.keep.ind <- function(x, ind.list, recalc=FALSE, mono.rm=FALSE, verbose=NULL){

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

  for (case in ind.list){
    if (!(case%in%indNames(x))){
      cat("  Warning: Listed individual",case,"not present in the dataset -- ignored\n")
      ind.list <- ind.list[!(ind.list==case)]
    }
  }
  if (length(ind.list) == 0) {
    stop("  Fatal Error: no individuals listed to keep!\n")
  }

# DO THE JOB

  if (verbose >= 2) {
    cat("    Deleteing all but the listed individuals", ind.list, "\n")
  }

  # Delete all but the listed individuals, recalculate relevant locus metadata and remove monomorphic loci
  
  # Remove rows flagged for deletion
    x <- x[x$ind.names%in%ind.list]
    
  # Monomorphic loci may have been created
    x@other$loc.metrics.flags$monomorphs == FALSE
    
  # Remove monomorphic loci
    if(mono.rm){
      if(verbose >= 2){cat("  Deleting monomorphic loc\n")}
      x <- gl.filter.monomorphs(x,verbose=0)
    } 
  # Check monomorphs have been removed
    if (x@other$loc.metrics.flags$monomorphs == FALSE){
      if (verbose >= 2){
        cat("  Warning: Resultant dataset may contain monomorphic loci\n")
      }  
    }
    
  # Recalculate statistics
    if (recalc) {
      x <- gl.recalc.metrics(x,verbose=0)
      if(verbose >= 2){cat("  Recalculating locus metrics\n")}
    } else {
      if(verbose >= 2){
        cat("  Locus metrics not recalculated\n")
        x <- utils.reset.flags(x,verbose=0)
      }
    }

# REPORT A SUMMARY
    
  if (verbose >= 3) {
    cat("  Summary of recoded dataset\n")
    cat(paste("    No. of loci:",nLoc(x),"\n"))
    cat(paste("    No. of individuals:", nInd(x),"\n"))
    cat(paste("    No. of populations:", nPop(x),"\n"))
  }
    
# ADD TO HISTORY
    nh <- length(x@other$history)
    x@other$history[[nh + 1]] <- match.call() 
    
# FLAG SCRIPT END

    if (verbose > 0) {
      cat("Completed: gl.keep.ind\n")
    }
  
    return(x)
}

