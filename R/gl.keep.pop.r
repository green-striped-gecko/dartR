#' Remove all but specified populations from a genelight \{adegenet\} object
#'
#' Individuals are assigned to populations based on the specimen metadata data file (csv) used with gl.read.dart(). 
#'
#' The script, having deleted the specified populations, optionally identifies resultant monomorphic loci or loci
#' with all values missing and deletes them (using gl.filter.monomorphs.r). The script also optionally
#' recalculates statistics made redundant by the deletion of individuals from the dataset.
#' 
#' The script returns a genlight object with the new population assignments and the recalculated locus metadata.
#'
#' @param x -- name of the genlight object containing SNP genotypes or a genind object containing presence/absence data [required]
#' @param pop.list -- a list of populations to be kept [required]
#' @param as.pop -- assign another metric to represent population [default NULL]
#' @param recalc -- Recalculate the locus metadata statistics [default FALSE]
#' @param mono.rm -- Remove monomorphic loci [default FALSE]
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2 or as specified using gl.set.verbosity]
#' @return A genlight object with the reduced data
#' @export
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#'  # SNP data
#'    gl2 <- gl.keep.pop(testset.gl, pop.list=c("EmsubRopeMata","EmvicVictJasp"))
#'    gl2 <- gl.keep.pop(testset.gl, pop.list=c("EmsubRopeMata","EmvicVictJasp"),mono.rm=TRUE,recalc=TRUE)
#'    gl2 <- gl.keep.pop(testset.gl, pop.list=c("Female"),as.pop="sex")
#'  # Tag P/A data  
#'    gs2 <- gl.keep.pop(testset.gs, pop.list=c("EmsubRopeMata","EmvicVictJasp"))
#'
#' @seealso \code{\link{gl.filter.monomorphs}} for when mono.rm=TRUE, \code{\link{gl.recalc.metrics}} for when recalc=TRUE
#' @seealso \code{\link{gl.drop.pop}} to drop rather than keep specified populations

gl.keep.pop <- function(x, pop.list, as.pop=NULL, recalc=FALSE, mono.rm=FALSE, verbose=NULL){

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
  
  # Population labels assigned?
  if(is.null(as.pop)){
    if (is.null(pop(x)) | is.na(length(pop(x))) | length(pop(x)) <= 0) {
      stop("Fatal Error: Population assignments not detected, run gl.compliance.check() and revisit population assignments\n")
    }
  }
    
  # Assign the new population list if as.pop is specified
  pop.hold <- pop(x)
  if (!is.null(as.pop)){    
    if(as.pop %in% names(x@other$ind.metrics)){
      pop(x) <- as.matrix(x@other$ind.metrics[as.pop])
      if (verbose >= 2) {cat("  Temporarily setting population assignments to",as.pop,"as specified by the as.pop parameter\n")}
    } else {
      stop("Fatal Error: individual metric assigned to 'pop' does not exist. Check names(gl@other$loc.metrics) and select again\n")
    }
  }

  if (verbose >= 2) {
    cat("  Checking for presence of nominated populations\n")
  }
  
  for (case in pop.list){
    if (!(case%in%popNames(x))){
      cat("  Warning: Listed population",case,"not present in the dataset -- ignored\n")
      pop.list <- pop.list[!(pop.list==case)]
    }
  }
  if (length(pop.list) == 0) {
    stop("  Fatal Error: no populations listed to keep!\n")
  }

# DO THE JOB

  if (verbose >= 2) {
    cat("  Retaining only populations", pop.list, "\n")
  }

  # Delete all but the listed populations, recalculate relevant locus metadata and remove monomorphic loci
  
  # Keep only rows flagged for retention
    x2 <- x[x$pop%in%pop.list]
    pop.hold <- pop.hold[x$pop%in%pop.list]
    x <- x2

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
      if (!is.null(as.pop)) {
        cat("  Summary of recoded dataset\n")
        cat(paste("    No. of loci:",nLoc(x),"\n"))
        cat(paste("    No. of individuals:", nInd(x),"\n"))
        cat(paste("    No. of levels of",as.pop,"remaining: ",nPop(x),"\n"))
        cat(paste("    No. of populations: ",nPop(pop.hold),"\n"))
      } else {
        cat("  Summary of recoded dataset\n")
        cat(paste("    No. of loci:",nLoc(x),"\n"))
        cat(paste("    No. of individuals:", nInd(x),"\n"))
        cat(paste("    No. of populations: ",nPop(x),"\n"))
      }  
    }
    
  # Reassign the initial population list if as.pop is specified
    
    if (!is.null(as.pop)){
      pop(x) <- pop.hold
      if (verbose >= 3) {cat("  Resetting population assignments to initial state\n")}
    }
    
# ADD TO HISTORY
  nh <- length(x@other$history)
  x@other$history[[nh + 1]] <- match.call() 
    
# FLAG SCRIPT END
    
  if (verbose >= 1) {
    cat("Completed:",funname,"\n")
  }
    
  return(x)
}
