#' Remove specified loci from a genelight \{adegenet\} object
#'
#' The script returns a genlight object with specified loci deleted.
#'
#' @param x -- name of the genlight object containing SNP genotypes or presence/absence data [required]
#' @param loc.list -- a list of loci to be deleted [required, if loc.range not specified]
#' @param first -- first of a range of loci to be deleted [required, if loc.list not specified]
#' @param last -- last of a range of loci to be deleted [if not specified, last locus in the dataset]
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2 or as specified using gl.set.verbosity]
#' @return A genlight object with the reduced data
#' @export
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#'    gl <- gl.drop.loc(testset.gl, loc.list=c("100051468|42-A/T", "100049816-51-A/G"))

gl.drop.loc <- function(x, loc.list=NULL, first=NULL, last=NULL, verbose=NULL){

# TRAP COMMAND, SET VERSION
  
  funname <- match.call()[[1]]
  build <- "Jacob"
  hold <- x
  
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
    stop("  Fatal Error: genlight object required!\n")
  }
  
  if (verbose >= 1){
    if (all(x@ploidy == 1)){
      cat("  Processing Presence/Absence (SilicoDArT) data\n")
    } else if (all(x@ploidy == 2)){
      cat("  Processing a SNP dataset\n")
    } else {
      stop ("Fatal Error: Ploidy must be universally 1 (fragment P/A data) or 2 (SNP data)")
    }
  }
  
# FUNCTION SPECIFIC ERROR CHECKING

  if (!is.null(loc.list) && !is.null(first)){
    flag <- 'both'
    if (verbose >= 2){
      cat("  Both a range of loci and a list of loci to keep has been specified\n")
    } 
  } else if (!is.null(loc.list)){
    flag <- 'list'
    if (verbose >= 2){
      cat("  List of loci to keep has been specified\n")
    } 
  } else if (!is.null(first)){
    flag <- 'range'
    if (verbose >= 2){
      cat("  Range of loci to keep has been specified\n")
    } 
  } else {
      stop("  Fatal Error: Need to specify either a range of loci to keep, or specific loci to keep\n")
  }
  
  if (flag=='both' || flag=='list'){
    for (case in loc.list){
      if (!(case%in%locNames(x))){
        cat("  Warning: Listed loci",case,"not present in the dataset -- ignored\n")
        loc.list <- loc.list[!(loc.list==case)]
      }
    }
  }

  if (flag=='range'){
    if (first <=0){
      cat("  Warning: Lower limit to range of loci cannot be less than 1, set to 1\n)")
      first <- 1
    }
    if (first > nLoc(x)){
      cat("  Warning: Upper limit to range of loci cannot be greater than the number of loci, set to",nLoc(x),"\n)")
      last <- nLoc(x)
    }
    if (first > last){
      cat("  Warning: Upper limit is smaller than lower limit, reversed\n")
      tmp <- first
      first <- last
      last <- tmp
    }
  }

# DO THE JOB

  if (verbose >= 2) {
    cat("    Deleteing the specified loci\n")
  }

  # Remove duplicated loci if specified

  if (!is.null(first) && !is.null(loc.list)){
    list.from.range <- locNames(x)[first:last]
    loc.list <- unique(c(loc.list,list.from.range))
  } else if (!is.null(first)) {
    loc.list <- locNames(x)[first:last]
  }
  if (length(loc.list) == 0) {
    cat("  Warning: no loci listed to delete! Genlight object returned unchanged\n")
    x2 <- x
  } else {
    # Remove loci flagged for deletion
    x2 <- x[,!x$loc.names%in%loc.list]
    x2@other$loc.metrics <- x@other$loc.metrics[!x$loc.names%in%loc.list,]
  }  

# REPORT A SUMMARY
    
  if (verbose >= 3) {
    cat("\n  Summary of recoded dataset\n")
    cat(paste("    Original No. of loci:",nLoc(hold),"\n"))
    cat(paste("    No. of loci deleted:",nLoc(hold)-nLoc(x2),"\n"))
    cat(paste("    No. of loci retained:",nLoc(x2),"\n"))
    cat(paste("    No. of individuals:", nInd(x2),"\n"))
    cat(paste("    No. of populations: ", length(levels(factor(pop(x2)))),"\n\n"))
  }

# ADD TO HISTORY    
    nh <- length(x2@other$history)
    x2@other$history[[nh + 1]] <- match.call()
  
# FLAG SCRIPT END

  if (verbose >= 1){  
     cat("Completed:",funname,"\n")
  }
    
  return(x2)
    
}

