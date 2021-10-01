#' @name gl.drop.loc
#' @title Remove specified loci from a genlight \{adegenet\} object
#' @description
#' The script returns a genlight object with specified loci deleted.
#'
#' @param x Name of the genlight object containing SNP genotypes or presence/absence data [required]
#' @param loc.list A list of loci to be deleted [required, if loc.range not specified]
#' @param first First of a range of loci to be deleted [required, if loc.list not specified]
#' @param last Last of a range of loci to be deleted [if not specified, last locus in the dataset]
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 
#' 3, progress and results summary; 5, full report [default 2 or as specified using gl.set.verbosity]
#' @return A genlight object with the reduced data
#' 
#' @author Custodian: Arthur Georges -- Post to \url{https://groups.google.com/d/forum/dartr}
#' 
#' @examples
#' # SNP data
#'   gl2 <- gl.drop.loc(testset.gl, loc.list=c("100051468|42-A/T", "100049816-51-A/G"),verbose=3)
#' # Tag P/A data
#'   gs2 <- gl.drop.loc(testset.gs, loc.list=c("20134188","19249144"),verbose=3)
#'   
#' @seealso \code{\link{gl.keep.loc}} to keep rather than drop specified loci
#' @export
#' 
gl.drop.loc <- function(x, 
                        loc.list = NULL, 
                        first = NULL, 
                        last = NULL, 
                        verbose = NULL){

# SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
# FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func=funname,build="Jody",v=verbose)
  
# CHECK DATATYPE 
  datatype <- utils.check.datatype(x,verbose=verbose)

# FUNCTION SPECIFIC ERROR CHECKING

  if (!is.null(loc.list) && !is.null(first)){
    flag <- 'both'
    if (verbose >= 2){
      cat(report("  Both a range of loci and a list of loci to keep has been specified\n"))
    } 
  } else if (!is.null(loc.list)){
    flag <- 'list'
    if (verbose >= 2){
      cat(report("  List of loci to drop has been specified\n"))
    } 
  } else if (!is.null(first)){
    flag <- 'range'
    if (verbose >= 2){
      cat(report("  Range of loci to drop has been specified\n"))
    } 
  } else {
      stop(error("Fatal Error: Need to specify either a range of loci to drop, or specific loci to drop\n"))
  }
  
  if (flag=='both' || flag=='list'){
    for (case in loc.list){
      if (!(case%in%locNames(x))){
        cat(warn("  Warning: Listed loci",case,"not present in the dataset -- ignored\n"))
        loc.list <- loc.list[!(loc.list==case)]
      }
    }
  }

  if (flag=='range'){
    if (first <=0){
      cat(warn("  Warning: Lower limit to range of loci cannot be less than 1, set to 1\n)"))
      first <- 1
    }
    if (first > nLoc(x)){
      cat(warn("  Warning: Upper limit to range of loci cannot be greater than the number of loci, set to",nLoc(x),"\n)"))
      last <- nLoc(x)
    }
    if (first > last){
      cat(warn("  Warning: Upper limit is smaller than lower limit, reversed\n"))
      tmp <- first
      first <- last
      last <- tmp
    }
  }

# DO THE JOB
  
  hold <- x

  if (verbose >= 2) {
    cat(report("  Deleting the specified loci\n"))
  }

  # Remove duplicated loci if specified

  if (!is.null(first) && !is.null(loc.list)){
    list.from.range <- locNames(x)[first:last]
    loc.list <- unique(c(loc.list,list.from.range))
  } else if (!is.null(first)) {
    loc.list <- locNames(x)[first:last]
  }
  if (length(loc.list) == 0) {
    cat(warn("  Warning: no loci listed to delete! Genlight object returned unchanged\n"))
    x2 <- x
  } else {
    # Remove loci flagged for deletion
    x2 <- x[,!x$loc.names%in%loc.list]
    x2@other$loc.metrics <- x@other$loc.metrics[!x$loc.names%in%loc.list,]
  }  

# REPORT A SUMMARY
    
  if (verbose >= 3) {
    cat("  Summary of recoded dataset\n")
    cat(paste("    Original No. of loci:",nLoc(hold),"\n"))
    cat(paste("    No. of loci deleted:",nLoc(hold)-nLoc(x2),"\n"))
    cat(paste("    No. of loci retained:",nLoc(x2),"\n"))
    # cat(paste("    No. of individuals:", nInd(x2),"\n"))
    # cat(paste("    No. of populations: ", nPop(x2),"\n"))
  }

# ADD TO HISTORY    
    nh <- length(x2@other$history)
    x2@other$history[[nh + 1]] <- match.call()
  
# FLAG SCRIPT END

  if (verbose >= 1){  
     cat(report("Completed:",funname,"\n"))
  }
    
  return(x2)
}

