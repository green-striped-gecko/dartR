#' Remove monomorphic loci, including those with all NAs
#'
#' This script deletes monomorphic loci from a genlight \{adegenet\} object
#'
#' A DArT dataset will not have monomorphic loci, but they can arise, along with loci that are scored all NA, when populations or individuals are deleted.
#' Retaining monomorphic loci unnecessarily increases the size of the dataset and will affect some calculations.
#' 
#' Note that for SNP data, NAs likely represent null alleles; in tag presence/absence data, NAs represent missing values (presence/absence could not 
#' be reliably scored)
#' 
#' @param x -- name of the input genlight object [required]
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @return A genlight object with monomorphic ( and all NA) loci removed
#' @import utils
#' @importFrom plyr count
#' @export
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' gl <- gl.filter.monomorphs(testset.gl, verbose=3)

gl.filter.monomorphs <- function (x, verbose=2) {

  # TIDY UP FILE SPECS
  
  build ='Jacob'
  funname <- match.call()[[1]]
  # Note does draw upon the monomorphs flag and will reset it to TRUE on completion
  
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
    stop("Fatal Error: genlight object required!")
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
  
  # mml <- !( colMeans(as.matrix(x), na.rm=TRUE)%%2 == 0) ; Code for readability

  hold <- x
  na.counter <- 0
  loc.list <- array(NA,nLoc(x))

  if (verbose >= 2){
    cat("Identifying monomorphic loci\n")
  }  
  
  # Tag presence/absence data
  if (data.type=="SilicoDArT"){
    matrix <- as.matrix(x)
    l.names <- locNames(x)
    for (i in 1:nLoc(x)){
      row <- matrix[,i] # Row for each locus
      if (all(row == 0, na.rm=TRUE) | all(row == 1, na.rm=TRUE) | all(is.na(row))){
        loc.list[i] <- l.names[i]
        if (all(is.na(row))){
          na.counter = na.counter + 1
        }
      }
    }                          
  } 
  
  # SNP data
  if (data.type=="SNP"){
    matrix <- as.matrix(x)
    l.names <- locNames(x)
    for (i in 1:nLoc(x)){
      row <- matrix[,i] # Row for each locus
      if (all(row == 0, na.rm=TRUE) | all(row == 2, na.rm=TRUE) | all(is.na(row))){
        loc.list[i] <- l.names[i]
        if (all(is.na(row))){
          na.counter = na.counter + 1
        }
      }
    }                          
  } 
  
  # Remove NAs from list of monomorphic loci and loci with all NAs
  loc.list <- loc.list[!is.na(loc.list)]
  
  # remove monomorphic loc and loci with all NAs
  if (verbose >= 2){
    cat("Removing monomorphic loci\n")
  } 
  x <- gl.drop.loc(x,loc.list=loc.list,verbose=0)
  
  # Report results
  if (verbose >= 3) {
    cat("  Original No. of loci:",nLoc(hold),"\n")
    cat("  Monomorphic loci:", nLoc(hold)-nLoc(x)-na.counter,"\n")
    cat("  Loci scored all NA:",na.counter,"\n")
    cat("  No. of loci deleted:",nLoc(hold)-nLoc(x),"\n")
    cat("  No. of loci retained:",nLoc(x),"\n")
    cat("  No. of individuals:",nInd(x),"\n")
    cat("  No. of populations:",nPop(x),"\n")
  }

# FLAG SCRIPT END

  # Add to history
  nh <- length(x@other$history)
  x@other$history[[nh + 1]] <- match.call()
  
  # Reset the flag
  x@other$loc.metrics.flags$monomorphs <- TRUE
  
  if (verbose >= 1){
    cat("Completed:",funname,"\n")
  }  
  
return (x)
  
}
