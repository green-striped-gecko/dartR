#' Remove specified loci from a genelight \{adegenet\} object
#'
#' This script deletes selected loci from the nominated dataset.
#' 
#' @param x -- name of the genlight object containing SNP genotypes or a genind object containing presence/absence data [required]
#' @param loc.list -- vector of loci names to be droped.
#' @param v -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @return A genlight object with the reduced data
#' @export
#' @author Arthur Georges (bugs? Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#'    gl <- gl.drop.loc(testset.gl, loc.list=c("100049687|12-A/G","100050106|50-G/A"))
#' 

gl.drop.loc <- function(x, loc.list, v=2){

# ERROR CHECKING
  
  if(class(x)!="genlight") {
    cat("Fatal Error: genlight object required!\n"); stop("Execution terminated\n")
  }
  if (length(loc.list) == 0) {
    cat("Fatal Error: list of loci to drop required!\n"); stop("Execution terminated\n")
  }
  test <- loc.list%in%locNames(x)
  if (!all(test,na.rm=FALSE)) {
    cat("Fatal Error: some of the listed loci are not present in the dataset!\n"); stop("Execution terminated\n")
  }
  if (v < 0 | v > 5){
    cat("    Warning: verbosity must be an integer between 0 [silent] and 5 [full report], set to 2\n")
    v <- 2
  }
  
# FLAG SCRIPT START
  
  if (v >= 1) {
    cat("Starting gl.drop.loc: Deleting selected loci\n")
  }

# REMOVE LOCI
  
  if (v >= 2) {
    cat("  Deleting selected loci", loc.list, "\n")
  }

  # Delete listed loci
  
  # Remove rows flagged for deletion
    index <- !locNames(x)%in%loc.list
    x <- x[,index]
    x@other$loc.metrics <- x@other$loc.metrics[index,]

# REPORT A SUMMARY
    
  if (v >= 3) {
    cat("Summary of recoded dataset\n")
    cat(paste("  No. of loci:",nLoc(x),"\n"))
    cat(paste("  No. of individuals:", nInd(x),"\n"))
    cat(paste("  No. of populations: ", length(levels(factor(pop(x)))),"\n"))
  }

# FLAG SCRIPT END
    
  if (v >= 1) {
      cat("Completed gl.drop.loc\n\n")
  }
    
  return <- x
  
}

