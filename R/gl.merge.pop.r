#' Merge two or more populations in a genelight \{adegenet\} object into one population
#'
#' Individuals are assigned to populations based on the specimen metadata data file (csv) used with gl.read.dart(). 
#'
#' This script assigns individuals from two nominated populations into a new single population. It can also be used
#' to rename populations.
#' 
#' The script returns a genlight object with the new population assignments.
#'
#' @param x -- name of the genlight object containing SNP genotypes or a genind object containing presence/absence data [required]
#' @param old -- a list of populations to be merged [required]
#' @param new -- name of the new population [required]
#' @param v -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @return A genlight object with the new population assignments
#' @export
#' @author Arthur Georges (glbugs@@aerg.canberra.edu.au)
#' @examples
#' \dontrun{
#'    gl <- gl.merge.pop(testset.gl, old=c("EmsubRopeMata","EmvicVictJasp"), new="Outgroup")
#' }

gl.merge.pop <- function(x, old=Null, new=Null, v=2) {

  if (v > 0) {
    if (length(old) == 1) {
      cat("Starting gl.merge.pop: Renaming a population\n")
    } else if (length(old) > 1) {
      cat("Starting gl.merge.pop: Merging a list of populations into one\n")
    } else {
      cat("Fatal Error: At least one old population label must be provided\n"); stop("Execution terminated\n")
    }
  }
  if (is.null(new)) {
    cat("Fatal Error: A new population label must be specified\n"); stop("Execution terminated\n")
  }
  if(class(x)!="genlight") {
    cat("Fatal Error: genlight object required for gl.keep.pop.r!\n"); stop("Execution terminated\n")
  }
  if (v > 1) {
    if (length(old) == 1) {
      cat("  Renaming",old,"as",new,"\n")
    } else {
      cat("  Merging",old,"into",new,"\n")
    } 
  }

# Merge or rename
  for (i in 1:length(old)) {
    #pop(x)[pop(x) == old[i]] <- new
    levels(pop(x))[levels(pop(x))==old[i]] <- new
  }

    if (v > 0) {
      cat("Completed gl.merge.pop\n\n")
    }
    
    return(x)
}

