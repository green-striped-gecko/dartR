#' Define a new population in a genelight \{adegenet\} object on the basis of specified individuals 
#'
#' The script reassigns existing individuals to a new population and removes their existing population assignment
#' 
#' The script returns a genlight object with the new population assignment.
#'
#' @param x -- name of the genlight object containing SNP genotypes [required]
#' @param ind.list -- a list of individuals to be assigned to the new population [required]
#' @param new -- name of the new population
#' @param v -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @return A genlight object with the redefined population structure
#' @export
#' @author Arthur Georges (glbugs@@aerg.canberra.edu.au)
#' @examples
#'    gl <- gl.define.pop(testset.gl, ind.list=c("AA019073","AA004859"), new="newguys")

gl.define.pop <- function(x, ind.list, new, v=2){

# ERROR CHECKING
  
  if(class(x)!="genlight") {
    cat("Fatal Error: genlight object required!\n"); stop("Execution terminated\n")
  }
  for (case in ind.list){
    if (!(case%in%indNames(x))){
      cat("Warning: Listed individual",case,"not present in the dataset -- ignored\n")
      ind.list <- ind.list[!(ind.list==case)]
    }
  }
  if (length(ind.list) == 0) {
    cat("Fatal Error: no individuals listed to assign to population",new,"\n"); stop("Execution terminated\n")
  }
  if (v < 0 | v > 5){
    cat("    Warning: verbosity must be an integer between 0 [silent] and 5 [full report], set to 2\n")
    v <- 2
  }
  
# FLAG SCRIPT START
  
  if (v >= 1) {
    cat("Starting gl.define.pop: Assigning individuals to a new population\n")
  }
  
# ASSIGN INDIVIDUALS
  
  if (v >= 2) {
    cat("Processing",class(x),"object\n")
    cat("  Assigned listed individuals", ind.list,"to new population",new, "\n")
  }

  tmp <- as.character(pop(x))
  for (case in ind.list){
    tmp[indNames(x) == case] <- new
  }
  pop(x) <- as.factor(tmp)
  
# REPORT A SUMMARY
    
  if (v >= 3) {
    cat("Summary of recoded dataset\n")
    cat(paste("  No. of loci:",nLoc(x),"\n"))
    cat(paste("  No. of individuals:", nInd(x),"\n"))
    cat(paste("  No. of populations: ", length(levels(factor(pop(x)))),"\n"))
  }

# FLAG SCRIPT END
    
    if (v >= 1) {
      cat("Completed gl.define.pop\n\n")
    }
    
    return <- x
}

