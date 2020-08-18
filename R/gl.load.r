#' Retrieves an object in compressed binary format earlier saved using gl.save.
#'
#' This is a wrapper for readRDS().
#'
#' The script retrieves the object from the current workspace and returns it by assignment.
#'
#' @param file -- name of the binary file from which to retrieve the object [required]
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2 or as specified using gl.set.verbosity]
#' @return the retrieved object
#' @export
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' @example
#' gl.save(testset.gl,"testset.rds")
#' gl.reloaded <- gl.load("testset.rds")

gl.load <- function(file, verbose=NULL){

# TRAP COMMAND, SET VERSION
  
  funname <- match.call()[[1]]
  build <- "Jacob"
  

# SET VERBOSITY
  
  if (is.null(verbose)){ 
    if(class(x)=="genlight" && !is.null(x@other$verbose)){ 
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
  
# DO THE JOB
  
  x <- readRDS(file)
  
# STANDARD ERROR CHECKING
  
  if(class(x)!="genlight") {
    if (all(x@ploidy == 1)){
      if (verbose >= 2){cat("  Loaded Presence/Absence (SilicoDArT) data\n")}
      data.type <- "SilicoDArT"
    } else if (all(x@ploidy == 2)){
      if (verbose >= 2){cat("  Loaded a SNP dataset\n")}
      data.type <- "SNP"
    }
  } else {
    if (verbose >= 2){cat("  Loaded an object of Class:",class(x),"\n")}
    data.type <- "other"
    }
  
# FLAG SCRIPT END

  if (verbose > 0) {
    cat("Completed:",funname,"\n")
  }
    
  return(x)
    
}
