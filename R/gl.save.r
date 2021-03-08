#' Save an object in compressed binary format for later rapid retrieval.
#'
#' This is a wrapper for saveRDS().
#'
#' The script saves the object to the current workspace and returns the input gl object.
#'
#' @param x -- name of the genlight object containing SNP genotypes [required]
#' @param file -- name of the file to receive the binary version of the object [required]
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2 or as specified using gl.set.verbosity]
#' @return the input object
#' @export
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' gl.save(testset.gl,file.path(tempdir(),"testset.rds"))
#' gl.reloaded <- gl.load(file.path(tempdir(),"testset.rds"))

gl.save <- function(x, file, verbose=NULL){

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
  
# STANDARD ERROR CHECKING
  
  if(class(x)!="genlight") {
    if (all(x@ploidy == 1)){
      if (verbose >= 2){cat("  Saving Presence/Absence (SilicoDArT) data\n")}
      data.type <- "SilicoDArT"
    } else if (all(x@ploidy == 2)){
      if (verbose >= 2){cat("  Saving a SNP dataset\n")}
      data.type <- "SNP"
    }
  } else {
    if (verbose >= 2){cat("  Saving an object of Class:",class(x),"\n")}
    data.type <- "other"
    }
  
# DO THE JOB

  saveRDS(x,file)
  
# ADD TO HISTORY
    nh <- length(x@other$history)
    x@other$history[[nh + 1]] <- match.call()  

# FLAG SCRIPT END

  if (verbose > 0) {
    cat("Completed:",funname,"\n")
  }
    
  return(invisible(2))
    
}
