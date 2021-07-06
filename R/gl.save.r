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

gl.save <- function(x, file, verbose=NULL){

  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func=funname,build="Jackson",v=verbose)
  
  # CHECK DATATYPE 
  datatype <- utils.check.datatype(x,verbose=0)
  
# DO THE JOB

  saveRDS(x,file)
  cat(report("  Saved object of type",datatype,"to",file,"\n\n"))
  
# FLAG SCRIPT END

  if (verbose > 0) {
    cat(report("Completed:",funname,"\n"))
  }
    
  invisible(x)
    
}
