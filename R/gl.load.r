#' Load an object from compressed binary format produced by gl.save().
#'
#' This is a wrapper for readRDS().
#'
#' The script loads the object from the current workspace and returns the gl object.
#'
#' @param file -- name of the file to receive the binary version of the object [required]
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2 or as specified using gl.set.verbosity]
#' @return the loaded object
#' @export
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' \dontrun{
#' gl <- gl.load("testset.rds")
#' }
#' @seealso \code{\link{gl.save}}
#' 
gl.load <- function(file, verbose=NULL){

  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(f=funname,build="Jackson",v=verbose)
  
  x <- readRDS(file)
  
# CHECK DATATYPE 
  datatype <- utils.check.datatype(x,verbose=verbose)
  cat(report("  Loaded object of type",datatype,"from",file,"\n\n"))
  
# FLAG SCRIPT END

  if (verbose > 0) {
    cat(report("Completed:",funname,"\n"))
  }
    
  invisible(x)
    
}
