#' @title gl.set.verbosity
#' @title Set the default verbosity level
#' @description
#' dartR functions have a verbosity parameter that sets the level of reporting during the execution of the function. The verbosity level,
#' set by parameter 'verbose' can be one of verbose 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary;
#' 5, full report. The default value for verbosity is stored in the r environment. This script sets the default value.
#' 
#' @param value -- set the default verbosity to be this value: 0, silent only fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @return verbosity value [set for all functions]
#' @export
#' @author Custodian: Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' gl <- gl.set.verbosity(value=2)

gl.set.verbosity <- function(value=2) {

  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func=funname,build="Jackson",v=value)
  
# SET GLOBAL VERBOSITY
  if (!is.null(value) & is.numeric(value) & value>=0 & value<=5)  {
    options("dartR_verbose"=value)
  }
  if(value>=2){
    cat(report("  Global verbosity set to:",value,"\n"))
  }
  
# FLAG SCRIPT END
  
  if (value >= 1) {
    cat("Completed:",funname,"\n")
  }
  
  invisible(NULL)
}
  
