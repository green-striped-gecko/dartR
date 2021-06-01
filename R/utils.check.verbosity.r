#' @name utils.check.verbosity
#' 
#' @title A utility script to set the  verbosity.
#'
#' @description 
#' The verbosity can be set in one of two ways -- (a) explicitly by the user by passing a value using the parameter
#' verbose= in a function, or (b) by setting the verbosity globally as part of the r environment.
#' 
#' @param verbosity -- user requested level of verbosity
#' @return The verbosity, in variable verbose
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' verbose <- utils.check.verbosity(NULL)

utils.check.verbosity <- function(verbose=verbose) {
  
# SET VERBOSITY
  hold <- verbose
  if (is.null(verbose)){ 
    verbose <- Sys.getenv("verbosity",unset=2)
  }

  if (verbose < 0 | verbose > 5){
    cat(warn("  Warning: Parameter 'verbose' must be an integer between 0 [silent] and 5 [full report], set to 2\n"))
    verbose <- 2
  }
  if(verbose>=2){
    if(is.null(hold)){
      cat(report("  Verbosity not specified, using global default\n"))
    }
  }

  return(verbose)
  
}
