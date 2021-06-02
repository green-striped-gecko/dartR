#' A function to check the current global verbosity 
#'
#' The verbosity can be set in one of two ways -- (a) explicitly by the user by passing a value using the parameter
#' verbose= in a function, or (b) by setting the verbosity globally as part of the r environment (gl.set.verbosity).
#' 
#' @param verbose -- user requested level of verbosity
#' @return The verbosity, in variable verbose
#' @export
#' @author Bernd Gruber (Post to \url{https://groups.google.com/d/forum/dartr})


gl.check.verbosity <- function(verbose=NULL) {
  
# SET VERBOSITY or GET it from global
  if (is.null(verbose)) verbose <- options()$dartR_verbose else {
  if (is.numeric(verbose) & verbose>=0 & verbose<=5) verbose <- verbose else verbose <- 2 
  }
    
  return(verbose)

}
