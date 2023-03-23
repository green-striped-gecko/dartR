#' Checks the current global verbosity
#'
#' The verbosity can be set in one of two ways -- (a) explicitly by the user by
#' passing a value using the parameter verbose in a function, or (b) by setting
#' the verbosity globally as part of the r environment (gl.set.verbosity).
#'
#' @param x User requested level of verbosity [default NULL].
#' @return The verbosity, in variable verbose
#' @examples 
#' gl.check.verbosity()
#' @export
#' @author Bernd Gruber (Post to \url{https://groups.google.com/d/forum/dartr})

gl.check.verbosity <- function(x = NULL) {
    # SET VERBOSITY or GET it from global
    if (is.null(x)) {
        if (is.null(options()$dartR_verbose)) {
            verbose <- 2
        } else {
            verbose <- options()$dartR_verbose
        }
    } else {
        if (is.numeric(x) & x >= 0 & x <= 5) {
            verbose <- x
        } else {
            cat(
                warn(
                    "Warning: Parameter verbose must be an integer in the range 
                    0 to 5, set to 2\n"
                )
            )
            verbose <- 2
        }
    }
    
    return(verbose)
    
}
