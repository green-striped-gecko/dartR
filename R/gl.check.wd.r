#' @name gl.check.wd
#' @title Checks the global working directory
#' @family environment
#' 
#' @description 
#' The working directory can be set in one of two ways -- (a) explicitly by the user by
#' passing a value using the parameter plot.dir in a function, or (b) by setting
#' the working directory globally as part of the r environment (gl.setwd). The default is in acccordance to CRAN set to tempdir().

#' @param wd path to the working directory [default: tempdir()].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].

#' @examples 
#' gl.check.wd()
#' 
#' @author Custodian: Bernd Gruber (Post to
#' \url{https://groups.google.com/d/forum/dartr})
#' 
#' @export
#' @return the working directory

gl.check.wd <- function(
    wd = NULL,
    verbose=NULL) {
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname,
                   build = "v.2023.2",
                   verbosity = verbose)
  # DO THE JOB
  # SET wd or GET it from global
  if (is.null(wd)) {
    # If wd is not provided, check if it's set in options()
    if (is.null(options()$dartR_wd)) {
      # If not set in options(), set wd to tempdir()
      wd <- tempdir()
    } else {
      # If set in options(), use that value for wd
      wd <- options()$dartR_wd
    }
  } else {
    # If wd is provided
    if (is.character(wd) & dir.exists(wd)) {
      # Check if it's a valid directory path
      wd <- wd
    } else {
      # If not a valid directory path, display a warning message and set wd to tempdir()
      cat(
        warn(
          "Warning: The path to the working directory does not exist! Set to tempdir().\n"
        )
      )
      wd <- tempdir()
    }
  }
  if(verbose >= 2){cat(report("  Working directory:",wd,"\n"))}
  # FLAG SCRIPT END
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }
    
    return(wd)
    
}
