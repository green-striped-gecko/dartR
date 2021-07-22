#' A utility script to flag the start of a script
#'
#' @param func name of the function that is starting [required].
#' @param build name of the build [NULL].
#' @param v verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @return calling function name
#' @author Arthur Georges -- Post to \url{https://groups.google.com/d/forum/dartr}

utils.flag.start <- function(func=NULL, build=NULL, v=2) {

  
 if(is.null(func)){
   stop(error("Fatal Error: The calling function must be specified.\n"))
 }
  if (v >= 1) {
    if (v == 5){
      if (!is.null(build)){
        cat(report("\n\nStarting", func, "\n[dartR vers.",packageVersion("dartR"),"Build =",build,"]\n"))
      } else {
        cat(report("\n\nStarting", func, "\n"))
      }  
    } else {
      cat(report("\n\nStarting", func, "\n"))

    }
  }
  
  invisible(func)
  
}