#' A utility script to flag the start of a script
#'
#' @param f name of the function that is starting [obtained from match.call].
#' @param build name of the build [NULL].
#' @param v verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @return function name
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' #out <- utils.function.start()

utils.flag.start <- function(func=funname, build=NULL, v=verbose) {
  
  if (v >= 1) {
    if (v == 5){
      if (!is.null(build)){
        cat(report("\n\nStarting", func, "\n[dartR vers.",packageVersion("dartR"),"Build =",build,"]\n\n"))
      } else {
        cat(report("\n\nStarting", func, "\n\n"))
      }  
    } else {
      cat(report("\n\nStarting", func, "\n\n"))
    }
  }
  
  invisible(func)
  
}