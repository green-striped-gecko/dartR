#' A utility script to flag the start of a script
#'
#' @param func Name of the function that is starting [required].
#' @param build Name of the build [default NULL].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#'  progress log; 3, progress and results summary; 5, full report [default 2].
#' @return calling function name
#' @author Custodian: Arthur Georges -- Post to
#'  \url{https://groups.google.com/d/forum/dartr}
#'  @export

utils.flag.start <- function(func = NULL,
                             build = NULL,
                             verbose = NULL) {
    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)
    
    if (is.null(func)) {
        stop(error("Fatal Error: The calling function must be specified.\n"))
    }
    if (verbose >= 1) {
        if (verbose == 5) {
            if (!is.null(build)) {
                cat(
                    report(
                        "Starting",
                        func,
                        "\n[dartR vers.",
                        packageVersion("dartR"),
                        "Build =",
                        build,
                        "]\n"
                    )
                )
            } else {
                cat(report("Starting", func, "\n"))
            }
        } else {
            cat(report("Starting", func, "\n"))
        }
    }
    invisible(func)
}
