#' Calculates basic statistics for each loci (Hs, Ho, Fis etc.)
#'
#' Based on function \code{\link[hierfstat]{basic.stats}}. Check ?basic.stats
#' for help.
#' @param x Name of the genlight object containing the SNP data [required].
#' @param digits Number of digits that should be returned [default 4].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#' @return Several tables and lists with all basic stats.
#' \code{\link[hierfstat]{basic.stats}} for details.
#' @author Bernd Gruber (bugs? Post to
#' \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' out <- gl.basic.stats(possums.gl[1:10,1:100])
#' @export

gl.basic.stats <- function(x,
                           digits = 4,
                           verbose = NULL) {
    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)
    
    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     build = "Jody",
                     verbosity = verbose)
    
    # CHECK DATATYPE
    datatype <- utils.check.datatype(x, verbose = verbose)
    
    # CHECK IF PACKAGES ARE INSTALLED
    pkg <- "hierfstat"
    if (!(requireNamespace(pkg, quietly = TRUE))) {
      cat(error(
        "Package",
        pkg,
        " needed for this function to work. Please install it.\n"
      ))
      return(-1)
    }
    
    # DO THE JOB
    
    out <-
        hierfstat::basic.stats(hierfstat::genind2hierfstat(gl2gi(x,
                                                                 verbose = 0)), 
                               digits = digits)
    
    # FLAG SCRIPT END
    
    if (verbose >= 1) {
        cat(report("Completed:", funname, "\n"))
    }
    
    # RETURN
    
    return(out)
}
