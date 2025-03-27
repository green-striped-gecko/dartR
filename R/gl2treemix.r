#' Converts a genlight object to a treemix input file
#'
#' The output file contains the SNP data in the format expected by treemix --
#' see the treemix manual. The file will be gzipped before in order to be
#' recognised by treemix. Plotting functions provided with treemix will need to
#'  be sourced from the treemix download page.
#'
#' @references Pickrell and Pritchard (2012). Inference of population splits and
#' mixtures from genome-wide allele frequency data. PLoS Genetics
#'  https://doi.org/10.1371/journal.pgen.1002967
#'
#' @param x Name of the genlight object containing the SNP data [required].
#' @param outfile File name of the output file (including gz extension)
#' [default 'treemix_input.gz'].
#' @param outpath Path where to save the output file
#'  [default tempdir(), mandated by CRAN]. Use outpath=getwd() when calling this
#'   function or set.tempdir <- getwd() elsewhere in your script to direct
#'   output files to your working directory.
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2 or as specified using gl.set.verbosity].
#' @return  returns no value (i.e. NULL)
#' @export
#' @author Custodian: Arthur Georges (Post to
#' \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' gl2treemix(testset.gl, outpath=tempdir())

gl2treemix <- function(x,
                       outfile = "treemix_input.gz",
                       outpath = tempdir(),
                       verbose = NULL) {
    outfilespec <- file.path(outpath, outfile)
    
    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)
    
    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     build = "Jody",
                     verbose = verbose)
    
    # CHECK DATATYPE
    datatype <- utils.check.datatype(x, verbose = verbose)
    
    # DO THE JOB
    
    freq <- gl.percent.freq(x, verbose = verbose)
    freq$ref <- freq$nobs * 2 - freq$sum
    freq$alt <- freq$sum
    freq$sum <- NULL
    freq$nobs <- NULL
    freq$nmissing <- NULL
    freq$frequency <- NULL
    freq$n <- NULL
    
    # Output the file
    outfile <- file.path(outpath, outfile)
    if (verbose >= 2) {
        cat(report(
            paste(
                "    Writing results to treemix input file",
                outfilespec,
                "\n"
            )
        ))
    }
    sink(gzfile(outfilespec))
    
    cat(unique(as.character(freq$popn)), "\n")
    k <- 1
    for (j in 1:nLoc(x)) {
        for (i in k:(k + nPop(x) - 1)) {
            cat(paste0(freq$ref[i], ",", freq$alt[i]), " ")
        }
        cat("\n")
        k <- k + nPop(x)
    }
    
    sink()
    if (verbose > 2) {
        cat(report(
            paste("    Records written to", outfilespec, ":", nInd(x), "\n")
        ))
    }
    
    # FLAG SCRIPT END
    
    if (verbose > 0) {
        cat(report("Completed:", funname, "\n"))
        cat(important("Output file has been gzipped for input to treemix\n"))
    }
    
    return(NULL)
    
}
