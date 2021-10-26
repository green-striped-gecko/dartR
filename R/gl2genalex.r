#' Convert a genlight object to format suitable for input to genalex
#'
#' The output csv file contains the snp data and other relevant lines suitable
#'  for genalex. This script is a wrapper for  \link[poppr]{genind2genalex}
#'  (package poppr).
#'
#' @references
#' Peakall, R. and Smouse P.E. (2012) GenAlEx 6.5: genetic analysis
#' in Excel. Population genetic software for teaching and research-an update.
#' Bioinformatics 28, 2537-2539.
#' http://bioinformatics.oxfordjournals.org/content/28/19/2537
#'
#' @param x Name of the genlight object containing the SNP data [required].
#' @param outfile File name of the output file (including extension)
#' [default 'genalex.csv'].
#' @param outpath Path where to save the output file [default tempdir()].
#' @param overwrite If FALSE and filename exists, then the file will not be
#' overwritten. Set this option to TRUE to overwrite the file [default FALSE].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end;
#' 2, progress log; 3, progress and results summary; 5, full report
#' [default 2 or as specified using gl.set.verbosity].
#' @return NULL
#' @export
#' @author Custodian: Luis Mijangos, Author: Katrin Hohwieler, wrapper Arthur
#' Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' \donttest{
#' gl2genalex(testset.gl, outfile='testset.csv')
#' }

gl2genalex <- function(x,
                       outfile = "genalex.csv",
                       outpath = tempdir(),
                       overwrite = FALSE,
                       verbose = NULL) {
    outfilespec <- file.path(outpath, outfile)
    
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
    pkg <- "poppr"
    if (!(requireNamespace(pkg, quietly = TRUE))) {
        stop(error(
            "Package",
            pkg,
            " needed for this function to work. Please install it."
        ))
    }
    
    # DO THE JOB
    
    gind <- gl2gi(x, verbose = 0)
    poppr::genind2genalex(
        gind,
        filename = outfilespec,
        sequence = TRUE,
        overwrite = overwrite
    )
    
    if (verbose > 2) {
        cat(report(paste(
            "    Records written to", outfile, ":", nInd(x), "\n"
        )))
    }
    
    # FLAG SCRIPT END
    
    if (verbose > 0) {
        cat(report("Completed:", funname, "\n"))
    }
    
    return(NULL)
    
}
