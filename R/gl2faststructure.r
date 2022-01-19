#' Exports DArT genlight object \{adegenet\} to faststructure format (to run
#'  faststructure elsewhere)
#'
#' Recodes in the quite specific faststructure format (e.g first six columns
#' need to be there, but are ignored...check faststructure documentation
#'  (if you find any :-( )))
#'
#' The script writes out the a file in faststructure format.
#' @param x Name of the genlight object containing the SNP data [required].
#' @param outfile File name of the output file (including extension)
#' [default "gl.str"].
#' @param outpath Path where to save the output file
#' [default tempdir(), mandated by CRAN]. Use outpath=getwd() or outpath='.'
#' when calling this function to direct output files to your working directory.
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2 or as specified using gl.set.verbosity].
#' @param probar Switch to show/hide progress bar [default FALSE].
#' @return NULL
#' @export
#' @author Bernd Gruber (Post to \url{https://groups.google.com/d/forum/dartr})
#' @importFrom utils getTxtProgressBar setTxtProgressBar txtProgressBar

gl2faststructure <- function(x,
                             outfile = "gl.str",
                             outpath = tempdir(),
                             probar = FALSE,
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
    
    # DO THE JOB
    
    x <- as.matrix(x)
    # add six dummy colums
    nc <- ncol(x) + 6
    if (probar)
        pb <- txtProgressBar(min = 0,
                             1,
                             style = 3,
                             initial = NA)
    zz <- file(outfilespec, "w")
    for (i in 1:nrow(x)) {
        dummy <- rbind(x[i, ], x[i, ])
        index <- colSums(dummy, na.rm = T) == 2
        dummy[, index] <- c(0, 2)
        dummy <- ifelse(is.na(dummy), -9, dummy)
        dummy <- ifelse(dummy == 0, 1, dummy)
        dummy <- cbind(i, i, i, i, i, i, dummy)
        write(
            t(dummy),
            file = outfilespec,
            sep = "\t",
            ncolumns = nc,
            append = TRUE
        )
        if (probar)
            setTxtProgressBar(pb, i / nrow(x))
    }
    close(zz)
    if (probar)
        close(pb)
   
    if (verbose >= 2) {
        cat(report(paste0("Saved faststructure file: ", outfilespec, "\n")))
    }
    if (verbose >= 3) {
        cat(report(paste(
            "Consists of",
            nrow(x),
            "individuals and ",
            ncol(x),
            "loci."
        )))
    }
    
    # FLAG SCRIPT END
    
    if (verbose > 0) {
        cat(report("Completed:", funname, "\n"))
    }
    
    return(NULL)
}
