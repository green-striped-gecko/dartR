#' Converts a genlight object to nexus format suitable for phylogenetic analysis
#'  by SNAPP (via BEAUti)
#'
#' The output nexus file contains the SNP data and relevant PAUP command lines
#' suitable for BEAUti.
#'
#' @references Bryant, D., Bouckaert, R., Felsenstein, J., Rosenberg, N.A. and
#' RoyChoudhury, A. (2012). Inferring species trees directly from biallelic
#' genetic markers: bypassing gene trees in a full coalescent analysis.
#'  Molecular Biology and Evolution 29:1917-1932.
#'
#' @param x Name of the genlight object containing the SNP data [required].
#' @param outfile File name of the output file (including extension)
#' [default "snapp.nex"].
#' @param outpath Path where to save the output file
#'  [default tempdir(), mandated by CRAN]. Use outpath=getwd() or outpath='.'
#'  when calling this function to direct output files to your working directory.
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2 or as specified using gl.set.verbosity].
#' @return  returns no value (i.e. NULL)
#' @export
#' @author Custodian: Arthur Georges (Post to
#' \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' gl2snapp(testset.gl)

gl2snapp <- function(x,
                     outfile = "snapp.nex",
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
    
    if (verbose >= 2) {
        cat(paste(
            report(
                "  Extacting SNP data and creating records for each individual\n"
            )
        ))
    }
    
    # Extract the reference base and the alternate base for each locus (excuse the contortion)
    m <- as.matrix(x)
    m[is.na(m)] <- "?"
    colnames(m) <- NULL
    df <- data.frame(m)
    df <- cbind(indNames(x), pop(x), df)
    indlabels <- df[, 1]
    poplabels <- df[, 2]
    
    # Create the snapp file
    if (verbose > 1) {
        cat(report(
            paste("  Writing results to nexus file", outfilespec, "\n")
        ))
    }
    
    sink(outfilespec)
    
    cat("#nexus\n")
    cat("BEGIN DATA;\n")
    cat(paste0("     dimensions ntax = ", nInd(x), " nchar = ", nLoc(x), " ;\n"))
    cat("     format datatype=integerdata missing=? symbols=\"012\";\n")
    cat("matrix\n")
    for (i in 1:nInd(x)) {
        cat(paste0(poplabels[i], "_", indlabels[i]))
        cat("  ")
        cat(m[i, ], sep = "")
        cat("\n")
    }
    cat(";\n")
    cat("end;\n\n")
    
    sink()
    
    if (verbose > 2) {
        cat(paste(
            report("    Records written to", outfilespec, ":", nInd(x), "\n")
        ))
    }
    
    # FLAG SCRIPT END
    
    if (verbose > 0) {
        cat(report("Completed:", funname, "\n"))
    }
    
    return(NULL)
    
}
