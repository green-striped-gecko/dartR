#' @name gl.write.csv
#' @title Writes out data from a genlight object to csv file
#' @description
#' This script writes to file the SNP genotypes with specimens as entities
#' (columns) and loci as attributes (rows). Each row has associated locus
#'  metadata. Each column, with header of specimen id, has population in the
#'  first row.
#'
#' The data coding differs from the DArT 1row format in that 0 = reference
#' homozygous, 2 = alternate homozygous, 1 = heterozygous, and NA = missing SNP
#' assignment.
#' @param x Name of the genlight object containing the SNP data [required].
#' @param outfile File name of the output file (including extension)
#' [default "outfile.csv"].
#' @param outpath Path where to save the output file
#' [default tempdir(), mandated by CRAN]. Use outpath=getwd() or outpath='.'
#' when calling this function to direct output files to your working directory.
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end;
#' 2, progress log; 3, progress and results summary; 5, full report
#' [default 2 or as specified using gl.set.verbosity].
#' @return Saves a genlight object to csv, returns NULL.
#' @export
#' @author Custodian: Arthur Georges -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' # SNP data
#'   gl.write.csv(testset.gl, outfile='SNP_1row.csv')
#' # Tag P/A data
#'   gl.write.csv(testset.gs, outfile='PA_1row.csv')

gl.write.csv <- function(x,
                         outfile = "outfile.csv",
                         outpath = tempdir(),
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
    
    outfilespec <- file.path(outpath, outfile)
    
    # DO THE JOB
    
    # Add individual names and population names to rows
    x1 <- cbind.data.frame(pop(x), as.matrix(x))
    # Transpose to have id as columns, loci as rows
    x1 <- t(x1)
    # Create two filler rows to bring number of data rows and number of 
    #locus.metadata rows together
    filler1 <- rep("*", length(x@other$loc.metrics[1, ]))
    # Bind the filler rows to the locus metadata
    x2 <- rbind(filler1, as.matrix(x@other$loc.metrics))
    # Bind the locus metadata to the data
    x3 <- cbind.data.frame(x2, x1)
    
    # Output
    if (verbose >= 2) {
        cat(report("  Writing records to", outfilespec, "\n"))
    }
    
    write.table(x3,
                file = outfilespec,
                sep = ",",
                row.names = FALSE)
    
    # FLAG SCRIPT END
    
    if (verbose > 0) {
        cat(report("Completed:", funname, "\n"))
    }
    
    return(NULL)
}
