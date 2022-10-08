#' @name gl.report.monomorphs
#' @title Reports monomorphic loci
#' @description
#' This script reports the number of monomorphic loci and those with all NAs in
#' a genlight \{adegenet\} object
#' @param x Name of the input genlight object [required].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default NULL, unless specified using gl.set.verbosity].
#' @details
#' A DArT dataset will not have monomorphic loci, but they can arise, along with
#' loci that are scored all NA, when populations or individuals are deleted.
#' Retaining monomorphic loci unnecessarily increases the size of the dataset
#' and will affect some calculations.
#'
#' Note that for SNP data, NAs likely represent null alleles; in tag
#' presence/absence data, NAs represent missing values (presence/absence could
#' not be reliably scored)
#' @return An unaltered genlight object
#' @rawNamespace import(adegenet, except = plot)
#' @author Custodian: Arthur Georges -- Post to
#'  \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' # SNP data
#'   gl.report.monomorphs(testset.gl)
#' # SilicoDArT data
#'   gl.report.monomorphs(testset.gs)
#' @seealso \code{\link{gl.filter.monomorphs}}
#' @family report functions
#' @export

gl.report.monomorphs <- function(x,
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
    
    # DO THE JOB
    
    hold <- x
    na.counter <- 0
    loc.list <- array(NA, nLoc(x))
    
    if (verbose >= 2) {
        cat(report("  Identifying monomorphic loci\n"))
    }
    # # Tag presence/absence data
    # if (datatype == "SilicoDArT") {
    #     mat <- as.matrix(x)
    #     lN <- locNames(x)
    #     for (i in 1:nLoc(x)) {
    #         row <- mat[, i]  # Row for each locus
    #         if (all(row == 0, na.rm = TRUE) |
    #             all(row == 1, na.rm = TRUE) | all(is.na(row))) {
    #             loc.list[i] <- lN[i]
    #             if (all(is.na(row))) {
    #                 na.counter <-na.counter + 1
    #             }
    #         }
    #     }
    # }
    # 
    # # SNP data
    # if (datatype == "SNP") {
    #     mat <- as.matrix(x)
    #     lN <- locNames(x)
    #     for (i in 1:nLoc(x)) {
    #         row <- mat[, i]  # Row for each locus
    #         if (all(row == 0, na.rm = TRUE) |
    #             all(row == 2, na.rm = TRUE) | all(is.na(row))) {
    #             loc.list[i] <- lN[i]
    #             if (all(is.na(row))) {
    #                 na.counter <-na.counter + 1
    #             }
    #         }
    #     }
    # }
    
    mono_tmp <- gl.alf(x)
    loc.list <- rownames(mono_tmp[which(mono_tmp$alf1==1 | 
                                          mono_tmp$alf1 == 0),])
    loc.list_NA <- rownames(mono_tmp[which(is.na(mono_tmp$alf1)),])
    
    # Remove NAs from list of monomorphic loci and loci with all NAs
    loc.list <- loc.list[!is.na(loc.list)]
    
    # remove monomorphic loci and loci with all NAs
    if (length(loc.list) > 0) {
        x <- gl.drop.loc(x, loc.list = loc.list, verbose = 0)
    }
    
    # PRINTING OUTPUTS Report results
    cat("\n  No. of loci:", nLoc(hold), "\n")
    cat("    Polymorphic loci:", nLoc(x), "\n")
    cat("    Monomorphic loci:", nLoc(hold) - nLoc(x), "\n")
    cat("    Loci scored all NA:", length(loc.list_NA), "\n")
    cat("  No. of individuals:", nInd(x), "\n")
    cat("  No. of populations:", nPop(x), "\n\n")
    
    # FLAG SCRIPT END
    
    if (verbose >= 1) {
        cat(report("Completed:", funname, "\n"))
    }
    
    # RETURN
    
    invisible(x)
    
}
