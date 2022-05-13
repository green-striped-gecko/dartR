#' A utility script to reset to FALSE (or TRUE) the locus metric flags after
#' some individuals or populations have been deleted.
#'
#' The locus metadata supplied by DArT has OneRatioRef, OneRatioSnp, PICRef,
#' PICSnp, and AvgPIC included, but the allelic composition will change when
#' some individuals are removed from the dataset and so the initial statistics
#' will no longer apply. This applies also to some variable calculated by dartR
#' (e.g. maf). This script resets the locus metrics flags to FALSE to indicate
#' that these statistics in the genlight object are no longer current. The
#' verbosity default is also set, and in the case of SilcoDArT, the flags PIC
#' and OneRatio are also set.
#'
#' If the locus metrics do not exist then they are added to the genlight object
#'  but not populated. If the locus metrics flags do not exist, then they are
#'  added to the genlight object and set to FALSE (or TRUE).
#'
#' @param x Name of the genlight object containing the SNP data or
#' tag presence/absence data (SilicoDArT) [required].
#' @param set Set the flags to TRUE or FALSE [default FALSE].
#' @param value Set the default verbosity for all functions, where verbosity is
#'  not specified [default 2].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report [default NULL].
#' @return The modified genlight object
#' @author Custodian: Luis Mijangos (Post to
#' \url{https://groups.google.com/d/forum/dartr})
#' @seealso \code{utils.recalc.metrics} for recalculating all metrics,
#' \code{utils.recalc.callrate} for recalculating CallRate,
#' \code{utils.recalc.freqhomref} for recalculating frequency of homozygous
#' reference, \code{utils.recalc.freqhomsnp} for recalculating frequency of
#' homozygous alternate, \code{utils.recalc.freqhet} for recalculating frequency
#' of heterozygotes, \code{gl.recalc.maf} for recalculating minor allele 
#' frequency, \code{gl.recalc.rdepth} for recalculating average read depth
#' @examples
#' #result <- utils.reset.flags(testset.gl)

utils.reset.flags <- function(x,
                              set = FALSE,
                              value = 2,
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
    
    # SCRIPT SPECIFIC ERROR TESTING
    
    if (value < 0 | value > 5) {
        cat(
            warn(
                "  Warning: Parameter 'value' must be an integer between 0 [silent] and 5 [full report], set to 2\n"
            )
        )
        value <- 2
    }
    
    # DO THE JOB
    if (datatype == "SNP") {
        if (verbose >= 2) {
            cat(
                report(
                    "  Resetting flags for AvgPIC, OneRatioRef, OneRatioSnp, PICRef, PICSnp, CallRate, maf, FreqHets, FreqHomRef, FreqHomSnp, monomorphs, OneRatio, PIC, allna to",
                    set,
                    "\n"
                )
            )
            cat(report(
                "  Resetting flags for OneRatio, PIC to FALSE\n"
            ))
        }
        
        # Check if the x@other$loc.metrics slot exists, if not, create as a dataframe
        if (is.null(x@other$loc.metrics)) {
            x@other$loc.metrics <- as.data.frame(array(NA, nLoc(x)))
        }
        # Check if the x@other$loc.metrics.flags slot exists, if not, create as a dataframe
        if (is.null(x@other$loc.metrics.flags)) {
            x@other$loc.metrics.flags <- as.data.frame(array(NA, 1))
        }
        
        # loc.metric should be a dataframe
        x@other$loc.metrics <- as.data.frame(x@other$loc.metrics)
        
        # AvgPIC
        if (is.null(x@other$loc.metrics$AvgPIC)) {
            x@other$loc.metrics$AvgPIC <- array(NA, nLoc(x))
            if (verbose >= 3) {
                cat(
                    report(
                        "  Locus metric AvgPIC does not exist, creating slot @other$loc.metrics$AvgPIC\n"
                    )
                )
            }
        }
        x@other$loc.metrics.flags$AvgPIC <- set
        
        # OneRatioRef
        if (is.null(x@other$loc.metrics$OneRatioRef)) {
            x@other$loc.metrics$OneRatioRef <- array(NA, nLoc(x))
            if (verbose >= 3) {
                cat(
                    report(
                        "  Locus metric OneRatioRef does not exist, creating slot @other$loc.metrics$OneRatioRef\n"
                    )
                )
            }
        }
        x@other$loc.metrics.flags$OneRatioRef <- set
        
        # OneRatioSnp
        if (is.null(x@other$loc.metrics$OneRatioSnp)) {
            x@other$loc.metrics$OneRatioSnp <- array(NA, nLoc(x))
            if (verbose >= 3) {
                cat(
                    report(
                        "  Locus metric OneRatioSnp does not exist, creating slot @other$loc.metrics$OneRatioSnp\n"
                    )
                )
            }
        }
        x@other$loc.metrics.flags$OneRatioSnp <- set
        
        # PICRef
        if (is.null(x@other$loc.metrics$PICRef)) {
            x@other$loc.metrics$PICRef <- array(NA, nLoc(x))
            if (verbose >= 3) {
                cat(
                    report(
                        "  Locus metric PICRef does not exist, creating slot @other$loc.metrics$PICRef\n"
                    )
                )
            }
        }
        x@other$loc.metrics.flags$PICRef <- set
        
        # PICSnp
        if (is.null(x@other$loc.metrics$PICSnp)) {
            x@other$loc.metrics$PICSnp <- array(NA, nLoc(x))
            if (verbose >= 3) {
                cat(
                    report(
                        "  Locus metric PICSnp does not exist, creating slot @other$loc.metrics$PICSnp\n"
                    )
                )
            }
        }
        x@other$loc.metrics.flags$PICSnp <- set
        
        # CallRate
        if (is.null(x@other$loc.metrics$CallRate)) {
            x@other$loc.metrics$CallRate <- array(NA, nLoc(x))
            if (verbose >= 3) {
                cat(
                    report(
                        "  Locus metric CallRate does not exist, creating slot @other$loc.metrics$CallRate\n"
                    )
                )
            }
        }
        x@other$loc.metrics.flags$CallRate <- set
        
        # FreqHomRef
        if (is.null(x@other$loc.metrics$FreqHomRef)) {
            x@other$loc.metrics$FreqHomRef <- array(NA, nLoc(x))
            if (verbose >= 3) {
                cat(
                    report(
                        "  Locus metric FreqHomRef does not exist, creating slot @other$loc.metrics$FreqHomRef\n"
                    )
                )
            }
        }
        x@other$loc.metrics.flags$FreqHomRef <- set
        
        # FreqHomSnp
        if (is.null(x@other$loc.metrics$FreqHomSnp)) {
            x@other$loc.metrics$FreqHomSnp <- array(NA, nLoc(x))
            if (verbose >= 3) {
                cat(
                    report(
                        "  Locus metric FreqHomSnp does not exist, creating slot @other$loc.metrics$FreqHomSnp\n"
                    )
                )
            }
        }
        x@other$loc.metrics.flags$FreqHomSnp <- set
        
        # FreqHets
        if (is.null(x@other$loc.metrics$FreqHets)) {
            x@other$loc.metrics$FreqHets <- array(NA, nLoc(x))
            if (verbose >= 3) {
                cat(
                    report(
                        "  Locus metric FreqHets does not exist, creating slot @other$loc.metrics$FreqHets\n"
                    )
                )
            }
        }
        x@other$loc.metrics.flags$FreqHets <- set
        
        # monomorphs
        if (is.null(x@other$loc.metrics$monomorphs)) {
            x@other$loc.metrics$monomorphs <- array(NA, nLoc(x))
            if (verbose >= 3) {
                cat(
                    report(
                        "  Locus metric monomorphs does not exist, creating slot @other$loc.metrics$monomorphs\n"
                    )
                )
            }
        }
        x@other$loc.metrics.flags$monomorphs <- set
        
        # maf
        if (is.null(x@other$loc.metrics$maf)) {
            x@other$loc.metrics$maf <- array(NA, nLoc(x))
            if (verbose >= 3) {
                cat(
                    report(
                        "  Locus metric maf does not exist, creating slot @other$loc.metrics$maf\n"
                    )
                )
            }
        }
        x@other$loc.metrics.flags$maf <- set
        
        # OneRatio
        if (is.null(x@other$loc.metrics$OneRatio)) {
            x@other$loc.metrics$OneRatio <- array(NA, nLoc(x))
            if (verbose >= 3) {
                cat(
                    report(
                        "  Locus metric OneRatio does not exist, creating slot @other$loc.metrics$OneRatio\n"
                    )
                )
            }
        }
        x@other$loc.metrics.flags$OneRatio <- FALSE
        
        # PIC
        if (is.null(x@other$loc.metrics$PIC)) {
            x@other$loc.metrics$PIC <- array(NA, nLoc(x))
            if (verbose >= 3) {
                cat(
                    report(
                        "  Locus metric PIC does not exist, creating slot @other$loc.metrics$PIC\n"
                    )
                )
            }
        }
        x@other$loc.metrics.flags$PIC <- FALSE
        
        # monomorphs if (is.null(x@other$loc.metrics$monomorphs)) { x@other$loc.metrics$monomorphs <- array(NA,nLoc(x)) if (verbose >= 3){
        # cat(' Locus metric monomorphs does not exist, creating slot @other$loc.metrics$monomorphs\n') } }
        x@other$loc.metrics.flags$monomorphs <- set
        
        # allna
        x@other$loc.metrics.flags$allna <- set
        
        
        # verbosity
        if (is.null(x@other$verbose)) {
            x@other$verbose <- 2
            if (verbose >= 3) {
                cat(
                    report(
                        "  Locus metric 'verbose' does not exist, creating slot @other$verbose, setting to default [2]\n"
                    )
                )
            }
        }
    }
    ##################################
    if (datatype == "SilicoDArT") {
        if (verbose >= 2) {
            cat(
                report(
                    "  Resetting flags for CallRate, PIC, OneRatio, monomorphs to",
                    set,
                    "\n"
                )
            )
            cat(
                report(
                    "  Setting SNP flags for AvgPIC, OneRatioRef, OneRatioSnp, PICRef, PICSnp, maf, FreqHets, FreqHomRef, FreqHomSnp to FALSE\n"
                )
            )
        }
        
        # Check if the x@other$loc.metrics slot exists, if not, create as a dataframe
        if (is.null(x@other$loc.metrics)) {
            x@other$loc.metrics <- as.data.frame(array(NA, nLoc(x)))
        }
        # Check if the x@other$loc.metrics.flags slot exists, if not, create as a dataframe
        if (is.null(x@other$loc.metrics.flags)) {
            x@other$loc.metrics.flags <- as.data.frame(array(NA, 1))
        }
        
        # AvgPIC
        if (is.null(x@other$loc.metrics$AvgPIC)) {
            x@other$loc.metrics$AvgPIC <- array(NA, nLoc(x))
            if (verbose >= 3) {
                cat(
                    report(
                        "  Locus metric AvgPIC does not exist, creating slot @other$loc.metrics$AvgPIC\n"
                    )
                )
            }
        }
        x@other$loc.metrics.flags$AvgPIC <- FALSE
        
        # OneRatioRef
        if (is.null(x@other$loc.metrics$OneRatioRef)) {
            x@other$loc.metrics$OneRatioRef <- array(NA, nLoc(x))
            if (verbose >= 3) {
                cat(
                    report(
                        "  Locus metric OneRatioRef does not exist, creating slot @other$loc.metrics$OneRatioRef\n"
                    )
                )
            }
        }
        x@other$loc.metrics.flags$OneRatioRef <- FALSE
        
        # OneRatioSnp
        if (is.null(x@other$loc.metrics$OneRatioSnp)) {
            x@other$loc.metrics$OneRatioSnp <- array(NA, nLoc(x))
            if (verbose >= 3) {
                cat(
                    report(
                        "  Locus metric OneRatioSnp does not exist, creating slot @other$loc.metrics$OneRatioSnp\n"
                    )
                )
            }
        }
        x@other$loc.metrics.flags$OneRatioSnp <- FALSE
        
        # PICRef
        if (is.null(x@other$loc.metrics$PICRef)) {
            x@other$loc.metrics$PICRef <- array(NA, nLoc(x))
            if (verbose >= 3) {
                cat(
                    report(
                        "  Locus metric PICRef does not exist, creating slot @other$loc.metrics$PICRef\n"
                    )
                )
            }
        }
        x@other$loc.metrics.flags$PICRef <- FALSE
        
        # PICSnp
        if (is.null(x@other$loc.metrics$PICSnp)) {
            x@other$loc.metrics$PICSnp <- array(NA, nLoc(x))
            if (verbose >= 3) {
                cat(
                    report(
                        "  Locus metric PICSnp does not exist, creating slot @other$loc.metrics$PICSnp\n"
                    )
                )
            }
        }
        x@other$loc.metrics.flags$PICSnp <- FALSE
        
        # CallRate
        if (is.null(x@other$loc.metrics$CallRate)) {
            x@other$loc.metrics$CallRate <- array(NA, nLoc(x))
            if (verbose >= 3) {
                cat(
                    report(
                        "  Locus metric CallRate does not exist, creating slot @other$loc.metrics$CallRate\n"
                    )
                )
            }
        }
        x@other$loc.metrics.flags$CallRate <- set
        
        # FreqHomRef
        if (is.null(x@other$loc.metrics$FreqHomRef)) {
            x@other$loc.metrics$FreqHomRef <- array(NA, nLoc(x))
            if (verbose >= 3) {
                cat(
                    report(
                        "  Locus metric FreqHomRef does not exist, creating slot @other$loc.metrics$FreqHomRef\n"
                    )
                )
            }
        }
        x@other$loc.metrics.flags$FreqHomRef <- FALSE
        
        # FreqHomSnp
        if (is.null(x@other$loc.metrics$FreqHomSnp)) {
            x@other$loc.metrics$FreqHomSnp <- array(NA, nLoc(x))
            if (verbose >= 3) {
                cat(
                    report(
                        "  Locus metric FreqHomSnp does not exist, creating slot @other$loc.metrics$FreqHomSnp\n"
                    )
                )
            }
        }
        x@other$loc.metrics.flags$FreqHomSnp <- FALSE
        
        # FreqHets
        if (is.null(x@other$loc.metrics$FreqHets)) {
            x@other$loc.metrics$FreqHets <- array(NA, nLoc(x))
            if (verbose >= 3) {
                cat(
                    report(
                        "  Locus metric FreqHets does not exist, creating slot @other$loc.metrics$FreqHets\n"
                    )
                )
            }
        }
        x@other$loc.metrics.flags$FreqHets <- FALSE
        
        # monomorphs if (is.null(x@other$loc.metrics$monomorphs)) { x@other$loc.metrics$monomorphs <- array(NA,nLoc(x)) if (verbose >= 3){
        # cat(' Locus metric monomorphs does not exist, creating slot @other$loc.metrics$monomorphs\n') } }
        x@other$loc.metrics.flags$monomorphs <- set
        
        # maf
        if (is.null(x@other$loc.metrics$maf)) {
            x@other$loc.metrics$maf <- array(NA, nLoc(x))
            if (verbose >= 3) {
                cat(
                    report(
                        "  Locus metric maf does not exist, creating slot @other$loc.metrics$maf\n"
                    )
                )
            }
        }
        x@other$loc.metrics.flags$maf <- FALSE
        
        # OneRatio
        if (is.null(x@other$loc.metrics$OneRatio)) {
            x@other$loc.metrics$OneRatio <- array(NA, nLoc(x))
            if (verbose >= 3) {
                cat(
                    report(
                        "  Locus metric OneRatio does not exist, creating slot @other$loc.metrics$OneRatio\n"
                    )
                )
            }
        }
        x@other$loc.metrics.flags$OneRatio <- set
        
        # PIC
        if (is.null(x@other$loc.metrics$PIC)) {
            x@other$loc.metrics$PIC <- array(NA, nLoc(x))
            if (verbose >= 3) {
                cat(
                    report(
                        "  Locus metric PIC does not exist, creating slot @other$loc.metrics$PIC\n"
                    )
                )
            }
        }
        x@other$loc.metrics.flags$PIC <- set
        
        # allna
        x@other$loc.metrics.flags$allna <- set
        
        
        # verbosity
        if (is.null(x@other$verbose)) {
            x@other$verbose <- 2
            if (verbose >= 3) {
                cat(
                    report(
                        "  Locus metric 'verbose' does not exist, creating slot @other$verbose, setting to default [2]\n"
                    )
                )
            }
        }
    }
    
    # ADD TO HISTORY not in utils functions
    
    # FLAG SCRIPT END
    
    if (verbose >= 1) {
        cat(report("Completed:", funname, "\n"))
    }
    
    return(x)
    
}
