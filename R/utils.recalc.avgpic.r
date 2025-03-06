#' A utility script to recalculate the OneRatioRef, OneRatioSnp, PICRef, PICSnp,
#'  and AvgPIC by locus after some individuals or populations have been deleted.
#'
#' The locus metadata supplied by DArT has OneRatioRef, OneRatioSnp, PICRef,
#'  PICSnp, and AvgPIC included, but the allelic composition will change when
#'  some individuals,or populations, are removed from the dataset and so the
#'  initial statistics will no longer apply. This script recalculates these
#'  statistics and places the recalculated values in the appropriate place in
#'  the genlight object.
#'
#' If the locus metadata OneRatioRef|Snp, PICRef|Snp and/or AvgPIC do not exist,
#'  the script creates and populates them.
#'
#' @param x Name of the genlight object containing the SNP data [required].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report [default 2].
#' @return The modified genlight object.
#' @author Custodian: Luis Mijangos (Post to
#'  \url{https://groups.google.com/d/forum/dartr})
#' @seealso \code{utils.recalc.metrics} for recalculating all metrics,
#' \code{utils.recalc.callrate} for recalculating CallRate,
#' \code{utils.recalc.freqhomref} for recalculating frequency of homozygous
#' reference, \code{utils.recalc.freqhomsnp} for recalculating frequency of
#' homozygous alternate, \code{utils.recalc.freqhet} for recalculating frequency
#'  of heterozygotes, \code{gl.recalc.maf} for recalculating minor allele
#'  frequency, \code{gl.recalc.rdepth} for recalculating average read depth
#' @examples
#' #out <- utils.recalc.avgpic(testset.gl)

utils.recalc.avgpic <- function(x,
                                verbose = NULL) {
    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)
    
    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     build = "Jody",
                     verbose = verbose)
    
    # CHECK DATATYPE
    datatype <- utils.check.datatype(x, verbose = verbose)
    
    # Check monomorphs have been removed up to date
    if (x@other$loc.metrics.flags$monomorphs == FALSE) {
        if (verbose >= 2) {
            cat(
                warn(
                    "  Warning: Dataset contains monomorphic loci which will be included in the ",
                    funname,
                    " calculations\n"
                )
            )
        }
    }
    
    # FUNCTION SPECIFIC ERROR CHECKING
    
    if (datatype == "SNP") {
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
    }
    
    if (datatype == "SilicoDArT") {
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
    }
    
    # DO THE JOB
    
    t <- as.matrix(x)
    
    if (datatype == "SNP") {
        if (verbose >= 2) {
            cat(
                report(
                    "  Recalculating OneRatioRef, OneRatioSnp, PICRef, PICSnp, AvgPIC\n"
                )
            )
        }
        
        c0 <- colSums(t == 0, na.rm = T)
        c1 <- colSums(t == 1, na.rm = T)
        c2 <- colSums(t == 2, na.rm = T)
        c <- (c0 + c1 + c2)
        x@other$loc.metrics$OneRatioRef <- (c0 + c1) / c
        x@other$loc.metrics.flags$OneRatioRef <- TRUE
        
        x@other$loc.metrics$OneRatioSnp <- (c1 + c2) / c
        x@other$loc.metrics.flags$OneRatioSnp <- TRUE
        
        OneRatioRef <- x@other$loc.metrics$OneRatioRef
        OneRatioSnp <- x@other$loc.metrics$OneRatioSnp
        ZeroRatioRef <- 1 - OneRatioRef
        ZeroRatioSnp <- 1 - OneRatioSnp
        x@other$loc.metrics$PICRef <-
            1 - ((OneRatioRef * OneRatioRef) + (ZeroRatioRef * ZeroRatioRef))
        x@other$loc.metrics.flags$PICRef <- TRUE
        
        x@other$loc.metrics$PICSnp <-
            1 - ((OneRatioSnp * OneRatioSnp) + (ZeroRatioSnp * ZeroRatioSnp))
        x@other$loc.metrics.flags$PICSnp <- TRUE
        
        x@other$loc.metrics$AvgPIC <-
            (x@other$loc.metrics$PICRef + x@other$loc.metrics$PICSnp) / 2
        x@other$loc.metrics.flags$AvgPIC <- TRUE
        
    }
    
    if (datatype == "SilicoDArT") {
        if (verbose >= 2) {
            cat(report("  Recalculating OneRatio, PIC\n"))
        }
        
        OneRatio <- colMeans(t == 1, na.rm = T)
        x@other$loc.metrics$OneRatio <- OneRatio
        x@other$loc.metrics.flags$OneRatio <- TRUE
        
        ZeroRatio <- 1 - OneRatio
        
        x@other$loc.metrics$PIC <-
            1 - ((OneRatio * OneRatio) + (ZeroRatio * ZeroRatio))
        x@other$loc.metrics.flags$PIC <- TRUE
        
    }
    
    # FLAG SCRIPT END
    
    if (verbose > 0) {
        cat(report("Completed:", funname, "\n"))
    }
    
    return(x)
}
