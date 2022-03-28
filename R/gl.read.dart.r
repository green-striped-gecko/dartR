#'@name gl.read.dart
#'
#'@title Imports DArT data into dartR and converts it into a genlight object
#'
#'@description
#'This function is a wrapper function that allows you to convert your DArT file
#'into a genlight object in one step. In previous versions you had to use
#'read.dart and then dart2genlight. In case you have individual metadata for
#'each individual/sample you can specify as before in the dart2genlight command
#'the file that combines the data.
#'
#'@param filename File containing the SNP data (csv file) [required].
#'@param ind.metafile File that contains additional information on individuals
#' [required].
#'@param covfilename Use ind.metafile parameter [depreciated, NULL].
#'@param nas A character specifying NAs [default '-'].
#'@param topskip A number specifying the number of rows to be skipped. If not
#'provided the number of rows to be skipped are 'guessed' by the number of rows
#' with '*' at the beginning [default NULL].
#'@param lastmetric Specifies the last non-genetic column (Default is 'RepAvg').
#' Be sure to check if that is true, otherwise the number of individuals will
#' not match. You can also specify the last column by a number
#'  [default 'RepAvg'].
#'@param service_row The row number in which the information of the DArT service
#'is contained [default 1].
#'@param plate_row The row number in which the information of the plate location
#' is contained [default 3].
#'@param recalc Force the recalculation of locus metrics, in case individuals
#'have been manually deleted from the input csv file [default TRUE].
#'@param mono.rm Force the removal of monomorphic loci (including all NAs), in
#'case individuals have been manually deleted from the input csv file
#' [default FALSE].
#'@param probar Show progress bar [default FALSE].
#'@param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#'progress log ; 3, progress and results summary; 5, full report
#' [default 2, or as set by gl.set.verbose()].
#'
#'@details
#' The dartR genlight object can then be fed into a number of initial screening,
#'  export and export functions provided by the package. For some of the
#'  functions it is necessary to have the metadata that was provided from DArT.
#'  Please check the vignette for more information. Additional information can
#'  also be found in the help documents for \code{\link{utils.read.dart}}.
#'
#'@return A genlight object that contains individual metrics
#'[if data were provided] and locus metrics [from a DArT report].
#'
#'@author Custodian: Bernd Gruber (Post to \url{https://groups.google.com/d/forum/dartr})
#'
#'@examples
#' dartfile <- system.file('extdata','testset_SNPs_2Row.csv', package='dartR')
#' metadata <- system.file('extdata','testset_metadata.csv', package='dartR')
#' gl <- gl.read.dart(dartfile, ind.metafile = metadata, probar=TRUE)
#'
#'@seealso \code{\link{utils.read.dart}}
#'
#'@family input data
#'
#'@export
#'

gl.read.dart <- function(filename,
                         ind.metafile = NULL,
                         recalc = TRUE,
                         mono.rm = FALSE,
                         nas = "-",
                         topskip = NULL,
                         lastmetric = "RepAvg",
                         covfilename = NULL,
                         service_row = 1,
                         plate_row = 3,
                         probar = FALSE,
                         verbose = NULL) {
    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)
    
    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     build = "Jody",
                     verbosity = verbose)
    
    if (verbose == 0) {
        probar = FALSE
    }
    
    # DO THE JOB
    
    # Deal with the redundant covfilename parameter
    if (is.null(ind.metafile)) {
        ind.metafile <- covfilename
    }
    
    dout <-
        utils.read.dart(
            filename = filename,
            nas = nas,
            topskip = topskip,
            lastmetric = lastmetric,
            service_row = service_row,
            plate_row = plate_row,
            verbose = verbose
        )
    glout <-
        utils.dart2genlight(
            dout,
            ind.metafile = ind.metafile,
            probar = probar,
            verbose = verbose
        )
    
    if (verbose >= 2) {
        cat(report(" Data read in. Please check carefully the output above\n"))
    }
    
    # Setting the recalc flags (TRUE=up-to-date, FALSE=no longer valid) for all locus metrics capable of being recalculated
    recalc.flags <-
        c(
            "AvgPIC",
            "OneRatioRef",
            "OneRatioSnp",
            "PICRef",
            "PICSnp",
            "CallRate",
            "maf",
            "FreqHets",
            "FreqHomRef",
            "FreqHomSnp",
            "monomorphs",
            "OneRatio",
            "PIC"
        )
    glout@other$loc.metrics.flags <-
        data.frame(matrix(TRUE, nrow = 1, ncol = length(recalc.flags)))
    names(glout@other$loc.metrics.flags) <- recalc.flags
    glout@other$verbose <- 2
    
    # Calculate locus metrics not provided by DArT Calculate Read Depth
    if (is.null(glout@other$loc.metrics$rdepth)) {
        if (verbose >= 2) {
            cat(report(
                "  Read depth calculated and added to the locus metrics\n"
            ))
        }
        glout@other$loc.metrics$rdepth <- array(NA, nLoc(glout))
        for (i in 1:nLoc(glout)) {
            called.ind <-
                round(nInd(glout) * glout@other$loc.metrics$CallRate[i],
                      0)
            ref.count <-
                called.ind * glout@other$loc.metrics$OneRatioRef[i]
            alt.count <-
                called.ind * glout@other$loc.metrics$OneRatioSnp[i]
            sum.count.ref <-
                ref.count * glout@other$loc.metrics$AvgCountRef[i]
            sum.count.alt <-
                alt.count * glout@other$loc.metrics$AvgCountSnp[i]
            glout@other$loc.metrics$rdepth[i] <-
                round((sum.count.alt + sum.count.ref) / called.ind, 1)
        }
    }
    
    # Calculate MAF
    if (is.null(glout@other$loc.metrics$maf)) {
        utils.recalc.maf(glout, verbose = 0)
        if (verbose >= 2) {
            cat(
                report(
                    "  Minor Allele Frequency (MAF) calculated and added to the locus metrics\n"
                )
            )
        }
    }
    
    # Calculate metrics provided by DArT, as a hedge against the user having deleted individuals from the input csv file
    if (recalc) {
        if (verbose >= 2) {
            cat(
                report(
                    "  Recalculating locus metrics provided by DArT (optionally specified)\n"
                )
            )
        }
        glout <- utils.recalc.avgpic(glout, verbose = 0)
        glout <- utils.recalc.callrate(glout, verbose = 0)
        glout <- utils.recalc.freqhets(glout, verbose = 0)
        glout <- utils.recalc.freqhomref(glout, verbose = 0)
        glout <- utils.recalc.freqhomsnp(glout, verbose = 0)
    }
    
    # Remove monomorphs, which should not be present, but might have been introduced it the user deleted individuals from the input csv
    # file
    
    glout@other$loc.metrics.flags$monomorphs <- FALSE
    if (mono.rm) {
        if (verbose >= 2) {
            cat(report(
                "  Deleting monomorphic loci (optionally requested)\n"
            ))
        }
        glout <- gl.filter.monomorphs(glout, verbose = 0)
    }
    
    # Set the SilicoDArT flags to FALSE
    glout@other$loc.metrics.flags$OneRatio <- FALSE
    glout@other$loc.metrics.flags$PIC <- FALSE
    
    # Provide a summary of the data
    if (verbose >= 3) {
        cat("\nSummary of the SNP dataset\n")
        cat("  No. of loci:", nLoc(glout), "\n")
        cat("  No. of individuals:", nInd(glout), "\n")
        cat("  No. of populations:", nPop(glout), "\n")
        if (!recalc) {
            cat(report(
                "  Locus metrics provided by DArT retained, not recalculated\n"
            ))
        }
        if (!mono.rm) {
            cat(report(
                "  Monomoporhic loci not deleted, assumed absent initially\n\n"
            ))
        }
    }
    
    # Create the history repository
    if (is.null(glout@other$history)) {
        glout@other$history <- list(match.call())
    }
    
    # FLAG SCRIPT END
    
    if (verbose > 0) {
        cat(report(paste("Completed:", funname, "\n")))
    }
    
    return(glout)
}
