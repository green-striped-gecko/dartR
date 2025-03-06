#' @name gl.read.dart
# Preliminaries -- Parameter specifications -------------- 
#'@title Imports DArT data into dartR and converts it into a dartR genlight object
#'
#'@description
#'This function is a wrapper function that allows you to convert your DArT file
#'into a genlight object of class dartR.
#'@details
#'The function will determine automatically if the data are in Diversity Arrays
#'one-row csv format or two-row csv format. 
#'
#'The first 
#'row of data is determined from the number of rows with an * in the first 
#'column. This can be alternatively specified with the topskip parameter.
#'
#'The DArT service code is added to the ind.metrics of the genlight object. 
#'The row containing the service code for each individual can be specified with 
#'the service.row parameter.
#'
#'#'The DArT plate well is added to the ind.metrics of the genlight object. 
#'The row containing the plate well for each individual can be specified with 
#'the plate.row parameter.
#'
#'If individuals have been deleted from the input file manually, then the locus
#'metrics supplied by DArT will no longer be correct and some loci may be
#'monomorphic. To accommodate this, set mono.rm and recalc to TRUE.
#'
#'@param filename File containing the SNP data (csv file) [required].
#'@param ind.metafile File that contains additional information on individuals
#' [required].
#'@param covfilename Deprecated, sse ind.metafile parameter [NULL].
#'@param nas A character specifying NAs [default '-'].
#'@param topskip A number specifying the number of initial rows to be skipped. [default NULL].
#'@param lastmetric Deprecated, specifies the last column of locus metadata. Can be 
#'specified as a column number [default NULL].
#'@param service.row The row number for the DArT service
#'is contained [default 1].
#'@param plate.row The row number the plate well [default 3].
#'@param recalc If TRUE, force the recalculation of locus metrics [default TRUE].
#'@param mono.rm If TRUE, force the removal of monomorphic loci (including all NAs.
#' [default FALSE].
#'@param probar Show progress bar [default FALSE].
#'@param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#'progress log ; 3, progress and results summary; 5, full report
#' [default 2, or as set by gl.set.verbose()].
#'
#'@return A dartR genlight object that contains individual and locus metrics
#'[if data were provided] and locus metrics [from a DArT report].
#'@export
#'
#'@family dartR-base
#'@author Custodian: Bernd Gruber (Post to \url{https://groups.google.com/d/forum/dartr})
#'
#'@examples
#' dartfile <- system.file('extdata','testset_SNPs_2Row.csv', package='dartR')
#' metadata <- system.file('extdata','testset_metadata.csv', package='dartR')
#' gl <- gl.read.dart(dartfile, ind.metafile = metadata, probar=TRUE)
#'
#'
# ------------------------
# Function
gl.read.dart <- function(filename,
                         ind.metafile = NULL,
                         recalc = TRUE,
                         mono.rm = FALSE,
                         nas = "-",
                         topskip = NULL,
                         lastmetric = NULL,
                         covfilename = NULL,
                         service.row = 1,
                         plate.row = 3,
                         probar = FALSE,
                         verbose = NULL) {
# Preliminaries -----------------
  
  # Function kindly provided by Andrew Kowalczyk
  
  getLastMarkerMetaDataField <- function(filepath){
    top <- read.csv(filepath,
                    header = FALSE,
                    nrows = 20,
                    stringsAsFactors = FALSE)
    
    last_metric <- top[last(which(top[,1]=="*"))+1, last(which(top[1,]=="*"))]  
    return(last_metric)
  }
  
  lastmetric <- getLastMarkerMetaDataField(filename)
  
    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)
    
    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     build = "v.2023.2",
                     verbose = verbose)
    
    if (verbose == 0) {
        probar <-FALSE
    }
    
    # DO THE JOB ----------------------
    
    # Deal with the depricated covfilename parameter
    if (is.null(ind.metafile)) {
        ind.metafile <- covfilename
    }
    
    # Read in the data
    dout <-
        utils.read.dart(
            filename = filename,
            nas = nas,
            topskip = topskip,
            lastmetric = lastmetric,
            service.row = service.row,
            plate.row = plate.row,
            verbose = verbose
        )
    
    glout <-
        utils.dart2genlight(
          dart = dout,
            ind.metafile = ind.metafile,
            probar = probar,
            verbose = verbose
        )
    
    if (verbose >= 2) {
        cat(report("  ",nrow(glout),"rows and",ncol(glout),"columns of data read\n"))
    }
    
    # Setting the recalc flags (TRUE=up-to-date, FALSE=no longer valid) for all
    # locus metrics capable of being recalculated
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
            "PIC",
            "allna"
        )
    glout@other$loc.metrics.flags <-
        data.frame(matrix(TRUE, nrow = 1, ncol = length(recalc.flags)))
    names(glout@other$loc.metrics.flags) <- recalc.flags
    glout@other$verbose <- 2
    
    # Calculate locus metrics not provided by DArT 
    #Calculate Read Depth
    # calculating "by hand" rather than getting them from DArT's report because
    # sometimes they are not reported
    # OneRatioRef	The proportion of samples for which the genotype score is "1", 
    # in the Reference allele row	
    # OneRatioSnp	The proportion of samples for which the genotype score is "1", 
    # in the SNP allele row	
    glout@other$loc.metrics$OneRatioRef  <- apply(as.matrix(glout),2,
                                                  function(y){
      (sum(!is.na(y[y==0])) + sum(!is.na(y[y==1])) )/ sum(!is.na(y))
    })
    
    glout@other$loc.metrics$OneRatioSnp  <- apply(as.matrix(glout),2,
                                                  function(y){
      (sum(!is.na(y[y==2])) + sum(!is.na(y[y==1])) )/ sum(!is.na(y))
    })

    if (is.null(glout@other$loc.metrics$rdepth)) {
        if (verbose >= 2) {
            cat(report(
                "  Read depth calculated and added to the locus metrics\n"
            ))
        }
        glout@other$loc.metrics$rdepth <- array(NA, nLoc(glout))
        for (i in 1:nLoc(glout)) {
            called.ind <-
                round(nInd(glout) * glout@other$loc.metrics$CallRate[i], 0)
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
                "  Deleting monomorphic loci and loci with all mising data (optionally requested)\n"
            ))
        }
        glout <- gl.filter.monomorphs(glout, verbose = 0)
        glout <- gl.filter.allna(glout, verbose = 0)
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
    
    # Create the history repository --------------------
    if (is.null(glout@other$history)) {
        glout@other$history <- list(match.call())
    }
    
    glout <- gl.compliance.check(glout)
    
    # FLAG SCRIPT END ---------------------
    if (verbose > 0) {
        cat(report(paste("Completed:", funname, "\n")))
    }
    # End Block -----------------
    
    return(glout)
}
