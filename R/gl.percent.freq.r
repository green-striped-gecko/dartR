#' Generate percentage allele frequencies by locus and population
#'
#' This is a support script, to take SNP data or SilicoDArT presence/absence
#' data grouped into populations in a genlight object \{adegenet\} and generate
#' a table of allele frequencies for each population and locus
#'
#' @param x Name of the genlight object containing the SNP or Tag P/A
#' (SilicoDArT) data [required].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2 or as specified using gl.set.verbosity].
#' @return A matrix with allele (SNP data) or presence/absence frequencies
#' (Tag P/A data) broken down by population and locus
#' @export
#' @importFrom plyr rbind.fill
#' @author Custodian: Arthur Georges (Post to
#' \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' m <-  gl.percent.freq(testset.gl)

gl.percent.freq <- function(x,
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
    
    # SCRIPT SPECIFIC ERROR CHECKING
    
    # Checking for and removing monomorphic loci
    if (!(x@other$loc.metrics.flags$monomorphs == TRUE)) {
        if (verbose >= 1) {
            cat(warn(
                "Warning: Monomorphic loci retained, used in calculations\n"
            ))
        }
    }
    # DO THE JOB
    x2 <- seppop(x)
    x2_list <- lapply(x2, as.matrix)
    
    if (datatype == "SilicoDArT") {
        if (verbose >= 2) {
            cat(
                report(
                    "Starting gl.percent.freq: Calculating Tag P/A frequencies for populations\n"
                )
            )
        }
        # Treat SilicoDArT as biallelic, no heterozygotes
        x2_list <- lapply(x2_list, function(x) {
            x[x == 1] <- 2
            return(x)
        })
        
    } else {
        if (verbose >= 1) {
            cat(
                report(
                    "Starting gl.percent.freq: Calculating allele frequencies for populations\n"
                )
            )
        }
        if (verbose >= 3) {
            cat(report("  This may take some time -- be patient\n"))
        }
    }
    
    loc_names <- lapply(x2_list, colnames)
    
    nmissing_temp <- lapply(x2_list, is.na)
    nmissing <- lapply(nmissing_temp, colSums)
    
    n_temp <- lapply(x2_list, nrow)
    n <- lapply(n_temp, rep, nLoc(x))
    
    nobs_temp <- lapply(nmissing, unname)
    nobs <- Map("-", n, nobs_temp)
    
    sum_res <- lapply(x2_list, colSums, na.rm = T)
    
    f <- lapply(x2_list, colMeans, na.rm = T)
    f <- lapply(f, "/", 2)
    f <- lapply(f, "*", 100)
    f <- lapply(f, "round", 2)
    
    m <-
        Map(cbind,
            names(sum_res),
            loc_names,
            sum_res,
            nobs,
            nmissing,
            f,
            n)
    m <- lapply(m, cbind, 1:nLoc(x))
    m <- lapply(m, as.data.frame)
    m <- rbind.fill(m)
    
    colnames(m) <-
        c("popn",
          "locus",
          "sum",
          "nobs",
          "nmissing",
          "frequency",
          "n",
          "loc_order")
    
    m$popn <- as.factor(m$popn)
    m$locus <- as.factor(m$locus)
    m$sum <- as.numeric(as.character(m$sum))
    m$nobs <- as.numeric(as.character(m$nobs))
    m$nmissing <- as.numeric(as.character(m$nmissing))
    m$frequency <- as.numeric(as.character(m$frequency))
    m$n <- as.numeric(as.character(m$n))
    m$loc_order <- as.numeric(as.character(m$loc_order))
    
    m <- m[order(m$loc_order, m$popn),]
    m <- m[,-ncol(m)]
    
    rownames(m) <- NULL
    
    # FLAG SCRIPT END
    
    if (verbose >= 1) {
        cat(report("Completed:", funname, "\n"))
    }
    
    return(m)
    
}
