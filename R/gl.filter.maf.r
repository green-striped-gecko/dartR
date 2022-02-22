#' @name gl.filter.maf
#' @title Filters loci on the basis of minor allele frequency (MAF) in a genlight
#'  {adegenet} object
#' @description
#' This script calculates the minor allele frequency for each locus and updates
#' the locus metadata for FreqHomRef, FreqHomSnp, FreqHets and MAF (if it
#' exists). It then uses the updated metadata for MAF to filter loci.
#'
#' @param x Name of the genlight object containing the SNP data [required].
#' @param threshold Threshold MAF -- loci with a MAF less than the threshold
#' will be removed [default 0.01].
#' @param by.pop Whether MAF should be calculated by population [default FALSE].
#' @param pop.limit Minimum number of populations in which MAF should be less 
#' than the threshold for a locus to be filtered out. Only used if by.pop=TRUE. 
#' The default value is half of the populations [default ceiling(nPop(x)/2)].
#' @param ind.limit Minimum number of individuals that a population should 
#' contain to calculate MAF. Only used if by.pop=TRUE [default 10].
#' @param recalc Recalculate the locus metadata statistics if any individuals
#' are deleted in the filtering [default FALSE].
#' @param plot.out Specify if histograms of call rate, before and after, are to
#' be produced [default TRUE].
#' @param plot_theme User specified theme for the plot [default theme_dartR()].
#' @param plot_colors List of two color names for the borders and fill of the
#' plots [default two_colors].
#' @param bins Number of bins to display in histograms [default 25].
#' @param save2tmp If TRUE, saves any ggplots and listings to the session
#'  temporary directory (tempdir) [default FALSE].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#'  [default 2, unless specified using gl.set.verbosity].
#' @return The reduced genlight dataset
#' @export
#' @author Custodian: Luis Mijangos -- Post to
#'  \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' result <- gl.filter.monomorphs(testset.gl)
#' result <- gl.filter.maf(result, threshold=0.05, verbose=3)

gl.filter.maf <- function(x,
                          threshold = 0.01,
                          by.pop = FALSE,
                          pop.limit = ceiling(nPop(x)/2),
                          ind.limit = 10,
                          recalc = FALSE,
                          plot.out = TRUE,
                          plot_theme = theme_dartR(),
                          plot_colors = two_colors,
                          bins = 25,
                          save2tmp = FALSE,
                          verbose = NULL) {
    hold <- x 
    
    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)
    
    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     build = "Jody",
                     verbosity = verbose)
    
    # CHECK DATATYPE
    datatype <- utils.check.datatype(x, verbose = verbose)
    
    # Work around a bug in adegenet if genlight object is created by subsetting
    if (nLoc(x) != nrow(x@other$loc.metrics)) {
        stop(
            error(
                "The number of rows in the loc.metrics table does not match the number of loci in your genlight object!"
            )
        )
    }
    
    # Set a population if none is specified (such as if the genlight object has been generated manually)
    if (is.null(pop(x)) |
        is.na(length(pop(x))) | length(pop(x)) <= 0) {
        if (verbose >= 2) {
            cat(
                report(
                    "  Population assignments not detected, individuals assigned to a single population labelled 'pop1'\n"
                )
            )
        }
        pop(x) <- array("pop1", dim = nInd(x))
        pop(x) <- as.factor(pop(x))
    }
    
    # Check for monomorphic loci
    tmp <- gl.filter.monomorphs(x, verbose = 0)
    if ((nLoc(tmp) < nLoc(x)) & verbose >= 2) {
        cat(warn("  Warning: genlight object contains monomorphic loci\n"))
    }
    
    # FUNCTION SPECIFIC ERROR CHECKING
    
    if (threshold > 0.5 | threshold <= 0) {
        cat(
            warn(
                "  Warning: threshold must be in the range (0,0.5], but usually small, set to 0.01\n"
            )
        )
        threshold <- 0.01
    }
    
    # DO THE JOB
    
    if(by.pop){
        if (verbose >= 3) {
            cat(report(
                "  Removing loci with MAF <",threshold, "in at least",pop.limit,"populations and recalculating FreqHoms and FreqHets\n"
            ))
        }
        x <- utils.recalc.maf(x, verbose = 0)
        pop.list <- seppop(x)
        # getting populations with more than ind.limit
        ind_per_pop <- which(unlist(lapply(pop.list, nInd))>=ind.limit)
        pop.list <- pop.list[ind_per_pop]
        # recalculating MAF by population
        pop.list <- lapply(pop.list,utils.recalc.maf,verbose=0)
        # getting loci with MAF < threshold
        loci.list <- lapply(pop.list,function(y){
             y$other$loc.metrics$maf <= threshold
            })
        # getting the loci in which MAF < threshold and in at least pop.limit
        # populations
        loci.list <- Reduce("+",loci.list)
        loci.list <- which(loci.list>=pop.limit)
        x2 <- x[, -loci.list]
        x2@other$loc.metrics <- x@other$loc.metrics[-loci.list,]
        x2 <- utils.recalc.maf(x2, verbose = 0)
    }else{
        # Recalculate the relevant loc.metrics
        if (verbose >= 3) {
            cat(report(
                "  Removing loci with MAF <",threshold, "over all the dataset and recalculating FreqHoms and FreqHets\n"
            ))
        }
        
        x <- utils.recalc.maf(x, verbose = 0)
        
        # Remove loci with NA count <= 1-threshold
        index <- x@other$loc.metrics$maf >= threshold
        x2 <- x[, index]
        x2@other$loc.metrics <- x@other$loc.metrics[index,]
        x2 <- utils.recalc.maf(x2, verbose = 0)
    }
    
    if(plot.out){
        maf <- NULL
    # Plot a histogram of MAF
    maf_pre <- data.frame(x@other$loc.metrics$maf)
    colnames(maf_pre) <- "maf"
    min <- min(maf_pre, threshold, na.rm = TRUE)
    min <- trunc(min * 100) / 100
 
    p1 <-
        ggplot(as.data.frame(maf_pre), aes(x = maf)) + 
        geom_histogram(bins = bins,color = plot_colors[1],fill = plot_colors[2]) +
        coord_cartesian(xlim = c(min, 1)) + 
        geom_vline(xintercept = threshold,color = "red",size = 1) + 
        xlab("Pre-filter SNP MAF\nOver all populations") + 
        ylab("Count") +
        plot_theme
    
    maf_post <- data.frame(x2@other$loc.metrics$maf)
    colnames(maf_post) <- "maf"
    min <- min(maf_post, threshold, na.rm = TRUE)
    min <- trunc(min * 100) / 100

    p2 <-
        ggplot(as.data.frame(maf_post), aes(x = maf)) + 
        geom_histogram(bins = bins,color = plot_colors[1],fill = plot_colors[2]) +
        coord_cartesian(xlim = c(min, 1)) + 
        geom_vline(xintercept = threshold,color = "red", size = 1) + 
        xlab("Post-filter SNP MAF\nOver all populations") + 
        ylab("Count") +
        plot_theme
    }
    
    if (recalc) {
        # Recalculate all metrics(flags reset in utils scripts)
        x2 <- gl.recalc.metrics(x2, verbose = verbose)
    } else {
        # Reset the flags as FALSE for all metrics except MAF (dealt with elsewhere)
        x2@other$loc.metrics.flags$AvgPIC <- FALSE
        x2@other$loc.metrics.flags$OneRatioRef <- FALSE
        x2@other$loc.metrics.flags$OneRatioSnp <- FALSE
        x2@other$loc.metrics.flags$PICRef <- FALSE
        x2@other$loc.metrics.flags$PICSnp <- FALSE
        x2@other$loc.metrics.flags$FreqHets <- FALSE
        x2@other$loc.metrics.flags$FreqHomRef <- FALSE
        x2@other$loc.metrics.flags$FreqHomSnp <- FALSE
        x2@other$loc.metrics.flags$CallRate <- FALSE
    }
    
    # REPORT A SUMMARY
    if (verbose >= 3) {
        cat("  Summary of filtered dataset\n")
        cat("  MAF for loci >", threshold, "\n")
        cat("  Initial number of loci:", nLoc(x), "\n")
        cat("  Number of loci deleted:", nLoc(x) - nLoc(x2), "\n")
        cat("  Final number of loci:", nLoc(x2), "\n")
    }
    
    # PRINTING OUTPUTS using package patchwork
    p3 <- (p1 / p2) + plot_layout(heights = c(1, 1))
    if (plot.out) {
        print(p3)
    }
    
    # SAVE INTERMEDIATES TO TEMPDIR
    if (save2tmp & plot.out) {
        # creating temp file names
        temp_plot <- tempfile(pattern = "Plot_")
        match_call <-
            paste0(names(match.call()),
                   "_",
                   as.character(match.call()),
                   collapse = "_")
        # saving to tempdir
        saveRDS(list(match_call, p3), file = temp_plot)
        if (verbose >= 2) {
            cat(report("  Saving ggplot(s) to the session tempfile\n"))
            cat(
                report(
                    "  NOTE: Retrieve output files from tempdir using gl.list.reports() and gl.print.reports()\n"
                )
            )
        }
    }
    
    # ADD TO HISTORY
    nh <- length(x2@other$history)
    x2@other$history[[nh + 1]] <- match.call()
    
    # FLAG SCRIPT END
    
    if (verbose > 0) {
        cat(report("Completed:", funname, "\n"))
    }
    
    return(x2)
}
