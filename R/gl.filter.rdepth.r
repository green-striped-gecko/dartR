#' @name gl.filter.rdepth
#' @title Filters loci based on counts of sequence tags scored at a locus (read
#'  depth)
#'
#' @description
#' SNP datasets generated by DArT report AvgCountRef and AvgCountSnp as counts
#' of sequence tags for the reference and alternate alleles respectively. These
#' can be used to back calculate Read Depth. Fragment presence/absence datasets
#' as provided by DArT (SilicoDArT) provide Average Read Depth and Standard
#'  Deviation of Read Depth as standard columns in their report.
#'
#' Filtering on Read Depth using the companion script gl.filter.rdepth can be on
#'  the basis of loci with exceptionally low counts,
#' or loci with exceptionally high counts.
#'
#' @param x Name of the genlight object containing the SNP or tag
#' presence/absence data [required].
#' @param lower Lower threshold value below which loci will be removed
#'  [default 5].
#' @param upper Upper threshold value above which loci will be removed
#'  [default 50].
#' @param plot.out Specify if plot is to be produced [default TRUE].
#' @param plot_theme Theme for the plot. See Details for options
#' [default theme_dartR()].
#' @param plot_colors List of two color names for the borders and fill of the
#'  plots [default two_colors].
#' @param save2tmp If TRUE, saves any ggplots and listings to the session
#' temporary directory (tempdir) [default FALSE].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#'  progress log; 3, progress and results summary; 5, full report
#'   [default 2, unless specified using gl.set.verbosity].
#'
#' @details
#'  For examples of themes, see:
#'  \itemize{
#'  \item \url{https://ggplot2.tidyverse.org/reference/ggtheme.html} and \item
#'  \url{https://yutannihilation.github.io/allYourFigureAreBelongToUs/ggthemes/}
#'  }
#' @return Returns a genlight object retaining loci with a Read Depth in the
#' range specified by the lower and upper threshold.
#' @author Custodian: Arthur Georges (Post to
#' \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' # SNP data
#'   gl.report.rdepth(testset.gl)
#'   result <- gl.filter.rdepth(testset.gl, lower=8, upper=50, verbose=3)
#' # Tag P/A data
#'   result <- gl.filter.rdepth(testset.gs, lower=8, upper=50, verbose=3)
#'
#' @seealso \code{\link{gl.filter.rdepth}}
#'
#' @family filters and filter reports
#' @import patchwork
#' @export

gl.filter.rdepth <-  function(x,
                              lower = 5,
                              upper = 50,
                              plot.out = TRUE,
                              plot_theme = theme_dartR(),
                              plot_colors = two_colors,
                              save2tmp = FALSE,
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
    
    # FUNCTION SPECIFIC ERROR CHECKING
    
    # DO THE JOB
    
    n0 <- nLoc(x)
    
    if (datatype == "SilicoDArT") {
        rdepth <- x@other$loc.metrics$AvgReadDepth
    } else if (datatype == "SNP") {
        rdepth <- x@other$loc.metrics$rdepth
    }
    
    # Remove SNP loci with rdepth < threshold
    
    if (verbose >= 2) {
        cat(report(
            "  Removing loci with rdepth <=",
            lower,
            "and >=",
            upper,
            "\n"
        ))
    }
    
    index <- (rdepth >= lower & rdepth <= upper)
    
    x2 <- x[, index]
    # Remove the corresponding records from the loci metadata
    x2@other$loc.metrics <- x@other$loc.metrics[index,]
    
    # PLOT HISTOGRAMS, BEFORE AFTER
    if (plot.out) {
        plotvar <- rdepth
        # min <- min(plotvar,lower) min <- trunc(min*100)/100
        max <- max(plotvar, upper, na.rm = TRUE)
        max <- ceiling(max / 10) * 10
        if (datatype == "SNP") {
            xlabel <- "Pre-filter SNP read depth"
        } else {
            xlabel <- "Pre-filter P/A read depth"
        }
        p1 <-
            ggplot(data.frame(plotvar), aes(x = plotvar)) + 
            geom_histogram(bins = 100,
                           color = plot_colors[1],
                           fill = plot_colors[2]) + 
            coord_cartesian(xlim = c(0, max)) + 
            geom_vline(xintercept = lower, color = "red", size = 1) +
            geom_vline(xintercept = upper, color = "red", size = 1) + 
            xlab(xlabel) +
            ylab("Count") + 
            plot_theme
        
        if (datatype == "SilicoDArT") {
            rdepth <- x2@other$loc.metrics$AvgReadDepth
        } else if (datatype == "SNP") {
            rdepth <- x2@other$loc.metrics$rdepth
        }
        plotvar <- rdepth
        # min <- min(plotvar,lower) min <- trunc(min*100)/100
        max <- max(plotvar, upper, na.rm = TRUE)
        max <- ceiling(max / 10) * 10
        if (datatype == "SNP") {
            xlabel <- "Post-filter SNP read depth"
        } else {
            xlabel <- "Post-filter P/A read depth"
        }
        p2 <-
            ggplot(data.frame(plotvar), aes(x = plotvar)) +
            geom_histogram(bins = 100,
                           color = plot_colors[1],
                           fill = plot_colors[2]) +
            coord_cartesian(xlim = c(0, max)) +
            geom_vline(xintercept = lower, color = "red",size = 1) + 
            geom_vline(xintercept = upper,color = "red",size = 1) + 
            xlab(xlabel) +
            ylab("Count") +
            plot_theme
        
        p3 <- (p1 / p2) + plot_layout(heights = c(1, 1))
        
        print(p3)
    }
    
    # REPORT A SUMMARY
    if (verbose >= 3) {
        cat("  Summary of filtered dataset\n")
        cat("    Initial no. of loci =", n0, "\n")
        # cat(paste(' read depth >=',lower,'and read depth <=',upper,'\n'))
        cat("    No. of loci deleted =", (n0 - nLoc(x2)), "\n")
        cat(paste("    No. of loci retained:", nLoc(x2), "\n"))
        cat(paste("    No. of individuals:", nInd(x2), "\n"))
        cat(paste("    No. of populations: ", length(levels(
            factor(pop(x2))
        )), "\n"))
    }
    
    # SAVE INTERMEDIATES TO TEMPDIR
    if (save2tmp & plot.out) {
        temp_plot <-
            tempfile(pattern = paste0(
                "Plot",
                paste0(
                    names(match.call()),
                    "_",
                    as.character(match.call()),
                    collapse = "_"
                ),
                "_"
            ))
        
        # saving to tempdir
        saveRDS(p3, file = temp_plot)
        if (verbose >= 2) {
            cat(report(
                "  Saving the plot in ggplot format to the session tempfile\n"
            ))
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
