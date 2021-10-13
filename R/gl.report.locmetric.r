#' @name gl.report.locmetric
#'
#' @title Report summary of the slot $other$loc.metrics
#'
#' @description 
#' This script uses any field with numeric values stored in $other$loc.metrics 
#' to produce summary statistics (mean, minimum, average, quantiles), histograms 
#' and boxplots to assist the decision of choosing thresholds for the filter 
#' function \code{\link{gl.filter.locmetric}}.
#'
#' @param x Name of the genlight object containing the SNP or presence/absence (SilicoDArT) data [required].
#' @param metric Name of the metric to be used for filtering [required].
#' @param plot.out Specify if plot is to be produced [default TRUE].
#' @param plot_theme Theme for the plot. See Details for options [default theme_dartR()].
#' @param plot_colors List of two color names for the borders and fill of the plots [default two_colors].
#' @param save2tmp If TRUE, saves any ggplots and listings to the session 
#' temporary directory (tempdir) [default FALSE].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default NULL, unless specified using gl.set.verbosity]
#'
#' @details 
#'The function \code{\link{gl.filter.locmetric}} will filter out the
#'  loci with a locmetric value below a specified threshold.
#'  
#'The fields that are included in dartR, and a short description, are found 
#'below. Optionally, the user can also set his/her own field by adding a vector
#' into $other$loc.metrics as shown in the example. You can check the names of 
#' all available loc.metrics via: names(gl$other$loc.metrics).
#'
#'\itemize{ 
#'\item SnpPosition - position (zero is position 1) in the sequence tag of the 
#'defined SNP variant base.
#'\item CallRate - proportion of samples for which the genotype call is 
#'non-missing (that is, not '-' ).
#'\item OneRatioRef - proportion of samples for which the genotype score is 0.
#'\item OneRatioSnp - proportion of samples for which the genotype score is 2.
#'\item FreqHomRef - proportion of samples homozygous for the Reference allele. 
#'\item FreqHomSnp - proportion of samples homozygous for the Alternate (SNP) 
#'allele.
#'\item FreqHets - proportion of samples which score as heterozygous, that is, 
#'scored as 1.
#'\item PICRef - polymorphism information content (PIC) for the Reference allele.
#'\item PICSnp - polymorphism information content (PIC) for the SNP.
#'\item AvgPIC - average of the polymorphism information content (PIC) of the
#' reference and SNP alleles.
#'\item AvgCountRef - sum of the tag read counts for all samples, divided by the
#' number of samples with non-zero tag read counts, for the Reference allele row.
#'\item AvgCountSnp - sum of the tag read counts for all samples, divided by the 
#'number of samples with non-zero tag read counts, for the Alternate (SNP) allele
#' row. 
#'\item RepAvg - proportion of technical replicate assay pairs for which the 
#'marker score is consistent.
#'\item rdepth - read depth.
#'}
#'
#'\strong{ Function's output }
#'
#' The minimum, maximum, mean and a tabulation of quantiles of the locmetric 
#' values against thresholds rate are provided. Output also includes a boxplot 
#' and a histogram.
#'  
#'  Quantiles are partitions of a finite set of values into q subsets of (nearly)
#'   equal sizes. In this function q = 20. Quantiles are useful measures because 
#'   they are less susceptible to long-tailed distributions and outliers.
#'  
#'  Plots and table were saved to the temporal directory (tempdir) and can be 
#'  accessed with the function \code{\link{gl.print.reports}} and listed with 
#'  the function \code{\link{gl.list.reports}}. Note that they can be accessed 
#'  only in the current R session because tempdir is cleared each time that the 
#'  R session is closed.
#'
#'  Examples of other themes that can be used can be consulted in \itemize{
#'  \item \url{https://ggplot2.tidyverse.org/reference/ggtheme.html} and \item
#'  \url{https://yutannihilation.github.io/allYourFigureAreBelongToUs/ggthemes/}
#'  }
#'
#' @return An unaltered genlight object
#'
#' @author Luis Mijangos (Post to \url{https://groups.google.com/d/forum/dartr})
#'
#' @examples
#' # adding dummy data
#' test <- testset.gl
#' test$other$loc.metrics$test <- 1:nLoc(test)
#' # SNP data
#' out <- gl.report.locmetric(test,metric='test')
#' 
#' # adding dummy data
#' test.gs <- testset.gs
#' test.gs$other$loc.metrics$test <- 1:nLoc(test.gs)
#' # Tag P/A data
#' out <- gl.report.locmetric(test.gs,metric='test')
#' 
#' @seealso \code{\link{gl.filter.locmetric}}, \code{\link{gl.list.reports}},
#'  \code{\link{gl.print.reports}}
#'  
#' @family filters and filter reports
#'
#' @export
#'

gl.report.locmetric <- function(x, 
                                metric, 
                                plot.out = TRUE,
                                plot_theme = theme_dartR(),
                                plot_colors = two_colors,
                                save2tmp = FALSE,
                                verbose = NULL) {
    # TRAP COMMAND
    
    funname <- match.call()[[1]]
    
    # SET VERBOSITY
    
    verbose <- gl.check.verbosity(verbose)
    
    # CHECKS DATATYPE 
    
    datatype <- utils.check.datatype(x, verbose=verbose)
    
    # FLAG SCRIPT START

    if (verbose >= 1) {
        if (verbose == 5) {
            cat(report("Starting", funname, "[ Build =", build, "]\n\n"))
        } else {
            cat(report("Starting", funname, "\n\n"))
        }
    }

    # FUNCTION SPECIFIC ERROR CHECKING

    # check whether the field exists in the genlight object
    if (!(metric %in% colnames(x$other$loc.metrics))) {
        stop(error("  Fatal Error: name of the metric not found\n"))
    }
    if (!is.numeric(unlist(x$other$loc.metrics[metric]))) {
        stop(error("  Fatal Error: metric is not numeric\n"))
    }

    # DO THE JOB

    field <- which(colnames(x@other$loc.metrics) == metric)

    # get title for plots
    if (all(x@ploidy == 2)) {
        title1 <- paste0("SNP data - ", metric, " by Locus")
    } else {
        title1 <- paste0("Fragment P/A data - ", metric, " by Locus")
    }

    metric_df <- data.frame(x$other$loc.metrics[field])
    colnames(metric_df) <- "field"

    p1 <- ggplot(metric_df, aes(y = field)) + 
        geom_boxplot(color = plot_colors[1], fill = plot_colors[2]) + 
        coord_flip() + 
        plot_theme + 
        xlim(range = c(-1, 1)) +
        ylab(metric) + 
        theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
        ggtitle(title1)

    p2 <- ggplot(metric_df, aes(x = field)) + 
        geom_histogram(bins = 50, color = plot_colors[1], fill = plot_colors[2]) + 
        xlab(metric) + 
        ylab("Count") + 
        plot_theme

    # Print out some statistics
    stats <- summary(metric_df)
    cat("  Reporting", metric, "by Locus\n")
    cat("  No. of loci =", nLoc(x), "\n")
    cat("  No. of individuals =", nInd(x), "\n")
    cat("    Minimum      : ", stats[1], "\n")
    cat("    1st quantile : ", stats[2], "\n")
    cat("    Median       : ", stats[3], "\n")
    cat("    Mean         : ", stats[4], "\n")
    cat("    3r quantile  : ", stats[5], "\n")
    cat("    Maximum      : ", stats[6], "\n\n")

    # Determine the loss of loci for a given threshold
    # using quantiles
    quantile_res <- quantile(metric_df$field, 
                             probs = seq(0, 1, 1/20))
    retained <- unlist(lapply(quantile_res, function(y) {
        res <- length(metric_df$field[metric_df$field>= 
                                            y])
    }))
    pc.retained <- round(retained * 100/nLoc(x), 
                         1)
    filtered <- nLoc(x) - retained
    pc.filtered <- 100 - pc.retained
    df <- data.frame(as.numeric(sub("%", "", names(quantile_res))), 
                     quantile_res, retained, pc.retained, filtered, 
                     pc.filtered)
    colnames(df) <- c("Quantile", "Threshold", 
                      "Retained", "Percent", "Filtered", "Percent")
    df <- df[order(-df$Quantile), ]
    df$Quantile <- paste0(df$Quantile, "%")
    rownames(df) <- NULL
    
    # PRINTING OUTPUTS
    if(plot.out){
        # using package patchwork
        p3 <- (p1/p2) + plot_layout(heights = c(1, 4))
        print(p3)
    }
    print(df)
    
    # SAVE INTERMEDIATES TO TEMPDIR             
    
    # creating temp file names
    if(save2tmp){
        if(plot.out){
            temp_plot <- tempfile(pattern = "Plot_")
            match_call <- paste0(names(match.call()),"_",as.character(match.call()),collapse = "_")
            # saving to tempdir
            saveRDS(list(match_call,p3), file = temp_plot)
            if(verbose>=2){
                cat(report("  Saving the ggplot to session tempfile\n"))
            }
        }
        temp_table <- tempfile(pattern = "Table_")
        saveRDS(list(match_call,df), file = temp_table)
        if(verbose>=2){
            cat(report("  Saving tabulation to session tempfile\n"))
            cat(report("  NOTE: Retrieve output files from tempdir using gl.list.reports() and gl.print.reports()\n"))
        }
    }
    
    # FLAG SCRIPT END
    
    if (verbose >= 1) {
        cat(report("Completed:", funname, "\n"))
    }
    
    # RETURN
    invisible(x)
    
}
