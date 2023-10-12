#' @name gl.report.maf
#' @title Reports minor allele frequency (MAF) for each locus in a SNP dataset
#' @description
#' This script provides summary histograms of MAF for each
#' population in the dataset and an overall histogram to assist the decision of
#' choosing thresholds for the filter function \code{\link{gl.filter.maf}}
#' @param x Name of the genlight object containing the SNP data [required].
#' @param maf.limit Show histograms MAF range <= maf.limit [default 0.5].
#' @param ind.limit Show histograms only for populations of size greater than
#' ind.limit [default 5].
#' @param plot.out Specify if plot is to be produced [default TRUE].
#' @param plot_theme Theme for the plot. See Details for options
#' [default theme_dartR()].
#' @param plot_colors_pop A color palette for population plots
#' [default discrete_palette].
#' @param plot_colors_all List of two color names for the borders and fill of
#' the overall plot [default two_colors].
#' @param bins Number of bins to display in histograms [default 25].
#' @param save2tmp If TRUE, saves any ggplots and listings to the session
#' temporary directory (tempdir) [default FALSE].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log ; 3, progress and results summary; 5, full report
#' [default NULL, unless specified using gl.set.verbosity].
#' @details
#'The function \code{\link{gl.filter.maf}} will filter out the
#'  loci with MAF below a specified threshold.
#'
#'\strong{ Function's output }
#'
#'  The minimum, maximum, mean and a tabulation of MAF quantiles against
#'  thresholds rate are provided. Output also includes a boxplot and a
#'  histogram.
#'
#'  This function reports the
#'  MAF for each of several quantiles. Quantiles are
#'  partitions of a finite set of values into q subsets of (nearly) equal sizes.
#'  In this function q = 20. Quantiles are useful measures because they are less
#'  susceptible to long-tailed distributions and outliers.
#'
#'  Plots and table are saved to the temporal directory (tempdir) and can be
#'  accessed with the function \code{\link{gl.print.reports}} and listed with
#'  the function \code{\link{gl.list.reports}}. Note that they can be accessed
#'  only in the current R session because tempdir is cleared each time that the
#'   R session is closed.
#'
#' Examples of other themes that can be used can be consulted in \itemize{
#'  \item \url{https://ggplot2.tidyverse.org/reference/ggtheme.html} and \item
#'  \url{https://yutannihilation.github.io/allYourFigureAreBelongToUs/ggthemes/}
#'  }
#' @return An unaltered genlight object
#' @author Custodian: Arthur Georges (Post to 
#' \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' gl <- gl.report.maf(platypus.gl)
#' @seealso \code{\link{gl.filter.maf}}, \code{\link{gl.list.reports}},
#'  \code{\link{gl.print.reports}}
#' @family report functions
#' @export

gl.report.maf <- function(x,
                          maf.limit = 0.5,
                          ind.limit = 5,
                          plot.out = TRUE,
                          plot_theme = theme_dartR(),
                          plot_colors_pop = discrete_palette,
                          plot_colors_all = two_colors,
                          bins = 25,
                          save2tmp = FALSE,
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
    
    # FUNCTION SPECIFIC ERROR CHECKING
    
    if (maf.limit > 0.5 | maf.limit <= 0) {
        cat(warn(
            "Warning: maf.limit must be in the range (0,0.5], set to 0.5\n"
        ))
        maf.limit <- 0.5
    }
    
    if (ind.limit <= 0) {
        cat(
            warn(
                "Warning: ind.limit must be an integer > 0 and less than population size, set to 5\n"
            )
        )
        ind.limit <- 5
    }
    
    # FLAG SCRIPT START
    
    if (verbose >= 1) {
        if (verbose == 5) {
            cat(report("Starting", funname, "[ Build =", build, "]\n\n"))
        } else {
            cat(report("Starting", funname, "\n\n"))
        }
    }
    
    # DO THE JOB
    
    pops_maf <- seppop(x)
    
    mafs_plots <- lapply(pops_maf, function(z) {
        z$other$loc.metrics <- as.data.frame(z$other$loc.metrics)
        z <- gl.filter.monomorphs(z, verbose = 0)
        z <- gl.recalc.metrics(z, verbose = 0)
        mafs_per_pop_temp <- z$other$loc.metrics$maf
        mafs_per_pop <-
            mafs_per_pop_temp[mafs_per_pop_temp < maf.limit]
        p_temp <-
            ggplot(as.data.frame(mafs_per_pop), aes(x = mafs_per_pop)) + 
            geom_histogram(bins = bins, color = "black", fill = plot_colors_pop(bins)) +
            xlab("Minor Allele Frequency") + 
            ylab("Count") + 
            xlim(0, maf.limit) + 
            plot_theme + 
            ggtitle(paste(popNames(z), "n =", nInd(z)))
        return(p_temp)
    })
    
    # Check for status -- any populations with ind > ind.limit; and is nPop > 1
    
    ind_per_pop <- unlist(lapply(pops_maf, nInd))
    
    test_pop <- as.data.frame(cbind(pop = names(ind_per_pop), ind_per_pop))
    test_pop$ind_per_pop <- as.numeric(test_pop$ind_per_pop)
    
    x2 <- x
    x2$other$loc.metrics <- as.data.frame(x2$other$loc.metrics)
    x2 <- gl.filter.monomorphs(x2, verbose = 0)
    x2 <- gl.recalc.metrics(x2, verbose = 0)
    maf <- data.frame(x2@other$loc.metrics$maf)
    colnames(maf) <- "maf"
    
    # Print out some statistics
    stats <- summary(x2@other$loc.metrics$maf)
    cat("  Reporting Minor Allele Frequency (MAF) by Locus\n")
    cat("  No. of loci =", nLoc(x), "\n")
    cat("  No. of individuals =", nInd(x), "\n")
    cat("    Minimum      : ", stats[1], "\n")
    cat("    1st quantile : ", stats[2], "\n")
    cat("    Median       : ", stats[3], "\n")
    cat("    Mean         : ", stats[4], "\n")
    cat("    3r quantile  : ", stats[5], "\n")
    cat("    Maximum      : ", stats[6], "\n")
    cat("    Missing Rate Overall: ", round(sum(is.na(as.matrix(
        x
    ))) / (nLoc(x) * nInd(x)), 2), "\n\n")
    
    # Determine the loss of loci for a given threshold using quantiles
    quantile_res <- quantile(maf$maf, probs = seq(0, 1, 1 / 20),type=1)
    retained <- unlist(lapply(quantile_res, function(y) {
        res <- length(maf$maf[maf$maf >= y])
    }))
    pc.retained <- round(retained * 100 / nLoc(x2), 1)
    filtered <- nLoc(x2) - retained
    pc.filtered <- 100 - pc.retained
    df <-
        data.frame(as.numeric(sub("%", "", names(quantile_res))),
                   quantile_res,
                   retained,
                   pc.retained,
                   filtered,
                   pc.filtered)
    colnames(df) <-
        c("Quantile",
          "Threshold",
          "Retained",
          "Percent",
          "Filtered",
          "Percent")
    df <- df[order(-df$Quantile), ]
    df$Quantile <- paste0(df$Quantile, "%")
    rownames(df) <- NULL
    
    # testing which populations comply with thresholds
    popn.hold <-
        test_pop[which(test_pop$ind_per_pop >= ind.limit), "pop"]
    mafs_plots_print <- mafs_plots[popn.hold]
    
    if (length(popn.hold) > 1) {
        title.str <- "Minor Allele Frequency\nOverall"
        
        p_all <-
            ggplot(as.data.frame(maf), aes(x = maf)) + 
            geom_histogram(bins = bins,color = plot_colors_all[1], fill = plot_colors_all[2]) +
            xlab("Minor Allele Frequency") + 
            ylab("Count") +
            xlim(0, maf.limit) + 
            plot_theme + 
            ggtitle(title.str)
        
        row_plots <- ceiling(length(popn.hold) / 3) + 1
        p3 <- p_all + mafs_plots_print + plot_layout(ncol = 3, nrow = row_plots)
    }
    
    if (length(popn.hold) == 0) {
        if (verbose >= 1) {
            cat(
                important(
                    "  No populations met minimum limits on number of individuals or loci, reporting for overall\n"
                )
            )
        }
        title.str <- "Minor Allele Frequency\nOverall"
        p3 <-
            ggplot(as.data.frame(maf), aes(x = maf)) + 
            geom_histogram(bins = bins,color = plot_colors_all[1],fill = plot_colors_all[2]) +
            xlab("Minor Allele Frequency") + 
            ylab("Count") + 
            xlim(0, maf.limit) + 
            plot_theme + 
            ggtitle(title.str)
    }
    
    if (length(popn.hold) == 1) {
        if (verbose >= 3) {
            cat(
                important(
                    "  Only one population met minimum limits on number of individuals or loci\n"
                )
            )
        }
        title.str <-
            paste("Minor Allele Frequency\n", popn.hold)
        p3 <-
            ggplot(as.data.frame(maf), aes(x = maf)) + 
            geom_histogram(bins = bins,color = plot_colors_all[1],fill = plot_colors_all[2]) +
            xlab("Minor Allele Frequency") + 
            ylab("Count") + xlim(0, maf.limit) + 
            plot_theme + 
            ggtitle(title.str)
    }
    
    if (nPop(x2) == 1) {
        if (verbose >= 1) {
            cat(important("  Only one population specified\n"))
        }
        title.str <-
            paste("Minor Allele Frequency\n", pop(x2)[1])
        p3 <-
            ggplot(as.data.frame(maf), aes(x = maf)) + 
            geom_histogram(bins = bins, color = plot_colors_all[1],fill = plot_colors_all[2]) +
            xlab("Minor Allele Frequency") +
            ylab("Count") + 
            xlim(0, maf.limit) + 
            plot_theme + 
            ggtitle(title.str)
    }
    
    # PRINTING OUTPUTS
    if (plot.out) {
        suppressWarnings(print(p3))
    }
    print(df)
    
    # SAVE INTERMEDIATES TO TEMPDIR
    
    # creating temp file names
    if (save2tmp) {
        if (plot.out) {
            temp_plot <- tempfile(pattern = "Plot_")
            match_call <-
                paste0(names(match.call()),
                       "_",
                       as.character(match.call()),
                       collapse = "_")
            # saving to tempdir
            saveRDS(list(match_call, p3), file = temp_plot)
            if (verbose >= 2) {
                cat(report("  Saving the ggplot to session tempfile\n"))
            }
        }
        temp_table <- tempfile(pattern = "Table_")
        saveRDS(list(match_call, df), file = temp_table)
        if (verbose >= 2) {
            cat(report("  Saving tabulation to session tempfile\n"))
            cat(
                report(
                    "  NOTE: Retrieve output files from tempdir using gl.list.reports() and gl.print.reports()\n"
                )
            )
        }
    }
    
    # FLAG SCRIPT END
    
    if (verbose >= 1) {
        cat(report("Completed:", funname, "\n"))
    }
    
    # RETURN
    invisible(x)
    
}
