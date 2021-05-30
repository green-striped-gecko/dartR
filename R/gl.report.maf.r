#'@name gl.report.maf
#'
#'@title Report minor allele frequency (MAF) for each locus in a SNP dataset
#'
#'@description This script provides summary histograms of MAF for each population in the dataset and an overall histogram as a basis for decisions on filtering.
#'
#'@param x Name of the genlight object containing the SNP data [required]
#'@param maf.limit Show histograms MAF range <= maf.limit [default 0.5]
#'@param ind.limit Show histograms only for populations of size greater than ind.limit [default 5]
#'@param loc.limit Show histograms only for populations with more than loc.limit polymorphic loci [default 30]
#'@param plot_theme Theme for the plot. See Details for options [default theme_dartR()].
#'@param plot_colours_pop A color palette for population plots [default discrete_palette].
#'@param plot_colours_all List of two color names for the borders and fill of the overall plot [default two_colors].
#'
#'@details 
#'The function \code{\link{gl.filter.maf}} will filter out the
#'  loci with MAF below a specified threshold.
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
#'  \strong{Plots and table are saved to the temporal directory (tempdir) and
#'  can be accessed with the function \code{\link{gl.access.report}}. Note that
#'  they can be accessed only in the current R session because tempdir is
#'  cleared each time that an R session is closed.}
#'
#' Examples of other themes that can be used can be consulted in \itemize{
#'  \item \url{https://ggplot2.tidyverse.org/reference/ggtheme.html} and \item
#'  \url{https://yutannihilation.github.io/allYourFigureAreBelongToUs/ggthemes/}
#'  }
#'
#'@return Returns a genlight object with the file names of plots and table that
#'  were saved in the tempdir stored in the slot other$history
#'
#'@author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#'
#'@examples
#' gl <- gl.report.maf(testset.gl)
#' # display of the slot history to identify which report the user wants to plot
#' gl.print.history(gl)
#' # reopen the plots from the last report
#' output <- gl.access.report(gl)
#'
#'@seealso \code{\link{gl.filter.maf}}, \code{\link{gl.access.report}},
#'  \code{\link{gl.print.history}}
#'
#'@export 
#'

gl.report.maf <- function(x, maf.limit = 0.5, ind.limit = 5, loc.limit = 30, plot_theme = theme_dartR(), plot_colours_pop = discrete_palette, 
    plot_colours_all = two_colors) {

    # TRAP COMMAND

    funname <- match.call()[[1]]

    # GENERAL ERROR CHECKING

    x <- utils.check.gl(x)

    # FUNCTION SPECIFIC ERROR CHECKING

    if (maf.limit > 0.5 | maf.limit <= 0) {
        cat(warn("Warning: maf.limit must be in the range (0,0.5], set to 0.5\n"))
        maf.limit <- 0.5
    }

    if (ind.limit <= 0) {
        cat(warn("Warning: ind.limit must be an integer > 0 and less than population size, set to 5\n"))
        ind.limit <- 5
    }

    if (loc.limit <= 1) {
        cat(warn("Warning: loc.limit must be an integer > 1 and less than the the total number of loci, set to 2\n"))
        loc.limit <- 2
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
        mafs_per_pop <- mafs_per_pop_temp[mafs_per_pop_temp < maf.limit]
        p_temp <- ggplot(as.data.frame(mafs_per_pop), aes(x = mafs_per_pop)) + geom_histogram(bins = 12, color = "black", 
            fill = plot_colours_pop(12)) + xlab("Minor Allele Frequency") + ylab("Count") + xlim(0, maf.limit) + plot_theme + 
            ggtitle(paste(popNames(z), "n =", nInd(z)))
        return(p_temp)
    })

    # Check for status -- any populations with loc > loc.limit; ind > ind.limit; and is nPop > 1
    locs_per_pop <- unlist(lapply(pops_maf, function(y) {
        y <- gl.filter.monomorphs(y, verbose = 0)
        return(nLoc(y))
    }))

    ind_per_pop <- unlist(lapply(pops_maf, nInd))

    test_pop <- as.data.frame(cbind(pop = names(locs_per_pop), locs_per_pop, ind_per_pop))
    test_pop$locs_per_pop <- as.numeric(test_pop$locs_per_pop)
    test_pop$ind_per_pop <- as.numeric(test_pop$ind_per_pop)
    
    x2 <- x
    x2$other$loc.metrics <- as.data.frame(x2$other$loc.metrics)
    x2 <- gl.filter.monomorphs(x2, verbose = 0)
    x2 <- gl.recalc.metrics(x2, verbose = 0)
    maf <- data.frame(x2@other$loc.metrics$maf)
    colnames(maf) <- "maf"

    # Print out some statistics
    cat("  Reporting Minor Allele Frequency (MAF) by Locus\n")
    cat("  No. of loci =", nLoc(x2), "\n")
    cat("  No. of individuals =", nInd(x2), "\n")
    cat("    Minimum MAF: ", round(min(x2@other$loc.metrics$maf), 2), "\n")
    cat("    Maximum MAF: ", round(max(x2@other$loc.metrics$maf), 2), "\n")
    cat("    Average MAF: ", round(mean(x2@other$loc.metrics$maf), 3), "\n")
    cat("    Missing Rate Overall: ", round(sum(is.na(as.matrix(x2)))/(nLoc(x2) * nInd(x2)), 2), "\n\n")

    # Determine the loss of loci for a given threshold using quantiles
    quantile_res <- quantile(maf$maf, probs = seq(0, 1, 1/20))
    retained <- unlist(lapply(quantile_res, function(y) {
        res <- length(maf$maf[maf$maf >= y])
    }))
    pc.retained <- round(retained * 100/nLoc(x2), 1)
    filtered <- nLoc(x2) - retained
    pc.filtered <- 100 - pc.retained
    df <- data.frame(as.numeric(sub("%", "", names(quantile_res))), quantile_res, retained, pc.retained, filtered, pc.filtered)
    colnames(df) <- c("Quantile", "Threshold", "Retained", "Percent", "Filtered", "Percent")
    df <- df[order(-df$Quantile), ]
    df$Quantile <- paste0(df$Quantile, "%")
    rownames(df) <- NULL

    # testing which populations comply with thresholds
    popn.hold <- test_pop[which(test_pop$locs_per_pop >= loc.limit & test_pop$ind_per_pop >= ind.limit), "pop"]
    mafs_plots_print <- mafs_plots[popn.hold]

    if (length(popn.hold) > 1) {
        title.str <- "Minor Allele Frequency\nOverall"
        
        p_all <- ggplot(as.data.frame(maf), aes(x = maf)) + geom_histogram(bins = 50, color = plot_colours_all[1], fill = plot_colours_all[2]) + 
          xlab("Minor Allele Frequency") + ylab("Count") + xlim(0, maf.limit) + plot_theme + ggtitle(title.str)
        
        row_plots <- ceiling(length(popn.hold)/3) + 1
        p3 <- p_all + mafs_plots_print + plot_layout(ncol = 3, nrow =row_plots)
    }

    if (length(popn.hold) == 0) {
        if (verbose >= 1) {
            cat(important("  No populations met minimum limits on number of individuals or loci, reporting for overall\n"))
        }
        title.str <- "Minor Allele Frequency\nOverall"
        p3 <- ggplot(as.data.frame(maf), aes(x = maf)) + geom_histogram(bins = 50, color = plot_colours_all[1], fill = plot_colours_all[2]) + 
          xlab("Minor Allele Frequency") + ylab("Count") + xlim(0, maf.limit) + plot_theme + ggtitle(title.str)
    }

    if (length(popn.hold) == 1) {
        if (verbose >= 3) {
            cat(important("  Only one population met minimum limits on number of individuals or loci\n"))
        }
        title.str <- paste("Minor Allele Frequency\n", popn.hold)
        p3 <- ggplot(as.data.frame(maf), aes(x = maf)) + geom_histogram(bins = 50, color = plot_colours_all[1], fill = plot_colours_all[2]) + 
          xlab("Minor Allele Frequency") + ylab("Count") + xlim(0, maf.limit) + plot_theme + ggtitle(title.str)
    }

    if (nPop(x2) == 1) {
        if (verbose >= 1) {
            cat(important("  Only one population specified\n"))
        }
        title.str <- paste("Minor Allele Frequency\n", pop(x2)[1])
        p3 <- ggplot(as.data.frame(maf), aes(x = maf)) + geom_histogram(bins = 50, color = plot_colours_all[1], fill = plot_colours_all[2]) + 
          xlab("Minor Allele Frequency") + ylab("Count") + xlim(0, maf.limit) + plot_theme + ggtitle(title.str)
    }
    
    # printing outputs
    suppressWarnings(print(p3))
    print(df)

    # creating temp file names
    temp_plot <- tempfile(pattern = "plot_")
    temp_table <- tempfile(pattern = "table_")

    # saving to tempdir
    saveRDS(p3, file = temp_plot)
    saveRDS(df, file = temp_table)

    # ADD TO HISTORY

    nh <- length(x@other$history)
    x@other$history[[nh + 1]] <- c(match.call(), temp_plot, temp_table)

    # FLAG SCRIPT END

    if (verbose >= 1) {
        cat(report("\n\nCompleted:", funname, "\n\n"))
    }

    cat(important(strwrap("Plots and table were saved to the temporal directory (tempdir) and can be accesed with the function gl.access.report(). Note that they can be accessed only in the current R session because tempdir is cleared each time that the R session is closed.\n\n")))

    invisible(x)
}
