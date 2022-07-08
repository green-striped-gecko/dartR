#' @name gl.report.parent.offspring
#' @title Identifies putative parent offspring within a population
#' @description
#' This script examines the frequency of pedigree inconsistent loci, that is,
#' those loci that are homozygotes in the parent for the reference allele, and
#' homozygous in the offspring for the alternate allele. This condition is not
#' consistent with any pedigree, regardless of the (unknown) genotype of the
#' other parent. The pedigree inconsistent loci are counted as an indication of
#' whether or not it is reasonable to propose the two individuals are in a
#' parent-offspring relationship.
#' @param x Name of the genlight object containing the SNP genotypes [required].
#' @param min.rdepth Minimum read depth to include in analysis [default 12].
#' @param min.reproducibility Minimum reproducibility to include in analysis
#' [default 1].
#' @param range Specifies the range to extend beyond the interquartile range for
#' delimiting outliers [default 1.5 interquartile ranges].
#' @param plot.out Creates a plot that shows the sex linked markers
#' [default TRUE].
#' @param plot_theme Theme for the plot. See Details for options
#'  [default theme_dartR()].
#' @param plot_colors List of two color names for the borders and fill of the
#'  plots [default two_colors].
#' @param save2tmp If TRUE, saves any ggplots and listings to the session
#' temporary directory (tempdir) [default FALSE].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#' @details
#' If two individuals are in a parent offspring relationship, the true number of
#' pedigree inconsistent loci should be zero, but SNP calling is not infallible.
#' Some loci will be miss-called. The problem thus becomes one of determining
#' if the two focal individuals have a count of pedigree inconsistent loci less
#' than would be expected of typical unrelated individuals. There are some quite
#' sophisticated software packages available to formally apply likelihoods to
#' the decision, but we use a simple outlier comparison.
#'
#' To reduce the frequency of miss-calls, and so emphasize the difference
#' between true parent-offspring pairs and unrelated pairs, the data can be
#' filtered on read depth.
#'
#' Typically minimum read depth is set to 5x, but you can examine the
#' distribution of read depths with the function \code{\link{gl.report.rdepth}}
#' and push this up with an acceptable loss of loci. 12x might be a good minimum
#' for this particular analysis. It is sensible also to push the minimum
#' reproducibility up to 1, if that does not result in an unacceptable loss of
#' loci. Reproducibility is stored in the slot \code{@other$loc.metrics$RepAvg}
#' and is defined as the proportion of technical replicate assay pairs for which
#' the marker score is consistent. You can examine the distribution of
#'  reproducibility with the function \code{\link{gl.report.reproducibility}}.
#'
#' Note that the null expectation is not well defined, and the power reduced, if
#' the population from which the putative parent-offspring pairs are drawn
#' contains many sibs. Note also that if an individual has been genotyped twice
#' in the dataset, the replicate pair will be assessed by this script as being
#' in a parent-offspring relationship.
#'
#' The function \code{\link{gl.filter.parent.offspring}} will filter out those
#' individuals in a parent offspring relationship.
#'
#' Note that if your dataset does not contain RepAvg or rdepth among the locus
#' metrics, the filters for reproducibility and read depth are no used.
#'
#'\strong{ Function's output }
#'
#'  Plots and table are saved to the temporal directory (tempdir) and can be
#'  accessed with the function \code{\link{gl.print.reports}} and listed with
#'  the function \code{\link{gl.list.reports}}. Note that they can be accessed
#'  only in the current R session because tempdir is cleared each time that the
#'   R session is closed.
#'
#'  Examples of other themes that can be used can be consulted in \itemize{
#'  \item \url{https://ggplot2.tidyverse.org/reference/ggtheme.html} and \item
#'  \url{https://yutannihilation.github.io/allYourFigureAreBelongToUs/ggthemes/}
#'  }
#' @return A set of individuals in parent-offspring relationship. NULL if no
#' parent-offspring relationships were found.
#' @author Custodian: Arthur Georges (Post to
#' \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' out <- gl.report.parent.offspring(testset.gl[1:10])
#' @seealso \code{\link{gl.list.reports}}, \code{\link{gl.report.rdepth}} ,
#'  \code{\link{gl.print.reports}},\code{\link{gl.report.reproducibility}},
#'  \code{\link{gl.filter.parent.offspring}}
#' @family reporting functions
#' @importFrom stats median IQR
#' @import patchwork
#' @export

gl.report.parent.offspring <- function(x,
                                       min.rdepth = 12,
                                       min.reproducibility = 1,
                                       range = 1.5,
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
    
    # DO THE JOB
    
    # Generate null expectation for pedigree inconsistency, and outliers
    if (verbose >= 2) {
        cat(
            report(
                "  Generating null expectation for distribution of counts of pedigree incompatibility\n"
            )
        )
    }
    # Assign individuals as populations
    pop(x) <- x$ind.names
    # Filter stringently on reproducibility to minimize miscalls
    if (is.null(x@other$loc.metrics$RepAvg)) {
        cat(
            warn(
                "  Dataset does not include RepAvg among the locus metrics, therefore the reproducibility filter was not used\n"
            )
        )
    } else {
        x <-
            gl.filter.reproducibility(x, threshold = min.reproducibility, verbose = 0)
    }
    # Filter stringently on read depth, to further minimize miscalls
    if (is.null(x@other$loc.metrics$rdepth)) {
        cat(
            warn(
                "  Dataset does not include rdepth among the locus metrics, therefore the read depth filter was not used\n"
            )
        )
    } else {
        x <- gl.filter.rdepth(x, lower = min.rdepth, verbose = 0)
    }
    
    # Preliminaries before for loops
    
    x2 <- as.matrix(x)
    split_vectors <- lapply(seq_len(nrow(x2)), function(i)
        x2[i,])
    names(split_vectors) <- popNames(x)
    
    fun <- function(x, y) {
        vect <- (x * 10) + y
        homalts <- sum(vect == 2 | vect == 20, na.rm = T)
    }
    
    count <- sapply(split_vectors, function(vect1) {
        sapply(split_vectors, function(vect2) {
            fun(vect1, vect2)
        })
    })
    
    counts <- count[lower.tri(count, diag = FALSE)]
    
    # Prepare for plotting
    
    if (verbose >= 2) {
        cat(
            report(
                "  Identifying outliers with lower than expected counts of pedigree inconsistencies\n"
            )
        )
    }
    title <-
        paste0("SNP data (DArTSeq)\nCounts of pedigree incompatible loci per pair")
    
    counts_plot <- as.data.frame(counts)
    
    # Boxplot
    p1 <-
        ggplot(counts_plot, aes(y = counts)) + geom_boxplot(color = plot_colors[1], fill = plot_colors[2]) + coord_flip() + plot_theme + xlim(range = c(-1,
                                                                                                                                                        1)) + ylim(min(counts), max(counts)) + ylab(" ") + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) + ggtitle(title)
    
    outliers_temp <- ggplot_build(p1)$data[[1]]$outliers[[1]]
    
    lower.extremes <-
        outliers_temp[outliers_temp < stats::median(counts)]
    if (length(lower.extremes) == 0) {
        outliers <- NULL
    } else {
        outliers <- data.frame(Outlier = lower.extremes)
    }
    
    # Ascertain the identity of the pairs
    if (verbose >= 2) {
        cat(report("  Identifying outlying pairs\n"))
    }
    if (length(lower.extremes) > 0) {
        tmp <- count
        tmp[lower.tri(tmp)] <-t(tmp)[lower.tri(tmp)]
        for (i in 1:length(outliers$Outlier)) {
            # Identify
            tmp2 <- tmp[tmp == outliers$Outlier[i]]
            outliers$ind1[i] <- popNames(x)[!is.na(tmp2)][1]
            outliers$ind2[i] <- popNames(x)[!is.na(tmp2)][2]
            # Z-scores
            zscore <-
                (mean(count, na.rm = TRUE) - outliers$Outlier[i]) / sd(count, na.rm = TRUE)
            outliers$zscore[i] <- round(zscore, 2)
            outliers$p[i] <-
                round(pnorm(
                    mean = mean(count, na.rm = TRUE),
                    sd = sd(count, na.rm = TRUE),
                    q = outliers$zscore[i]
                ), 4)
        }
        # ordering by number of outliers
        outliers <- outliers[order(outliers, decreasing = T),]
        # removing duplicated values
        outliers <- outliers[!duplicated(outliers),]
        # removing NAs
        outliers <- outliers[complete.cases(outliers),]
    }
    
    # Extract the quantile threshold
    iqr <- stats::IQR(counts, na.rm = TRUE)
    qth <- quantile(counts, 0.25, na.rm = TRUE)
    cutoff <- qth - iqr * range
    
    # Histogram
    p2 <-
        ggplot(counts_plot, aes(x = counts)) + geom_histogram(bins = 50,
                                                              color = plot_colors[1],
                                                              fill = plot_colors[2]) + geom_vline(xintercept = cutoff,
                                                                                                  color = "red",
                                                                                                  size = 1) + coord_cartesian(xlim = c(min(counts), max(counts))) + xlab("No. Pedigree incompatible") + ylab("Count") +
        plot_theme
    
    # Output the outlier loci
    if (length(lower.extremes) == 0) {
        cat(important("  No outliers detected\n"))
    }
    if (length(lower.extremes) > 0) {
        outliers <- outliers[order(outliers$Outlier),]
        if (verbose >= 3) {
            print(outliers)
        }
    }
    
    df <- outliers
    # PRINTING OUTPUTS
    if (plot.out) {
        # using package patchwork
        p3 <- (p1 / p2) + plot_layout(heights = c(1, 4))
        print(p3)
    }
    
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
