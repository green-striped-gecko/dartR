#' @name gl.report.hamming
#' @title Calculates the pairwise Hamming distance between DArT trimmed DNA
#' sequences
#'
#' @description Hamming distance is calculated as the number of base differences
#' between two sequences which can be expressed as a count or a proportion.
#' Typically, it is calculated between two sequences of equal length. In the
#' context of DArT trimmed sequences, which differ in length but which are
#' anchored to the left by the restriction enzyme recognition sequence, it is
#' sensible to compare the two trimmed sequences starting from immediately after
#' the common recognition sequence and terminating at the last base of the
#' shorter sequence.
#'
#' @param x Name of the genlight object containing the SNP data [required].
#' @param rs Number of bases in the restriction enzyme recognition sequence
#' [default 5].
#' @param threshold Minimum acceptable base pair difference for display on the
#' boxplot and histogram [default 3].
#' @param taglength Typical length of the sequence tags [default 69].
#' @param plot.out Specify if plot is to be produced [default TRUE].
#' @param plot_theme Theme for the plot. See Details for options
#' [default theme_dartR()].
#' @param plot_colors List of two color names for the borders and fill of the
#' plots [default two_colors].
#' @param probar If TRUE, then a progress bar is displayed on long loops
#' [default TRUE].
#' @param save2tmp If TRUE, saves any ggplots and listings to the session
#' temporary directory (tempdir) [default FALSE].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#'
#' @details The function \code{\link{gl.filter.hamming}} will filter out one of
#' two loci if their Hamming distance is less than a specified percentage
#'
#' Hamming distance can be computed by exploiting the fact that the dot product
#' of two binary vectors x and (1-y) counts the corresponding elements that are
#' different between x and y. This approach can also be used for vectors that
#' contain more than two possible values at each position (e.g. A, C, T or G).
#'
#' If a pair of DNA sequences are of differing length, the longer is truncated.
#'
#' The algorithm is that of Johann de Jong
#' \url{https://johanndejong.wordpress.com/2015/10/02/faster-hamming-distance-in-r-2/}
#' as implemented in \code{\link{utils.hamming}}
#'
#'  Plots and table are saved to the session's temporary directory (tempdir)
#'
#'  Examples of other themes that can be used can be consulted in \itemize{
#'  \item \url{https://ggplot2.tidyverse.org/reference/ggtheme.html} and \item
#'  \url{https://yutannihilation.github.io/allYourFigureAreBelongToUs/ggthemes/}
#'  }
#'
#' @return Returns unaltered genlight object
#' @author Custodian: Arthur Georges -- Post to 
#' \url{https://groups.google.com/d/forum/dartr}
#'
#' @examples
#' gl.report.hamming(testset.gl[,1:100])
#' gl.report.hamming(testset.gs[,1:100])
#'
#' @seealso \code{\link{gl.filter.hamming}}
#'
#' @family report functions
#' @importFrom stats sd
#' @import patchwork
#' @export

gl.report.hamming <- function(x,
                              rs = 5,
                              threshold = 3,
                              taglength = 69,
                              plot.out = TRUE,
                              plot_theme = theme_dartR(),
                              plot_colors = two_colors,
                              probar = FALSE,
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
    
    if (length(x@other$loc.metrics$TrimmedSequence) == 0) {
        stop(error("Fatal Error: Data must include Trimmed Sequences\n"))
    }
    
    if (rs < 0 | rs > taglength) {
        stop(
            error(
                "Fatal Error: Length of restriction enzyme recognition sequence
                must be greater than zero, and less that the maximum length of a
                sequence tag; usually it is less than 9\n"
            )
        )
    }
    
    if (nLoc(x) == 1) {
        stop(error("Fatal Error: Data must include more than one locus\n"))
    }
    
    # DO THE JOB
    
    s <- as.character(x@other$loc.metrics$TrimmedSequence)
    tld <- threshold / (taglength - rs)
    
    if (probar) {
        pb <-
            txtProgressBar(
                min = 0,
                max = 1,
                style = 3,
                initial = 0,
                label = "Working ...."
            )
        getTxtProgressBar(pb)
    }
    
    if (verbose >= 3) {
        cat(
            report(
                "  Hamming distance ranges from zero (sequence identity) to 1 
                (no bases shared at any position)\n"
            )
        )
    }
    if (verbose >= 2) {
        cat(
            report(
                "  Calculating pairwise Hamming distances between trimmed 
                Reference sequence tags\n"
            )
        )
    }
    
    count <-0
    nL <- nLoc(x)
    d <- rep(NA, (((nL - 1) * nL) / 2))
    
    for (i in 1:(nL - 1)) {
        for (j in ((i + 1):nL)) {
            count <- count + 1
            d[count] <- utils.hamming(s[i], s[j], r = rs)
        }
        if (probar) {
            setTxtProgressBar(pb, i / (nL - 1))
        }
    }
    
    # get title for plots
    if (datatype == "SNP") {
        title <-
            paste0("SNP data (DArTSeq)\nPairwise Hamming Distance between 
                   sequence tags")
    } else {
        title <-
            paste0(
                "Fragment P/A data (SilicoDArT)\nPairwise Hamming Distance 
                between sequence tags"
            )
    }
    
    if (verbose >= 2) {
        if (plot.out) {
            cat(
                report(
                    "  Plotting boxplot and histogram of Hamming distance, 
                    showing a threshold of",
                    threshold,
                    "bp [HD",
                    round(tld, 2),
                    "]\n"
                )
            )
        }
    }
    
    # Boxplot
    p1 <-
        ggplot(as.data.frame(d), aes(y = d)) +
      geom_boxplot(color = plot_colors[1], fill = plot_colors[2]) + 
      geom_hline(yintercept = tld,color = "red", size = 1) + 
      coord_flip() + 
      plot_theme + 
      xlim(range = c(-1, 1)) + 
      ylim(0, 1) +
      ylab(" ") + 
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
      ggtitle(title)
    
    # Histogram
    p2 <-
        ggplot(as.data.frame(d), aes(x = d)) + 
      geom_histogram(bins = 50, color = plot_colors[1],fill = plot_colors[2]) +
      geom_vline(xintercept = tld,color = "red",size = 1) + 
      coord_cartesian(xlim = c(0, 1)) +
      xlab("Hamming distance") +
      ylab("Count") + 
      annotate(geom = "text", 
               x = tld + 0.2, 
               y = max(graphics::hist(d, breaks = seq(0, 1, by = 1 / 50),
                                      plot = FALSE)$counts) * 0.75,
               label = paste("Threshold of\n", threshold,
                             "bp [HD", round(tld, 2), "]")) + 
      plot_theme
    
    cat("    No. of loci =", nLoc(x), "\n")
    cat("    No. of individuals =", nInd(x), "\n")
    cat("    Minimum Hamming distance: ", round(min(d), 2), "\n")
    cat("    Maximum Hamming distance: ", round(max(d), 2), "\n")
    cat(paste0(
        "    Mean Hamming Distance ",
        round(mean(d), 2),
        "+/-",
        round(sd(d), 3),
        " SD\n"
    ))
    n.outliers <- sum(d <= (threshold / (taglength - rs)))
    cat(
        "    No. of pairs with Hamming Distance less than or equal to",
        threshold,
        "base pairs:",
        n.outliers,
        "\n\n"
    )
    
    # Determine the loss of loci for a given threshold using quantiles
    nl <- nLoc(x)
    quantile_res <- quantile(d, probs = seq(0, 1, 1 / 20),type=1)
    retained <- unlist(lapply(quantile_res, function(y) {
        res <- sum(d >= y)
    }))
    pc.retained <- round(retained * 100 / ((((
        nL - 1
    ) * nL) / 2)), 1)
    filtered <- ((((nL - 1) * nL) / 2)) - retained
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
    
    # PRINTING OUTPUTS
    if (plot.out) {
        # using package patchwork
        p3 <- (p1 / p2) + plot_layout(heights = c(1, 4))
        print(p3)
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
                    "  NOTE: Retrieve output files from tempdir using 
                    gl.list.reports() and gl.print.reports()\n"
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
