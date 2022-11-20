#' @name gl.report.bases
#'
#' @title Reports summary of base pair frequencies
#'
#' @description
#' This script calculates the frequencies of the four DNA nucleotide bases:
#' adenine (A), cytosine (C), 'guanine (G) and thymine (T), and the frequency of
#' transitions (Ts) and transversions (Tv) in a DArT genlight object.
#'
#' @param x Name of the genlight object containing the SNP or presence/absence
#' (SilicoDArT) data [required].
#' @param plot.out If TRUE, histograms of base composition are produced
#' [default TRUE].
#' @param plot_theme Theme for the plot. See Details for options
#' [default theme_dartR()].
#' @param plot_colors List of two color names for the borders and fill of the
#'  plots [default two_colors].
#' @param save2tmp If TRUE, saves any ggplots and listings to the session
#' temporary directory (tempdir) [default FALSE]
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#'  [default NULL, unless specified using gl.set.verbosity]
#'
#' @details 
#' The script checks first if trimmed sequences are included in the
#' locus metadata (@@other$loc.metrics$TrimmedSequence), and if so, tallies up
#' the numbers of A, T, G and C bases. Only the reference state at the SNP locus
#' is counted. Counts of transitions (Ts) and transversions (Tv) assume that
#' there is no directionality, that is C->T is the same as T->C, because the
#' reference state is arbitrary.
#'
#' For presence/absence data (SilicoDArT), it is not possible to count
#' transversions or transitions or transversions/transitions ratio because the
#'  SNP data is not available, only a single sequence tag.
#'
#' Examples of other themes that can be used can be consulted in \itemize{
#'  \item \url{https://ggplot2.tidyverse.org/reference/ggtheme.html} and \item
#'  \url{https://yutannihilation.github.io/allYourFigureAreBelongToUs/ggthemes/}
#'  }
#'
#' @return The unchanged genlight object
#' @author Custodian: Arthur Georges -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#'
#' @examples
#' # SNP data
#'   out <- gl.report.bases(testset.gl)
#'   #' # Tag P/A data
#'   out <- gl.report.bases(testset.gs)
#'
#' @family report functions
#' @import stringr
#' @import patchwork
#' @export


gl.report.bases <- function(x,
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
    
    if (!any(names(x@other$loc.metrics) == "TrimmedSequence")) {
        stop(error(
            "  Fatal Error: Dataset does not include variable 
            TrimmedSequence!\n"
        ))
    }
    
    # DO THE JOB
    
    # Count up the number of bases, and the number of each of ATGC, and other
    if (verbose >= 2) {
        cat(report("  Counting the bases\n"))
    }
    
    A <- sum(stringr::str_count(x@other$loc.metrics$TrimmedSequence, "A"))
    G <- sum(stringr::str_count(x@other$loc.metrics$TrimmedSequence, "G"))
    C <- sum(stringr::str_count(x@other$loc.metrics$TrimmedSequence, "C"))
    T <- sum(stringr::str_count(x@other$loc.metrics$TrimmedSequence, "T"))
    total <- sum(stringr::str_length(x@other$loc.metrics$TrimmedSequence))
    total.ATGC <- sum(A, G, C, T)
    if (verbose >= 2) {
        if (total != total.ATGC) {
            cat(warn("  Warning: Codes other than A, T, G and C present\n"))
        }
    }
    other <- total - total.ATGC
    other <- other * 100 / total
    A <- A * 100 / total
    G <- G * 100 / total
    T <- T * 100 / total
    C <- C * 100 / total
    
    # Calculate the fragment lengths
    mn <- mean(stringr::str_length(x@other$loc.metrics$TrimmedSequence))
    mx <- max(stringr::str_length(x@other$loc.metrics$TrimmedSequence))
    mi <- min(stringr::str_length(x@other$loc.metrics$TrimmedSequence))
    
    if (datatype == "SNP") {
        # Extract the SNPs
        matrix <- stringr::str_split_fixed(x@other$loc.metrics$SNP, ":", 2)
        state.change <- matrix[, 2]
        
        if (verbose >= 2) {
            cat(report("  Counting Transitions and Transversions\n"))
        }
        
        # Sum the transitions and transversions
        tv <-
            sum(str_count(state.change, "A>C")) + 
          sum(stringr::str_count(state.change, "C>A")) +
          sum(stringr::str_count(state.change, "G>T")) +
            sum(stringr::str_count(state.change, "T>G")) + 
          sum(stringr::str_count(state.change, "A>T")) + 
          sum(stringr::str_count(state.change, "T>A")) + 
          sum(stringr::str_count(state.change, "G>C")) + 
          sum(stringr::str_count(state.change, "C>G"))
        
        ts <-
            sum(stringr::str_count(state.change, "A>G")) + 
          sum(stringr::str_count(state.change, "G>A")) + 
          sum(stringr::str_count(state.change, "C>T")) + 
          sum(stringr::str_count(state.change, "T>C"))
        
        if (verbose >= 2) {
            if (ts + tv != length(x@other$loc.metrics$TrimmedSequence)) {
                cat(
                    warn(
                        "  Warning: Sum of transitions plus transversions does 
                        not equal number of loci.\n"
                    )
                )
            }
        }
        ts <- ts * 100 / length(x@other$loc.metrics$TrimmedSequence)
        tv <- tv * 100 / length(x@other$loc.metrics$TrimmedSequence)
        ratio <- ts / tv
    }
    
    # PRINTING OUTPUTS
    cat(paste(
        "  Average trimmed sequence length:",
        round(mn, digits = 1),
        "(",
        mi,
        "to",
        mx,
        ")"
    ),
    "\n")
    cat(paste(
        "  Total number of trimmed sequences:",
        length(x@other$loc.metrics$TrimmedSequence)
    ), "\n")
    cat("  Base frequencies (%)\n")
    cat(paste("    A:", round(A, 2)), "\n")
    cat(paste("    G:", round(G, 2)), "\n")
    cat(paste("    T:", round(T, 2)), "\n")
    cat(paste("    C:", round(C, 2)), "\n\n")
    
    if (datatype == "SilicoDArT") {
        if (verbose >= 2) {
            cat(
                important(
                    "  Tag P/A data (SilicoDArT), transition/transversions 
                    cannot be calculated\n"
                )
            )
        }
        tv <- NA
        ts <- NA
    } else {
        cat(paste("  Transitions  :", round(ts, 2), "\n"))
        cat(paste("  Transversions:", round(tv, 2), "\n"))
        cat(paste("  tv/ts ratio:", round(ratio, 4), "\n\n"))
    }
    
    if (plot.out) {
        if (datatype == "SNP") {
            title <- paste0("SNP: Base Frequencies")
        } else {
            title <- paste0("Tag P/A: Base Frequencies")
        }
        
        bases <- c("A", "C", "T", "G")
        freq <- round(c(A, C, T, G), 1)
        df <- data.frame(bases = bases, freq = freq)
        
        p1 <-
            ggplot(data = df, aes(x = bases, y = freq)) +
          geom_bar(stat="identity",color=plot_colors[1],fill=plot_colors[2]) + 
          xlab("Bases") +
            ylab("Percent Frequency") + 
          ggtitle(title) +
          plot_theme
        
        if (datatype == "SNP") {
            bases <- c("Ts", "Tv")
            freq <- round(c(ts, tv), 1)
            df2 <- data.frame(bases = bases, freq = freq)
            
            p2 <-
                ggplot(data = df2, aes(x = bases, y = freq)) +
          geom_bar(stat="identity",color=plot_colors[1],fill=plot_colors[2]) +
                xlab("Mutation Type") + 
              ylab("Percent Frequency") + 
              ggtitle(paste("SNP: Ts/Tv Rates [ratio =",round(ratio,2),"]")) +
              plot_theme
            
            p3 <- (p1 / p2)  # Using package patchwork
        } else {
            p3 <- p1
        }
        print(p3)
    }
    
    # Create return list
    if (verbose >= 2) {
        cat(
            report(
                "  Returning a list containing
[[1]] $freq -- the table of base frequencies and transition/transversion ratios;
[[2]] $plotbases -- ggplot bargraph of base frequencies;
[[3]] $plottstv -- ggplot bargraph of transitions and transversions."
            )
        )
    }
    
    out <-
        c(round(A, 2),
          round(G, 2),
          round(T, 2),
          round(C, 2),
          round(tv, 2),
          round(ts, 2))
    names(out) <- c("A", "G", "T", "C", "tv", "ts")
    
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
