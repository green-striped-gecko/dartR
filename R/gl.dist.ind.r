#' @name gl.dist.ind
#' @title Calculates a distance matrix for individuals defined in a genlight object
#' @description
#' This script calculates various distances between individuals based on allele
#'  frequencies or presence-absence data 
#' @details
#' The distance measure for SNP genotypes can be one of:
#' \itemize{
#'  \item Euclidean Distance
#'  \item Scaled Euclidean Distance
#'  \item Simple Mismatch Distance
#'  \item Absolute Mismatch Distance
#'  \item Czekanowski (Manhattan) Distance
#'  }
#'
#' The distance measure for Sequence Tag Presence/Absence data (binary) can be one of:
#' \itemize{
#'  \item Euclidean Distance
#'  \item Scaled Euclidean Distance
#'  \item Simple Matching Distance
#'  \item Jaccard Distance
#'  \item Bray-Curtis Distance
#'  }
#'
#' Refer to the dartR Technical Note on Distances in Genetics.
#'
#' @param x Name of the genlight containing the SNP genotypes or presence-absence data [required].
#' @param method Specify distance measure [SNP: Euclidean; P/A: SMD].
#' @param scale If TRUE, the distances are scaled to fall in the range [0,1] [default TRUE]
#' @param swap If TRUE and working with presence-absence data, then presence 
#' (no disrupting mutation) is scored as 0 and absence (presence of a disrupting 
#' mutation) is scored as 1 [default FALSE].
#' @param output Specify the format and class of the object to be returned, 
#' 'dist' for a object of class dist, 'matrix' for an object of class matrix [default "dist"].
#' @param plot.out If TRUE, display a histogram and a boxplot of the genetic distances [TRUE].
#' @param plot_theme User specified theme [default theme_dartR].
#' @param plot_colors Vector with two color names for the borders and fill [default two_colors].
#' @param save2tmp If TRUE, saves any ggplots to the session temporary directory [default FALSE].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log ; 3, progress and results summary; 5, full report
#' [default 2 or as specified using gl.set.verbosity].
#' @return An object of class 'matrix' or dist' giving distances between individuals
#' @importFrom ape dist.gene
#' @importFrom stats dist
#' @export
#' @author Author(s): Arthur Georges. Custodian: Arthur Georges -- Post to #' \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' D <- gl.dist.ind(testset.gl, method='euclidean')
#' D <- gl.dist.ind(testset.gs, method='Jaccard',swap=TRUE)

gl.dist.ind <- function(x,
                        method = NULL,
                        scale = FALSE,
                        swap=FALSE,
                        output="dist",
                        plot.out = TRUE,
                        plot_theme = theme_dartR(),
                        plot_colors = two_colors,
                        save2tmp = FALSE,
                        verbose = NULL) {
    
    # CHECK IF PACKAGES ARE INSTALLED
    pkg <- "rrBLUP"
    if (!(requireNamespace(pkg, quietly = TRUE))) {
        stop(error(
            "Package ",
            pkg,
            " needed for this function to work. Please install it."
        ))
    }
    
    pkg <- "poppr"
    if (!(requireNamespace(pkg, quietly = TRUE))) {
        stop(error(
            "Package ",
            pkg,
            " needed for this function to work. Please install it."
        ))
    }
    
    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)
    
    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     build = "Jody",
                     verbosity = verbose)
    
    # CHECK DATATYPE
    datatype <-
        utils.check.datatype(x,
                             accept = c("SNP", "SilicoDArT"),
                             verbose = verbose)
    
    # FUNCTION SPECIFIC ERROR CHECKING
    
    if (is.null(method) && datatype == "SNP") {
        method <- "Euclidean"
    }
    if (is.null(method) && datatype == "SilicoDArT") {
        method <- "Simple"
    }
    method <- tolower(method)
    
    if (!(
        method %in% c(
            "euclidean",
            "locus.count",
            "allele.count",
            "relatedness",
            "simple",
            "jaccard",
            "bray-curtis",
            "sorenson",
            "czekanowski",
            "phi"
        )
    )) {
        cat(
            warn(
                " Warning: Method not in the list of options, set to Euclidean for SNP data; Simple Matching for Tag P/A data\n"
            )
        )
        if (datatype == "SNP") {
            method <- "euclidean"
        }
        if (datatype == "SilicoDArT") {
            method <- "simple"
        }
    }
    
    # DO THE JOB
    
    if (datatype == "SNP") {
        # Calculate euclidean distance using dist {adegenet}
        if (method == "euclidean") {
            dd <- stats::dist(x)
            if (verbose >= 2) {
                cat(report(
                    "  Calculating Euclidean Distances between individuals\n"
                ))
            }
        }
        # Calculate the number of loci that are different between individuals using dist.gene {ape}
        if (method == "locus.count") {
            dd <-
                ape::dist.gene(
                    as.matrix(x),
                    method = "pairwise",
                    pairwise.deletion = FALSE,
                    variance = FALSE
                )
            if (verbose >= 2) {
                cat(report(
                    "  Calculating number of loci for which individuals differ\n"
                ))
            }
        }
        
        # Calculate the number of allelic differences between individuals using dist.gene {poppr}
        if (method == "allele.count") {
            dd <- poppr::diss.dist(gl2gi(x), percent = FALSE, mat = FALSE)
            if (verbose >= 2) {
                cat(
                    report(
                        "  Calculating number of allelic differences between individuals\n"
                    )
                )
            }
        }
        
        # Calculate the genetic relatedness G matrix
        if (method == "relatedness") {
            dd <- rrBLUP::A.mat(as.matrix(x) - 1)
            if (verbose >= 2) {
                cat(report(
                    "  Calculating relatedness among individuals (G matrix)\n"
                ))
            }
        }
        dd <- as.dist(dd)
        
        # # Revert to original order ord <- rank(pop(x)) mat <- as.matrix(dd)[ord, ord] dd <- as.dist(mat)
        mat <- as.matrix(dd)
    }
    
    if (datatype == "SilicoDArT") {
        if (method == "euclidean" && scale==FALSE) {
            if (verbose >= 2) {
                cat(report("  Calculating the Unscaled Euclidean Distances\n"))
            }
        }
        if (method == "euclidean" && scale==TRUE) {
            if (verbose >= 2) {
                cat(report("  Calculating the Scaled Euclidean Distances\n"))
            }
        }
        if (method == "simple") {
            if (verbose >= 2) {
                cat(report("  Calculating the Simple Matching Distances\n"))
            }
        }
        if (method == "jaccard") {
            if (verbose >= 2) {
                cat(report("  Calculating distances based on the Jaccard Coefficient\n"))
            }
        }
        if (method == "bray-curtis") {
            if (verbose >= 2) {
                cat(report(
                    "  Calculating the Bray-Curtis Distance\n"
                ))
            }
        }
        # if (method == "phi") {
        #     if (verbose >= 2) {
        #         cat(report(
        #             "  Calculating the Pearson Phi Index (= Binary correlation\n"
        #         ))
        #     }
        # }
        mat <- utils.dist.binary(x, 
                                 method = method, 
                                 swap=swap, 
                                 output="matrix", 
                                 scale = scale, 
                                 verbose = 0)
        dd <- as.dist(mat)
    }
    
    # PLOT
    if (plot.out) {
        if (datatype == "SNP") {
            title_plot <-
                paste0("SNP data (DArTSeq)\nInter-individual ",
                       method,
                       " distance")
        } else {
            if(method=="euclidean" && scale == TRUE){
                title_plot <-
                paste0(
                    "Presence[1]/Absence[0] data (SilicoDArT)\nInter-individual scaled ",
                    method,
                    " distance"
                )
            } else {
                if(swap==TRUE){
                    title_plot <- paste0(
                        "Presence[0]/Absence[1] data (SilicoDArT swapped)\nInter-individual ",
                        method, " distance")
                } else {
                    title_plot <- paste0(
                        "Presence[1]/Absence[0] data (SilicoDArT)\nInter-individual ",
                        method, " distance")
                }
            }
        }
        values <- NULL
        df_plot <- data.frame(values = as.vector(mat))
        # colnames(df_plot) <- 'values'
        
        # Boxplot
        p1 <-
            ggplot(df_plot, aes(y = values)) + 
            geom_boxplot(color = plot_colors[1], 
            fill = plot_colors[2]) + 
            coord_flip() + 
            plot_theme + 
            xlim(range = c(-1,1)) + 
            ylim(min(df_plot$values, na.rm = TRUE),
            max(df_plot$values, na.rm = TRUE)) + 
            ylab(" ") + 
            theme(axis.text.y = element_blank(), 
                  axis.ticks.y = element_blank()) + 
            ggtitle(title_plot)
        
        # Histogram
        p2 <-
            ggplot(df_plot, aes(x = values)) + 
            geom_histogram(bins = 100,
                           color = plot_colors[1],
                           fill = plot_colors[2]) + 
            xlim(min(df_plot$values, na.rm = TRUE), 
                 max(df_plot$values, na.rm = TRUE)) + 
            xlab("Distance") + 
            ylab("Count") + 
            plot_theme
        
        # PRINTING OUTPUTS
        if (plot.out) {
            # using package patchwork
            p3 <- (p1 / p2) + plot_layout(heights = c(1, 4))
            suppressWarnings(print(p3))
        }
    }
    
    # SUMMARY Print out some statistics
    if (verbose >= 3) {
        cat("  Reporting inter-individual distances\n")
        cat("  Distance measure:", method, "\n")
        cat("    No. of populations =", nPop(x), "\n")
        cat("    Average no. of individuals per population =",
            round(nInd(x) / nPop(x), 1),
            "\n")
        cat("    No. of loci =", nLoc(x), "\n")
        cat("    Minimum Distance: ", round(min(dd), 2), "\n")
        cat("    Maximum Distance: ", round(max(dd), 2), "\n")
        cat("    Average Distance: ", round(mean(dd), 3), "\n\n")
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
        saveRDS(list(match_call, dd), file = temp_table)
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
    
    if(output=="matrix"){
        if(verbose >= 2){cat(report("  Returning a square matrix\n"))}
        return(mat)
        if (verbose > 0) {
            cat(report("Completed:", funname, "\n"))
        }
    } else {
        if(verbose >= 2){cat(report("  Returning a stat::dist object\n"))}
        return(dd)
        if (verbose > 0) {
            cat(report("Completed:", funname, "\n"))
        }
    }
}
