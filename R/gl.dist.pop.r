#' @name gl.dist.pop
#' @title Calculates a distance matrix for populations with SNP genotypes in a
#' genlight object
#' @description
#' This script calculates various distances between populations based on allele
#' frequencies (SNP genotypes) or frequency of presences in presence-absence data 
#' (Euclidean and Fixed-diff distances only). 
#' @details
#' The distance measure can be one of 'euclidean', 'fixed-diff', 'Reynolds',
#' 'nei' and 'chord'. Refer to the documentation of functions
#'   described in the the dartR Distance Analysis tutorial for algorithms
#'   and definitions.
#'
#' @param x Name of the genlight containing the SNP genotypes [required].
#' @param method Specify distance measure [default euclidean].
#' @param plot.out If TRUE, display a histogram of the genetic distances, and a
#'  whisker plot [default TRUE].
#' @param scale If TRUE and method='Euclidean', the distance will be scaled to 
#'  fall in the range [0,1] [default FALSE].
#' @param output Specify the format and class of the object to be returned, 
#' dist for a object of class dist, matrix for an object of class matrix [default "dist"].
#' @param plot_theme User specified theme [default theme_dartR()].
#' @param plot_colors Vector with two color names for the borders and fill
#' [default two_colors].
#' @param save2tmp If TRUE, saves any ggplots and listings to the session
#'  temporary directory (tempdir) [default FALSE].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#'  progress log ; 3, progress and results summary; 5, full report
#'   [default 2 or as specified using gl.set.verbosity].
#' @return An object of class 'dist' giving distances between populations
#' @export
#' @author author(s): Arthur Georges. Custodian: Arthur Georges -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' # SNP genotypes
#' D <- gl.dist.pop(testset.gl, method='euclidean')
#' D <- gl.dist.pop(testset.gl, method='euclidean',scale=TRUE)
#' D <- gl.dist.pop(testset.gl, method='nei')
#' D <- gl.dist.pop(testset.gl, method='reynolds')
#' D <- gl.dist.pop(testset.gl, method='chord')
#' D <- gl.dist.pop(testset.gl, method='fixed-diff')
#' #Presence-Absence data
#' D <- gl.dist.pop(testset.gs, method='euclidean')

gl.dist.pop <- function(x,
                        method = "euclidean",
                        plot.out = TRUE,
                        scale = FALSE,
                        output="dist",
                        p = NULL,
                        plot_theme = theme_dartR(),
                        plot_colors = two_colors,
                        save2tmp = FALSE,
                        verbose = NULL) {

    # CHECK IF PACKAGES ARE INSTALLED
    pkg <- "reshape2"
    if (!(requireNamespace(pkg, quietly = TRUE))) {
        stop(error(
            "Package",
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
        utils.check.datatype(x, accept = c("SNP","SilicoDArT"), verbose = verbose)
    
    # FUNCTION SPECIFIC ERROR CHECKING
    
    # DO THE JOB
    
    available_methods <-
        c(
            "euclidean",
            "nei",
            "reynolds",
            "chord",
            "fixed-diff"
        )

    if (!(method %in% available_methods)) {
        cat(error("Fatal Error: Specified distance method is not among those 
                available.\n"))
            stop("Specify one of ",paste(available_methods, 
                collapse = ", ")," or fixed-diff.\n")
    }
    # hard.min.p <- 0.25

    nI <- nInd(x)
    nL <- nLoc(x)
    nP <- nPop(x)
    dd <- array(NA, c(nPop(x), nPop(x)))
    
    # Calculate distances 
    # if (method %in% distmethod) {
        if (verbose >= 2) {
            cat(report(paste(
                "  Calculating distances: ", method, "\n"
            )))
            cat(report(
                "  Refer to the dartR Distance Analysis tutorial for algorithms\n"
            ))
        }
    
    # Calculate allele frequencies for each population and locus
    f <- gl.percent.freq(x, verbose = 0)
    # Select only pop, locus, frequency columns
    f <- f[, c(1, 2, 6)]
    # Convert to a pop x locus matrix
    f <- reshape2::dcast(f, popn ~ locus, value.var = "frequency")
    # Reassign names to the populations, and convert from percentages to proportions
    row.names(f) = f[, 1]
    f <- f[,-c(1)]
    p <- f / 100

# For both DArTseq and SilicoDArT
    if (method == "euclidean") {
        for (i in (1:(nP - 1))) {
            for (j in ((i + 1):nP)) {
                row1 <- p[i,]
                row2 <- p[j,]
                sq <- (row1-row2)**2
                sq <- sq[!is.na(sq)]
                L <- length(sq)
                if(scale==TRUE){
                    dd[j,i] <- sqrt(sum(sq)/L)
                } else {
                    dd[j,i] <- sqrt(sum(sq))
                }
            }
        }
    }
    
# For DArTseq only
    if (method == "reynolds") {
        if(datatype=="SilicoDArT"){
            stop(error("Fatal Error: Reynolds Distance is not available 
                       for presence-absence data\n"))
        }
        for (i in (1:(nP - 1))) {
            for (j in ((i + 1):nP)) {
                # Pull the loci for individuals i and j
                prow1 <- p[i,]
                prow2 <- p[j,]
                # Delete the pairwise missing
                tmp <- prow1+prow2
                prow1 <- prow1[!is.na(tmp)]
                prow2 <- prow2[!is.na(tmp)]
                # Squares
                psq <- (prow1-prow2)**2
                # Repeat for q
                qrow1 <- 1-prow1
                qrow2 <- 1-prow2
                qsq <- (qrow1-qrow2)**2
                # Cross products
                p12 <- prow1*prow2
                q12 <- qrow1*qrow2
                # Non-missing loci
                L <- length(psq)
                dd[j,i] <- sqrt(sum(psq+qsq)/(2*sum(1-p12-q12)))
            }
        }
    }
    
    if (method == "nei") {
        if(datatype=="SilicoDArT"){
            stop(error("Fatal Error: Nei Standard Distance is not available 
                       for presence-absence data\n"))
        }
        for (i in (1:(nP - 1))) {
            for (j in ((i + 1):nP)) {
                # Pull the loci for individuals i and j
                prow1 <- p[i,]
                prow2 <- p[j,]
                # Delete the pairwise missing
                tmp <- prow1+prow2
                prow1 <- prow1[!is.na(tmp)]
                prow2 <- prow2[!is.na(tmp)]
                # Squares
                p1sq <- prow1*prow1
                p2sq <- prow2*prow2
                # Repeat for q=1-p
                qrow1 <- 1-prow1
                qrow2 <- 1-prow2
                q1sq <- qrow1*qrow1
                q2sq <- qrow2*qrow2
                # Cross products
                p12 <- prow1*prow2
                q12 <- qrow1*qrow2
                # Number of non-missing loci
                L <- length(p12)  
                
                dd[j,i] <- -log(sqrt(sum(p12+q12))/(sqrt(sum(p1sq+q1sq))*sqrt(sum(p2sq+q2sq))))
            }
        }
    }
    
    if (method == "chord") {
        if(datatype=="SilicoDArT"){
            stop(error("Fatal Error: Czfordi-Edwards Chord Distance is not available 
                       for presence-absence data\n"))
        }
        for (i in (1:(nP - 1))) {
            for (j in ((i + 1):nP)) {
                # Pull the loci for individuals i and j
                prow1 <- p[i,]
                prow2 <- p[j,]
                # Delete the pairwise missing
                tmp <- prow1+prow2
                prow1 <- prow1[!is.na(tmp)]
                prow2 <- prow2[!is.na(tmp)]
                # create proportions for allele 2
                qrow1 <- 1-prow1
                qrow2 <- 1-prow2
                # Cross products
                p12 <- prow1*prow2
                q12 <- qrow1*qrow2
                # Non-missing Loci
                L <- length(p12)

                dd[j,i] <- (2/pi)*sqrt(2*(1 - (sum(sqrt(p12))/L + sum(sqrt(q12)/L))))
            }
        }
    }
    if (method == "fixed-diff") {
        dd <- gl.fixed.diff(x, verbose = 0)[[3]]/100
        if (verbose >= 2) {
            cat(report("  Calculating proportion of fixed differences\n"))
            cat(
                warn(
                    "Note: this distance may be non-metric, and so should be considered a dissimilarity measure\n"
                )
            )
        }
    }

    # # Revert to original order ord <- rank(popNames(x)) mat <- as.matrix(dd)[ord, ord] dd <- as.dist(mat)

    # PLOT Plot Box-Whisker plot
    
    if (plot.out) {
        if (datatype == "SNP") {
            title_plot <- paste0("SNP data\nUsing ", method, " distance")
        } else {
            title_plot <-
                paste0("Tag P/A data (SilicoDArT)\nUsing ",
                       method,
                       " distance")
        }
        values <- NULL
        df_plot <- data.frame(values = as.vector(dd))
        
        # Boxplot
        p1 <- ggplot(df_plot, aes(y = values)) +
        geom_boxplot(color = plot_colors[1], fill = plot_colors[2]) +
        coord_flip()  +
        plot_theme  +
        xlim(range = c(-1,1)) + 
        ylim(min(df_plot$values, na.rm = TRUE),max(df_plot$values, na.rm = TRUE)) + 
        ylab(" ") + 
        theme(axis.text.y = element_blank(),axis.ticks.y = element_blank()) + 
        ggtitle(title_plot)
        
        # Histogram
        p2 <- ggplot(df_plot, aes(x = values)) +
        geom_histogram(bins = 20,color = plot_colors[1], fill = plot_colors[2]) +
        xlim(min(df_plot$values, na.rm = TRUE), max(df_plot$values, na.rm = TRUE)) +
        xlab("Distance") +
        ylab("Count") +
        plot_theme
    }
    
    # SUMMARY Print out some statistics
    if (verbose >= 3) {
        cat("  Reporting inter-population distances\n")
        cat("  Distance measure:", method, "\n")
        cat("    No. of populations =", nPop(x), "\n")
        cat("    Average no. of individuals per population =",
            nInd(x) / nPop(x),
            "\n")
        cat("    No. of loci =", nLoc(x), "\n")
        cat("    Minimum Distance: ", round(min(dd), 2), "\n")
        cat("    Maximum Distance: ", round(max(dd), 2), "\n")
        cat("    Average Distance: ", round(mean(dd), 3), "\n")
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
    
    # PRINTING OUTPUTS
    if (plot.out) {
        # using package patchwork
        p3 <- (p1 / p2) + plot_layout(heights = c(1, 4))
        suppressWarnings(print(p3))
    }
    
    if(output=="dist"){
        dd <- as.dist(dd)
        if(verbose >= 2){cat(report("  Returning a stats::dist object\n"))}
    } else {
        dd <- as.matrix(dd)
        if(verbose >= 2){cat(report("  Returning a square matrix object\n"))}
    }
    
    # FLAG SCRIPT END
    
    if (verbose > 0) {
        cat(report("Completed:", funname, "\n"))
    }
   
    return(dd)
}
