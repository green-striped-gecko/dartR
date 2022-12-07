#' @name gl.grm
#' @title Calculates an identity by descent matrix
#' @description
#' This function calculates the mean probability of identity by state (IBS)
#' across loci that would result from all the possible crosses of the
#' individuals analyzed. IBD is calculated by an additive relationship matrix
#' approach developed by Endelman and Jannink (2012) as implemented in the
#' function \link[rrBLUP]{A.mat} (package rrBLUP).
#'
#' @param x Name of the genlight object containing the SNP data [required].
#' @param plotheatmap A switch if a heatmap should be shown [default TRUE].
#' @param palette_discrete A discrete palette for the color of populations or a
#' list with as many colors as there are populations in the dataset
#'  [default discrete_palette].
#' @param palette_convergent A convergent palette for the IBD values
#'  [default convergent_palette].
#' @param legendx x coordinates for the legend[default 0].
#' @param legendy y coordinates for the legend[default 1].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#'  progress log ; 3, progress and results summary; 5, full report
#'  [default 2 or as specified using gl.set.verbosity].
#' @param ... Parameters passed to function A.mat from package rrBLUP.
#'
#' @details
#' Two or more alleles are identical by descent (IBD) if they are identical
#' copies of the same ancestral allele in a base population. The additive
#' relationship matrix is a theoretical framework for estimating a relationship
#' matrix that is consistent with an approach to estimate the probability that
#' the alleles at a random locus are identical in state (IBS).
#'
#' This function also plots a heatmap, and a dendrogram, of IBD values where
#' each diagonal element has a mean that equals 1+f, where f is the inbreeding
#' coefficient (i.e. the probability that the two alleles at a randomly chosen
#' locus are IBD from the base population). As this probability lies between 0
#'  and 1, the diagonal elements range from 1 to 2. Because the inbreeding
#'  coefficients are expressed relative to the current population, the mean of
#'  the off-diagonal elements is -(1+f)/n, where n is the number of loci.
#'  Individual names are shown in the margins of the heatmap and colors
#'  represent different populations.
#'
#' @return An identity by descent matrix
#' @author Custodian: Arthur Georges -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' @references \itemize{
#' \item Endelman, J. B. (2011). Ridge regression and other kernels for genomic
#'  selection with r package rrblup. The Plant Genome 4, 250.
#' \item Endelman, J. B. , Jannink, J.-L. (2012). Shrinkage estimation of the
#' realized relationship matrix. G3: Genes, Genomics, Genetics 2, 1405.
#' }
#' @examples
#' gl.grm(platypus.gl[1:10,1:100])
#'
#' @seealso \code{\link{gl.grm.network}}
#' @family inbreeding functions
#' @export

gl.grm <- function(x,
                   plotheatmap = TRUE,
                   palette_discrete = discrete_palette,
                   palette_convergent = convergent_palette,
                   legendx = 0,
                   legendy = 0.5,
                   verbose = NULL,
                   ...) {
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
    
    # check if package is installed
    pkg <- "rrBLUP"
    if (!(requireNamespace(pkg, quietly = TRUE))) {
      cat(error(
        "Package",
        pkg,
        " needed for this function to work. Please install it.\n"
      ))
      return(-1)
    }
    
    pkg <- "gplots"
    if (!(requireNamespace(pkg, quietly = TRUE))) {
      cat(error(
        "Package",
        pkg,
        " needed for this function to work. Please install it.\n"
      ))
      return(-1)
    }
    
    # Set a population if none is specified (such as if the genlight object has been generated manually)
    if (is.null(pop(x)) |
        is.na(length(pop(x))) | length(pop(x)) <= 0) {
        if (verbose >= 2) {
            cat(
                important(
                    "  Population assignments not detected, individuals assigned to a single population labelled 'pop1'\n"
                )
            )
        }
        pop(x) <- array("pop1", dim = nInd(x))
        pop(x) <- as.factor(pop(x))
    }
    
    # DO THE JOB
    
    # assigning colors to populations
    if (is(palette_discrete, "function")) {
        colors_pops <- palette_discrete(length(levels(pop(x))))
    }
    
    if (!is(palette_discrete, "function")) {
        colors_pops <- palette_discrete
    }
    
    names(colors_pops) <- as.character(levels(x$pop))
    
    # calculating the realized additive relationship matrix
    
    G <- rrBLUP::A.mat(as.matrix(x) - 1, ...)
    
    df_colors_temp_1 <-
        as.data.frame(cbind(indNames(x), as.character(pop(x)), 1:nInd(x)))
    colnames(df_colors_temp_1) <- c("ind", "pop", "order")
    df_colors_temp_2 <-
        as.data.frame(cbind(names(colors_pops), colors_pops))
    colnames(df_colors_temp_2) <- c("pop", "color")
    df_colors <-
        merge(df_colors_temp_1, df_colors_temp_2, by = "pop")
    df_colors$order <- as.numeric(df_colors$order)
    df_colors <- df_colors[order(df_colors$order), ]
    df_colors_2 <- df_colors[, c("pop", "color")]
    df_colors_2 <- unique(df_colors_2)
    
    if (plotheatmap == T) {
        # plotting heatmap
        par(mar = c(1, 1, 1, 1))
        gplots::heatmap.2(
            G,
            col = palette_convergent(255),
            dendrogram = "column",
            ColSideColors = df_colors$color,
            RowSideColors = df_colors$color,
            trace = "none",
            density.info = "none",
            scale = "none",
            main = "Probability of identity by descent"
        )
        legend(
            legendx,
            legendy,
            legend = df_colors_2$pop,
            fill = df_colors_2$color,
            cex = 0.75,
            title = "Populations"
        )
    }
    
    # FLAG SCRIPT END
    
    if (verbose >= 1) {
        cat(report("Completed:", funname, "\n"))
    }
    
    # RETURN
    invisible(G)
    
}
