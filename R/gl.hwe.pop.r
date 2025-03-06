#' @name gl.hwe.pop
#' @title Performs Hardy-Weinberg tests over loci and populations
#' @description
#' Hardy-Weinberg tests are performed for each loci in each of the populations
#' as defined by the pop slot in a genlight object.
#'
#' @param x A genlight object with a population defined
#' [pop(x) does not return NULL].
#' @param alpha_val Level of significance for testing [default 0.05].
#' @param HWformat Switch if data should be returned in HWformat (counts of
#' Genotypes to be used in package \code{HardyWeinberg})
#' @param plot.out If TRUE, returns a plot object compatible with ggplot,
#' otherwise returns a dataframe [default TRUE].
#' @param plot_theme User specified theme [default theme_dartR()].
#' @param plot_colors Vector with two color names for the borders and fill
#' [default two_colors].
#'  [default discrete_palette].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log ; 3, progress and results summary; 5, full report
#' [default NULL, unless specified using gl.set.verbosity].
#'
#' @details
#' This function employs the \code{HardyWeinberg} package, which needs to be
#' installed. The function that is used is
#' \code{\link[HardyWeinberg]{HWExactStats}}, but there are several other great
#' functions implemented in the package regarding HWE. Therefore, this function
#' can return the data in the format expected by the HWE package expects, via
#' \code{HWformat=TRUE} and then use this to run other functions of the package.
#'
#' This functions performs a HWE test for every population (rows) and loci
#' (columns) and returns a true false matrix. True is reported if the p-value of
#' an HWE-test for a particular loci and population was below the specified
#' threshold (alpha_val, default=0.05). The thinking behind this approach is
#' that loci that are not in HWE in several populations have most likely to be
#' treated (e.g. filtered if loci under selection are of interest). If plot=TRUE
#'  a barplot on the loci and the sum of deviation over all population is
#'  returned. Loci that deviate in the majority of populations can be identified
#'   via colSums on the resulting matrix.
#'
#' Plot themes can be obtained from \itemize{
#'  \item \url{https://ggplot2.tidyverse.org/reference/ggtheme.html} and \item
#'  \url{https://yutannihilation.github.io/allYourFigureAreBelongToUs/ggthemes/}
#'  }
#'
#' Resultant ggplots and the tabulation are saved to the session's temporary
#' directory.
#'
#' @return The function returns a list with up to three components:
#' \itemize{
#'  \item 'HWE' is the matrix over loci and populations
#'  \item 'plot' is a plot (ggplot) which shows the significant results
#'  for population and loci (can be amended further using ggplot syntax)
#'  \item 'HWEformat=TRUE' the 'HWformat' entails SNP data for each population
#'  in 'HardyWeinberg'-format to be used with other functions of the package
#'  (e.g \code{\link[HardyWeinberg]{HWPerm}} or
#'  \code{\link[HardyWeinberg]{HWExactPrevious}}).
#'  }
#' @author Custodian: Arthur Georges -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' out <- gl.hwe.pop(bandicoot.gl[,1:33], alpha_val=0.05, plot.out=TRUE, HWformat=FALSE)
#' @export

gl.hwe.pop <-  function(x,
                        alpha_val = 0.05,
                        plot.out = TRUE,
                        plot_theme = theme_dartR(),
                        plot_colors = c("gray90", "deeppink"),
                        HWformat = FALSE,
                        verbose = NULL) {
    
    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)
    
    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     build = "Jackson",
                     verbose = verbose)
    
    # CHECK DATATYPE
    datatype <- utils.check.datatype(x, verbose = verbose)
    
    # FUNCTION SPECIFIC ERROR CHECKING check if packages is installed
    pkg <- "HardyWeinberg"
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
        pop(x) <- array("pop1", dim = nLoc(x))
        pop(x) <- as.factor(pop(x))
    }
    
    # DO THE JOB
    
    pops <- adegenet::seppop(x)
    
    out <- lapply(pops, function(x) {
        pp <- as.matrix(x)
        xx <-
            data.frame(
                AA = colSums(pp == 0, na.rm = T),
                AB = colSums(pp == 1, na.rm = T),
                BB = colSums(pp == 2, na.rm = T)
            )
        HardyWeinberg::HWExactStats(xx) < alpha_val
    })
    
    # convert to matrix
    output <- matrix(unlist(out), ncol = nLoc(x), byrow = TRUE)
    rownames(output) <- names(pops)
    colnames(output) <- locNames(x)
    
    p1 <- NULL
    
    if (plot.out) {
        Var1 <- Var2 <- value <- NA
        longData <- reshape2::melt((as.matrix(output)) * 1)
        
        p1 <-
            ggplot(longData, aes(x = Var2, y = Var1, )) + geom_raster(aes(fill = plot_colors[value + 1])) + scale_fill_identity() + labs(x = "Loci",
                                                                                                                                         y = "Populations",
                                                                                                                                         title = "HWE over populations and loci") + plot_theme + theme(
                                                                                                                                             axis.title.x = element_blank(),
                                                                                                                                             axis.text.x = element_blank(),
                                                                                                                                             axis.ticks.x = element_blank(),
                                                                                                                                             panel.border = element_blank(),
                                                                                                                                             panel.grid.major = element_blank(),
                                                                                                                                             panel.grid.minor = element_blank()
                                                                                                                                         )
        # PRINTING OUTPUTS
        print(p1)
    }
    
    out <- NULL
    
    if (HWformat) {
        out <- lapply(pops, function(x) {
            pp <- as.matrix(x)
            xx <-
                data.frame(
                    AA = colSums(pp == 0, na.rm = T),
                    AB = colSums(pp == 1, na.rm = T),
                    BB = colSums(pp == 2, na.rm = T)
                )
            rownames(xx) <- locNames(x)
            xx
        })
    }
    
    # FLAG SCRIPT END
    
    if (verbose >= 1) {
        cat(report("Completed:", funname, "\n\n"))
    }
    
    # RETURN
    invisible(list(
        HWE = output,
        plot = p1,
        HWformat = out
    ))
    
}
