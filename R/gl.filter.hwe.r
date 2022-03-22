#' @name gl.filter.hwe
#' @title Filters loci that show significant departure from Hardy-Weinberg
#'  Equilibrium
#' @description
#' This function filters out loci showing significant departure from H-W
#' proportions based on observed frequencies of reference homozygotes,
#'  heterozygotes and alternate homozygotes.
#'
#' Loci are filtered out if they show HWE departure in any one population.
#'
#' @param x Name of the genlight object containing the SNP data [required].
#' @param subset Way to group individuals to perform H-W tests. Either a vector
#' with population names, 'each', 'all' (see details) [default 'each'].
#' @param method_sig Method for determining statistical significance: 'ChiSquare'
#' or 'Exact' [default 'Exact'].
#' @param multi_comp Whether to adjust p-values for multiple comparisons
#' [default FALSE].
#' @param multi_comp_method Method to adjust p-values for multiple comparisons:
#' 'holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr'
#' (see details) [default 'fdr'].
#' @param alpha_val Level of significance for testing [default 0.05].
#' @param pvalue_type Type of p-value to be used in the Exact method.
#' Either 'dost','selome','midp' (see details) [default 'midp'].
#' @param cc_val The continuity correction applied to the ChiSquare test
#'  [default 0.5].
#' @param min_sample_size Minimum number of individuals per population in which
#' perform H-W tests [default 5].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#' @details
#'  There are several factors that can cause deviations from Hardy-Weinberg
#'  proportions including: mutation, finite population size, selection,
#'  population structure, age structure, assortative mating, sex linkage,
#'  nonrandom sampling and genotyping errors. Therefore, testing for
#'  Hardy-Weinberg proportions should be a process that involves a careful
#'  evaluation of the results, a good place to start is Waples (2015).
#'
#'  Note that tests for H-W proportions are only valid if there is no population
#'  substructure (assuming random mating) and have sufficient power only when
#'  there is sufficient sample size (n individuals > 15).
#'
#' Populations can be defined in three ways:
#' \itemize{
#' \item Merging all populations in the dataset using subset = 'all'.
#' \item Within each population separately using: subset = 'each'.
#' \item Within selected populations using for example: subset = c('pop1','pop2').
#' }
#'
#' Two different statistical methods to test for deviations from Hardy Weinberg
#' proportions:
#' \itemize{
#' \item The classical chi-square test (method_sig='ChiSquare') based on the
#' function \code{\link[HardyWeinberg]{HWChisq}} of the R package HardyWeinberg.
#' By default a continuity correction is applied (cc_val=0.5). The
#' continuity correction can be turned off (by specifying cc_val=0), for example
#' in cases of extreme allele frequencies in which the continuity correction can
#' lead to excessive type 1 error rates.
#' \item The exact test (method_sig='Exact') based on the exact calculations
#' contained in the function \code{\link[HardyWeinberg]{HWExactStats}} of the R
#' package HardyWeinberg, and described in Wigginton  et al. (2005). The exact
#' test is recommended in most cases (Wigginton  et al., 2005).
#' Three different methods to estimate p-values (pvalue_type) in the Exact test
#' can be used:
#' \itemize{
#' \item 'dost' p-value is computed as twice the tail area of a one-sided test.
#' \item 'selome' p-value is computed as the sum of the probabilities of all
#' samples less or equally likely as the current sample.
#' \item 'midp', p-value is computed as half the probability of the current
#' sample + the probabilities of all samples that are more extreme.
#' }
#' The standard exact p-value is overly conservative, in particular
#' for small minor allele frequencies. The mid p-value ameliorates this problem
#' by bringing the rejection rate closer to the nominal level, at the price of
#' occasionally exceeding the nominal level (Graffelman & Moreno, 2013).
#' }
#'
#' Correction for multiple tests can be applied using the following methods
#' based on the function \code{\link[stats]{p.adjust}}:
#' \itemize{
#' \item 'holm' is also known as the sequential Bonferroni technique (Rice, 1989).
#' This method has a greater statistical power than the standard Bonferroni test,
#' however this method becomes very stringent when many tests are performed and
#' many real deviations from the null hypothesis can go undetected (Waples, 2015).
#' \item 'hochberg' based on Hochberg, 1988.
#' \item 'hommel' based on Hommel, 1988. This method is more powerful than
#' Hochberg's, but the difference is usually small.
#' \item 'bonferroni' in which p-values are multiplied by the number of tests.
#' This method is very stringent and therefore has reduced power to detect
#' multiple departures from the null hypothesis.
#' \item 'BH' based on Benjamini & Hochberg, 1995.
#' \item 'BY' based on Benjamini & Yekutieli, 2001.
#' }
#'
#' The first four methods are designed to give strong control of the family-wise
#' error rate. The last two methods control the false discovery rate (FDR),
#' the expected proportion of false discoveries among the rejected hypotheses.
#' The false discovery rate is a less stringent condition than the family-wise
#' error rate, so these methods are more powerful than the others, especially
#' when number of tests is large.
#' The number of tests on which the adjustment for multiple comparisons is
#' the number of populations times the number of loci.
#'
#' @return A genlight object with the loci departing significantly from H-W
#' proportions removed.
#' @author Custodian: Luis Minjangos -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' result <- gl.filter.hwe(x = bandicoot.gl)
#' @references
#' \itemize{
#'  \item Benjamini, Y., and Yekutieli, D. (2001). The control of the false
#'  discovery rate in multiple testing under dependency. Annals of Statistics,
#'  29, 1165–1188.
#' \item Graffelman, J. (2015). Exploring Diallelic Genetic Markers: The Hardy
#' Weinberg Package. Journal of Statistical Software 64:1-23.
#' \item Graffelman, J. & Morales-Camarena, J. (2008). Graphical tests for
#' Hardy-Weinberg equilibrium based on the ternary plot. Human Heredity 65:77-84.
#' \item Graffelman, J., & Moreno, V. (2013). The mid p-value in exact tests for
#' Hardy-Weinberg equilibrium. Statistical applications in genetics and molecular
#' biology, 12(4), 433-448.
#' \item Hochberg, Y. (1988). A sharper Bonferroni procedure for multiple tests
#'  of significance. Biometrika, 75, 800–803.
#' \item Hommel, G. (1988). A stagewise rejective multiple test procedure based
#'  on a modified Bonferroni test. Biometrika, 75, 383–386.
#' \item Rice, W. R. (1989). Analyzing tables of statistical tests. Evolution,
#'  43(1), 223-225.
#' \item Waples, R. S. (2015). Testing for Hardy–Weinberg proportions: have we
#' lost the plot?. Journal of heredity, 106(1), 1-19.
#' \item Wigginton, J.E., Cutler, D.J., & Abecasis, G.R. (2005). A Note on Exact
#' Tests of Hardy-Weinberg Equilibrium. American Journal of Human Genetics
#' 76:887-893.
#' }
#' @seealso \code{\link{gl.report.hwe}}
#' @family filters/filter reports
#' @export

gl.filter.hwe <- function(x,
                          subset = "each",
                          method_sig = "Exact",
                          multi_comp = FALSE,
                          multi_comp_method = "BY",
                          alpha_val = 0.05,
                          pvalue_type = "midp",
                          cc_val = 0.5,
                          min_sample_size = 5,
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
    
    # FUNCTION SPECIFIC ERROR CHECKING check if packages are installed
    pkg <- "HardyWeinberg"
    if (!(requireNamespace(pkg, quietly = TRUE))) {
        stop(error(
            "Package",
            pkg,
            " needed for this function to work. Please install it."
        ))
    }
    
    if (datatype == "SilicoDArT") {
        cat(error("  Detected Presence/Absence (SilicoDArT) data\n"))
        stop(
            error(
                "Cannot calculate HWE from fragment presence/absence data. Please provide a SNP dataset.\n"
            )
        )
    }
    
    if (alpha_val < 0 | alpha_val > 1) {
        cat(
            warn(
                "    Warning: level of significance per locus alpha must be an integer between 0 and 1, set to 0.05\n"
            )
        )
        alpha_val <- 0.05
    }
    
    # DO THE JOB
    hold <- x
    
    #### Interpret options for subset all
    if (subset[1] == "all") {
        if (verbose >= 2) {
            cat(report("  Pooling all populations for HWE calculations\n"))
        }
        if (verbose >= 3) {
            cat(
                warn(
                    "  Warning: Significance of tests may indicate heterogeneity among populations\n\n"
                )
            )
        }
        # assigning the same population to all individuals
        pop(x) <- array("pop", dim = nInd(x))
        pop(x) <- as.factor(pop(x))
    }
    
    ########### for subset each
    if (subset[1] == "each") {
        if (verbose >= 2) {
            cat(report("  Analysing each population separately\n"))
        }
    }
    
    ########### for subset selected populations
    if (subset[1] != "each" & subset[1] != "all") {
        # check whether the populations exist in the dataset
        pops_hwe_temp <- pop(x) %in% subset
        pops_hwe <- sum(pops_hwe_temp[pops_hwe_temp == TRUE])
        # if the populations are not in the dataset
        if (pops_hwe == 0) {
            stop(
                error(
                    "Fatal Error: subset parameter must be \"each\", \"all\", or a list of populations existing in the dataset\n"
                )
            )
        }
        # subsetting the populations
        x <- x[pop(x) %in% subset]
        # assigning the same population to all individuals
        pop(x) <- array("pop", dim = nInd(x))
        pop(x) <- as.factor(pop(x))
        if (verbose >= 2) {
            cat(report(
                paste(
                    "  Pooling populations",
                    paste(subset, collapse = " "),
                    "together for HWE calculations\n"
                )
            ))
        }
        if (verbose >= 3) {
            cat(
                warn(
                    "  Warning: Significance of tests may indicate heterogeneity among populations\n\n"
                )
            )
        }
    }
    
    poplist_temp <- seppop(x)
    # filtering monomorphs
    poplist <-
        lapply(poplist_temp, gl.filter.monomorphs, verbose = 0)
    
    # testing whether populations have heteromorphic loci
    monomorphic_pops_temp <- unlist(lapply(poplist, nLoc))
    monomorphic_pops <-
        monomorphic_pops_temp[which(monomorphic_pops_temp == 0)]
    
    if (length(monomorphic_pops) > 0) {
        if (verbose >= 2) {
            cat(
                warn(
                    " Warning: No heteromorphic loci in population",
                    names(monomorphic_pops),
                    "... skipped\n"
                )
            )
            # removing pops that do not have heteromorphic loci
            pops_to_remove <-
                which(names(poplist) %in% names(monomorphic_pops))
            poplist <- poplist[-pops_to_remove]
        }
    }
    
    # testing whether populations have small sample size
    n_ind_pops_temp <- unlist(lapply(poplist, nInd))
    n_ind_pops <-
        n_ind_pops_temp[which(n_ind_pops_temp <= min_sample_size)]
    
    if (length(n_ind_pops) > 0) {
        if (verbose >= 2) {
            cat(
                warn(
                    " Warning: population",
                    names(n_ind_pops),
                    "has less than",
                    min_sample_size,
                    "individuals... skipped\n"
                )
            )
            # removing pops that have low sample size
            pops_to_remove_2 <-
                which(names(poplist) %in% names(n_ind_pops))
            poplist <- poplist[-pops_to_remove_2]
        }
    }
    
    if (length(poplist) < 1) {
        stop(
            error(
                "No populations left after removing populations with low sample size and populations with monomorphic loci"
            )
        )
    }
    
    result <- as.data.frame(matrix(nrow = 1, ncol = 10))
    colnames(result) <-
        c(
            "Population",
            "Locus",
            "Hom_1",
            "Het",
            "Hom_2",
            "N",
            "Prob",
            "Sig",
            "Prob.adj",
            "Sig.adj"
        )
    
    for (i in poplist) {
        mat_HWE_temp <- t(as.matrix(i))
        mat_HWE <- matrix(nrow = nLoc(i), ncol = 3)
        colnames(mat_HWE) <- c("AA", "AB", "BB")
        mat_HWE[, "AA"] <- apply(mat_HWE_temp, 1, function(y) {
            length(y[which(y == 0)])
        })
        mat_HWE[, "AB"] <- apply(mat_HWE_temp, 1, function(y) {
            length(y[which(y == 1)])
        })
        mat_HWE[, "BB"] <- apply(mat_HWE_temp, 1, function(y) {
            length(y[which(y == 2)])
        })
        
        if (method_sig == "ChiSquare") {
            p.values <- apply(mat_HWE, 1, function(x) {
                HardyWeinberg::HWChisq(x, verbose = F)$pval
            })
        }
        
        if (method_sig == "Exact") {
            p.values <-
                HardyWeinberg::HWExactStats(mat_HWE, pvaluetype = pvalue_type)
        }
        
        total <- rowSums(mat_HWE, na.rm = T)
        
        sig2 <- rep(NA, length(p.values))
        p.values_adj <- rep(NA, length(p.values))
        bonsig2 <- rep(NA, length(p.values))
        
        # Assemble results into a dataframe
        result_temp <-
            cbind.data.frame(
                popNames(i),
                locNames(i),
                mat_HWE,
                total,
                p.values,
                sig2,
                p.values_adj,
                bonsig2,
                stringsAsFactors = FALSE
            )
        names(result_temp) <-
            c(
                "Population",
                "Locus",
                "Hom_1",
                "Het",
                "Hom_2",
                "N",
                "Prob",
                "Sig",
                "Prob.adj",
                "Sig.adj"
            )
        
        result <-
            rbind.data.frame(result, result_temp, stringsAsFactors = FALSE)
    }
    result <- result[-1, ]
    
    if (multi_comp == TRUE) {
        result$Prob.adj <-
            stats::p.adjust(result$Prob, method = multi_comp_method)
    }
    
    result[which(result$Prob < alpha_val), "Sig"] <- "sig"
    result[which(result$Prob > alpha_val), "Sig"] <- "no_sig"
    result[which(result$Prob.adj < alpha_val), "Sig.adj"] <- "sig"
    result[which(result$Prob.adj > alpha_val), "Sig.adj"] <-
        "no_sig"
    
    df <- result
    #### Report the results
    if (multi_comp == F) {
        df <- df[which(df$Prob <= alpha_val), ]
    }
    if (multi_comp == T) {
        df <- df[which(df$Prob.adj <= alpha_val), ]
    }
    
    failed.loci <- as.character(unique(df$Locus))
    
    if (verbose >= 2) {
        cat("  Loci examined:", nLoc(hold), "\n")
    }
    
    index <- !locNames(hold) %in% failed.loci
    hold <- hold[, index]
    hold@other$loc.metrics <- hold@other$loc.metrics[index, ]
    
    #### Report the results
    if (verbose >= 2) {
        if (multi_comp == T) {
            cat(
                "  Deleted",
                length(failed.loci),
                "loci with significant departure from HWE, after correction for multiple tests using the",
                multi_comp_method,
                "method at experiment-wide alpha =",
                alpha_val,
                "\n"
            )
        }
        if (multi_comp == F) {
            cat(
                "  Deleted",
                length(failed.loci),
                "loci with significant departure from HWE at alpha =",
                alpha_val,
                "applied locus by locus\n"
            )
        }
        cat("  Loci retained:", nLoc(hold), "\n\n")
        cat(
            important(
                "    Adjustment of p-values for multiple comparisons vary with sample size\n"
            )
        )
    }
    
    # ADD TO HISTORY
    nh <- length(hold@other$history)
    hold@other$history[[nh + 1]] <- match.call()
    
    # FLAG SCRIPT END
    if (verbose >= 1) {
        cat(report("Completed:", funname, "\n\n"))
    }
    
    # RETURN
    invisible(hold)
    
}
