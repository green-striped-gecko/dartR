#' @name gl.assign.pa
#' @title Eliminates populations as possible source populations for an individual
#'  of unknown provenance, using private alleles
#' @description
#' This script eliminates from consideration as putative source populations,
#' those populations for which the individual has too many private alleles. The
#' populations that remain are putative source populations, subject to further
#' consideration.
#'
#' The algorithm identifies those target populations for which the individual
#' has no private alleles or for which the number of private alleles does not
#' exceed a user specified threshold.
#'
#' An excessive count of private alleles is an indication that the unknown does
#'  not belong to a target population (provided that the sample size is
#'  adequate, say >=10).
#'
#' @param x Name of the input genlight object [required].
#' @param unknown SpecimenID label (indName) of the focal individual whose
#' provenance is unknown [required].
#' @param nmin Minimum sample size for a target population to be included in the
#' analysis [default 10].
#' @param threshold Populations to retain for consideration; those for which the
#'  focal individual has less than or equal to threshold loci with private
#'  alleles [default 0].
#' @param n.best If given a value, dictates the best n=n.best populations to
#' retain for consideration (or more if their are ties) based on private alleles
#'  [default NULL].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#'  [default 2 or as specified using gl.set.verbosity].
#' @return A genlight object containing the focal individual (assigned to
#' population 'unknown') and populations for which the focal individual is not
#' distinctive (number of loci with private alleles less than or equal to the
#' threshold). If no such populations, the genlight object contains only data
#' for the unknown individual.
#'
#' @importFrom stats dnorm qnorm
#' @export
#'
#' @author Custodian: Arthur Georges -- Post to
#'  \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' # Test run with a focal individual from the Macleay River (EmmacMaclGeor)
#'   x <- gl.assign.pa(testset.gl, unknown='UC_00146', nmin=10, threshold=1,verbose=3)
#'
#' @seealso \code{\link{gl.assign.pca}}

gl.assign.pa <- function(x,
                         unknown,
                         nmin = 10,
                         threshold = 0,
                         n.best = NULL,
                         verbose = NULL) {
    
    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)
    
    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     build = "Jackson",
                     verbosity = verbose)
    
    # CHECK DATATYPE
    datatype <- utils.check.datatype(x, verbose = verbose)
    
    # FUNCTION SPECIFIC ERROR CHECKING
    
    test <- unknown %in% indNames(x)
    if (!all(test, na.rm = FALSE)) {
        stop(
            error(
                "Fatal Error: nominated focal individual (of unknown provenance) is not present in the dataset!\n"
            )
        )
    }
    if (nmin <= 0 & verbose >= 1) {
        cat(
            warn(
                "  Warning: the minimum size of the target population must be greater than zero, set to 10\n"
            )
        )
        nmin <- 10
    }
    if (threshold < 0 & verbose >= 1) {
        cat(
            warn(
                "  Warning: the threshold for private alleles must be non-negative, set to 0\n"
            )
        )
        threshold <- 0
    }
    
    if (!is.null(n.best)) {
        if (n.best < 1 & verbose >= 1) {
            cat(
                warn(
                    "  Warning: the n.best parameter for retention of best match populations must be a positive integer, set to NULL\n"
                )
            )
            n.best <- NULL
        }
    }
    
    # DO THE JOB
    
    # Set a hard recommended minimum population size
    hard.min <- 10
    if (verbose >= 1 & (nmin < hard.min)) {
        cat(warn(
            "  Warning: The specified minimum sample size is less than 10 individuals\n"
        ))
        cat(
            warn(
                "    Risk of alleles present in the unknown being missed during sampling of populations with sample sizes less than 10\n"
            )
        )
    }
    
    # Assign the unknown individual to population 'unknown'
    vec <- as.vector(pop(x))
    vec[indNames(x) == unknown] <- "unknown"
    pop(x) <-
        as.factor(vec)  # Note, population containing the unknown has been reduced in size by 1
    
    # Remove loci scored as NA for the unknown
    a <- x[pop(x) == "unknown",]
    b <-
        data.frame(as.matrix(a))  # Fuck, it changed the locus names, replaced hyphens with periods
    names(b) <- locNames(a)  # Change them back
    c <- names(b)[is.na(b)]
    if (length(c) > 0) {
        x <- gl.drop.loc(x, loc.list = c, verbose = 0)
    }
    
    # Split the genlight object into one containing the unknown and one containing the remaining populations
    unknowns <- gl.keep.pop(x, pop.list = "unknown", verbose = 0)
    knowns <- gl.drop.pop(x, pop.list = "unknown", verbose = 0)
    
    # Remove all known populations with less than nmin individuals
    pop.keep <- levels(pop(knowns))[table(pop(knowns)) >= nmin]
    pop.toss <- levels(pop(knowns))[table(pop(knowns)) < nmin]
    if (verbose >= 2) {
        cat(report(
            "  Discarding",
            length(pop.toss),
            "populations with sample size <",
            nmin,
            ":\n"
        ))
        if (verbose >= 3) {
            cat(paste(pop.toss, collapse = ", "), "\n")
        }
    }
    
    # knowns <- knowns[pop(knowns) %in% pop.keep]
    knowns <- gl.keep.pop(knowns, pop.list = pop.keep, verbose = 0)
    
    # Warn of populations retained with sizes less than the hard wired minimum
    pop.warn <- levels(pop(knowns))[table(pop(knowns)) < hard.min]
    if (length(pop.warn >= 1)) {
        if (verbose >= 1) {
            cat(
                warn(
                    "  Warning: Some retained populations have sample sizes less than",
                    hard.min,
                    ":",
                    pop.warn,
                    "\n"
                )
            )
        }
        if (verbose >= 1) {
            cat(warn(
                "    Substantial risk of private alleles arising as sampling error\n\n"
            ))
        }
    }
    
    # Report number of focal individuals (1) and number of target populations
    n <- length(pop(unknowns))
    N <- length(levels(pop(knowns)))
    if (n != 1) {
        if (verbose >= 0) {
            stop(
                error(
                    "Fatal Error: Number of unknown focal individuals > 1; population label 'unknown' already in use\n"
                )
            )
        }
    }
    if (verbose >= 2) {
        cat(
            report(
                "  Assigning",
                n,
                "unknown individual",
                unknown,
                "to",
                N,
                "target populations using",
                nLoc(x),
                "loci\n"
            )
        )
    }
    
    # CALCULATE NUMBER OF LOCI WITH PRIVATE ALLELES
    
    # Genotype of the unknown individual
    unknown.ind <- as.matrix(unknowns)
    # for each population
    pop.list <- seppop(knowns)
    count <- rep(0, length(pop.list))
    # For each population
    for (i in 1:length(pop.list)) {
        gen <- as.matrix(pop.list[[i]])
        # For each locus
        for (j in 1:nLoc(pop.list[[i]])) {
            unknown.gen <- unknown.ind[j]
            pop.gen <- gen[, j]
            if (length(pop.gen) >= nmin) {
                # Where the unknown focal individual is homozygous reference allele [aa]
                if (unknown.gen == 0) {
                    # If all individuals in the target population are bb [not aa or ab] then focal individual has private allele [a]
                    if (all(pop.gen == 2, na.rm = TRUE)) {
                        count[i] <- count[i] + 1
                    }
                }
                # Where the unknown focal individual is homozygous for the alternate allele [bb]
                if (unknown.gen == 2) {
                    # If all individuals in the target population are aa [not bb or ab] then focal individual has private allele [b]
                    if (all(pop.gen == 0, na.rm = TRUE)) {
                        count[i] <- count[i] + 1
                    }
                }
                # Where the unknown focal individual is heterozgous [ab]
                if (unknown.gen == 1) {
                    # If all individuals in the target population are aa, then [b] is private or if bb, then [a] is private
                    if ((all(pop.gen == 0, na.rm = TRUE)) ||
                        (all(pop.gen == 2, na.rm = TRUE))) {
                        count[i] <- count[i] + 1
                    }
                }
            }
        }
    }
    
    # Print out results
    
    if (verbose >= 3) {
        cat("  Table showing populations against number of loci with private alleles\n")
    }
    counter <- 1
    retain <- NULL
    for (m in levels(as.factor(count))) {
        if (as.numeric(as.character(m)) > threshold) {
            if (verbose >= 3) {
                cat(paste0("  >", threshold, "---"))
            }
        }
        if (verbose >= 3) {
            cat("  ", m, levels(pop(knowns))[count == m], "\n")
        }
        if (!is.null(n.best)) {
            if (counter < n.best) {
                retain <- c(retain, levels(pop(knowns))[count == m])
                counter <-
                    counter + length(levels(pop(knowns))[count == m])
            }
        }
    }
    
    # Save the data in a new gl object
    
    if (is.null(n.best)) {
        gl <-
            gl.keep.pop(
                x,
                pop.list = c(levels(pop(knowns))[count <= threshold], "unknown"),
                mono.rm = TRUE,
                verbose = 0
            )
    } else {
        gl <-
            gl.keep.pop(
                x,
                pop.list = c(retain, "unknown"),
                mono.rm = TRUE,
                verbose = 0
            )
    }
    
    # Check that there is more than one population to assign (excluding 'unknown')
    if (is.null(n.best)) {
        if (verbose >= 2) {
            if (nPop(gl) == 1) {
                # Taking into account the unknown as a population
                if (verbose >= 2) {
                    cat(
                        report(
                            "  There are no populations retained for assignment.",
                            "  The unknown may not belong to one of the target populations.",
                            "  Returning genlight object for the unknown individual only\n"
                        )
                    )
                }
            } else {
                if (verbose >= 2) {
                    cat(
                        report(
                            "  Identified and retained",
                            nPop(gl) - 1,
                            "putative source populations for",
                            unknown,
                            "based on the specified threshold\n"
                        )
                    )
                    cat(report("  Monomorphic loci removed\n"))
                }
            }
        }
    } else {
        if (verbose >= 2) {
            if (length(retain) == 0) {
                if (verbose >= 2) {
                    cat(
                        report(
                            "  There are no populations retained for assignment.",
                            "  The unknown may not belong to one of the target populations.",
                            "  Returning genlight object for the unknown individual only\n"
                        )
                    )
                }
            } else {
                if (verbose >= 2) {
                    cat(
                        report(
                            "  Identified and retained the top",
                            length(retain),
                            "best listed putative source populations for",
                            unknown,
                            "\n"
                        )
                    )
                    cat(report("  Monomorphic loci removed\n"))
                }
            }
        }
    }
    
    # FLAG SCRIPT END
    
    if (verbose > 0) {
        cat(report("Completed:", funname, "\n"))
    }
    
    return(gl)
}
