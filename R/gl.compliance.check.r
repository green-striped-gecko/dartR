#' @name gl.compliance.check
#' @title Checks a genlight object to see if it complies with dartR
#'  expectations and amends it to comply if necessary
#' @description
#' This function will check to see that the genlight object conforms to
#' expectation in regard to dartR requirements (see details), and if it does
#' not, will rectify it.
#' @param x Name of the input genlight object [required].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log ; 3, progress and results summary; 5, full report
#' [default 2 or as specified using gl.set.verbosity].
#' @details
#' A genlight object used by dartR has a number of requirements that allow
#' functions within the package to operate correctly. The genlight object
#' comprises:
#' \enumerate{
#' \item The SNP genotypes or Tag Presence/Absence data (SilicoDArT);
#' \item An associated dataframe (gl@other$loc.metrics) containing the locus
#' metrics (e.g. Call Rate, Repeatability, etc);
#' \item An associated dataframe (gl@other$ind.metrics) containing the
#' individual/sample metrics (e.g. sex, latitude (=lat), longitude(=lon), etc);
#' \item A specimen identity field (indNames(gl)) with the unique labels applied
#' to each individual/sample;
#' \item A population assignment (popNames) for each individual/specimen;
#' \item Flags that indicate whether or not calculable locus metrics have been
#' updated.
#' }
#' @return A genlight object that conforms to the expectations of dartR
#' @export
#' @author Custodian: Luis Mijangos -- Post to
#'  \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' x <- gl.compliance.check(testset.gl)
#' x <- gl.compliance.check(testset.gs)

gl.compliance.check <- function(x,
                                verbose = NULL) {
    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)
    
    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     build = "Jody",
                     verbosity = verbose)
    
    # DO THE JOB
    
    # if ploidy is null
    if (is.null(x@ploidy)) {
        # if diploid
        if (unique(ploidy(x)) == 2) {
            ploidy(x) <- rep(2, nInd(x))
            # if haploid
        } else if (unique(ploidy(x)) == 1) {
            ploidy(x) <- rep(1, nInd(x))
        } else {
            stop(error(
                "Ploidy cannot be determined, please check your input file"
            ))
        }
    }
    
    # CHECKS DATATYPE
    datatype <- utils.check.datatype(x, verbose = verbose)
    
    # if loci have no name
    if (is.null(locNames(x))) {
        locNames(x) <- paste0("Loc", 1:nLoc(x))
    }
    
    # Check that the data exist, and that they are restricted to the
    # appropriate values
    
    if (datatype == "SNP") {
        mat <- as.matrix(x)
        scores <- c(0, 1, 2, NA)
        if (verbose >= 2) {
            cat(report("  Checking coding of SNPs\n"))
        }
        if (max(mat) %in% scores) {
            if (verbose >= 1) {
                cat(report("    SNP data scored NA, 0, 1 or 2 confirmed\n"))
            }
        } else {
            if (verbose >= 1) {
                cat(
                    error(
                        "    Error: SNP data must be scored NA, 0 or 1 or 2, 
                        revisit data input\n"
                    )
                )
            }
        }
    } else {
        mat <- as.matrix(x)
        scores <- c(0, 1, NA)
        if (verbose >= 2) {
            cat(report("  Checking coding of Tag P/A data\n"))
        }
        if (max(mat) %in% scores) {
            if (verbose >= 1) {
                cat(
                    report(
                        "    Tag P/A data (SilicoDArT) scored 1, 0 (present or 
                        absent) confirmed\n"
                    )
                )
            }
        } else {
            if (verbose >= 1) {
                cat(
                    error(
               "    Error: Tag P/A data (SilicoDArT) must be scored NA
              for missing, 0 for absent or 1 for present, revisit data input\n"
                    )
                )
            }
        }
    }
    
    # Check for the locus metrics, and create if they do not exist.
    # Check for the locus metrics flags, and create if they do not exist.
    # Check for the verbosity flag, and create if it does not exist.
    
    if (verbose >= 2) {
        cat(report("  Checking locus metrics and flags\n"))
    }
    x <- utils.reset.flags(x, set = FALSE, verbose = 0)
    
    # Calculate locus metrics
    if (verbose >= 2) {
        cat(report("  Recalculating locus metrics\n"))
    }
    x <- gl.recalc.metrics(x, verbose = 0)
    
    # Check for monomorphic loci
    if (verbose >= 2) {
        cat(report("  Checking for monomorphic loci\n"))
    }
    x2 <- gl.filter.monomorphs(x, verbose = 0)
    if (nLoc(x2) == nLoc(x)) {
        if (verbose >= 1) {
            cat(report("    No monomorphic loci detected\n"))
        }
        x@other$loc.metrics.flags$monomorphs <- TRUE
    } else {
        if (verbose >= 1) {
            cat(warn("    Dataset contains monomorphic loci\n"))
        }
        x@other$loc.metrics.flags$monomorphs <- FALSE
    }
    
    # Check for loci with all NAs
    if (verbose >= 2) {
        cat(report("  Checking for loci with all missing data\n"))
    }
    x2 <- gl.filter.allna(x, verbose = 0)
    if (nLoc(x2) == nLoc(x)) {
        if (verbose >= 1) {
            cat(report("    No loci with all missing data detected\n"))
        }
        x@other$loc.metrics.flags$allna <- TRUE
    } else {
        if (verbose >= 1) {
            cat(warn("    Dataset contains loci with all missing dat\n"))
        }
        x@other$loc.metrics.flags$allna <- FALSE
    }
    
    # Check that the number of values in the loc.metrics dataframe is the same
    # as the number of loci
    if (nLoc(x) != nrow(x@other$loc.metrics)) {
        cat(
            warn(
                "  The number of rows in the loc.metrics table does not match 
                the number of loci! This is potentially a major problem if there
                is a mismatch of the loci with the metadata. Trace back to 
                identify the cause.\n"
            )
        )
    }
    
    # check that individual names are unique, and if not, add underscore and
    # letters
    if (verbose >= 2) {
        cat(report("  Checking whether individual names are unique.\n"))
    }
    
    if (any(duplicated(indNames(x)))) {
        if (verbose >= 2) {
            cat(
                warn(
                    "  Individual names are not unique. Appending an extra 
                    number to make them unique.\n"
                )
            )
        }
        indNames(x) <- make.unique(indNames(x), sep = '_')
    }
    
    # Check that the individual metrics exist, and if not, create the df
    
    if (verbose >= 2) {
        cat(report("  Checking for individual metrics\n"))
    }
    if (is.null(x@other$ind.metrics)) {
        if (verbose >= 1) {
            cat(warn("  Warning: Creating a slot for individual metrics\n"))
        }
        x@other$ind.metrics <- as.data.frame(matrix(nrow= nInd(x),ncol = 1))
        colnames(x@other$ind.metrics) <- "id"
        x@other$ind.metrics$id <- indNames(x)
    } else {
        if (verbose >= 1) {
            cat(report("    Individual metrics confirmed\n"))
        }
    }
    
    # convert the ind.metric slot into a dataframe
    x@other$ind.metrics <- as.data.frame(x@other$ind.metrics)

    # Check that the population variable exists, and if it does not, create it
    # with a single population 'pop1'
    
    if (verbose >= 2) {
        cat(report("  Checking for population assignments\n"))
    }
    if (is.null(pop(x)) |
        is.na(length(pop(x))) | length(pop(x)) <= 0) {
        if (verbose >= 1) {
            cat(
                warn(
                    "  Population assignments not detected, individuals assigned
                    to a single population labelled 'pop1'\n"
                )
            )
        }
        pop(x) <- array("pop1", dim = nInd(x))
        pop(x) <- as.factor(pop(x))
    } else {
        if (verbose >= 1) {
            cat(report("    Population assignments confirmed\n"))
        }
    }
    
    # check if coordinates are in the right place and not misspell
    if (!is.null(x@other$latlong)) {
        x@other$latlon <- x@other$latlong
    }
    
    if (!is.null(x@other$latlon$long)) {
        x@other$latlon$lon <- x@other$latlon$long
    }
    
    # remove misspelt columns if they exist...
    x@other$latlong <- NULL
    x@other$latlon$long <- NULL
    if (verbose >= 2) {
        cat(report(
            "  Spelling of coordinates checked and changed if necessary to 
            lat/lon\n"
        ))
    }
    
    # ADD TO HISTORY
    nh <- length(x@other$history)
    x@other$history[[nh + 1]] <- match.call()
    
    # FLAG SCRIPT END
    
    if (verbose > 0) {
        cat(report("Completed:", funname, "\n"))
    }
    
    invisible(x)
}

