#' @name gl.assign.pca
#' @title Assign an individual of unknown provenance to population based on PCA
#' @description
#' This script assigns an individual of unknown provenance to one or more target
#'  populations based on its proximity to each population defined by a
#'  confidence ellipse in ordinated space of two dimensions.
#'
#' The following process is followed:
#' \enumerate{
#' \item The space defined by the loci is ordinated to yield a series of
#' orthogonal axes (independent), and the top two dimensions are considered.
#' Populations for which the unknown lies outside the specified confidence
#' limits are no longer removed from the dataset.
#' }
#' @details
#' There are three considerations to assignment. First, consider only those
#' populations for which the unknown has no private alleles. Private alleles are
#' an indication that the unknown does not belong to a target population
#' (provided that the sample size is adequate, say >=10). This can be evaluated
#'  with gl.assign.pa().
#'
#' A next step is to consider the PCoA plot for populations where no private
#' alleles have been detected and the position of the unknown in relation to the
#' confidence ellipses as is plotted by this script. Note, this plot is
#' considering only the top two dimensions of the ordination, and so an unknown
#' lying outside the confidence ellipse can be unambiguously interpreted as it lying outside
#' the confidence envelope. However, if the unknown lies inside the confidence
#' ellipse in two dimensions, then it may still lie outside the confidence
#' envelope in deeper dimensions. This second step is good for eliminating populations from
#' consideration, but does not provide confidence in assignment.
#'
#' The third step is to consider the assignment probabilities, using the script
#' gl.assign.mahalanobis(). This approach calculates the squared Generalised 
#' Linear Distance (Mahalanobis distance) of the unknown from the centroid 
#' for each population, and calculates the probability associated with its 
#' quantile under the zero truncated normal distribution. This index takes 
#' into account position of the unknown in relation to the confidence envelope 
#' in all selected dimensions of the ordination. 
#'
#' Each of these approaches provides evidence, none are 100% definitive. They
#'  need to be interpreted cautiously. They are best applied sequentially.
#'  
#' In deciding the assignment, the script considers an individual to be an
#' outlier with respect to a particular population at alpha = 0.001 as default.
#'
#' @param x Name of the input genlight object [required].
#' @param unknown Identity label of the focal individual whose provenance is
#' unknown [required].
#' @param plevel Probability level for bounding ellipses in the PCoA plot
#' [default 0.999].
#' @param plot.out If TRUE, plot the 2D PCA showing the position 
#' of the unknown [default TRUE]
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2 or as specified using gl.set.verbosity].
#'
#' @return A genlight object containing only those populations that are
#' putative source populations for the unknown individual. 
#'
#' @importFrom stats dnorm qnorm
#' @export
#'
#' @author Custodian: Arthur Georges -- Post to \url{https://groups.google.com/d/forum/dartr}
#' @examples 
#' \dontrun{
#' #Test run with a focal individual from the Macleay River (EmmacMaclGeor) 
#' test <- gl.assign.pa(testset.gl, unknown='UC_00146', nmin=10, threshold=1,
#' verbose=3) 
#' test_2 <- gl.assign.pca(test, unknown='UC_00146', plevel=0.95, verbose=3)
#' }

gl.assign.pca <- function(x,
                          unknown,
                          plevel = 0.999,
                          plot.out=TRUE,
                          verbose = NULL) {
    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)
    if(verbose==0){
        plot.out <- FALSE
    }

    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     build = "Josh",
                     verbosity = verbose)
    
    
    
    # check if package is installed
    pkg <- "SIBER"
    if (!(requireNamespace(pkg, quietly = TRUE))) {
        stop(error(
            "Package",
            pkg,
            "needed for this function to work. Please install it."
        ))
    }
    requireNamespace(pkg, quietly = TRUE)

    # CHECK DATATYPE
    datatype <- utils.check.datatype(x, verbose = 0)
    if (nPop(x) < 2) {
        stop(
            error(
                "Fatal Error: Only one population, including the unknown, no putative source\n"
            )
        )
    }

    # FUNCTION SPECIFIC ERROR CHECKING

    if (!(unknown %in% indNames(x))) {
        stop(
            error(
                "Fatal Error: Unknown must be listed among the individuals in the genlight object!\n"
            )
        )
    }
    if (plevel > 1 || plevel < 0) {
        cat(warn(
            "  Warning: Value of plevel must be between 0 and 1, set to 0.95\n"
        ))
        plevel <- 0.95
    }

    if (nLoc(x) < nPop(x)) {
        stop(error(
            "Fatal Error: Number of loci less than number of populations\n"
        ))
    }

    # DO THE JOB
    vec <- as.vector(pop(x))
    vec[indNames(x) == unknown] <- "unknown"
    pop(x) <- as.factor(vec)
    
    # Ordinate a reduced space of 2 dimensions
    if (verbose >= 2){
        cat(report("  Calculating a PCA to represent the unknown in the context of putative sources\n"))
    }
    pcoa <- gl.pcoa(x, nfactors = 2, verbose = 0)
    # Plot
    if (verbose >= 0 && plot.out==TRUE) {
        suppressWarnings(suppressMessages(gl.pcoa.plot(
            pcoa, x, ellipse = TRUE, plevel = plevel, verbose=0
        )))  # Because the unknown pop throws an ellipse error
    }
    # Combine Pop names and pca scores
    df <- data.frame(pcoa$scores)
    df <- cbind(as.character(pop(x)), df)
    names(df) <- c("pop", "x", "y")
    # Determine if the unknown lies within the confidence ellipses specified by plevel
    result <- data.frame()
    count <- 0
    for (i in popNames(x)) {
        if (i == "unknown" | nInd(x[pop = i]) <= 1) {
            next
        }
        count <- count + 1
        A <- pcoa$scores[df$pop == i,]
        mu <- colMeans(A)
        sigma <- stats::cov(A)
        testset <- rbind(pcoa$scores[df$pop == "unknown",], A)
        transform <- SIBER::pointsToEllipsoid(testset, sigma, mu)
        inside.or.out <- SIBER::ellipseInOut(transform, p = plevel)
        result[count, 1] <- i
        result[count, 2] <- inside.or.out[1]
    }
    names(result) <- c("pop", "hit")
    nhits <- length(result$pop[result$hit])
    nohits <- length(result$pop[!result$hit])
    if(verbose >= 2){
        cat(report("  Eliminating populations for which the unknown is outside their confidence envelope\n"))
    }
    if (verbose >= 3) {
        if (nhits > 0) {
            cat("  Putative source populations:",
                paste(result[result$hit == TRUE, "pop"], collapse = ", "),
                "\n")
        } else {
            cat("  No putative source populations identified\n")
        }
        if (all(result$hit == TRUE)) {
            cat("  No populations eliminated from consideration\n")
        } else {
            cat(
                "  Populations eliminated from consideration:",
                paste(result[result$hit == FALSE, "pop"], collapse = ", "),
                "\n"
            )
        }
    }
    if (all(result$hit == TRUE)) {
        x2 <- x
    } else {
        x2 <-
            gl.drop.pop(x, pop.list = result[result$hit == FALSE, "pop"], verbose = 0)
    }
    if(verbose >= 2){
        cat(report("  Returning a genlight object with remaining putative source populations plus the unknown\n"))
    }
    
    # FLAG SCRIPT END
    
    if (verbose > 0) {
        cat(report("Completed:", funname, "\n"))
    }
    
    return(x2)
}
