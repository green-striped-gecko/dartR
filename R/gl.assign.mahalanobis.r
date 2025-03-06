#' @name gl.assign.mahalanobis
#' @title Assign an individual of unknown provenance to population based on 
#' Mahalanobis Distance
#' @description
#' This script assigns an individual of unknown provenance to one or more target
#' populations based on the unknown individual's proximity to population 
#' centroids; proximity is estimated using Mahalanobis Distance. 
#'
#' The following process is followed:
#' \enumerate{
#' \item An ordination is undertaken on the populations to again yield a
#' series of orthogonal (independent) axes.
#' \item A workable subset of dimensions is chosen, that specified, or
#' equal to the number of dimensions with substantive eigenvalues, whichever is
#' the smaller.
#' \item The Mahalobalis Distance is calculated for the unknown against each
#' population and probability of membership of each population is calculated.
#' The assignment probabilities are listed in support of a decision.
#' }
#' @details
#' There are three considerations to assignment. First, consider only those
#' populations for which the unknown has no private alleles. Private alleles are
#' an indication that the unknown does not belong to a target population
#' (provided that the sample size is adequate, say >=10). This can be evaluated
#'  with gl.assign.pa().
#'
#' A next step is to consider the PCoA plot for populations where no private
#' alleles have been detected. The position of the unknown in relation to the
#' confidence ellipses is plotted by this script as a basis for narrowing down
#' the list of putative source populations. This can be evaluated with 
#' gl.assign.pca().
#' 
#' The third step (delivered by this script) is to consider the assignment 
#' probabilities based on the squared Generalised Linear Distance 
#' (Mahalanobis distance) of the unknown from the centroid for each population, 
#' then to consider the probability associated with its quantile using the 
#' Chisquare approximation. In effect, this index takes into account position 
#' of the unknown in relation to the confidence envelope in all selected 
#' dimensions of the ordination. The larger the assignment probability, 
#' the greater the confidence in the assignment. 
#' 
#' If dim.limit is set to 2, to correspond with the dimensions used in
#' gl.assign.pa(), then the output provides a ranking of the final set
#' of putative source populations.
#' 
#' If dim.limit is set to be > 2, then this script provides a basis for
#' further narrowing the set of putative populations.If the unknown individual
#' is an extreme outlier, say at less than 0.001 probability of population 
#' membership (0.999 confidence envelope), then the associated population 
#' can be eliminated from further consideration.
#' 
#' Warning: gl.assign.mahal() treats each specified dimension equally, without
#' regard to the percentage variation explained after ordination. If the 
#' unknown is an outlier in a lower dimension with an explanatory variance of,
#' say, 0.1%, the putative population will be eliminated. Use only substantive
#' dimensions from the ordination.
#'
#' Each of these above approaches provides evidence, none are 100% definitive. 
#' They need to be interpreted cautiously.
#' 
#' In deciding the assignment, the script considers an individual to be an
#' outlier with respect to a particular population at alpha = 0.001 as default
#'
#' @param x Name of the input genlight object [required].
#' @param unknown Identity label of the focal individual whose provenance is
#' unknown [required].
#' @param dim.limit Maximum number of dimensions to consider for the
#' confidence ellipses [default 2]
#' @param plevel Probability level for bounding ellipses
#' [default 0.999].
#' @param plot.out If TRUE, produces a plot showing the position of the 
#' unknown in relation to putative source populations [default TRUE]
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2 or as specified using gl.set.verbosity].
#'
#' @return A data frame with the results of the assignment analysis. 
#'
#' @importFrom stats dnorm qnorm
#' @export
#'
#' @author Custodian: Arthur Georges --
#' Post to \url{https://groups.google.com/d/forum/dartr}
#'
#' @examples 
#' \dontrun{
#' #Test run with a focal individual from the Macleay River (EmmacMaclGeor) 
#' test <- gl.assign.pa(testset.gl, unknown='UC_01044', nmin=10, threshold=1,
#' verbose=3) 
#' test_2  <- gl.assign.pca(test, unknown='UC_01044', plevel=0.95, verbose=3)
#' df <- gl.assign.mahalanobis(test_2, unknown='UC_01044', verbose=3)
#' }

gl.assign.mahalanobis <- function(x,
                                  dim.limit=2,
                                  plevel=0.999,
                                  plot.out=TRUE,
                                  unknown,
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
                     verbose = verbose)
    
    # CHECK DATATYPE
    datatype <- utils.check.datatype(x, verbose = 0)
    
    # FUNCTION SPECIFIC ERROR CHECKING
    
    if (!(unknown %in% indNames(x))) {
        stop(
            error(
                "Fatal Error: Unknown must be listed among the individuals in 
                the genlight object!\n"
            )
        )
    }
    
    if (plevel > 1 || plevel < 0) {
        cat(warn(
            "  Warning: Value of plevel must be between 0 and 1, set to 0.999\n"
        ))
        plevel <- 0.999
    }
    
    if (nLoc(x) < nPop(x)) {
        stop(error(
            "Fatal Error: Number of loci less than number of populations"
        ))
    }
    
    # DO THE JOB
    if (nPop(x) < 2) {
        if(verbose >= 2){
            cat(warn("  Only one population, including the unknown, no putative 
                     source identified.\n"))
        }
        # Add a row
        df[i,1] <- unknown
        df[i,2] <- NA
        df[i,3] <- NA
        df[i,4] <- NA
        df[i,5] <- NA
        df[i,6] <- NA
        # FLAG SCRIPT END
        if (verbose > 0) {cat(report("Completed:", funname, "\n"))}
        return(df)
        
    } else {
        
    vec <- as.vector(pop(x))
    vec[indNames(x) == unknown] <- "unknown"
    pop(x) <- as.factor(vec)
    
    # Run the pcoa 
    if (nInd(x) < 2) {
        df <- NULL
    } else {
        pcoa <- gl.pcoa(x,nfactors=dim.limit,verbose=0)
        if(plot.out==TRUE){
          if(verbose >= 2){
            cat(report("  Plotting the unknown against putative populations\n"))
          }
            suppressWarnings(suppressMessages(gl.pcoa.plot(pcoa,
                                                           x,
                                                           ellipse=TRUE,
                                                           plevel=plevel,
                                                           verbose=0)))
        }    

        # Determine the number of dimensions for confidence envelope (the 
        # ordination and dimension reduction) From the eigenvalue
        # distribution
        s <- sum(pcoa$eig)
        e <- round(pcoa$eig * 100 / s, 1)
        e.sign <- e[e > mean(e,na.rm=TRUE)]
        first.est <- length(e.sign)
        # From the number of populations, including the unknown 
        # sec.est <- nPop(x)
        
        # cat(' Number of populations, including the unknown:',sec.est,'\n')
        if (verbose >= 2) {
            cat(
                report(
                    "  Number of dimensions with substantial eigenvalues:",
                    first.est,
                    ". Hardwired limit",
                    dim.limit,
                    "\n"
                )
            )
            cat(report("    Selecting the smallest of the two\n"))
        }
        dim <- min(first.est, dim.limit)
        if (verbose >= 2) {
            cat(report("    Dimension of confidence envelope set at", dim,
                       "\n"))
        }
        pcoa$scores <- pcoa$scores[, 1:dim]
        
        # Add population names to the scores
        c <-
            data.frame(cbind(pcoa$scores, as.character(pop(x))), 
                       stringsAsFactors = FALSE)
        colnames(c)[dim + 1] <- "pop"
        
        # Create a set of data without the unknown
        clouds <- c[c[, "pop"] != "unknown",]
        Unknown <- c[c[, "pop"] == "unknown",]

        # For each population
        p <- as.factor(unique(clouds[, "pop"]))
        for (i in 1:length(levels(p))) {
            # Pull out population i
            m <- clouds[clouds[, "pop"] == levels(p)[i],]
            # Discard the population labels
            m <- m[, 1:dim]
            hold <- row.names(m)
            row.names(m) <- NULL
            Unknown <- Unknown[1:dim]
            row.names(Unknown) <- NULL
            # Convert to numeric, and reformat as a matrix
            n <- as.matrix(sapply(m, as.numeric))
            Unknown <- as.numeric(Unknown)
            # Calculate statistics for the data without the unknown
            means <- colMeans(n)
            covariance <- stats::cov(n)
            # Add back in the unknown
            all <- rbind(n,Unknown)
            # Calculate Mahalanobis Distances
            D <- stats::mahalanobis(all, means, covariance, toll=1e-20)
#            wtD <- WMDB::wmahalanobis(all,means,covariance,weight=e)
            names(D) <- c(hold,"unknown")
            # Calculate the associated probabilities
            pval <- (pchisq(D, df=length(D)-1, lower.tail=FALSE))
            # Is the result non-significant, then assign=yes
            if (pval["unknown"] >= 1-plevel) {
                assign <- "yes"
            } else {
                assign <- "no"
            }
            # Create a dataframe to hold the results
            if(i==1){
              df <- data.frame(unknown=NA,pop=NA,MahalD=NA,pval=NA,critval=NA,
                               assign="NA")
 #             df <- data.frame(unknown=NA,pop=NA,MahalD=NA,wtMahalD=NA,pval=NA,
              # critval=NA,assign="NA")
            }
            # # Add a row
            # df[i,1] <- unknown
            # df[i,2] <- levels(p)[i]
            # df[i,3] <- D["unknown"]
            # df[i,4] <- wtD["unknown"]
            # df[i,5] <- pval["unknown"]
            # df[i,6] <- 1 - plevel
            # df[i,7] <- assign
            # Add a row
            df[i,1] <- unknown
            df[i,2] <- levels(p)[i]
            df[i,3] <- D["unknown"]
            df[i,4] <- pval["unknown"]
            df[i,5] <- 1 - plevel
            df[i,6] <- assign
        }
        # Order the dataframe in descending order on pval
        df <- df[order(df$pval,decreasing=TRUE),]
        # Display
        if (verbose >= 3) {
            print(df)
        }
        # Extract the best match, and report
        best <- as.character(df[df$assign == "yes",][1,2])
       
        if (verbose >= 3) {
            cat("  Best assignment is the population with the larges probability
                of assignment, in this case",
                    best,
                    "\n"
                )
        }
    }
    if(verbose >= 2){
      cat(report("  Returning a dataframe with the Mahalanobis Distances\n"))
    }
    
    # FLAG SCRIPT END
    
    if (verbose > 0) {
        cat(report("Completed:", funname, "\n"))
    }
  }
    
    return(df)
}
