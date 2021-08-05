#'@name gl.pcoa
#'
#'@title Ordination applied to genotypes in a genlight object (PCA), in an fd object, or to a distance matrix (PCoA)
#'
#'@description This function takes the genotypes for individuals and undertakes a Pearson Principal Component analysis (PCA) on SNP or Tag P/A (SilicoDArT) 
#' data; it undertakes a Gower Principal Coordinate analysis (PCoA) if supplied with a distance matrix. Technically, any distance matrix can 
#' be represented in an ordinated space using PCoA.
#'
#' @param x Name of the genlight object or fd object containing the SNP data, or a distance matrix of type dist [required].
#' @param nfactors Number of axes to retain in the output of factor scores [default 5].
#' @param correction Method applied to correct for negative eigenvalues, either 'lingoes' or 'cailliez' [Default NULL].
#' @param mono.rm If TRUE, remove monomorphic loci [default TRUE].
#' @param parallel TRUE if parallel processing is required (does fail under Windows) [default FALSE].
#' @param n.cores Number of cores to use if parallel processing is requested [default 16].
#' @param plot_theme Theme for the plot. See Details for options [default theme_dartR()].
#' @param plot_colours List of two color names for the borders and fill of the plot [default two_colors].
#' @param save2tmp If TRUE, saves any ggplots and listings to the session temporary directory (tempdir) [default FALSE]
#' @param verbose verbose= 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2 or as specified using gl.set.verbosity].
#'
#'@details 
#'The function is essentially a wrapper for glPca {adegenet} or pcoa \{ape\} 
#'with default settings apart from those specified as parameters in this function.
#'\strong{ Sources of stress in the visual representation }
#'
#' While, technically, any distance matrix can be represented in an ordinated space, the representation will not typically be exact.There are three 
#' major sources of stress in a reduced-representation of distances or dissimilarities among entities using PCA or PCoA. By far the greatest
#' source comes from the decision to select only the top two or three axes from the ordinated set of axes derived from the PCA or PCoA. The representation of
#' the entities such a heavily reduced space will not faithfully represent the distances in the input distance matrix simply because of the loss of information
#' in deeper informative dimensions. For this reason, it is not sensible to be too precious about managing the other two sources of stress in
#' the visual representation.
#' 
#' The measure of distance between entities in a PCA is the Pearson Correlation Coefficent, essentially a standardized Euclidean distance. This is both a 
#' metric distance and a Euclidean distance. In PCoA, the second source of stress is the choice of distance measure or dissimilarity measure. While any 
#' distance or dissimilarity matrix can be represented in an ordinated space, the distances between entities can be faithfully represented 
#' in that space (that is, without stress) only if the distances are metric. Furthermore, for distances between entities to be faithfully 
#' represented in a rigid Cartesian space, the distance measure needs to be Euclidean. If this is not the case, 
#' the distances between the entities in the ordinated visualized space will not exactly represent the distances in the input matrix 
#' (stress will be non-zero). This source of stress will be evident as negative eigenvalues in the deeper dimensions. 
#' 
#' A third source of stress arises from having a sparse dataset, one with missing values. This affects both PCA and PCoA. If the original data matrix 
#' is not fully populated, that is, if there are missing values, then even a Euclidean distance matrix will not necessarily be 'positive definite'. 
#' It follows that some of the eigenvalues may be negative, even though the distance metric is Euclidean. This issue is exacerbated when the number 
#' of loci greatly exceeds the number of individuals, as is typically the case when working with SNP data. The impact of missing values can be minimized 
#' by stringently filtering on Call Rate, albeit with loss of data. An alternative is given in a paper 'Honey, I shrunk the sample covariance matrix' 
#' and more recently by Ledoit and Wolf (2018), but their approach has not been implemented here. 
#' 
#' The good news is that, unless the sum of the negative eigenvalues, arising from a non-Euclidean distance measure or from missing values, approaches those 
#' of the final PCA or PCoA axes to be displayed, the distortion is probably of no practical consequence and certainly not comparable to the stress arising from
#' selecting only two or three final dimensions out of several informative dimensions for the visual representation.
#' 
#'\strong{ Function's output }
#'
#' Two diagnostic plots are produced. The first is a Scree Plot, showing the percentage variation explained by each of the PCA or PCoA axes, for those axes that 
#' explain more than the original variables (loci) on average. That is, only informative axes are displayed. The scree plot informs the number of dimensions
#' to be retained in the visual summaries. As a rule of thumb, axes with more than 10% of variation explained should be included.
#' 
#' The second graph shows the distribution of eigenvalues for the remaining uninformative (noise) axes, including those with negative eigenvalues. 
#' Action is recommended (verbose >= 2) if the negative eigenvalues are dominant, their sum approaching in magnitude the eigenvalues for axes selected for 
#' the final visual solution. 
#' 
#' Output is a glPca object conforming to adegenet::glPca but with only the following retained.
#'\itemize{ 
#'\item  $call - The call that generated the PCA/PCoA
#'\item  $eig - Eigenvalues -- All eigenvalues (positive, null, negative).
#'\item  $scores - Scores (coefficients) for each individual
#'\item  $loadings - Loadings of each SNP for each principal component  
#'    }
#' 
#' Plots and table were saved to the temporal directory (tempdir) and can be accessed with the function \code{\link{gl.print.reports}} and listed with the function \code{\link{gl.list.reports}}. Note that they can be accessed only in the current R session because tempdir is cleared each time that the R session is closed.
#'   
#'  Examples of other themes that can be used can be consulted in \itemize{
#'  \item \url{https://ggplot2.tidyverse.org/reference/ggtheme.html} and \item
#'  \url{https://yutannihilation.github.io/allYourFigureAreBelongToUs/ggthemes/}
#'  }
#' 
#' PCA was developed by Pearson (1901) and Hotelling (1933), whilst the best modern reference is Jolliffe (2002). PCoA was developed by Gower (1966) while the
#' best modern reference is Legendre & Legendre (1998).
#' 
#'@return An object of class pcoa containing the eigenvalues and factor scores
#'@author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#'
#'@examples
#' fd <- gl.fixed.diff(testset.gl)
#' fd <- gl.collapse(fd)
#' pca <- gl.pcoa(fd)
#' gl.pcoa.plot(pca,fd$gl)
#'

#'@references
#'\itemize{ 
#'\item Cailliez, F. (1983) The analytical solution of the additive constant problem. Psychometrika, 48, 305-308.
#'\item Gower, J. C. (1966) Some distance properties of latent root and vector methods used in multivariate analysis. Biometrika, 53, 325-338.
#'\item Hotelling, H., 1933. Analysis of a complex of statistical variables into Principal Components. Journal of Educational Psychology 24:417-441, 498-520.
#'\item Jolliffe, I. (2002) Principal Component Analysis. 2nd Edition, Springer, New York. 
#'\item Ledoit, O. and Wolf, M. (2018). Analytical nonlinear shrinkage of large-dimensional covariance matrices. University of Zurich, Department of Economics, Working Paper No. 264, Revised version. Available at SSRN: https://ssrn.com/abstract=3047302 or http://dx.doi.org/10.2139/ssrn.3047302 
#'\item Legendre, P. and Legendre, L. (1998). Numerical Ecology, Volume 24, 2nd Edition. Elsevier Science, NY.
#'\item Lingoes, J. C. (1971) Some boundary conditions for a monotone analysis of symmetric matrices. Psychometrika, 36, 195-203.
#'\item Pearson, K. (1901). On lines and planes of closest fit to systems of points in space. Philosophical Magazine. Series 6, vol. 2, no. 11, pp. 559-572.
#' }
#' 
#'@seealso \code{\link{gl.pcoa.plot}}
#'@family data exploration functions
#'@importFrom ape pcoa
#'@export 

gl.pcoa <- function(x, 
                    nfactors = 5, 
                    correction = NULL, 
                    mono.rm = TRUE, 
                    parallel = FALSE, 
                    n.cores = 16, 
                    plot_theme = theme_dartR(),
                    plot_colours = two_colors, 
                    save2tmp = FALSE,
                    verbose = 2) {

# SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
# FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func=funname,build="Jackson",v=verbose)
  
# CHECK DATATYPE 
  datatype <- utils.check.datatype(x,verbose=verbose)
  if(datatype=="fd"){
    datatype <- utils.check.datatype(x$gl,verbose=0)
    x <- x$gl
  }

# SCRIPT SPECIFIC ERROR CHECKING

    if (mono.rm == TRUE & (datatype=="SNP" | datatype=="SilicoDArT")) {
        x <- gl.filter.monomorphs(x,verbose=0)
    }

    if (is.null(correction)) {
        correction <- "none"
    } else {
        correction <- tolower(correction)
        if (correction != "lingoes" && correction != "cailliez") {
            if (verbose >= 2) {
                cat(warn("  Warning: Correction if specified needs to be lingoes or cailliez, set to the default 'None'"))
            }
            correction <- "none"
        }
    }

    # DO THE JOB

    ######## DISTANCE ANALYSIS

    if (datatype == "dist"){
            D <- x
            
            # Calculate the pcoa
            if (verbose >= 2) {
                if (correction == "none") {
                  cat(report("  Performing a PCoA, individuals as entities, no correction applied\n"))
                  title <- "PCoA on Distance Matrix (no correction)\nScree Plot (informative axes only)"
                } else {
                  cat(report("  Performing a PCoA, individuals as entities,", correction, "correction for stress (-ive eigenvalues) applied\n"))
                  title <- paste0("PCoA on Distance Matrix(", correction, ")\nScree Plot (informative axes only)")
                }
            }
            pco <- ape::pcoa(D, correction = correction, rn = labels(D))

            # Extract relevant variables

            if (correction == "none") {
                eig.raw <- pco$values$Eigenvalues
            } else {
                eig.raw <- pco$values$Corr_eig
            }

            # Identify the number of axes with explanatory value greater than the original variables on average
            eig.raw.pos <- eig.raw[eig.raw >= 0]
            eig.raw.pos.pc <- eig.raw.pos * 100/sum(eig.raw.pos)
            eig.top <- eig.raw.pos[eig.raw.pos > mean(eig.raw.pos)]
            eig.top.pc <- round(eig.top * 100/sum(eig.raw.pos), 1)
            eig.raw.noise <- eig.raw[eig.raw <= mean(eig.raw)]

            if (any(eig.raw < 0)) {
                if (verbose >= 2) {
                  problem <- (-sum(eig.raw[eig.raw < 0])/mean(eig.raw[1:3])) * 100
                  cat(warn("  Warning: Some eigenvalues negative -- sum to", round(problem, 2), "% of the mean eigenvalue for PCoA axes 1-3\n"))
                  cat(report("    Tolerable negative eigenvalues should sum to much less than the eigenvalues of displayed PCoA axes (say, less than 20%)\n"))
                  if (problem > 20) {
                    cat(report("    If the stress (negative eigenvalues) is considered a problem, and you might reasonably choose to ignore it, you have the following options:\n"))
                    cat(report("    (a) Apply more stringent filtering on Call Rate before generating the distance matrix; or\n"))
                    cat(report("    (b) Select an alternate distance measure, preferably a metric distance or better still, a Euclidean distance, if you have not already; or\n"))
                    cat(report("    (c) Apply a transformation (correction) to eliminate the negative eigenvalues. If this was already done, try another correction; or\n"))
                    cat(report("    (d) Interpret the visual representation of the ordination with caution, seeking corroborating evidence.\n"))
                  }
                }
            }

            # Provide a summary
            if (verbose >= 3) {
                if (correction == "lingoes" | correction == "cailliez") {
                  cat(" Correction", correction, "applied to remove negative eigenvalues\n")
                  cat(paste("  Uncorrected ordination yielded", length(eig.top), "informative dimensions from", nInd(x) - 1, "original dimensions\n"))
                } else {
                  cat(paste("  Ordination yielded", length(eig.top), "informative dimensions from", dim(as.matrix(D))[1], "original dimensions\n"))
                }
                cat(paste("    PCoA Axis 1 explains", round(eig.raw.pos.pc[1], 1), "% of the total variance\n"))
                cat(paste("    PCoA Axis 1 and 2 combined explain", round(eig.raw.pos.pc[1] + eig.raw.pos.pc[2], 1), "% of the total variance\n"))
                cat(paste("    PCoA Axis 1-3 combined explain", round(eig.raw.pos.pc[1] + eig.raw.pos.pc[2] + eig.raw.pos.pc[3], 1), "% of the total variance\n"))
            }
            # Construct a universal output file
            p.object <- list()
            p.object$scores <- pco$vectors[, 1:nfactors]
            p.object$eig <- pco$values$Eigenvalues
            p.object$loadings <- pco$vectors.cor[, 1:nfactors]
            p.object$call <- match.call()

        }  ######## END DISTANCE DATA

    ######## SNP or P/A DATA, PCA

    if (datatype == "SNP" || datatype == "SilicoDArT"){

            if (verbose >= 2) {
                if (datatype == "SNP") {
                  cat(report("  Performing a PCA, individuals as entities, loci as attributes, SNP genotype as state\n\n"))
                  title <- "PCA on SNP Genotypes\nScree Plot (informative axes only)"
                }
                if (datatype == "SilicoDArT") {
                  cat(report("  Performing a PCA, individuals as entities, loci as attributes, Tag P/A as state\n\n"))
                  title <- "PCA on Tag P/A Data\nScree Plot (informative axes only)"
                }
            }
            pca <- glPca(x, nf = nfactors, parallel = parallel, n.cores = n.cores)

            # Identify the number of axes with explanatory value greater than the original variables on average
            eig.raw <- pca$eig

            eig.raw.pos <- eig.raw[eig.raw >= 0]
            eig.raw.pos.pc <- eig.raw.pos * 100/sum(eig.raw.pos)
            eig.top <- eig.raw.pos[eig.raw.pos >= mean(eig.raw.pos)]
            eig.top.pc <- round(eig.top * 100/sum(eig.raw.pos), 1)
            eig.raw.noise <- eig.raw[eig.raw <= mean(eig.raw)]

            if (any(eig.raw < 0)) {
                if (verbose >= 2) {
                  problem <- (-sum(eig.raw[eig.raw < 0])/mean(eig.raw[1:3])) * 100
                  cat(warn("  Warning: Some eigenvalues negative -- sum to", round(problem, 2), "% of the mean eigenvalue for PCA axes 1-3\n"))
                  cat(report("    Tolerable negative eigenvalues should sum to much less than the eigenvalues of displayed PCA axes (say, less than 20%)\n"))
                  if (problem > 20) {
                    cat(report("    If the stress (negative eigenvalues) is considered a problem, and you might reasonably choose to ignore it, you have the following options:\n"))
                    cat(report("    (a) Apply more stringent filtering on Call Rate and repeat the PCA; or\n"))
                    cat(report("    (b) Undertake a PCoA with an appropriate distance measure and a transformation (correction) to eliminate the negative eigenvalues; or\n"))
                    cat(report("    (c) Interperate the visual representation of the ordination with caution, seeking corroborating evidence.\n"))
                  }
                }
            }

            e <- pca$eig[pca$eig > sum(pca$eig/length(pca$eig))]
            e <- round(e * 100/sum(pca$eig), 1)
            if (verbose >= 3) {
                cat(paste("  Ordination yielded", length(e), "informative dimensions from", nInd(x) - 1, "original dimensions\n"))
                cat(paste("    PCA Axis 1 explains", e[1], "% of the total variance\n"))
                cat(paste("    PCA Axis 1 and 2 combined explain", e[1] + e[2], "% of the total variance\n"))
                cat(paste("    PCA Axis 1-3 combined explain", e[1] + e[2] + e[3], "% of the total variance\n"))
            }
            # Construct a universal output file
            p.object <- list()
            p.object$scores <- pca$scores
            p.object$eig <- pca$eig
            p.object$loadings <- pca$loadings
            p.object$call <- match.call()

        }  #### END SNP ANALYSIS
    
      # PLOT THE DIAGNOSTICS
    
    # Plot Scree plot
    #avoid no visible binding probl
    eigenvalue <- percent <-NULL
    
    df <- data.frame(eigenvalue=seq(1:length(eig.top.pc)), percent=eig.top.pc)
    if (datatype == "SNP") {
        xlab <- paste("PCA Axis")
    } else {
        xlab <- paste("PCoA Axis")
    }
    ylab <- paste("Percentage Contribution")
    
    p1 <- ggplot(df,aes(x= eigenvalue, y=percent )) +
            geom_line( color=plot_colours[2]) +
            geom_point(color=plot_colours[1],size=4 ) +
            geom_hline(yintercept = 10, color="blue" ) +
            plot_theme +
            xlab(xlab) +
            ylab(ylab) +
            ggtitle(title)

    if (any(eig.raw < 0)) {
        main <- "Noise Axes -- Warning: some eigenvalues < 0"
    } else {
        main <- "Noise Axes -- all eigenvalues positive"
    }
    
    p2 <- ggplot(as.data.frame(eig.raw.noise),aes(x= eig.raw.noise)) +
        geom_histogram(bins = 50, color = plot_colours[1], fill = plot_colours[2]) +
        geom_vline(xintercept = 0, color="blue" ) +
        plot_theme +
        xlab("Eigenvalue") +
        ylab("Count") +
        ggtitle(main)
    
    # printing outputs
    p3 <- (p1/p2)
    print(p3)

    # SAVE INTERMEDIATES TO TEMPDIR 
    if(save2tmp){
      # creating temp file names
      temp_plot <- tempfile(pattern = "Plot_")
      match_call <- paste0(names(match.call()),"_",as.character(match.call()),collapse = "_")
      # saving to tempdir
      saveRDS(list(match_call,p3), file = temp_plot)
      if(verbose>=2){
        cat(report("  Saving ggplot(s) to the session tempfile\n"))
        cat(report("  NOTE: Retrieve output files from tempdir using gl.list.reports() and gl.print.reports()\n"))
      }
    }

    # FLAG SCRIPT END

    if (verbose > 0) {
        cat(report("Completed:", funname, "\n\n"))
    }

     class(p.object) <- "glPca"
     invisible(p.object)
}

