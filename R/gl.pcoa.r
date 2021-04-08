#' Ordination applied to genotypes in a genlight object (PCA), in an fd object, or to a distance matrix (PCoA)
#'
#' This function takes the genotypes for individuals and undertakes a Pearson Principal Component analysis (PCA) on SNP or Tag P/A (SilicoDArT) 
#' data; it undertakes a Gower Principal Coordinate analysis (PCoA) if supplied with a distance matrix. Technically, any distance matrix can 
#' be represented in an ordinated space using PCoA.
#' 
#' The function is essentially a wrapper for glPca {adegenet} or pcoa \{ape\} with default settings apart from those specified as parameters in this 
#' function.
#' 
#' While, technically, any distance matrix can be represented in an ordinated space, the representation will not typically be exact.There are three 
#' major sources of stress in a reduced-reprentation of distances or dissimilarities among entities using PCA or PCoA. By far the greatest
#' source comes from the decision to select only the top two or three axes from the ordinated set of axes derived from the PCA or PCoA. The representation of
#' the entities such a heavily reduced space will not faithfully represent the distances in the input distance matrix simply because of the loss of information
#' in deeper informative dimensions. For this reason, it is not sensible to be too precious about managing the other two sources of stress in
#' the visual representation.
#' 
#' The measure of distance between entities in a PCA is the Pearson Correlation Coefficent, essentially a standardized Euclidean distance. This is both a 
#' metric distance and a Euclidean distance. In PCoA, the second source of stress is the choice of distance measure or dissimilarity measure. While any 
#' distance or dissimilarity matrix can be represented in an ordinated space, the distances between entities can befaithfully represented 
#' in that space (that is, without stress) only if the distances are metric. Furthermore, for distances between entities to be faithfully 
#' represented in a rigid Cartesian space, the distance measure needs to be Euclidean. If this is not the case, 
#' the distances between the entities in the ordinated visualized space will not exactly represent the distances in the input matrix 
#' (stress will be non-zero). This source of stress will be evident as negative eigenvalues in the deeper dimensions. 
#' 
#' A third source of stress arises from having a sparse dataset, one with missing values. This affects both PCA and PCoA. If the original data matrix 
#' is not fully populated, that is, if there are missing values, then even a Euclidean distance matrix will not necessarily be 'positive definite'. 
#' It follows that some of the eigenvalues may be negative, even though the distance metric is Euclidean. This issue is exacerbated when the number 
#' of loci greatly exceeds the number of individuals, as is typically the case when working with SNP data. The impact of missing values can be minimized 
#' by stringently filtering on Call Rate, albeit with loss of data. An alternative is given in a paper "Honey, I shrunk the sample covariance matrix" 
#' and more recently by Ledoit and Wolf (2018), but their approach has not been implemented here. 
#' 
#' The good news is that, unless the sum of the negative eigenvalues, arising from a non-Euclidean distance measure or from missing values, approaches those 
#' of the final PCA or PCoA axes to be displayed, the distortion is probably of no practical consequence and certainly not comparable to the stress arising from
#' selecting only two or three final dimensions out of several informative dimensions for the visual representation.
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
#'  $call
#'    The call that generated the PCA/PCoA
#'  $eig
#'    Eigenvalues	-- All eigenvalues (positive, null, negative).
#'  $scores 
#'    Scores (coefficients) for each individual
#'  $loadings
#'    Loadings of each SNP for each principal component  
#' 
#' PCA was developed by Pearson (1901) and Hotelling (1933), whilst the best modern reference is Jolliffe (2002). PCoA was developed by Gower (1966) while the
#' best modern reference is Legendre & Legendre (1998).
#' 
#' @references
#' Cailliez, F. (1983) The analytical solution of the additive constant problem. Psychometrika, 48, 305-308.
#' Gower, J. C. (1966) Some distance properties of latent root and vector methods used in multivariate analysis. Biometrika, 53, 325-338.
#' Hotelling, H., 1933. Analysis of a complex of statistical variables into Principal Components. Journal of Educational Psychology 24:417-441, 498-520.
#' Jolliffe, I. (2002) Principal Component Analysis. 2nd Edition, Springer, New York. 
#' Ledoit, O. and Wolf, M. (2018). Analytical nonlinear shrinkage of large-dimensional covariance matrices. University of Zurich, Department of Economics, Working Paper No. 264, Revised version. Available at SSRN: https://ssrn.com/abstract=3047302 or http://dx.doi.org/10.2139/ssrn.3047302 
#' Legendre, P. and Legendre, L. (1998). Numerical Ecology, Volume 24, 2nd Edition. Elsevier Science, NY.
#' Lingoes, J. C. (1971) Some boundary conditions for a monotone analysis of symmetric matrices. Psychometrika, 36, 195-203.
#' Pearson, K. (1901). On lines and planes of closest fit to systems of points in space. Philosophical Magazine. Series 6, vol. 2, no. 11, pp. 559-572.
#'  
#' @param x -- name of the genlight object or fd object containing the SNP data, or a distance matrix of type dist [required]
#' @param nfactors -- number of axes to retain in the output of factor scores.
#' @param correction Method applied to correct for negative eigenvalues, either 'lingoes' or 'cailliez' [Default NULL]
#' @param parallel TRUE if parallel processing is required (does fail under Windows) [default FALSE]
#' @param n.cores Number of cores to use if parallel processing is requested [default 16]
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2 or as specified using gl.set.verbosity]
#' @return An object of class pcoa containing the eigenvalues and factor scores
#' @importFrom ape pcoa
#' @export
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' fd <- gl.fixed.diff(testset.gl)
#' fd <- gl.collapse(fd)
#' pca <- gl.pcoa(fd)
#' gl.pcoa.plot(pca,fd)

gl.pcoa <- function(x, nfactors=5, correction=NULL, parallel=FALSE, n.cores=16, verbose=NULL) {

# TRAP COMMAND, SET VERSION
  
  funname <- match.call()[[1]]
  build <- "Jacob"
  
# SET VERBOSITY
  
  if(class(x)=="genlight"){
    if (is.null(verbose)){ 
      if(!is.null(x@other$verbose)){ 
        verbose <- x@other$verbose
      } else { 
        verbose <- 2
      }
    }
  }
  if(class(x)=="fd"){
    x <- x$gl
    if (is.null(verbose)){
      verbose <- 2
    }  
  }
  if(class(x)=="dist"){
    if (is.null(verbose)){
      verbose <- 2
    }  
  }
  
  if (verbose < 0 | verbose > 5){
    cat(paste("  Warning: Parameter 'verbose' must be an integer between 0 [silent] and 5 [full report], set to 2\n"))
    verbose <- 2
  }
  
# FLAG SCRIPT START
  
  if (verbose >= 1){
    if(verbose==5){
      cat("Starting",funname,"[ Build =",build,"]\n")
    } else {
      cat("Starting",funname,"\n")
    }
  }
  
# STANDARD ERROR CHECKING
  
  if(class(x)=="genlight") {
    
    if (all(x@ploidy == 1)){
      if (verbose >= 2){cat("  Processing  Presence/Absence (SilicoDArT) data\n")}
      data.type <- "SilicoDArT"
    } else if (all(x@ploidy == 2)){
      if (verbose >= 2){cat("  Processing a SNP dataset\n")}
      data.type <- "SNP"
    } else {
      stop("Fatal Error: Ploidy must be universally 1 (fragment P/A data) or 2 (SNP data)")
    }
  } else if(class(x)=='dist'){
    if (verbose >= 2){cat("  Processing a Distance Matrix, D \n")}
    data.type <- "dist"
  } else if(class(x)=='fd'){
    if (verbose >= 2){cat("  Processing a genlight object after a fixed difference analysis, x$fd \n")}
    data.type <- "SNP"
  } else {
    stop("Fatal Error: Expecting either a genlight object (SNP or SilicoDArT) or a distance matrix\n")
  }
  
# SCRIPT SPECIFIC ERROR CHECKING
  
  if (is.null(correction)) {
    correction <- "none"
  } else {
    correction <- tolower(correction)
    if (correction != "lingoes" && correction != "cailliez"){
      if(verbose >= 2){cat("  Warning: Correction if specified needs to be lingoes or cailliez, set to the default 'None'")}
      correction <- "none"
    }
  }  
  
# DO THE JOB
  
######## DISTANCE ANALYSIS
  
if (data.type == "dist"){
  
  D <- x

  # Calculate the pcoa
  if (verbose >= 2) {
    if (correction=='none'){
      cat("  Performing a PCoA, individuals as entities, no correction applied\n")
      title <- "PCoA on Distance Matrix (no correction)\nScree Plot (informative axes only)"
    } else {
      cat("  Performing a PCoA, individuals as entities,",correction,"correction for stress (-ive eigenvalues) applied\n")
      title <- paste0("PCoA on Distance Matrix(",correction,")\nScree Plot (informative axes only)")
    }
  }
    pco <- ape::pcoa(D,correction=correction,rn=labels(D))
    
  # Extract relevant variables
    
    if (correction=='none'){
      eig.raw <- pco$values$Eigenvalues
    } else {
      eig.raw <- pco$values$Corr_eig
    }
  
  # Identify the number of axes with explanatory value greater than the original variables on average
    eig.raw.pos <- eig.raw[eig.raw >= 0] 
    eig.raw.pos.pc <- eig.raw.pos*100/sum(eig.raw.pos)
    eig.top <- eig.raw.pos[eig.raw.pos > mean(eig.raw.pos)]
    eig.top.pc <- round(eig.top*100/sum(eig.raw.pos),1)
    eig.raw.noise <- eig.raw[eig.raw <= mean(eig.raw)]

    if (any(eig.raw < 0)){
      if(verbose >= 2){
        problem <- (-sum(eig.raw[eig.raw<0])/mean(eig.raw[1:3]))*100
          cat("  Warning: Some eigenvalues negative -- sum to",round(problem,2),"% of the mean eigenvalue for PCoA axes 1-3\n")
          cat("    Tolerable negative eigenvalues should sum to much less than the eigenvalues of displayed PCoA axes (say, less than 20%)\n")
          if (problem > 20){
            cat("    If the stress (negative eigenvalues) is considered a problem, and you might reasonably choose to ignore it, you have the following options:\n")
            cat("    (a) Apply more stringent filtering on Call Rate before generating the distance matrix; or\n")
            cat("    (b) Select an alternate distance measure, preferably a metric distance or better still, a Euclidean distance, if you have not already; or\n")
            cat("    (c) Apply a transformation (correction) to eliminate the negative eigenvalues. If this was already done, try another correction; or\n")
            cat("    (d) Interpret the visual representation of the ordination with caution, seeking corroborating evidence.\n")
          }
      }  
    }
    
    # Provide a summary  
    if (verbose >=3) {
      if(correction == "lingoes" | correction == "cailliez"){
        cat(" Correction",correction,"applied to remove negative eigenvalues\n")
        cat(paste("  Uncorrected ordination yielded",length(eig.top),"informative dimensions from",nInd(x)-1,"original dimensions\n"))
      } else {
        cat(paste("  Ordination yielded",length(eig.top),"informative dimensions from",dim(as.matrix(D))[1],"original dimensions\n"))
      }  
      cat(paste("    PCoA Axis 1 explains",round(eig.raw.pos.pc[1],1),"% of the total variance\n"))
      cat(paste("    PCoA Axis 1 and 2 combined explain",round(eig.raw.pos.pc[1]+eig.raw.pos.pc[2],1),"% of the total variance\n"))
      cat(paste("    PCoA Axis 1-3 combined explain",round(eig.raw.pos.pc[1]+eig.raw.pos.pc[2]+eig.raw.pos.pc[3],1),"% of the total variance\n"))
    }
    # Construct a universal output file
    p.object <- list()
    p.object$scores <- pco$vectors[,1:nfactors]
    p.object$eig <- pco$values$Eigenvalues
    p.object$loadings <- pco$vectors.cor[,1:nfactors]
    p.object$call <- match.call()
    
} ######## END DISTANCE DATA
  
######## SNP or P/A DATA, PCA  
  
if (data.type == "SNP" || data.type == "SilicoDArT"){  
  
  if (verbose >= 2) {
    if (data.type == "SNP"){
      cat("  Performing a PCA, individuals as entities, loci as attributes, SNP genotype as state\n")
      title <- "PCA on SNP Genotypes\nScree Plot (informative axes only)"
    }
    if (data.type == "SilicoDArT"){
      cat("  Performing a PCA, individuals as entities, loci as attributes, Tag P/A as state\n")
      title <- "PCA on Tag P/A Data\nScree Plot (informative axes only)"
    }
  }
  pca <- glPca(x, nf=nfactors, parallel=parallel, n.cores=n.cores)
  
  # Identify the number of axes with explanatory value greater than the original variables on average
  eig.raw <- pca$eig
  
  eig.raw.pos <- eig.raw[eig.raw >= 0] 
  eig.raw.pos.pc <- eig.raw.pos*100/sum(eig.raw.pos)
  eig.top <- eig.raw.pos[eig.raw.pos > mean(eig.raw.pos)]
  eig.top.pc <- round(eig.top*100/sum(eig.raw.pos),1)
  eig.raw.noise <- eig.raw[eig.raw <= mean(eig.raw)]
  
  if (any(eig.raw < 0)){
    if(verbose >= 2){
      problem <- (-sum(eig.raw[eig.raw<0])/mean(eig.raw[1:3]))*100
      cat("  Warning: Some eigenvalues negative -- sum to",round(problem,2),"% of the mean eigenvalue for PCA axes 1-3\n")
      cat("    Tolerable negative eigenvalues should sum to much less than the eigenvalues of displayed PCA axes (say, less than 20%)\n")
      if (problem > 20){
        cat("    If the stress (negative eigenvalues) is considered a problem, and you might reasonably choose to ignore it, you have the following options:\n")
        cat("    (a) Apply more stringent filtering on Call Rate and repeat the PCA; or\n")
        cat("    (b) Undertake a PCoA with an appropriate distance measure and a transformation (correction) to eliminate the negative eigenvalues; or\n")
        cat("    (c) Interperate the visual representation of the ordination with caution, seeking corroborating evidence.\n")
      }
    }  
  }
  
  e <- pca$eig[pca$eig > sum(pca$eig/length(pca$eig))]
  e <- round(e*100/sum(pca$eig),1)
  if (verbose >=3) {
    cat(paste("  Ordination yielded",length(e),"informative dimensions from",nInd(x)-1,"original dimensions\n"))
    cat(paste("    PCA Axis 1 explains",e[1],"% of the total variance\n"))
    cat(paste("    PCA Axis 1 and 2 combined explain",e[1]+e[2],"% of the total variance\n"))
    cat(paste("    PCA Axis 1-3 combined explain",e[1]+e[2]+e[3],"% of the total variance\n"))
  }
  # Construct a universal output file
    p.object <- list()
    p.object$scores <- pca$scores
    p.object$eig <- pca$eig
    p.object$loadings <- pca$loadings
    p.object$call <- match.call()
    
}  #### END SNP ANALYSIS
  
# PLOT THE DIAGNOSTICS
  
    # Save the prior settings for mfrow, oma, mai and pty, and reassign
    op <- par(mfrow = c(2, 1), oma=c(1,1,1,1), mai=c(0.5,0.5,0.5,0.5),pty="m")
    
    # Set margins for first plot
    par(mai=c(1,1,0.5,0.5))
    
    # Plot Scree plot

    m <- cbind(seq(1:length(eig.top.pc)),eig.top.pc)
    df <- data.frame(m)
    colnames(df) <- c("eigenvalue","percent")
    if (data.type == "SNP"){
      xlab <- paste("PCA Axis")
    } else {
      xlab <- paste("PCoA Axis")
    }  
    ylab <- paste("Percentage Contribution")
    plot(df$eigenvalue,
         df$percent,col='red', 
         main=title, 
         type="o",
         xlab=xlab,
         ylab=ylab)
    abline(h=10,col="blue")

    # Set margins for second plot
    par(mai=c(1,1,0.5,0.5))

      if (any(eig.raw < 0)){
          main <- "Noise Axes -- Warning: some eigenvalues < 0"
      } else {
          main <- "Noise Axes -- all eigenvalues positive"
      }  
      hist(eig.raw.noise, 
           main=main, 
           xlab="Eigenvalue", 
           col="red",
           breaks=100)
      abline(v=0,col="blue",lwd=2)

# FLAG SCRIPT END

  if (verbose > 0) {
    cat("Completed:",funname,"\n")
  }
    
  # Reset the par options    
    par(op)

    class(p.object) <- "glPca"
    return(p.object)
}

