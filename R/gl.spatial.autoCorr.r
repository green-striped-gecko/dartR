#' @name gl.spatial.autoCorr
#' @title Spatial autocorrelation following Smouse and Pekall 1999
#' 
#' @description  Global spatial autocorrelation is a multivariate approach combining all loci
#' into a single analysis. The autocorrelation coefficient r is calculated for
#' each pairwise genetic distance pairs for all specified distance classes. For
#' more information see Smouse and Peakall 1999, Peakall et a. 2003 and Smouse
#' et al. 2008.
#' 
#' @details The code of this function execute behind the scene a modified 
#' version of \code{spautocorr} from the package 
#' \code{PopGenReport}, and it also includes 
#' bootstraps to compute the 95% confidence intervals around the r 
#' estimates, the one-tail test, and the correction factor described by 
#' Peakall et al 2003.
#' 
#' If \code{bins} is of length = 1 it is interpreted as the number of (even) 
#' bins to use. In this case the starting point is always the minimum 
#' value in the distance matrix, and the last is the maximum. If it is a numeric
#' vector, it is interpreted as the braking points. In this case, the first has 
#' to be the lowest value, and the last has to be the highest. The number of 
#' bins will be  \code{length(bins) - 1}.
#' 
#' The permutation constructs the 95% confidence intervals around the null 
#' hypothesis of no spatial structure (this is a two-tail test). The same data
#' are also used to calculate the probability of the one-tail test (See reference 
#' below for details).
#' 
#' @param GD a matrix of individual pairwise genetic distances. In principle 
#' any other squared distance matrix can be used. see example
#' @param GGD A geographical distance matrix, based on the coordinates of
#' individuals. This is typically an Euclidean distance but it can be 
#' any meaningful (geographical) distance metrics.
#' @param bins The number of bins for the distance classes 
#' (i.e. \code{length(bins) == 1)} or a vectors with the break points. See details.
#' @param reps The number to be used for permutation and bootstrap analyses
#' @param permutation Whether permutation calculations for the null hypothesis of 
#' no spatial structure should be carried out
#' @param bootstrap Whether bootstrap calculations to compute the 95% confidence intervals 
#' around r should be carried out
#' @inheritParams gl.diagnostic.hwe
#' @return Returns a data frame with the following columns:
#' \enumerate{
#' \item Bin  The distance classes
#' \item N The number of pairwise comparisons within each distance class
#' \item r.uc The uncorrected autocorrelation coefficient
#' \item Correction the correction 
#' \item r The corrected autocorrelation coefficient
#' \item L.r The corrected autocorrelation coefficient lower limit 
#' (if \code{bootstap = TRUE})
#' \item U.r The corrected autocorrelation coefficient upper limit 
#' (if \code{bootstap = TRUE})
#' \item L.r.null.uc The uncorrectedlower limit for the null hypothesis of no spatial 
#' autocorrelation (if \code{permutation = TRUE}) 
#' \item U.r.null.uc  The uncorrected upper limit for the null hypothesis of no spatial 
#' autocorrelation (if \code{permutation = TRUE})
#' \item L.r.null The corrected lower limit for the null hypothesis of no spatial 
#' autocorrelation (if \code{permutation = TRUE})
#' \item U.r.null The corrected upper limit for the null hypothesis of no spatial 
#' autocorrelation (if \code{permutation = TRUE})
#' \item p.one.tail The p value of the one tail statistical test
#' }
#' 
#' @author Carlo Pacioni & Bernd Gruber
#' @references \itemize{
#' \item Smouse PE, Peakall R. 1999. Spatial autocorrelation analysis of
#' individual multiallele and multilocus genetic structure. Heredity 82:
#' 561-573.
#' 
#' \item Double, MC, et al. 2005. Dispersal, philopatry and infidelity: dissecting
#' local genetic structure in superb fairy-wrens (Malurus cyaneus). Evolution
#' 59, 625-635.
#' 
#' \item Peakall, R, et al. 2003. Spatial autocorrelation analysis offers new
#' insights into gene flow in the Australian bush rat, Rattus fuscipes.
#' Evolution 57, 1182-1195.
#' 
#' \item Smouse, PE, et al. 2008. A heterogeneity test for fine-scale genetic
#' structure. Molecular Ecology 17, 3389-3400.
#' 
#' \item Gonzales, E, et al. 2010. The impact of landscape disturbance on spatial
#' genetic structure in the Guanacaste tree, Enterolobium
#' cyclocarpum(Fabaceae). Journal of Heredity 101, 133-143.
#' 
#' \item Beck, N, et al. 2008. Social constraint and an absence of sex-biased
#' dispersal drive fine-scale genetic structure in white-winged choughs.
#' Molecular Ecology 17, 4346-4358.
#' }
#' @examples
#'  # Select one pop only
#' plat_Tent <- gl.keep.ind(platypus.gl, 
#' ind.list = platypus.gl@ind.names[pop(platypus.gl) == "TENTERFIELD"], 
#' mono.rm = TRUE)
#' # Compute a simple distance matrix and reverse it so that correlated values 
#' # indicated more similar individuals as we are used to see plots in GenAleEx
#' gd <- 1 - as.matrix(gl.dist.ind(plat_Tent, method = "Simple", plot.out = FALSE))
#' # Replace the diagonal with zeros
#' diag(gd) <- 0
#' # Compute the geographical distance matrix
#' ggd <- as.matrix(dist(plat_Tent@other$latlon))
#' 
#' # Spatial autocorrelation
#' spa <- gl.spatial.autoCorr(gd, ggd, bins = 3, reps = 100, 
#' permutation = TRUE, bootstrap = TRUE)
#' 
#' # Alternatively, use Smouse distance from PopGneReport package
#' smouse.plat<-as.matrix(PopGenReport::gd.smouse(gl2gi(plat_Tent), verbose=T))
#' spa.sm <- gl.spatial.autoCorr(smouse.plat, ggd, bins = 3, reps = 100, 
#' permutation = TRUE, bootstrap = TRUE)
#' @import ggplot2
#' @export


gl.spatial.autoCorr <- function(GD, GGD, bins=1, reps=100, 
                                permutation=FALSE, 
                                bootstrap=FALSE, 
                                plot_theme = theme_classic(),
                                save2tmp = FALSE,
                                verbose=NULL) {
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname,
                   build = "Jackson",
                   verbosity = verbose)
  
  # CHECK DATATYPE
  if(!is.matrix(GD))
    if(class(GD) == "dist") GD <- as.matrix(GD)
  else
    stop(error("  The argument 'GD' should be a matrix"))
  
  if(!is.matrix(GGD))
    if(class(GGD) == "dist") GGD <- as.matrix(GGD)
  else
    stop(error("  The argument 'GGD' should be a matrix"))
  
  if(is.character(all.equal(dim(GD), dim(GGD))))
    stop(error("  The arguments GD and GDD should have identical dimensions"))
  
  if(!is.numeric(bins)) 
    stop(error("  The argument 'bins' should be a numeric vector"))
  
  # DO THE JOB
  sample.size <- nrow(GGD)
  crt <- 1/(sample.size - 1) # correction
  nbins <- if(length(bins) == 1) bins else length(bins) - 1
  splist<- utils.spautocor(GD, GGD, permutation=FALSE, bins=bins)
  
  if(permutation) {
    bssplist <- replicate(reps, utils.spautocor(GD, GGD, permutation=TRUE, bins=bins))
    
    #convert the output into a matrix
    bs <-matrix(unlist(bssplist), nrow=reps, ncol=nbins, byrow=TRUE)
    bs.l <- apply(bs,2, quantile, probs=0.025, na.rm=TRUE)
    bs.u <- apply(bs,2, quantile, probs=0.975, na.rm=TRUE)
    
    p.one.tail <- sapply(seq_along(splist$r.uc), function(i, r.rc, r, crt=crt) {
      if(is.na(r[i])) NA 
      else
        if (r[i] >= 0) sum(r.rc[, i] >= r[i])/length(r.rc[, i]) 
        else 
        sum(r.rc[, i] <= r[i])/length(r.rc[, i])
    }, r=splist$r.uc + crt,  r.rc=bs + crt)
    
  }
  if(bootstrap) {
    errors <- replicate(reps, utils.spautocor(GD, GGD, bootstrap=TRUE, bins=bins))
    errors <- matrix(unlist(errors), nrow=reps, ncol=nbins, byrow=TRUE)
    err.l <- apply(errors, 2, quantile, probs=0.025, na.rm=TRUE)
    err.u <- apply(errors, 2, quantile, probs=0.975, na.rm=TRUE)
  }
  
  res <- cbind(splist, Correction=crt, r=splist$r.uc + crt)
  if(bootstrap) res <- cbind(res, L.r=err.l + crt, U.r=err.u + crt)
  if(permutation) res <- cbind(res, L.r.null.uc=bs.l, U.r.null.uc = bs.u,  
                               L.r.null=bs.l + crt,  U.r.null=bs.u + crt, 
                               p.one.tail=p.one.tail)
  
  p <- ggplot(res, aes(Bin, r)) + geom_line() + geom_point() + 
    geom_hline(yintercept=0, col="black") +
    scale_x_continuous(sec.axis=sec_axis(trans = ~., breaks = res$Bin, labels = res$N)) +
    xlab("Distance class") +
    plot_theme
  
  if(bootstrap) p <- p + geom_errorbar(aes(ymin=L.r, ymax=U.r)) 
  if(permutation) p <- p + geom_line(aes(y=L.r.null), col="Red", linetype="dashed") + 
    geom_point(aes(y=L.r.null), col="Red") +
    geom_line(aes(y=U.r.null), col="Red", linetype="dashed") +
    geom_point(aes(y=U.r.null), col="Red")
  
  print(p)
  
  # SAVE INTERMEDIATES TO TEMPDIR
    # creating temp file names
  if (save2tmp) {
    temp_plot <- tempfile(pattern = "Plot_")
    match_call <-
      paste0(names(match.call()),
             "_",
             as.character(match.call()),
             collapse = "_")
    # saving to tempdir
    saveRDS(list(match_call, p), file = temp_plot)
    
    if (verbose >= 2) {
      cat(report("  Saving the ggplot to session tempfile\n"))
      cat(
        report(
          "  NOTE: Retrieve output files from tempdir using gl.list.reports() and gl.print.reports()\n"
        )
      )
    }
    
  }
    
  
  return(res)
}