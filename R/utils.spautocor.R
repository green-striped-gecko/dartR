#' @name utils.spautocorr
#' @title Spatial autocorrelation coefficient calculations
#'
#' @description  Carries out calculation for spatial autocorrelation coefficient
#' starting from a genetic and geogaphic distance matrix.
#'
#' @details The code of this function is based one \code{spautocorr} from the package
#' \code{PopGenReport}, which has been modified to fix a few bugs (as of
#' \code{PopGenReport v 3.0.4} and allow calculations of bootstraps estimates.
#'
#' See details from \code{gl.spatial.autoCorr} for a detailed explanation.
#' @param GD Genetic distance matrix.
#' @param GGD Geographic distance matrix.
#'
#' @inheritParams gl.spatial.autoCorr
#' @return Returns a data frame with the following columns:
#' \enumerate{
#' \item Bin  The distance classes
#' \item N The number of pairwise comparisons within each distance class
#' \item r.uc The uncorrected autocorrelation coefficient
#' }
#' if both \code{bootstap} and \code{permutation} are \code{FALSE} otherwise only
#' \code{r} estimates are returned
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
#' # See gl.spatial.autoCorr
#' @seealso \code{\link{gl.spatial.autoCorr}}
#' @export

utils.spautocor <- function(GD,
                            GGD,
                            permutation = FALSE,
                            bootstrap = FALSE,
                            bins = 10,
                            reps) {
  gd <- GD
  ed <- GGD
  
  if (permutation == TRUE) {
    gdsample <- sample(nrow(gd))
    gd <- gd[gdsample, gdsample]
  }
  
  cdmat <- function(gd) {
    dimen <- nrow(gd)
    sgd <- sum(gd, na.rm = TRUE)
    cscd <- matrix(colSums(gd, na.rm = TRUE), dimen, dimen)
    rscd <-
      matrix(rowSums(gd, na.rm = TRUE), dimen, dimen, byrow = TRUE)
    cd <- 0.5 * (-gd + 1 / dimen * (cscd + rscd) - 1 / dimen ^ 2 * (sgd))
    cd
  }
  cd <- cdmat(gd)
  #remove upper triangel to speed things up....
  ed[upper.tri(ed)] <- NA
  diag(ed) <- NA
  
  if (length(bins) == 1) {
    steps <- seq_len(bins) * signif(diff(range(ed, na.rm = TRUE)) / bins, 4)
    steps <- c(min(ed, na.rm = TRUE), steps)
    if (steps[length(steps)] < max(ed, na.rm = TRUE))
      steps[length(steps)] <- max(ed, na.rm = TRUE)
  } else {
    steps <- bins
  }
  compute.r <- function(d, steps, ed, cd, bootstrap = FALSE) {
    if (length(steps) == d) {
      index <- which(ed <= steps[d] & ed >= steps[d - 1], arr.ind = TRUE)
    } else {
      index <- which(ed < steps[d] & ed >= steps[d - 1], arr.ind = TRUE)
    }
    
    N <- length(index) / 2
    distance <- steps[d]
    
    if (bootstrap == TRUE) {
      if (reps > choose(N + N - 1, N)) {
        # The number of possible combination with replacement are C(n+k-1, k)
        return(data.frame(
          Bin = distance,
          N = N,
          r.uc = NA
        ))
      } else {
        if (nrow(index) > 1) {
          row.sample <- sample(nrow(index), replace = TRUE)
          index <- index[row.sample,]
        }
      }
      
    }
    
    cx <- sum(cd[index])
    cxii <- sum(diag(cd)[index[, 1]])
    cxjj <- sum(diag(cd)[index[, 2]])
    r <-  2 * cx / (cxii + cxjj)
    
    
    
    return(data.frame(Bin = distance, N = N, r.uc = r))
  }
  l.res <-
    lapply(
      2:length(steps),
      compute.r,
      steps = steps,
      ed = ed,
      cd = cd,
      bootstrap = bootstrap
    )
  res <- do.call(rbind, l.res)
  if (permutation == FALSE & bootstrap == FALSE)
    return(res)
  else
    return(res[, "r.uc"])
  
}