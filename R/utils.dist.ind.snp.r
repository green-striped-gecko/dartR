#' @name utils.dist.ind.snp
#' @title Calculates a distance matrix for individuals defined in a dartR
#' genlight object using SNP data (DArTseq)
#' @description
#' This script calculates various distances between individuals based on 
#' SNP genotypes.
#' @details
#' The distance measure can be one of:
#'  \itemize{
#'   \item Euclidean -- Euclidean Distance applied to Cartesian coordinates defined
#'   by the loci, scored as 0, 1 or 2. 
#'  \item Simple -- simple mismatch, 0 where no alleles are shared, 1 where one
#'  allele is shared, 2 where both alleles are shared. 
#'  \item Absolute -- absolute mismatch, 0 where no alleles are shared, 1 where
#'  one or both alleles are shared.
#'  \item Czekanowski (or Manhattan) calculates the city block metric distance
#'  by summing the scores on each axis (locus).
#'  }
#'
#' @param x Name of the genlight containing the genotypes [required].
#' @param method Specify distance measure [default Euclidean].
#' @param scale If TRUE and method='Euclidean', the distance will be scaled to 
#'  fall in the range [0,1] [default FALSE].
#' @param output Specify the format and class of the object to be returned, 
#' dist for a object of class dist, matrix for an object of class matrix [default "dist"].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#'  progress log; 3, progress and results summary; 5, full report [default 2].
#' @return An object of class 'dist' or 'matrix' giving distances between individuals
#' @export
#' @author Author(s): Arthur Georges. Custodian: Arthur Georges -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' D <- utils.dist.ind.snp(testset.gl, method='Manhattan')
#' D <- utils.dist.ind.snp(testset.gl, method='Simple')
#' D <- utils.dist.ind.snp(testset.gl, method='Euclidean',scale=TRUE)
#' 

utils.dist.ind.snp <- function(x,
                              method = "Euclidean",
                              scale=FALSE,
                              output="dist",
                              verbose = NULL) {
    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)
    
    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     build = "Jody",
                     verbosity = verbose)
    
    # CHECK DATATYPE
    datatype <-
        utils.check.datatype(x, accept = "SNP", verbose = verbose)
    
    # SCRIPT SPECIFIC ERROR CHECKING
    
    method <- tolower(method)
    
    # FUNCTION SPECIFIC ERROR CHECKING
    
    if (!(method %in% c(
        "euclidean",
        "simple",
        "absolute",
        "czekanowski",
        "manhattan"
        ))) {
        if (verbose >= 2) {
            cat(warn(
                "  Warning: Method not in the list of options, set to Euclidean\n"
            ))
        }
        method <- "euclidean"
    }
    
    if(scale==TRUE && !(method == "euclidean")){
        cat(warn("  Warning: parameter scale only applies to Euclidean Distance, ignored\n"))
    }
    
    # DO THE JOB
    
    mat <- as.matrix(x)

    dd <- array(NA, c(nInd(x), nInd(x)))
    nI <- nInd(x)
    nL <- nLoc(x)
    
    if (verbose >= 2) {
        if(method=="euclidean"){
            if(scale==TRUE){
                cat(report("  Calculating the scaled distance matrix --", method, "\n"))
            } else {
                cat(report("  Calculating the unscaled distance matrix --", method, "\n"))
            }
        } else {
            cat(report("  Calculating the distance matrix --", method, "\n"))
        }
    }
    for (i in (1:(nI - 1))) {
        for (j in ((i + 1):nI)) {
            row1 <- mat[i,]
            row2 <- mat[j,]
            
            if (method == "euclidean") {
              sq <- (row1-row2)**2
              sq <- sq[!is.na(sq)]
              L <- length(sq)
                if(scale==TRUE){
                    dd[j,i] <- sqrt(sum(sq)/L)
                } else {
                    dd[j,i] <- sqrt(sum(sq))
                }
            } else if (method == "simple") {
              row <- array(1,dim=nL)
              row[((row1 + row2) == 4)] <- 2
              row[((row1 + row2) == 0)] <- 0
              row[is.na(row1 + row2)] <- NA
              row <- row[!is.na(row)]
              L <- length(row)
              dd[j,i] <- 1 - sum(row)/(2*L)
            } else if (method == "absolute") {
              row <- array(1,dim=nL)
              row[((row1 + row2) == 0)] <- 0
              row[is.na(row1 + row2)] <- NA
              row <- row[!is.na(row)]
              L <- length(row)
              dd[j,i] <- 1 - sum(row)/(L)
            } else if (method == "manhattan" || method == "czekanowski") {
              modq <- abs(row1-row2)
              modq <- modq[!is.na(modq)]
              L <- length(modq)
              dd[j,i] <- sum(modq)/(2*L)
            } else {
                # Programming error
                stop(error("Fatal Error: Notify dartR development team\n"))
            }
        }
        dd[i, i] <- 0
        dd[i,j] <- dd[j,i]
    }

    if(output=="dist"){
      dd <- as.dist(dd)
      if(verbose >= 2){cat(report("  Returning a stats::dist object\n"))}
    } else {
        if(verbose >= 2){cat(report("  Returning a square matrix object\n"))}
    }
    
    # FLAG SCRIPT END
    
    if (verbose > 0) {
        cat(report("Completed:", funname, "\n"))
    }
    
    return(dd)
}
