#' @title utils.dist.binary
#' @title Calculate a distance matrix for individuals defined in an \{adegenet\} genlight object using binary P/A data (SilicoDArT)
#' @description 
#' This script calculates various distances between individuals based on Tag Presence/Absence data.  
#' @details  
#' The distance measure can be one of
#'  
#'  simple -- simple matching, both 1 or both 0 = 0; one 1 and the other 0 = 1. Presence and absence equally weighted.
#'  Jaccard -- ignores matching 0, both 1 = 0; one 1 and the other 0 = 1. Absences could be for different reasons.
#'  Dice -- both 0 = 0; both 1 = 2; one 1 and the other 0 = 1. Absences could be for different reasons. Sometimes called the Czekanowski or Sorensen distance.
#'  Phi -- binary analogue of the Pearson Correlation coefficient.
#'  
#'  One might choose to disregard or downweight absences in comparison with presences because the homology of absences is less clear (mutation at one or
#'  the other, or both restriction sites). Your call.
#'  
#' @param x -- name of the genlight containing the SNP genotypes [required]
#' @param method -- Specify distance measure [simple]
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [2]
#' @return An object of class 'dist' giving distances between individuals
#' @export
#' @author Custodian: Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' D <- utils.dist.binary(testset.gs, method="Jaccard")

utils.dist.binary <- function(x, method="simple", verbose=NULL) {
  
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func=funname,build="Jackson",v=verbose)
  
  # CHECK DATATYPE 
  datatype <- utils.check.datatype(x,accept="SNP",verbose=verbose)
  
# SCRIPT SPECIFIC ERROR CHECKING
  
  method <- tolower(method)
  
# FUNCTION SPECIFIC ERROR CHECKING
  
  if (!(method %in% c("simple", "jaccard", "dice", "sorenson", "czekanowski", "phi"))){
    if(verbose >= 2){cat(warn(" Warning: Method not in the list of options, set to simple matching\n"))}
    method <- 'simple'
  }
  
# DO THE JOB
  
  mat <- as.matrix(x)
  
  dd <- array(NA,c(nInd(x),nInd(x)))
  #dd[1:10,1:10]
  nI <- nInd(x)

  if(verbose >= 2)cat(report("  Calculating the distance matrix --",method,"\n"))  
  for (i in (1:(nI-1))) {
  for (j in ((i+1):nI)){
    row1 <- mat[i,]
    row2 <- mat[j,]
    #row1[1:10]
    #row2[1:10]
    a11 <- (row1+row2)==2
    a10 <- ((row1+row2)==1)*row1
    a01 <- ((row1+row2)==1)*row2
    a00 <- (row1+row2)==0
    a <- sum(a11==1,na.rm=TRUE)
    b <- sum(a01==1,na.rm=TRUE)
    c <- sum(a10==1,na.rm=TRUE)
    d <- sum(a00==1,na.rm=TRUE)
    #a;b;c;d
    if (method == 'simple'){
      dd[j,i] <- 1 - (a+d)/(a+b+c+d)
    } else if (method == 'jaccard'){  
      dd[j,i] <- 1 - a/(a+b+c)
    } else if (method == 'dice' || method == 'sorenson' || method == 'czekanowski'){
      dd[j,i] <- 1 -2*a/(2*a+b+c)
    } else { 
      # method == phi
      dd[j,i] <- 1 -((a*d - b*c)/sqrt((a + b)*(a + c)*(d + b)*(d + c)))
    } 
  }
    dd[i,i] <- 0
  }
#dd[1:10,1:10]
  if(verbose >= 2){ cat(report("  Converting to a distance object\n"))}
  dd <- as.dist(dd)

# FLAG SCRIPT END

  if (verbose > 0) {
    cat("Completed:",funname,"\n")
  }

  return(dd)
}
