#' Report monomorphic loci, including those with all NAs
#'
#' This script reports the number of monomorphic loci from a SNP or tag presence/absence dataset
#'
#' A DArT dataset will not have monomorphic loci, but they can arise when populations or individuals are deleted.
#' Retaining monomorphic loci may unnecessarily increases the size of the dataset and will affect some calculations.
#'
#' @param x -- name of the input genlight object [required]
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @return NULL
#' @import utils
#' @export
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' gl2 <- gl.report.monomorphs(testset.gl)

gl.report.monomorphs <- function (x, probar=FALSE, verbose=2) {
  
  # TIDY UP FILE SPECS
  
  build ='Jacob'
  funname <- match.call()[[1]]
  # Note does not draw upon or modify the loc.metrics.flags
  
  # FLAG SCRIPT START
  
  if (verbose < 0 | verbose > 5){
    cat("  Warning: Parameter 'verbose' must be an integer between 0 [silent] and 5 [full report], set to 2\n")
    verbose <- 2
  }
  
  cat("Starting",funname,"[ Build =",build,"]\n")
  
  # STANDARD ERROR CHECKING
  
  if(class(x)!="genlight") {
    stop("Fatal Error: genlight object required!\n")
  }
  
  if (all(x@ploidy == 1)){
    cat("  Processing Presence/Absence (SilicoDArT) data\n")
  } else if (all(x@ploidy == 2)){
    cat("  Processing a SNP dataset\n")
  } else {
    stop("Fatal Error: Ploidy must be universally 1 (fragment P/A data) or 2 (SNP data)!\n")
  }
  
# DO THE JOB

  #if (verbose >= 2){cat("  Identifying monomorphic loci\n")}
# Create vectors to hold test results
  # homozygote reference
  a <- vector(mode="logical", length=nLoc(x))
  for (i in 1:nLoc(x)) {a[i] <- FALSE}
  # homozygote alternate
  b <- vector(mode="logical", length=nLoc(x))
  for (i in 1:nLoc(x)) {b[i] <- FALSE}
  # heterozygote 
  c <- vector(mode="logical", length=nLoc(x))
  for (i in 1:nLoc(x)) {c[i] <- FALSE}
  # NA
  #if (verbose >= 2){cat("  Identifying loci scored NA over all individuals\n")}
  d <- vector(mode="logical", length=nLoc(x))
  for (i in 1:nLoc(x)) {d[i] <- FALSE}
  
# Identify polymorphic, monomorphic and 'all na' loci
  # Set a,b,c <- TRUE if monomorphic, d <- TRUE if all NAs
  xmat <-as.matrix(x)
  for (i in (1:nLoc(x))) {
    if (all(is.na(xmat[,i]))) {
      d[i] <- TRUE
      a[i] <- FALSE
      b[i] <- FALSE
      c[i] <- FALSE
    } else {
      d[i] <- FALSE
      a[i] <- all(xmat[,i]==0,na.rm=TRUE)
      b[i] <- all(xmat[,i]==2,na.rm=TRUE)
      c[i] <- all(xmat[,i]==1,na.rm=TRUE)
    }
    ##cat(xmat[,i],a[i],b[i],c[i],d[i],"\n")
  }
  # Calculate the number of monomorphic loci with values
  monom <- sum(a,na.rm=TRUE) + sum(b,na.rm=TRUE)
  allna <- sum(d,na.rm=TRUE)
  polym <- nLoc(x) - (monom + allna)
    cat("  Breakdown of", nLoc(x), "loci\n")
    if(all(x@ploidy)==2){
      cat("    Monomorphic loci:", monom,"\n    Polymorphic loci:", polym, "\n    Loci with no scores (all NA):" , allna ,"\n")
    } else {
      cat("    Loci scored all present or all absent:", monom,"\n    Loci with some individals scored present, some absent:", polym, "\n    Loci with no scores (all NA):" , allna ,"\n")
    }

# FLAG SCRIPT END

    cat("Completed:",funname,"\n")

return(NULL)

}

