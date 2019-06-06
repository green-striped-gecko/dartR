#' Report monomorphic loci, including those with all NAs
#'
#' This script reports the number of monomorphic loci from a genlight \{adegenet\} object
#'
#' A DArT dataset will not have monomorphic loci, but they can arise when populations or individuals are deleted.
#' Retaining monomorphic loci may unnecessarily increases the size of the dataset.
#'
#' @param x -- name of the input genlight object [required]
#' @param probar -- if TRUE, a progress bar will be displayed for long loops [default = TRUE]
#' @param verbose level of verbosity. verbose=0 is silent, verbose=1 returns more detailed output during conversion.
#' @return NULL
#' @import utils
#' @export
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' gl2 <- gl.report.monomorphs(testset.gl)

# Last amended 3-Feb-19

gl.report.monomorphs <- function (x, probar=FALSE, verbose=0) {
  
# TIDY UP FILE SPECS

  funname <- match.call()[[1]]

# FLAG SCRIPT START

  if (verbose > 0) {
    cat("Starting",funname,"\n")
  }

# STANDARD ERROR CHECKING
  
  if(class(x)!="genlight") {
    cat("  Fatal Error: genlight object required!\n"); stop("Execution terminated\n")
  }

  # Set a population if none is specified (such as if the genlight object has been generated manually)
    if (is.null(pop(x)) | is.na(length(pop(x))) | length(pop(x)) <= 0) {
      if (verbose >= 1){ cat("  Population assignments not detected, individuals assigned to a single population labelled 'pop1'\n")}
      pop(x) <- array("pop1",dim = nInd(x))
      pop(x) <- as.factor(pop(x))
    }

# DO THE JOB

  cat("Identifying monomorphic loci\n")
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
  d <- vector(mode="logical", length=nLoc(x))
  for (i in 1:nLoc(x)) {d[i] <- FALSE}
  
# Set up the progress counter
  if(probar) {
    pb <- txtProgressBar(min=0, max=1, style=3, initial=0, label="Working ....")
    getTxtProgressBar(pb)
  }
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
    if (probar) {setTxtProgressBar(pb, i/nLoc(x))}
  }
  # Calculate the number of monomorphic loci with values
  monom <- sum(a,na.rm=TRUE) + sum(b,na.rm=TRUE)
  allna <- sum(d,na.rm=TRUE)
  polym <- nLoc(x) - (monom + allna)
  cat("\nBreakdown of", nLoc(x), "loci\n")
  cat("  Monomorphic loci:", monom,"\n  Polymorphic loci:", polym, "\n  Loci with no scores (all NA):" , allna ,"\n")

# FLAG SCRIPT END

    cat("Completed:",funname,"\n")

return(NULL)

}

