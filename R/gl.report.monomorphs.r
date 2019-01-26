#' Report monomorphic loci, including those with all NAs
#'
#' This script reports the number of monomorphic loci from a genlight \{adegenet\} object
#'
#' A DArT dataset will not have monomorphic loci, but they can arise when populations or individuals are deleted.
#' Retaining monomorphic loci may unnecessarily increases the size of the dataset.
#'
#' @param x -- name of the input genlight object [required]
#' @param probar -- if TRUE, a progress bar will be displayed for long loops [default = TRUE]
#' @return NULL
#' @import adegenet plyr utils
#' @export
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' gl2 <- gl.report.monomorphs(testset.gl)

gl.report.monomorphs <- function (x, probar=FALSE) {
  
# ERROR CHECKING
  
  if(class(x)!="genlight") {
    cat("Fatal Error: genlight object required for gl.report.repavg!\n"); stop()
  }
  # Work around a bug in adegenet if genlight object is created by subsetting
  x@other$loc.metrics <- x@other$loc.metrics[1:nLoc(x),]
  
# FLAG SCRIPT START
    cat("Starting gl.report.monomorphs: Reporting frequency of monomorphic loci\n")

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

    cat("gl.report.monomorphs Completed\n")

return(NULL)

}

