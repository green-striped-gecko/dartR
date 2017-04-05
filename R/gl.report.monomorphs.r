#' Report monomorphic loci, including those with all NAs
#'
#' This script reports the number of monomorphic loci from a genlight \{adegenet\} object
#'
#' A DArT dataset will not have monomorphic loci, but they can arise when populations or individuals are deleted.
#' Retaining monomorphic loci unnecessarily increases the size of the dataset.
#'
#' @param gl -- name of the input genlight object [required]
#' @return A report on loci, polymorphic, monomorphic, all NAs
#' @import adegenet plyr utils
#' @export
#' @author Arthur Georges (glbugs@@aerg.canberra.edu.au)
#' @examples
#' gl <- gl.report.monomorphs(testset.gl)

gl.report.monomorphs <- function (gl) {
x <- gl

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
  pb <- txtProgressBar(min=0, max=1, style=3, initial=0, label="Working ....")
  getTxtProgressBar(pb)
# Identify polymorphic, monomorphic and 'all na' loci
  # Set a,b,c,d <- TRUE if monomorphic, or if all NAs
  xmat <-as.matrix(x)
  for (i in (1:nLoc(x))) {
    a[i] <- all(xmat[,i]==0,na.rm=TRUE)
    b[i] <- all(xmat[,i]==2,na.rm=TRUE)
    c[i] <- all(xmat[,i]==1,na.rm=TRUE)
    if (all(is.na(xmat[,i]))) {d[i] <- TRUE}
    setTxtProgressBar(pb, i/nLoc(x))
  }
  polym <- nLoc(x) - sum(a) - sum(b) - sum(d) - sum(c)
  cat("\nBreakdown of", nLoc(x), "loci\n")
  cat("  Polymorphic loci:", polym, "\n  Monomorphic loci:", sum(a)+sum(b)+sum(c), "\n  Loci with no scores (all NA):" , sum(d) ,"\n")

return <- x

}



