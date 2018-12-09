#' Report monomorphic loci, including those with all NAs
#'
#' This script reports the number of monomorphic loci from a genlight \{adegenet\} object
#'
#' A DArT dataset will not have monomorphic loci, but they can arise when populations or individuals are deleted.
#' Retaining monomorphic loci unnecessarily increases the size of the dataset.
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
  
  if (v < 0 | v > 5){
    cat("    Warning: verbosity must be an integer between 0 [silent] and 5 [full report], set to 2\n")
    v <- 2
  }
  
# FLAG SCRIPT START
  if (v >= 1) {
    cat("Starting gl.report.monomorphs: Reporting frequency of monomorphic loci\n")
  }
 
  cat("Identifying monomorphic loci\n")
# Create vectors to hold test results
  # homozygote reference
  a <- vector(mode="logical", length=nLoc(x))
  for (i in 1:nLoc(x)) {a[i] <- NA}
  # homozygote alternate
  b <- vector(mode="logical", length=nLoc(x))
  for (i in 1:nLoc(x)) {b[i] <- NA}
  # heterozygote 
  c <- vector(mode="logical", length=nLoc(x))
  for (i in 1:nLoc(x)) {c[i] <- NA}
  # NA
  d <- vector(mode="logical", length=nLoc(x))
  for (i in 1:nLoc(x)) {d[i] <- NA}
  
# Set up the progress counter
  if(probar) {
    pb <- txtProgressBar(min=0, max=1, style=3, initial=0, label="Working ....")
    getTxtProgressBar(pb)
  }
# Identify polymorphic, monomorphic and 'all na' loci
  # Set a,b,c,d <- TRUE if monomorphic, or if all NAs
  xmat <-as.matrix(x)
  for (i in (1:nLoc(x))) {
    if (all(is.na(xmat[,i]))) {
      d[i] <- TRUE
      a[i] <- FALSE
      b[i] <- FALSE
      c[i] <- FALSE
    } else {
      a[i] <- all(xmat[,i]==0,na.rm=TRUE)
      b[i] <- all(xmat[,i]==2,na.rm=TRUE)
      c[i] <- all(xmat[,i]==1,na.rm=TRUE)
      d[i] <- FALSE
    }
    ##cat(xmat[,i],a[i],b[i],c[i],d[i],"\n")
    if (probar) {setTxtProgressBar(pb, i/nLoc(x))}
  }
  s1 <- sum(a,na.rm=TRUE) + sum(b,na.rm=TRUE) + sum(c,na.rm=TRUE)
  s2 <- s1 - sum(d,na.rm=TRUE)
  polym <- nLoc(x) - s2
  cat("\nBreakdown of", nLoc(x), "loci\n")
  cat("  Monomorphic loci:", s1,"\n  Polymorphic loci:", polym, "\n  Loci with no scores (all NA):" , sum(d) ,"\n")

# FLAG SCRIPT END
  if (v >= 1) {
    cat("gl.report.monomorphs Completed\n")
  }
  
return(NULL)

}

