#' Report number of private alleles possessed by an individual of unknown provenance
#'
#' This script calculates the number of private alleles possessed by individuals of unknown
#' provinence when compared to a series of populations.
#'
#' The genlight object should identify individuals to a number of known populations with individuals of unknown
#' provinence assigned to population "unknown".
#'
#' @param gl -- name of the input genlight object [required]
#' @return A table with unknown individuals as rows; populations as columns; number of loci with private alleles in the body
#' @import 
#' @export
#' @author Arthur Georges (glbugs@@aerg.canberra.edu.au)
#' @examples
#' gl <- gl.report.pa(testset.gl)

gl.report.pa <- function (gl) {
x <- gl

x <- testset.gl
pop(x[indNames(x)=="UC_00161",])<-"unknown"

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
  # NA
  index <- vector(mode="logical", length=nLoc(x))
  for (i in 1:nLoc(x)) {index[i] <- NA}
  
# Set up the progress counter
  pb <- txtProgressBar(min=0, max=1, style=3, initial=0, label="Working ....")
  getTxtProgressBar(pb)
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
    setTxtProgressBar(pb, i/nLoc(x))
  }
  s1 <- sum(a,na.rm=TRUE) + sum(b,na.rm=TRUE) + sum(c,na.rm=TRUE)
  s2 <- s1 + sum(d,na.rm=TRUE)
  polym <- nLoc(x) - s2
  cat("\nBreakdown of", nLoc(x), "loci\n")
  cat("  Polymorphic loci:", polym, "retained\n  Monomorphic loci:", s1, "deleted\n  Loci with no scores (all NA):" , sum(d) ,"deleted\n")

# Write the polymorphic loci to a new genlight object
#  cat("Deleting monomorphic loci and loci with no scores\n")
  for (i in (1:nLoc(x))) {
    index[i] <- !(a[i] == TRUE || b[i] == TRUE || c[i] == TRUE || d[i] == TRUE)
  }
  x <- x[,(index)]
  x@other$loc.metrics <- x@other$loc.metrics[(index),]

return <- x

}

