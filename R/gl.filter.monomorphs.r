#' Remove monomorphic loci, including those with all NAs
#'
#' This script deletes monomorphic loci from a genlight \{adegenet\} object
#'
#' A DArT dataset will not have monomorphic loci, but they can arise when populations are deleted by assignment or by using
#' the delete option in gl.pop.recode(). Retaining monomorphic loci unnecessarily increases the size of the dataset.
#'
#' @param gl -- name of the input genlight object [required]
#' @return A genlight object with monomorphic loci removed
#' @import adegenet plyr utils
#' @export
#' @author Arthur Georges (glbugs@@aerg.canberra.edu.au)
#' @examples
#' gl <- gl.filter.monomorphs(testset.gl)

gl.filter.monomorphs <- function (gl) {
x <- gl

  cat("Identifying monomorphic loci\n")
# Create a vector to hold test results
  a <- vector(mode="logical", length=nLoc(x))
  for (i in 1:nLoc(x)) {a[i] <- NA}
# Set up the progress counter
  pb <- txtProgressBar(min=0, max=1, style=3, initial=0, label="Working ....")
  getTxtProgressBar(pb)
# Identify polymorphic, monomorphic and 'all na' loci
  # Set a <- TRUE if monomorphic, or if all NAs
  xmat <-as.matrix(x)
  for (i in (1:nLoc(x))) {
    a[i] <- all(xmat[,i]==0,na.rm=TRUE) || all(xmat[,i]==2,na.rm=TRUE)
    if (all(is.na(xmat[,i]))) {a[i] <- NA}
    setTxtProgressBar(pb, i/nLoc(x))
  }
# Count the number of monomorphic loci (TRUE), polymorphic loci (FALSE) and loci with no scores (all.na)
  counts <- count(a)
  cat("\nPolymorphic loci:", counts[1,2], "\nMonomorphic loci:", counts[2,2], "\nLoci with no scores (all NA):" , counts[3,2] ,"\n")
#Treat all na loci as monomorphic
  # TRUE if monomorphic or all na
  a[is.na(a)] <- TRUE
# Write the polymorphic loci to a new genlight object
#  cat("Deleting monomorphic loci and loci with no scores\n")
  x <- x[,(a==FALSE)]
  x@other$loc.metrics <- x@other$loc.metrics[(a==FALSE),]

# REPORT A SUMMARY
#  cat("\nSummary of filtered dataset\n")
#  cat(paste("  No. of loci:",nLoc(x),"\n"))
#  cat(paste("  No. of individuals:", nInd(x),"\n"))
#  cat(paste("  No. of populations: ", length(levels(factor(pop(x)))),"\n"))

return <- x

}



