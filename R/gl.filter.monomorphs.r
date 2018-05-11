#' Remove monomorphic loci, including those with all NAs
#'
#' This script deletes monomorphic loci from a genlight \{adegenet\} object
#'
#' A DArT dataset will not have monomorphic loci, but they can arise when populations are deleted by assignment or by using
#' the delete option in gl.pop.recode(). Retaining monomorphic loci unnecessarily increases the size of the dataset.
#' @param x -- name of the input genlight object [required]
#' @param v -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @param pb -- display progress bar [FALSE]
#' @return A genlight object with monomorphic loci removed
#' @import adegenet plyr utils
#' @export
#' @author Arthur Georges (glbugs@@aerg.canberra.edu.au)
#' @examples
#' gl <- gl.filter.monomorphs(testset.gl)

# Last edit:25-Apr-18

gl.filter.monomorphs <- function (x, v=2, pb=FALSE) {

  if (v > 0) {
    cat("Starting gl.filter.monomorphs: Deleting monomorphic loci\n")
  }

# Create a vector to hold test results
  a <- vector(mode="logical", length=nLoc(x))
  for (i in 1:nLoc(x)) {a[i] <- NA}
# Set up the progress counter
  if (v > 1 && pb == TRUE) {
    progress <- txtProgressBar(min=0, max=1, style=3, initial=0, label="Working ....")
    getTxtProgressBar(progress)
  }
# Identify polymorphic, monomorphic and 'all na' loci
  # Set a <- TRUE if monomorphic, or if all NAs
  xmat <-as.matrix(x)
  for (i in (1:nLoc(x))) {
    a[i] <- all(xmat[,i]==0,na.rm=TRUE) || all(xmat[,i]==2,na.rm=TRUE)
    if (all(is.na(xmat[,i]))) {a[i] <- NA}
    if (v > 1 && pb == TRUE) {setTxtProgressBar(progress, i/nLoc(x))}
  }
# Count the number of monomorphic loci (TRUE), polymorphic loci (FALSE) and loci with no scores (all.na)
  counts <- plyr::count(a)
  if (v > 2) {
    if (pb) {cat("\n")}
    cat("  Polymorphic loci:", counts[1,2], "\n")
    cat("  Monomorphic loci:", counts[2,2], "\n")
    if (is.na(counts[3,2])) {counts[3,2] <- 0}
    cat("  Loci with no scores (all NA):" , counts[3,2] ,"\n")
  }

  #Treat all na loci as monomorphic
  a[is.na(a)] <- TRUE
# Write the polymorphic loci to a new genlight object
  if (v > 1) {cat("  Deleting monomorphic loci and loci with all NA scores\n")}

  x <- x[,(a==FALSE)]
  x@other$loc.metrics <- x@other$loc.metrics[(a==FALSE),]
  
  if (v > 0) {
    cat("Completed gl.filter.monomorphs\n\n")
  }
  
return (x)
}
