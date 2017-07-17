#' Remove monomorphic loci, including those with all NAs
#'
#' This script deletes monomorphic loci from a genlight \{adegenet\} object
#'
#' A DArT dataset will not have monomorphic loci, but they can arise when populations are deleted by assignment or by using
#' the delete option in gl.pop.recode(). Retaining monomorphic loci unnecessarily increases the size of the dataset.
#'

#' @param gl -- name of the input genlight object [required]
#' @param v -- verbosty: 0, silent; 1, brief, 2; verbose if TRUE, silent if FALSE [default 1]
#' @return A genlight object with monomorphic loci removed
#' @import adegenet plyr utils
#' @export
#' @author Arthur Georges (glbugs@@aerg.canberra.edu.au)
#' @examples
#' gl2 <- gl.filter.monomorphs(testset.gl)

gl.filter.monomorphs <- function (gl, v=1) {
  
   if (v==1) {cat("Identifying monomorphic loci\n")}
# Create a vector to hold test results
  a <- vector(mode="logical", length=nLoc(gl))
  for (i in 1:nLoc(gl)) {a[i] <- NA}
# Set up the progress counter

  if (v==1) {pb <- txtProgressBar(min=0, max=1, style=3, initial=0, label="Working ....")}
  if (v==1) {getTxtProgressBar(pb)}
# Identify polymorphic, monomorphic and 'all na' loci
  # Set a <- TRUE if monomorphic, or if all NAs
  xmat <-as.matrix(gl)
  for (i in (1:nLoc(gl))) {
    a[i] <- all(xmat[,i]==0,na.rm=TRUE) || all(xmat[,i]==2,na.rm=TRUE)
    if (all(is.na(xmat[,i]))) {a[i] <- NA}
    if (v==1) {setTxtProgressBar(pb, i/nLoc(gl))}

  }
# Count the number of monomorphic loci (TRUE), polymorphic loci (FALSE) and loci with no scores (all.na)
  counts <- count(a)
  if (v==1) {
    cat("\nPolymorphic loci:", counts[1,2], "\nMonomorphic loci:", counts[2,2], "\nLoci with no scores (all NA):" , counts[3,2] ,"\n")
  }
    #Treat all na loci as monomorphic
  # TRUE if monomorphic or all na
  a[is.na(a)] <- TRUE
# Write the polymorphic loci to a new genlight object
#  cat("Deleting monomorphic loci and loci with no scores\n")
  gl <- gl[,(a==FALSE)]
  gl@other$loc.metrics <- gl@other$loc.metrics[(a==FALSE),]

return (gl)

}



