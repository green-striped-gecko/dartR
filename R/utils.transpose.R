# this is a function to transpose a genlight object
utils.transpose <- function(x,
                            parallel = FALSE) {
  hold <- x
  x@gen <- matrix2gen(t(as.matrix(x)), parallel = parallel)
  x@n.loc <- nInd(hold)
  indNames(x) <- locNames(hold)
  locNames(x) <- indNames(hold)
  # This is just a dummy vector to comply with the attributes of a genlight object
  alleles(x) <-
    paste(rep("A", nInd(hold)), rep("A", nInd(hold)), sep = "/")
  ploidy(x) <- unique(ploidy(hold))
  pop(x) <- rep("NA", nLoc(hold))
  x@other$loc.metrics <- hold@other$ind.metrics
  x@other$ind.metrics <- hold@other$loc.metrics
  return(x)
}