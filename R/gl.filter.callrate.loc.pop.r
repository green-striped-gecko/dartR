#' Filter loci in a genlight adegenet object based on call rate and population level. 
#'
#' This funciton is a convenience function that filters callrate for the set threshold for each population and then returns a genlight object where each loci of every population passes the threshold of the call rate filter. For more details on call rates see gl.filter.callrate.
#' @param  x - genlight object to be filtered
#' @param threshold threshold on call rate for loci for each population (see gl.filter.callrate for details.)
#' @return The filtered genlight or genind object
#' @export
#' @examples
#' gg <- gl.filter.callrate.loc.pop(testset.gl, 0.9)

gl.filter.callrate.loc.pop <- function(x, threshold){
  
  pops <- seppop(x) 
  ll <- lapply(pops, function(x) locNames(gl.filter.callrate(x, method = "loc", threshold = threshold, verbose = 0)))
  locall <- Reduce(intersect, ll)
  index <- which(locNames(x) %in% locall)
  x <- x[ , locall]
  x@other$loc.metrics <- x@other$loc.metrics[locall,]
  #add to history
  nh <- length(x@other$history)
  x@other$history[[nh + 1]] <- match.call()
  return(x)
}