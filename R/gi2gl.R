#' Converts a genind object to genlight object
#' 
#' @param gi -- a genind object
#' @param parallel -- switch to deactivate parallel version. Default set to TRUE. Only for testing purpose. 
#' @return A genlight object, with all slots filled.
#' @export
#' @author Bernd Gruber (bugs? Post to \url{https://groups.google.com/d/forum/dartr})
#' @details Be aware due to ambiguity which one is the reference allele a combination of gi2gl(gl2gi(gl)) does not return an identical object (but in terms of analysis this conversions are equivalent)


gi2gl <- function(gi, parallel=TRUE)
{
  locna <- gi@loc.n.all
  ccc<-1
  for (i in 2:length(locna)) 
  {
    if (locna[i-1]==1)  ccc[i] <- ccc[i-1]+1 else ccc[i]<- ccc[i-1]+2
  }
  gl <-new("genlight", gi@tab[,ccc], pop = pop(gi), other=gi@other, ploidy=2, loc.names=locNames(gi), ind.names=indNames(gi), parallel=parallel)
  return(gl)
}

