#' Recalculate locus metrics after the deletion of individuals or populations
#'
#' The script recalculates statistics made redundant by the deletion of individuals from the dataset and
#' returns the genlight object with the recalculated locus metadata.
#'
#' @param gl -- name of the genlight object containing SNP genotypes [required]
#' @param v -- verbosity: 0, silent; 1, brief; 2, verbose [default 1]
#' @return A genlight object with the revised locus metadata
#' @export
#' @author Arthur Georges (glbugs@@aerg.canberra.edu.au)
#' @examples
#' \dontrun{
#'    gl <- gl.keep.pop(testset.gl, pop.list=c("EmsubRopeMata","EmvicVictJasp"))
#'    gl <- gl.recalc.metrics(gl,v=2)
#' }
#' @seealso \code{\link{gl.filter.monomorphs}}
#' 

gl.recalc.metrics <- function(gl, v=1){

  # Recalculate statistics
    x2 <- dartR:::utils.recalc.avgpic(gl,v=v)
    x2 <- dartR:::utils.recalc.callrate(x2,v=v)
    x2 <- dartR:::utils.recalc.freqhets(x2,v=v)
    x2 <- dartR:::utils.recalc.freqhomref(x2,v=v)
    x2 <- dartR:::utils.recalc.freqhomsnp(x2,v=v)
    
  return(x2)  
}