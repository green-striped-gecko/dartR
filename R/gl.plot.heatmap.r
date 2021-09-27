#' @title gl.plot.heatmap
#' @title Represent a distance matrix as a heatmap 
#' @description 
#' The script plots a heat map to represent the distances in the distance or dissimilarity matrix
#'
#' @param D Name of the distance matrix or class fd object [required]
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 
#' 3, progress and results summary; 5, full report [default 2 or as specified using gl.set.verbosity]
#' @return NULL
#' @export
#' @author Custodian: Luis Mijangos -- Post to \url{https://groups.google.com/d/forum/dartr}) as a wrapper for pheatmap by Raivo Kolde.
#' @examples
#'    gl <- testset.gl[1:10,]
#'    D <- dist(as.matrix(gl),upper=TRUE,diag=TRUE)
#'    gl.plot.heatmap(D)
#'    D2 <- gl.dist.pop(testset.gl)
#'    gl.plot.heatmap(D2)
#'    D3 <- gl.fixed.diff(testset.gl)
#'    gl.plot.heatmap(D3)

gl.plot.heatmap <- function(D,verbose=NULL){

# CHECK IF PACKAGES ARE INSTALLED
  pkg <- "RColorBrewer"
  if (!(requireNamespace(pkg, quietly = TRUE))) {
    stop("Package ",pkg," needed for this function to work. Please install it.") 
    } 
  
  pkg <- "pheatmap"
  if (!(requireNamespace(pkg, quietly = TRUE))) {
    stop("Package ",pkg," needed for this function to work. Please install it.")
    } 
  
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func=funname,build="Jody",v=verbose)
  
  # CHECK DATATYPE 
  datatype <- utils.check.datatype(D,accept=c("dist","fd"),verbose=verbose)

# DO THE JOB
  
  if (class(D)=="dist"){
    pheatmap::pheatmap(as.matrix(D),color=RColorBrewer::brewer.pal(n = 8, name = 'Blues'))
  }
  if (class(D)=="fd"){
    pheatmap::pheatmap(as.matrix(D$fd),color=RColorBrewer::brewer.pal(n = 8, name = 'Blues'))
  }
  
# FLAG SCRIPT END

  if (verbose > 0) {
    cat(report("Completed:",funname,"\n"))
  }
  
  return(NULL)

}

