#' @title gl.plot.heatmap
#' @title Represent a distance matrix as a heatmap 
#' @description 
#' The script plots a heat map to represent the distances in the distance or 
#' dissimilarity matrix. This function is a wrapper for 
#' \link[gplots]{heatmap.2} (package gplots).
#'
#' @param D Name of the distance matrix or class fd object [required].
#' @param palette_divergent A divergent palette for the distance values
#'  [default diverging_palette].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 
#' 3, progress and results summary; 5, full report [default 2 or as specified using gl.set.verbosity]
#' @param ... Parameters passed to function \link[gplots]{heatmap.2} (package gplots)
#' @return NULL
#' @export
#' @author Custodian: Luis Mijangos -- Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#'    gl <- testset.gl[1:10,]
#'    D <- dist(as.matrix(gl),upper=TRUE,diag=TRUE)
#'    gl.plot.heatmap(D)
#'    D2 <- gl.dist.pop(testset.gl)
#'    gl.plot.heatmap(D2)
#'    D3 <- gl.fixed.diff(testset.gl)
#'    gl.plot.heatmap(D3)

gl.plot.heatmap <- function(D,
                            palette_divergent = diverging_palette,
                            verbose = NULL,
                            ...){

  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func=funname,build="Jody",v=verbose)
  
  # CHECK DATATYPE 
  datatype <- utils.check.datatype(D,accept=c("dist","fd"),verbose=verbose)

# DO THE JOB
  
  if (datatype=="dist"){
    p3 <- gplots::heatmap.2(as.matrix(D),
                            col = palette_divergent(255),
                            dendrogram = "column",
                            trace = "none",...)
  }
  if (datatype=="fd"){
    p3 <- gplots::heatmap.2(as.matrix(D$fd),
                            col = palette_divergent(255),
                            dendrogram = "column",
                            trace = "none",...)
  }
  
# FLAG SCRIPT END

  if (verbose > 0) {
    cat(report("Completed:",funname,"\n"))
  }
  
  invisible(p3)

}

