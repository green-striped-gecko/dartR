#' Represent a distance matrix as a heatmap 
#'
#' The script plots a heat map to represent the distances in the distance or dissimilarity matrix
#'
#' @param D -- name of the distance matrix or class fd object [required]
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @return NULL
#' @export
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr}) as a wrapper for pheatmap by Raivo Kolde.
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
    stop("Package",pkg," needed for this function to work. Please install it.") } 
  
  pkg <- "pheatmap"
  if (!(requireNamespace(pkg, quietly = TRUE))) {
    stop("Package",pkg," needed for this function to work. Please install it.") } 
  
  
# TRAP COMMAND, SET VERSION
  
  funname <- match.call()[[1]]
  build <- "Jacob"
  
# SET VERBOSITY
  
  if (is.null(verbose)){ 
          verbose <- 2
      } 
  
  if (verbose < 0 | verbose > 5){
    cat(paste("  Warning: Parameter 'verbose' must be an integer between 0 [silent] and 5 [full report], set to 2\n"))
    verbose <- 2
  }
  
# FLAG SCRIPT START
  
  if (verbose >= 1){
    if(verbose==5){
      cat("Starting",funname,"[ Build =",build,"]\n")
    } else {
      cat("Starting",funname,"\n")
    }
  }

# FUNCTION SPECIFIC ERROR CHECKING

  if(class(D)!="dist" & class(D)!="matrix" & class(D)!="fd") {
    cat("  Fatal Error: distance matrix of class 'dist' or class 'fd' required!\n"); stop("Execution terminated\n")
  }
  if (class(D)=="dist") {
    if (max(D)==0){cat("    Warning: matrix contains no nonzero distances\n")}
  } 
  if (class(D)=="fd") {
    if (max(D$fd)==0){cat("    Warning: matrix contains no nonzero distances\n")}
  } 

# DO THE JOB
  
  if (class(D)=="dist"){
    pheatmap::pheatmap(as.matrix(D),color=RColorBrewer::brewer.pal(n = 8, name = 'Blues'))
  }
  if (class(D)=="fd"){
    pheatmap::pheatmap(as.matrix(D$fd),color=RColorBrewer::brewer.pal(n = 8, name = 'Blues'))
  }
  
# FLAG SCRIPT END

  if (verbose > 0) {
    cat("Completed:",funname,"\n")
  }
  
  return(NULL)

}

