#' Represent a distance matrix as a heatmap 
#'
#' The script plots a heat map to represent the distances in the distance or dissimilarity matrix
#'
#' @param dst -- name of the distance matrix [required]
#' @param ncolors -- number of colors to display [default 5]
#' @param labels -- if TRUE, and the number of rows is <= 20, labels are added to the heatmap [default = TRUE]
#' @param values -- if TRUE, and the number of rows is <= 20, distances are added to the body of the heatmap [default = TRUE]
#' @param rank -- if TRUE, then the distance matrix will be reordered to group like with like, otherwise order will be displayed as given [default FALSE]
#' @param v -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @import graphics
#' @importFrom stats dist
#' @importFrom grDevices heat.colors
#' @return NULL
#' @export
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#'    gl <- testset.gl[1:10,]
#'    d <- dist(as.matrix((gl)))
#'    gl.dist.heatmap(d)
#'    gl.dist.heatmap(d, ncolors=10, rank=TRUE)

gl.dist.heatmap <- function(dst, ncolors=5, labels=TRUE, values=TRUE, rank=FALSE, v=2){
  
# ERROR CHECKING
  
  if(class(dst)!="dist") {
    cat("Fatal Error: distance matrix required!\n"); stop("Execution terminated\n")
  }

  if (v < 0 | v > 5){
    cat("    Warning: verbosity must be an integer between 0 [silent] and 5 [full report], set to 2\n")
    v <- 2
  }
  
  if (ncolors < 0){
    cat("    Warning: ncolors must be a positive integer, set to 5\n")
    ncolors <- 5
  }
  
  if (max(dst) == 0) {
    cat("    Warning: matrix contains no nonzero distances\n")
  } 
  
# FLAG SCRIPT START
  
  if (v >= 1) {
    cat("Starting gl.dist.heatmap: Displaying distance matrix\n")
  }
  
# DO THE JOB
  
  # Convert the distance matrix to a numeric matrix
  x <- as.matrix(dst)
  #x <- x[1:10,1:10]
  dim <- ncol(x)
  
  # If the matrix is to be ordered on rank
  if(rank){
    x <- x[order(rowMeans(x),decreasing=TRUE),order(colMeans(x),decreasing=TRUE)]
  }
  
  # Check if labels and values can be plotted
  if (dim > 20){
    if (labels){
      cat("    Warning: too many labels to display (more than 20)\n")
      labels=FALSE
    }
    if (values){
      cat("    Warning: too many cells to display values within (more than 20x20)\n")
      values=FALSE
    }
  }
  
  # Hold the raw values
  vals <- x
  
  # Scale the values to fall between 0 and 1
  if (max(x) > 0){
    x <- 1-x/max(x)
  }  
  
  # Invert the matrix so the diagonal runs top left to bottom right
  x <- apply(x, 2, rev) 
  vals <- apply(vals,2,rev) 

  # Plot the heat map  
  par(pty="s")
  image(1:dim, 1:dim, x, axes = FALSE, xlab="", ylab="", col=heat.colors(ncolors))

  # Add the labels
  if (labels) {
    axis(1, 1:dim, row.names(x), cex.axis=(10/dim), las=3)
    axis(2, 1:dim, colnames(x), cex.axis = (10/dim), las=1)
  } 

  # Add the values
  if (values) {
    text(expand.grid(1:dim, 1:dim), sprintf("%0.1f", vals), cex=10/dim)
  }  

  # FLAG SCRIPT END
  
  if (v >= 1) {
    cat("Completed gl.dist.heatmap\n\n")
  }
  
  return()

}
