#' Plotting genlight object as a smear plot (loci by individuals color coded for scores of 0, 1, 2 and NA)
#' 
#' This function is based on the glPlot function from adegenet. It adds the option to put labels on the individuals and scales them accordingly.
#' If there are too many individuals, it is best to use labels=FALSE.
#'  
#' For arguments please refer to the original adegenet function ?glPlot. 
#' 
#' @param x -- a genlight object [required]
#' @param labels -- if TRUE, individual labels are added
#' @param indlabels -- labels for individuals [default = first 8 letters from indNames]
#' @param col -- optional color vector (see ?glPlot) [default NULL]
#' @param legend -- if TRUE, a legend will be added [default = TRUE]
#' @param posi -- position of the legend [default = "bottomleft"]
#' @param bg -- background color of the legend [default transparent white]
#' @param v -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @param ... --- additional arguments passed to glPlot function.
#' @export
#'@examples 
#'gl.plot(foxes.gl[1:30,])
#'gl.plot(foxes.gl[1:10,],labels=TRUE)

gl.plot <- function (x, labels=FALSE, indlabels=indNames(x), col=NULL, legend=TRUE, posi="bottomleft", bg=rgb(1,1,1,.5), v=2,...) 
{
  
# ERROR CHECKING
  
  if(class(x)!="genlight") {
    cat("Fatal Error: genlight object required for gl.report.callrate!\n"); stop()
  }
  if (v < 0 | v > 5){
    cat("    Warning: verbosity must be an integer between 0 [silent] and 5 [full report], set to 2\n")
    v <- 2
  }
  
  # FLAG SCRIPT START
  
  if (v >= 1) {
    cat("Starting gl.plot: Smear plot\n")
  }
  
  if (!labels){
    glPlot(x)
  } else {
    X <- t(as.matrix(x))
    X <- X[, ncol(X):1]
  
    if (is.null(indlabels)) indlabes <- pretty(1:nInd(x),5)
    if (is.null(col)) {
      myCol <- colorRampPalette(c("royalblue3", "firebrick1"))(max(X, 
                                                                 na.rm = TRUE) + 1)
    }
    else {
      myCol <- col
    }
    image(x = 1:nLoc(x), y = 1:nInd(x), z = X, xlab = "SNP index", yaxt = "n", col = myCol, ylab="",...)
    axis(side = 2, at = nInd(x):1 , labels = indlabels, las=2)
    if (legend) {
      legend(posi, fill = myCol, legend = 0:max(X, na.rm = TRUE), 
           horiz = TRUE, bg = bg, title = "Number of 2nd allele")
    }
  }  
  
# FLAG SCRIPT END
  if (v >= 1) {
    cat("\ngl.plot Completed\n")
  }
  
  return(invisible())
}


