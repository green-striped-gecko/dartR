#' Plotting genlight object
#' 
#' This function is based on the glPlot function from adegenet. It simply aims to put labels on the individuals and scales them accordingly (if there are not too many). 
#' For arguments please refer to the original adegenet function ?glPlot. 
#' 
#' @param x -- a genlight object
#' @param indlabels -- labels for individuals. if not provided labels are taken from the first 8 letters from indNames.
#' @param col -- optional color vector (see ?glPlot)
#' @param legend -- a logical indicating whether a legend should be added
#' @param posi -- position of the legend
#' @param bg -- background color of the legend [default is transparent white]

#' @param ... --- additional arguments passed to glPlot function.
#' @export
#'@examples 
#'gl.plot(foxes.gl[1:30,])
gl.plot <- function (x, indlabels=indNames(x), col=NULL, legend=TRUE, posi="bottomleft", bg=rgb(1,1,1,.5),...) 
{
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
  return(invisible())
}


