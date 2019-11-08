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
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default NULL]
#' @param ... --- additional arguments passed to glPlot function.
#' @export
#'@examples 
#'gl.plot(bandicoot.gl[1:30,])
#'gl.plot(bandicoot.gl[1:30,])
#'gl.plot(bandicoot.gl[1:10,],labels=TRUE)

gl.plot <- function (x, labels=FALSE, indlabels=indNames(x), col=NULL, legend=TRUE, posi="bottomleft", bg=rgb(1,1,1,.5), verbose=NULL,...)
{
  
# TRAP COMMAND, SET VERSION

  funname <- match.call()[[1]]
  build <- "Jacob"
  
# SET VERBOSITY
  
  if (is.null(verbose)){
    verbose <- x@other$verbose
  }
  if (verbose < 0 | verbose > 5){
      cat(paste("  Warning: Parameter 'verbose' must be an integer between 0 [silent] and 5 [full report], set to default",x@other$verbose,"\n"))
      verbose <- x@other$verbose
  }

# FLAG SCRIPT START

  if (verbose >= 1){
    cat("Starting",funname,"[ Build =",build,"]\n")
  }

# STANDARD ERROR CHECKING
  
  if(class(x)!="genlight") {
    cat("  Fatal Error: genlight object required!\n"); stop("Execution terminated\n")
  }

  # Set a population if none is specified (such as if the genlight object has been generated manually)
    if (is.null(pop(x)) | is.na(length(pop(x))) | length(pop(x)) <= 0) {
      if (verbose >= 2){ cat("  Population assignments not detected, individuals assigned to a single population labelled 'pop1'\n")}
      pop(x) <- array("pop1",dim = nInd(x))
      pop(x) <- as.factor(pop(x))
    }

  # Check for monomorphic loci
    tmp <- gl.filter.monomorphs(x, verbose=0)
    if ((nLoc(tmp) < nLoc(x)) & verbose >= 2) {cat("  Warning: genlight object contains monomorphic loci\n")}

# DO THE JOB
  
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

  if (verbose > 0) {
    cat("Completed:",funname,"\n")
  }
  
  return(invisible())
}


