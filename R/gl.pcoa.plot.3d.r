#' 3D interactive plot of the results of a PCoA ordination
#'
#' This script takes output from the ordination undertaken using gl.pcoa() and plots the individuals in 3D space. The visualisation
#' can be rotated, zoomed in and zoomed out with the mouse to examine the structure.
#'
#' The factor scores are taken from the output of gl.pcoa(), an object of class glPca, and the population assignments from the original data file
#' and plots the specimens in a 3D plot.
#'
#' Axes can be specified from the ordination, provided they are within the range of the nfactors value provided to gl.pcoa().
#' 
#' @param glPca -- name of the glPca object containing the factor scores and eigenvalues [required]
#' @param x -- name of the genlight object or fd object from which the PCoA was generated [required]
#' @param xaxis -- identify the x axis from those available in the ordination (xaxis <= nfactors) [default 1]
#' @param yaxis -- identify the y axis from those available in the ordination (yaxis <= nfactors) [default 2]
#' @param zaxis -- identify the z axis from those available in the ordination (zaxis <= nfactors) [default 3]
#' @param radius -- size of the points [default 8]
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2 or as specified using gl.set.verbosity]
#' @return NULL, plots an interactive 3D plot of the ordination in a separate window
#' @export
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' pcoa <- gl.pcoa(possums.gl, nfactor=5)
#' gl.pcoa.plot.3d(pcoa, possums.gl, xaxis=1, yaxis=2, zaxis=3)

gl.pcoa.plot.3d <- function(glPca, 
                            x, 
                            xaxis=1,
                            yaxis=2,
                            zaxis=3,
                            radius=8,
                            verbose=NULL) {
# CHECK IF PACKAGES ARE INSTALLED
  pkg <- "dplyr"
  if (!(requireNamespace(pkg, quietly = TRUE))){
    stop("Package ",pkg," needed for this function to work. Please install it.")
    }
  pkg <- "plotly"
  if (!(requireNamespace(pkg, quietly = TRUE))){
    stop("Package ",pkg," needed for this function to work. Please install it.")
  }else{

  
  # TRAP COMMAND, SET VERSION
  
  funname <- match.call()[[1]]
  build <- "Jacob"
  
# SET VERBOSITY
  
  if(class(x)=="genlight"){
    if (is.null(verbose)){ 
      if(!is.null(x@other$verbose)){ 
        verbose <- x@other$verbose
      } else { 
        verbose <- 2
      }
    }
  }
  if(class(x)=="fd"){
    x <- x$gl
    if (is.null(verbose)){
      verbose <- 2
    }  
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
  
# STANDARD ERROR CHECKING
  
  if(class(x)!="genlight") {
    cat("  Fatal Error: genlight object required!\n"); stop("Execution terminated\n")
  }

  # Set a population if none is specified (such as if the genlight object has been generated manually)
    if (is.null(pop(x)) | is.na(length(pop(x))) | length(pop(x)) <= 0) {
      if (verbose >= 2){ cat("  Population assignments not detected, individuals assigned to a single population labelled 'pop1'\n")}
      pop(x) <- array("pop1",dim = nLoc(x))
      pop(x) <- as.factor(pop(x))
    }

# FUNCTION SPECIFIC ERROR CHECKING

# DO THE JOB
  
  # function to replicate defaults colors of ggplot
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  
  # assigning colors to populations
  colors_pop <- gg_color_hue(length(levels(pop(x))))
  names(colors_pop) <- as.character(levels(x$pop))
  
  # extracting PCs
  PCOA_scores <- as.data.frame(glPca$scores)
  PCOA_scores$pop <- x$pop
  
  # PCA 3D
  print(
  plot_ly(PCOA_scores,x=~PCOA_scores[,xaxis],y=~PCOA_scores[,yaxis],z=~PCOA_scores[,zaxis],
          marker = list(size = radius),colors = colors_pop)%>% 
    add_markers(color=~pop)%>%
    layout(legend=list(title=list(text='<b> Populations <b>')),
           scene = list(xaxis = list(title = paste0('PC ',xaxis),titlefont = list(size = 25)),
                        yaxis = list(title = paste0('PC ',yaxis),titlefont = list(size = 25)),
                        zaxis = list(title = paste0('PC ',zaxis),titlefont = list(size = 25))))
  )

# FLAG SCRIPT END

  if (verbose > 0) {
    cat("Completed:",funname,"\n")
  }

  return(NULL)
  }
}
