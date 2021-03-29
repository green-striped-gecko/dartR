#' Isolation by distance
#' 
#' This functions performs an isolation by distance analysis based on a mantel test and also produces an isolation by distance plot. If a genlight object with coordinates is provided) then a Euclidean and genetic distance matrix are calculated (currently. Currently only pairwise Fst between population is implemented. Coordinates are expected as lat long and converted to Google Earth Mercator projection. If coordinates are already projected, set projected=TRUE. If such an object is provided an isolation by distance analysis and plot is performed on log(Euclidean distance) against population based pairwise Fst/1-Fst (see  Rousseau's distance measure. Genetics April 1, 1997 vol. 145 no. 4 1219-1228)
#' You can provide also your own genetic and Euclidean distance matrix. The function is based on the code provided by the adegenet tutorial (\url{http://adegenet.r-forge.r-project.org/files/tutorial-basics.pdf}), using the functions  \link[vegan]{mantel} (package vegan), \link[StAMPP]{stamppFst} (package StAMPP) and Mercator in package dismo.
#' 
#' @importFrom vegan mantel
#' @importFrom MASS kde2d
#' @importFrom grDevices colorRampPalette
#' @importFrom graphics abline title points
#' @importFrom stats as.dist lm
#' @importFrom StAMPP stamppFst
#' @param x genlight object. If provided a standard analysis on Fst/1-Fst and log(distance) is performed
#' @param Dgen genetic distance matrix if no genlight object with coordinates is provided
#' @param Dgeo Euclidean distance matrix if no genlight object is provided
#' @param projected Switch to indicate that coordinates are already projected (not in lat long) and therefore no projection is carried out. Default is FALSE, so it is assumed coordinates are in lat/longs.
#' @param permutations number of permutations in the mantel test
#' @param plot should an isolation by distance plot be returned. Default is plot=TRUE

#' @return returns a list of the following components: Dgen (the genetic distance matrix), Dgeo (the Euclidean distance matrix), mantel (the statistics of the mantel test)
#' @export
#' @author Bernd Gruber (bugs? Post to \url{https://groups.google.com/d/forum/dartr})
#' @seealso \link[vegan]{mantel}, \link[StAMPP]{stamppFst}
#' @references 
#' Rousset (1997) Genetic Differentiation and Estimation of Gene Flow from F-Statistics Under Isolation by Distancenetics 145(4), 1219-1228.
#' @examples 
#' \donttest{
#' ibd <- gl.ibd(bandicoot.gl)
#' ibd <- gl.ibd(bandicoot.gl,plot = FALSE)
#' }


gl.ibd <- 
  function(x, 
           Dgen=NULL, 
           Dgeo=NULL, 
           projected=FALSE, 
           permutations=999, 
           plot=TRUE) {

    # TRAP COMMAND, SET VERSION
    #funname <- match.call()[[1]]
    #build <- "Jacob"
    
    # ERROR CHECKING
    #x <- utils.check.gl(x,verbose)
    #verbose <- x@other$verbose
    
    
  # CHECK IF PACKAGES ARE INSTALLED
  if (!(requireNamespace("dismo", quietly = TRUE))) {
    stop("Package dismo needed for this function to work. Please install it.") } else {
  
  
if (!is.null(Dgen) & !is.null(Dgeo)) cat("Analysis performed on provided genetic and Euclidean distance matrices.")


if (class(x)=="genlight") 
{
  cat("Standard analysis performed on the genlight object. Mantel test and plot will be Fst/1-Fst versus log(distance)\n")
if (nrow(x@other$latlong)!=nInd(x)) stop("Cannot find coordinates for each individual in slot @other$latlong")
#rename long to lon if necessary
if  (sum(match(names(x@other$latlong),"long"), na.rm=T)==1) x@other$latlong$lon <- x@other$latlong$long

  
#project coordinates into Mercator () needs lon/lat order
if (!projected) {
  xy <- dismo::Mercator(x@other$latlong[,c("lon","lat")]) 
  cat("Coordinates transformed to Mercator (google) projection to calculate distances in meters.\n")
  } else 
  {
  xy=x@other$latlong[,c("lon","lat")]
  cat("Coordinates not transformed. Distances calculated on the provided coordinates.")
  }

pop.xy <- apply(xy, 2, function(a) tapply(a, pop(x), mean, na.rm=T) )
Dgeo <- dist(pop.xy)
Dgeo <- log(Dgeo)
Dgen <- as.dist(StAMPP::stamppFst(x, nboots=1))
Dgen <- Dgen/(1-Dgen)

### order both matrices to be alphabetically as levels in genlight

ordering <- levels(pop(x))
Dgen <- as.dist(as.matrix(Dgen)[ordering, ordering])
Dgeo <- as.dist(as.matrix(Dgeo)[ordering, ordering])
}

miss=FALSE
if (sum(is.na(Dgen)>0)) {
  miss=TRUE
  cat("There are missing values in the genetic distance matrix. No kernel distance plot is possible.\n")
  }
if (sum(is.na(Dgeo)>0)) {
  miss=TRUE
  cat("There are missing values in the genetic distance matrix. No kernel distance plot is possible.\n")
}


manteltest <- vegan::mantel(Dgen, Dgeo, na.rm=TRUE, permutations = 999)
print(manteltest)

if (plot) 
  {
  if (!miss) {
  # from adegenet tutorial
  dens <- MASS::kde2d(as.numeric(Dgeo),as.numeric(Dgen), n=300)
  myPal <- colorRampPalette(c("white","blue","gold", "orange", "red"))
  plot(Dgeo, Dgen, pch=20,cex=0.8)
  image(dens, col=transp(myPal(300),.7), add=TRUE)
  points(Dgeo, Dgen, pch=20,cex=0.8)
  abline(lm(as.numeric(Dgen)~as.numeric(Dgeo)))
  title("Isolation by distance")
  } else {
    plot(Dgeo, Dgen)
    abline(lm(as.numeric(Dgen)~as.numeric(Dgeo)))
    title("Isolation by distance")
  }
  
}

out <- list(Dgen=Dgen, Dgeo=Dgeo, mantel=manteltest)
return(out)
    }

}

