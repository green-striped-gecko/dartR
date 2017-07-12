#' Isolation by distance
#' 
#' This functions performs an isolation by distance analysis based on a mantel test and also produces an isolation by distance plot. If a genlight object with coordinates is provided) then a Euclidean and genetic distance matrix are calculated (currently. Currently only pairwise Fst between population is implemented. Coordinates are expected as lat long and converted to Google Earth Mercator projection. If coordinates are already projected, set projected=TRUE. If such an object is provided an isolation by distance analysis and plot is performed on log(Euclidean distance) against population based pairwise Fst/1-Fst (see  Rousseau's distance measure. Genetics April 1, 1997 vol. 145 no. 4 1219-1228)
#' You can provide also your own genetic and Euclidean distance matrix. The function is based on the code provided by the adegenet tutorial \url{adegenet.r-forge.r-project.org/files/tutorial-basics.pdf}, using functions from the \link[vegan]{mantel}, \link[StAMPP]{stamppFst} and \link[dismo]{Mercator} packages.
#' 
#' @importFrom vegan mantel
#' @importFrom MASS kde2d
#' @importFrom grDevices colorRampPalette
#' @importFrom graphics abline title
#' @importFrom stats as.dist lm
#' @importFrom StAMPP stamppFst
#' @importFrom dismo Mercator
#' @param gl genlight object. If provided a standard analysis on Fst/1-Fst and log(distance) is performed
#' @param Dgen genetic distance matrix if no genlight object with coordinates is provided
#' @param Dgeo Euclidean distance matrix if no genlight object is provided
#' @param projected Switch to indicate that coordinates are already projected (not in lat long) and therefore no projection is carried out. Default is FALSE, so it is assumed coordinates are in lat/longs.
#' @param permutations number of permutations in the mantel test
#' @param plot should an isolation by distance plot be returned. Default is plot=TRUE

#' @return returns a list of the following components: Dgen (the genetic distance matrix), Dgeo (the Euclidean distance matrix), mantel (the statistics of the mantel test)
#' @export
#' @author Bernd Gruber (glbugs@@aerg.canberra.edu.au)
#' @seealso \link[vegan]{mantel}, \link[StAMPP]{stamppFst}, \link[dismo]{Mercator}
#' @references 
#' Rousset (1997) Genetic Differentiation and Estimation of Gene Flow from F-Statistics Under Isolation by Distancenetics 145(4), 1219-1228.
#' @examples 
#' \dontrun{
#' gl <- gl.ibd(gl=testset.gl)
#' }

gl.ibd <- function(gl=NULL, Dgen=NULL, Dgeo=NULL, projected=FALSE, permutations=999, plot=TRUE) {

if (!is.null(Dgen) & !is.null(Dgeo)) cat("Analysis performed on provided genetic and Euclidean distance matrices.")


if (class(gl)=="genlight") 
{
  cat("Standard analysis performed on the genlight object. Mantel test and plot will be Fst/1-Fst versus log(distance)\n")
if (nrow(gl@other$latlong)!=nInd(gl)) stop("Cannot find coordinates for each individual in slot @other$latlong")


#project coordinates into Mercator () needs lon/lat order
if (!projected) {
  xy <- Mercator(gl@other$latlong[,c("lon","lat")]) 
  cat("Coordinates transformed to Mercator (google) projection to calculate distances in meters.\n")
  } else 
  {
  xy=gl@other$latlong[,c("lon","lat")]
  cat("Coordinates not transformed. Distances calculated on the provided coordinates.")
  }

pop.xy <- apply(xy, 2, function(a) tapply(a, pop(gl), mean, na.rm=T) )
Dgeo <- dist(pop.xy)
Dgeo <- log(Dgeo)
Dgen <- as.dist(stamppFst(gl, nboots=1))
Dgen <- Dgen/(1-Dgen)

### order both matrices to be alphabetically as levels in genlight

ordering <- levels(pop(gl))
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


manteltest <- mantel(Dgen, Dgeo, na.rm=TRUE, permutations = 999)
print(manteltest)

if (plot) 
  {
  if (!miss) {
  # from adegenet tutorial
  dens <- kde2d(Dgeo,Dgen, n=300)
  myPal <- colorRampPalette(c("white","blue","gold", "orange", "red"))
  plot(Dgeo, Dgen, pch=16,cex=0.8)
  image(dens, col=transp(myPal(300),.7), add=TRUE)
  points(Dgeo, Dgen, pch=16,cex=0.8)
   abline(lm(Dgen~Dgeo))
  title("Isolation by distance")
  } else {
    plot(Dgeo, Dgen)
    abline(lm(Dgen~Dgeo))
    title("Isolation by distance")
  }
  
}

out <- list(Dgen=Dgen, Dgeo=Dgeo, mantel=manteltest)
return(out)
}


