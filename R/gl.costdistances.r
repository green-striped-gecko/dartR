#'Calculates cost distances for a given landscape (resistance matrix)
#'
#'@param landscape a raster object coding the resistance of the landscape
#'@param locs coordinates of the subpopulations. If a genlight object is provided coordinates are taken from @other$latlon and centers for population (pop(gl)) are calculated. In case you want to calculate costdistances between individuals redefine pop(gl) via: \code{pop(gl)<- indNames(gl)}.
#'@param method defines the type of cost distance, types are "least-cost", "rSPDistance" or "commute" (Circuitscape type).
#'@param NN number of next neighbours recommendation is 8
#'@return a costdistance matrix between all pairs of locs
#'@description calculates a cost distance matrix, to be used with run.popgensim
#' @export
#' @examples
#' \dontrun{
#' data(possums.gl)
#' library(raster)  #needed for that example
#' landscape.sim <- readRDS(system.file("extdata","landscape.sim.rdata", package="dartR"))
#' #calculate mean centers of individuals per population
#' xy <- apply(possums.gl@other$xy, 2, function(x) tapply(x, pop(possums.gl), mean))
#' cd <- gl.costdistances(landscape.sim, xy, method="leastcost", NN=8)
#' round(cd,3)
#' }

gl.costdistances <- function(landscape, locs, method, NN)
{
# CHECK IF PACKAGES ARE INSTALLED
  pkg <- "gdistance"
  if (!(requireNamespace(pkg, quietly = TRUE))) {
    stop("Package",pkg," needed for this function to work. Please install it.") } 
  

  if (class(locs)=="genlight")
  {
    if (is.null(locs@other$latlon)) stop("no locations were provided in the genlight object [@other$latlon].\n")
    if (is.null(pop(locs))) 
      {
      cat("No population definition provided, hence I will calculate costdistances between individuals\n")
      pop(locs) <- indNames(locs)
      if (is.null(pop(locs))) pop(locs)<- 1:nInd(locs)
    }
      locs <- apply(locs@other$latlon,2,function(x) tapply(x, pop(locs), mean))
   } else locs <- as.matrix(locs)
   
  fric.mat <- gdistance::transition(landscape,function(x) 1/x[2],NN)
  #set distances to meters  if no projected already
  fric.mat@crs@projargs<- "+proj=merc +units=m"
  fric.mat.cor <- gdistance::geoCorrection(fric.mat)
  if (method=="leastcost") cd.mat <-gdistance::costDistance(fric.mat.cor, locs, locs)
  if (method=="rSPDistance") cd.mat <- gdistance::rSPDistance(fric.mat.cor, locs, locs, theta=1)
  if (method=="commute") cd.mat <-as.matrix(gdistance::commuteDistance(fric.mat.cor,locs))
  colnames(cd.mat) <- row.names(locs)
  rownames(cd.mat) <- row.names(locs)
  return (cd.mat)
}
