#'Calculates cost distances for a given landscape (resistance matrix)
#'
#'@param landscape a raster object coding the resistance of the landscape
#'@param locs coordinates of the subpopulations. If a genlight object is provided coordinates are taken from @other$latlong and centers for population (pop(gl)) are calculated. In case you want to calculate costdistances between individuals redefine pop(gl) via: \code{pop(gl)<- indNames(gl)}.
#'@param method defines the type of cost distance, types are "least-cost", "rSPDistance" or "commute (Circuitscape type)"
#'@param NN number of next neighbours recommendation is 8
#'@return a costdistance matrix between all pairs of locs
#'@description calculates a cost distance matrix, to be used with run.popgensim
#'@importFrom gdistance costDistance rSPDistance commuteDistance
#' @export

gl.costdistances <- function(landscape, locs, method, NN)
{
  if (is(locs,"genlight"))
  {
    if (is.null(locs@other$latlong)) stop("no locations were provided in the genlight object [@other$latlong].\n")
    if (is.null(pop(locs))) 
      {
      cat("No population definition provided, hence I will calculate costdistances between individuals\n")
      pop(locs) <- indNames(locs)
      if (is.null(pop(locs))) pop(locs)<- 1:nInd(locs)
    }
      locs <- apply(locs@other$latlong,2,function(x) tapply(x, pop(locs), mean))
   } else locs <- as.matrix(locs)
   
  fric.mat <- transition(landscape,function(x) 1/x[2],NN)
  #set distances to meters  if no projected already
  fric.mat@crs@projargs<- "+proj=merc +units=m"
  fric.mat.cor <- geoCorrection(fric.mat)
  if (method=="leastcost") cd.mat <-costDistance(fric.mat.cor, locs, locs)
  if (method=="rSPDistance") cd.mat <- rSPDistance(fric.mat.cor, locs, locs, theta=1)
  if (method=="commute") cd.mat <-as.matrix(commuteDistance(fric.mat.cor,locs))
  colnames(cd.mat) <- row.names(locs)
  rownames(cd.mat) <- row.names(locs)
  return (cd.mat)
}
