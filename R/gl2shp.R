#' Convert genlight objects to ESRI shapefiles
#'
#' This function exports cordinates in a genlight object to a point shape file (including also individual meta data if available). Coordinates are provided under gl@other$latlong and assumed to be in WGS84 coordinates if no proj4 string is provided. 
#' @param gl -- genlight containing lat longs  [required]
#' @param proj4 -- proj4string of data set. If not provided WGS84 is taken as default. (see spatialreference.org for other projections)
#' @param outfile -- name (path) of the output shape file
#' @param outpath -- path of the output file. Default is to tempdir(). If to be saved in the current working directory change to "."
#' @param v -- verbosity: if v=0 no output, v=1 reports name and path of output file. default 1
#' @export
#' @importFrom rgdal writeOGR
#' @importFrom sp SpatialPointsDataFrame coordinates<- CRS proj4string<-
#' @author Bernd Guber (glbugs@@aerg.canberra.edu.au)
#' @examples
#' \dontrun{
#' gl2shp(testset.gl)
#'}

gl2shp <- function(gl, proj4="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs",  outfile="gl", outpath=tempdir(), v=1)
{
  if (is.null(gl@other$latlong)) stop("No coordinates provided in slot: gl@other$latlong")
  if (nrow(gl@other$latlong)!=nInd(gl)) stop("Number of coordinates provided is different from the number of individuals in the data set.")
  #check if names are given as lat long instead lat lon
  if  (sum(match(names(gl@other$latlong),"long"), na.rm=T)==1) gl@other$latlong$lon <- gl@other$latlong$long
  
  
  glpoints <- gl@other$ind.metrics
  glpoints$lat <- gl@other$latlong$lat
  glpoints$lon <- gl@other$latlong$lon
  
  toremove <- which(!complete.cases(gl@other$latlong))
  if (v==1) cat(paste("Removed", length(toremove),"individual(s) due to missing coordinates.\n" ))
  if (v==1) cat(paste("Removed: ", indNames(gl)[toremove],"\n"))
  glpoints <- glpoints[complete.cases(glpoints),]
  
  
  glpoints$id <- 1:nrow(glpoints)
  
  coordinates(glpoints) <- c("lon","lat")
  
  #create allsites point shp files
  spdf = SpatialPointsDataFrame(glpoints, data.frame(glpoints))
  proj4string(spdf) <- CRS(proj4)
  #if (!is.null(reproj4)) spdf <- project(spdf, proj = reproj4, inv = TRUE)
  writeOGR(spdf, dsn=outpath, layer=outfile, driver="ESRI Shapefile", overwrite_layer=TRUE)

  if (v==1)  cat(paste("Shapefile saved as:", paste(outfile,".shp", sep=""),"\nin folder:",outpath))
}
