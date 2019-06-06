#' Convert genlight objects to ESRI shapefiles or kml files
#'
#' This function exports cordinates in a genlight object to a point shape file (including also individual meta data if available).
#' Coordinates are provided under x@other$latlong and assumed to be in WGS84 coordinates, if not proj4 string is provided.
#' @param x -- name of the genlight object containing the SNP data and location data, lat longs  [required]
#' @param type -- type of output "kml" or "shp" [default 'shp']
#' @param proj4 -- proj4string of data set (see spatialreference.org for projections) [default WGS84]
#' @param outfile -- name (path) of the output shape file [default 'gl']. shp extension is added automatically.
#' @param outpath -- path where to save the output file [default tempdir(), mandated by CRAN]. Use outpath=getwd() or outpath="." when calling this function to direct output files to your working directory.
#' @param verbose -- specify the level of verbosity: 0, silent, fatal errors only; 1, flag function begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @export
#' @importFrom rgdal writeOGR
#' @importFrom sp SpatialPointsDataFrame coordinates<- CRS proj4string<-
#' @importFrom stats complete.cases
#' @author Bernd Guber (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' gl2shp(testset.gl)

gl2shp <- function(x, type ="shp", proj4="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs",  outfile="gl", outpath=tempdir(), verbose=2){

# TIDY UP FILE SPECS

  outfilespec <- file.path(outpath, outfile)
  funname <- match.call()[[1]]

# FLAG SCRIPT START

  if (verbose < 0 | verbose > 5){
    cat("  Warning: Parameter 'verbose' must be an integer between 0 [silent] and 5 [full report], set to 2\n")
    verbose <- 2
  }

  if (verbose > 0) {
    cat("Starting",funname,"\n")
  }

# STANDARD ERROR CHECKING
  
  if(class(x)!="genlight") {
    cat("  Fatal Error: genlight object required!\n"); stop("Execution terminated\n")
  }

  # Work around a bug in adegenet if genlight object is created by subsetting
      if (nLoc(x)!=nrow(x@other$loc.metrics)) { stop("The number of rows in the loc.metrics table does not match the number of loci in your genlight object!")  }

  # Set a population if none is specified (such as if the genlight object has been generated manually)
    if (is.null(pop(x)) | is.na(length(pop(x))) | length(pop(x)) <= 0) {
      if (verbose >= 2){ cat("  Population assignments not detected, individuals assigned to a single population labelled 'pop1'\n")}
      pop(x) <- array("pop1",dim = nInd(x))
      pop(x) <- as.factor(pop(x))
    }

  # Check for monomorphic loci
    tmp <- gl.filter.monomorphs(x, verbose=0)
    if ((nLoc(tmp) < nLoc(x)) & verbose >= 2) {cat("  Warning: genlight object contains monomorphic loci\n")}

# FUNCTION SPECIFIC ERROR CHECKING

  if (is.null(x@other$latlong)){stop("Fatal Error: No coordinates provided in slot: gl@other$latlong\n")}

  if (nrow(x@other$latlong)!=nInd(x)){stop("Fatal Error: Number of coordinates provided is different from the number of individuals in the data set\n")}

  #check if names are given as lat long instead lat lon
  if (sum(match(names(x@other$latlong),"long"), na.rm=T)==1) {
    cat("  Warning: Names given as lat long, instead of lat lon. Rectifying\n")
    x@other$latlong$lon <- x@other$latlong$long
  }  
  
# DO THE JOB
  
  glpoints <- x@other$ind.metrics
  glpoints$lat <- x@other$latlong$lat
  glpoints$lon <- x@other$latlong$lon
  
  toremove <- which(!complete.cases(x@other$latlong))
  if (verbose >= 2) {cat(paste("Removed", length(toremove),"individual(s) due to missing coordinates.\n" ))}
  if (verbose >= 2 & (length(toremove) > 0)) {cat(paste("Removed: ", indNames(x)[toremove],"\n"))}
  glpoints <- glpoints[complete.cases(glpoints),]
  
  glpoints$id <- 1:nrow(glpoints)
  
  sp::coordinates(glpoints) <- c("lon","lat")
  
#create allsites point shp files
  spdf = SpatialPointsDataFrame(glpoints, data.frame(glpoints))
  proj4string(spdf) <- CRS(proj4)
  # if (!is.null(reproj4)) spdf <- project(spdf, proj = reproj4, inv = TRUE)
  if (type=="shp") writeOGR(spdf, dsn=outpath, layer=outfile, driver="ESRI Shapefile", overwrite_layer=TRUE)

  if (type=="kml") writeOGR(spdf,  driver = 'KML', dsn = paste0(file.path(outpath,outfile),".kml"), layer=outfile, overwrite_layer=TRUE)

  if (verbose >= 2)  cat(paste("Shapefile saved as:", paste0(outfile,".",type),"\nin folder:",outpath,"\n"))

# FLAG SCRIPT END

  if (verbose > 0) {
    cat("Completed:",funname,"\n")
  }

  return(NULL)
}
