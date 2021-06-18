#' Convert genlight objects to ESRI shapefiles or kml files
#'
#' This function exports cordinates in a genlight object to a point shape file (including also individual meta data if available).
#' Coordinates are provided under x@other$latlon and assumed to be in WGS84 coordinates, if not proj4 string is provided.
#' @param x -- name of the genlight object containing the SNP data and location data, lat longs  [required]
#' @param type -- type of output "kml" or "shp" [default 'shp']
#' @param proj4 -- proj4string of data set (see spatialreference.org for projections) [default WGS84]
#' @param outfile -- name (path) of the output shape file [default 'gl']. shp extension is added automatically.
#' @param outpath -- path where to save the output file [default tempdir(), mandated by CRAN]. Use outpath=getwd() or outpath="." when calling this function to direct output files to your working directory.
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2 or as specified using gl.set.verbosity]
#' @export
#' @importFrom sp SpatialPointsDataFrame coordinates<- CRS proj4string<-
#' @importFrom stats complete.cases
#' @author Bernd Guber (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' gl2shp(testset.gl)

gl2shp <- function(x, type ="shp", proj4="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs",  outfile="gl", outpath=tempdir(), verbose=NULL){

# CHECK IF PACKAGES ARE INSTALLED
  pkg <- "rgdal"
  if (!(requireNamespace(pkg, quietly = TRUE))) {
    stop("Package",pkg," needed for this function to work. Please   install it.") } else {
  
  
# TRAP COMMAND, SET VERSION
  
  funname <- match.call()[[1]]
  build <- "Jacob"
  outfilespec <- file.path(outpath, outfile)
  
# SET VERBOSITY
  
  if (is.null(verbose)){ 
    if(!is.null(x@other$verbose)){ 
      verbose <- x@other$verbose
    } else { 
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
    stop("  Fatal Error: genlight object required!\n")
  }

  if (verbose >= 2){
    if (all(x@ploidy == 1)){
      stop("Fatal Error: Detected Presence/Absence (SilicoDArT) data. Please provide a SNP dataset\n")
    } else if (all(x@ploidy == 2)){
      cat("  Processing a SNP dataset\n")
    } else {
      stop("Fatal Error: Ploidy must be universally 1 (fragment P/A data) or 2 (SNP data)")
    }
  }

# FUNCTION SPECIFIC ERROR CHECKING

  if (is.null(x@other$latlon)){stop("Fatal Error: No coordinates provided in slot: gl@other$latlon\n")}

  if (nrow(x@other$latlon)!=nInd(x)){stop("Fatal Error: Number of coordinates provided is different from the number of individuals in the data set\n")}

  #check if names are given as lat long instead lat lon
  if (sum(match(names(x@other$latlon),"long"), na.rm=T)==1) {
    cat("  Warning: Names given as lat long, instead of lat lon. Rectifying\n")
    x@other$latlon$lon <- x@other$latlon$long
    x@other$latlon$long <- NULL
  }  
  
# DO THE JOB
  
  glpoints <- x@other$ind.metrics
  glpoints$lat <- x@other$latlon$lat
  glpoints$lon <- x@other$latlon$lon
  
  toremove <- which(!complete.cases(x@other$latlon))
  if (verbose >= 2) {cat(paste("Removed", length(toremove),"individual(s) due to missing coordinates.\n" ))}
  if (verbose >= 2 & (length(toremove) > 0)) {cat(paste("Removed: ", indNames(x)[toremove],"\n"))}
  glpoints <- glpoints[complete.cases(glpoints),]
  
  glpoints$id <- 1:nrow(glpoints)
  
  sp::coordinates(glpoints) <- c("lon","lat")
  
#create allsites point shp files
  spdf = SpatialPointsDataFrame(glpoints, data.frame(glpoints))
  proj4string(spdf) <- CRS(proj4)
  # if (!is.null(reproj4)) spdf <- project(spdf, proj = reproj4, inv = TRUE)
  if (type=="shp") rgdal::writeOGR(spdf, dsn=outpath, layer=outfile, driver="ESRI Shapefile", overwrite_layer=TRUE)

  if (type=="kml") rgdal::writeOGR(spdf,  driver = 'KML', dsn = paste0(file.path(outpath,outfile),".kml"), layer=outfile, overwrite_layer=TRUE)

  if (verbose >= 2)  cat(paste("Shapefile saved as:", paste0(outfile,".",type),"\nin folder:",outpath,"\n"))

# FLAG SCRIPT END

  if (verbose > 0) {
    cat("Completed:",funname,"\n")
  }

  return(NULL)
}
}