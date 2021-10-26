#' Convert genlight objects to ESRI shapefiles or kml files
#'
#' This function exports coordinates in a genlight object to a point shape file
#'  (including also individual meta data if available).
#' Coordinates are provided under x@other$latlon and assumed to be in WGS84
#' coordinates, if not proj4 string is provided.
#' @param x Name of the genlight object containing the SNP data and location
#' data, lat longs [required].
#' @param type Type of output 'kml' or 'shp' [default 'shp'].
#' @param proj4 Proj4string of data set (see spatialreference.org for
#' projections) [default WGS84].
#' @param outfile Name (path) of the output shape file [default 'gl']. shp
#'  extension is added automatically.
#' @param outpath Path where to save the output file
#' [default tempdir(), mandated by CRAN]. Use outpath=getwd() or outpath='.'
#' when calling this function to direct output files to your working directory.
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2 or as specified using gl.set.verbosity].
#' @export
#' @importFrom sp SpatialPointsDataFrame coordinates<- CRS proj4string<-
#' @importFrom stats complete.cases
#' @author Bernd Guber (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' gl2shp(testset.gl)

gl2shp <- function(x,
                   type = "shp",
                   proj4 = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs",
                   outfile = "gl",
                   outpath = tempdir(),
                   verbose = NULL) {
    outfilespec <- file.path(outpath, outfile)
    
    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)
    
    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     build = "Jody",
                     verbosity = verbose)
    
    # CHECK DATATYPE
    datatype <- utils.check.datatype(x, verbose = verbose)
    
    # FUNCTION SPECIFIC ERROR CHECKING
    
    # CHECK IF PACKAGES ARE INSTALLED
    pkg <- "rgdal"
    if (!(requireNamespace(pkg, quietly = TRUE))) {
        stop(error(
            "Package",
            pkg,
            " needed for this function to work. Please   install it."
        ))
    } else {
        if (is.null(x@other$latlon)) {
            stop(error(
                "Fatal Error: No coordinates provided in slot: gl@other$latlon\n"
            ))
        }
        
        if (nrow(x@other$latlon) != nInd(x)) {
            stop(
                error(
                    "Fatal Error: Number of coordinates provided is different from the number of individuals in the data set\n"
                )
            )
        }
        
        # check if names are given as lat long instead lat lon
        if (sum(match(names(x@other$latlon), "long"), na.rm = T) == 1) {
            cat(warn(
                "  Warning: Names given as lat long, instead of lat lon. Rectifying\n"
            ))
            x@other$latlon$lon <- x@other$latlon$long
            x@other$latlon$long <- NULL
        }
        
        # DO THE JOB
        
        glpoints <- x@other$ind.metrics
        glpoints$lat <- x@other$latlon$lat
        glpoints$lon <- x@other$latlon$lon
        
        toremove <- which(!complete.cases(x@other$latlon))
        if (verbose >= 2) {
            cat(warn(
                paste(
                    "Removed",
                    length(toremove),
                    "individual(s) due to missing coordinates.\n"
                )
            ))
        }
        if (verbose >= 2 & (length(toremove) > 0)) {
            cat(warn(paste(
                "Removed: ", indNames(x)[toremove], "\n"
            )))
        }
        glpoints <- glpoints[complete.cases(glpoints),]
        
        glpoints$id <- 1:nrow(glpoints)
        
        sp::coordinates(glpoints) <- c("lon", "lat")
        
        # create all sites point shp files
        spdf = SpatialPointsDataFrame(glpoints, data.frame(glpoints))
        proj4string(spdf) <- CRS(proj4)
        # if (!is.null(reproj4)) spdf <- project(spdf, proj = reproj4, inv = TRUE)
        if (type == "shp")
            rgdal::writeOGR(
                spdf,
                dsn = outpath,
                layer = outfile,
                driver = "ESRI Shapefile",
                overwrite_layer = TRUE
            )
        
        if (type == "kml")
            rgdal::writeOGR(
                spdf,
                driver = "KML",
                dsn = paste0(file.path(outpath, outfile), ".kml"),
                layer = outfile,
                overwrite_layer = TRUE
            )
        
        if (verbose >= 2)
            cat(report(
                paste(
                    "Shapefile saved as:",
                    paste0(outfile, ".", type),
                    "\nin folder:",
                    outpath,
                    "\n"
                )
            ))
        
        # FLAG SCRIPT END
        
        if (verbose > 0) {
            cat(report("Completed:", funname, "\n"))
        }
        
        return(NULL)
    }
}
