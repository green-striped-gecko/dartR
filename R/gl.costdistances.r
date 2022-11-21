#' Calculates cost distances for a given landscape (resistance matrix)
#'
#' @param landscape A raster object coding the resistance of the landscape
#' [required].
#' @param locs Coordinates of the subpopulations. If a genlight object is
#' provided coordinates are taken from @other$latlon and centers for population
#' (pop(gl)) are calculated. In case you want to calculate costdistances between
#' individuals redefine pop(gl) via: \code{pop(gl)<- indNames(gl)} [required].
#' @param method Defines the type of cost distance, types are 'leastcost',
#' 'rSPDistance' or 'commute' (Circuitscape type) [required].
#' @param NN Number of next neighbours recommendation is 8 [required].
#' @return A costdistance matrix between all pairs of locs.
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log ; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#' @description
#' Calculates a cost distance matrix, to be used with run.popgensim.
#' @export
#' @examples
#' \dontrun{
#' data(possums.gl)
#' library(raster)  #needed for that example
#' landscape.sim <- readRDS(system.file('extdata','landscape.sim.rdata', 
#' package='dartR'))
#' #calculate mean centers of individuals per population
#' xy <- apply(possums.gl@other$xy, 2, function(x) tapply(x, pop(possums.gl),
#'  mean))
#' cd <- gl.costdistances(landscape.sim, xy, method='leastcost', NN=8)
#' round(cd,3)
#' }

gl.costdistances <- function(landscape,
                             locs,
                             method,
                             NN,
                             verbose = NULL) {
    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)
    
    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     build = "Jody",
                     verbosity = verbose)
    
    # CHECK IF PACKAGES ARE INSTALLED
    pkg <- "gdistance"
    if (!(requireNamespace(pkg, quietly = TRUE))) {
      cat(error(
        "Package",
        pkg,
        " needed for this function to work. Please install it.\n"
      ))
      return(-1)
    }
    
    if (is(locs, "genlight")) {
        if (is.null(locs@other$latlon)) {
            stop(
                error(
                    "No locations were provided in the genlight object
                    [@other$latlon].\n"
                )
            )
        }
        if (is.null(pop(locs))) {
            cat(
                warn(
                    "No population definition provided, hence I will calculate 
                    costdistances between individuals\n"
                )
            )
            pop(locs) <- indNames(locs)
            if (is.null(pop(locs))) {
                pop(locs) <- 1:nInd(locs)
            }
        }
        locs <-
            apply(locs@other$latlon, 2, function(x)
                tapply(x, pop(locs), mean))
    } else {
        locs <- as.matrix(locs)
    }
    
    fric.mat <-
        gdistance::transition(landscape, function(x)
            1 / x[2], NN)
    # set distances to meters if no projected already
    fric.mat@crs@projargs <- "+proj=merc +units=m"
    fric.mat.cor <- gdistance::geoCorrection(fric.mat)
    if (method == "leastcost")
        cd.mat <- gdistance::costDistance(fric.mat.cor, locs, locs)
    if (method == "rSPDistance")
        cd.mat <-
        gdistance::rSPDistance(fric.mat.cor, locs, locs, theta = 1)
    if (method == "commute")
        cd.mat <-
        as.matrix(gdistance::commuteDistance(fric.mat.cor, locs))
    colnames(cd.mat) <- row.names(locs)
    rownames(cd.mat) <- row.names(locs)
    
    # FLAG SCRIPT END
    
    if (verbose > 0) {
        cat(report("Completed:", funname, "\n"))
    }
    
    return(cd.mat)
}
