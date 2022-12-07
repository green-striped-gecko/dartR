#' Creates an interactive map (based on latlon) from a genlight object
#'
#' @param x A genlight object (including coordinates within the latlon slot) 
#' [required].
#' @param matrix A distance matrix between populations or individuals. The
#' matrix is visualised as lines between individuals/populations. If matrix is
#' asymmetric two lines with arrows are plotted [default NULL].
#' @param standard If a matrix is provided line width will be standardised to be
#' between 1 to 10, if set to true, otherwise taken as given [default TRUE].
#' @param symmetric If a symmetric matrix is provided only one line is drawn
#' based on the lower triangle of the matrix. If set to false arrows indicating
#' the direction are used instead [default TRUE].
#' @param pop.labels Population labels at the center of the individuals of
#'  populations [default TRUE].
#' @param pop.labels.cex Size of population labels [default 12].
#' @param ind.circles Should individuals plotted as circles [default TRUE].
#' @param ind.circle.cols Colors of circles. Colors can be provided as usual by 
#' names (e.g. "black") and are re-cycled. So a color c("blue","red") colors 
#' individuals alternatively between blue and red using the genlight object
#'  order of individuals. For transparency see parameter 
#'  ind.circle.transparency. Defaults to rainbow colors by population  if not
#'   provided. If you want to have your own colors for each population, check
#'    the platypus.gl example below.
#' @param ind.circle.cex (size or circles in pixels ) [default 10].
#' @param ind.circle.transparency Transparency of circles between 0=invisible 
#' and 1=no transparency. Defaults to 0.8.
#' @param palette_links Color palette for the links in case a matrix is provided
#'  [default NULL].
#' @param leg_title Legend's title for the links in case a matrix is provided
#'  [default NULL].
#' @param provider Passed to leaflet [default "Esri.NatGeoWorldMap"].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#' @return plots a map
#' @importFrom methods is
#' @export
#' @details 
#' A wrapper around the \pkg{leaflet} package. For possible background 
#' maps check as specified via the provider:
#' \url{http://leaflet-extras.github.io/leaflet-providers/preview/index.html}
#' 
#' The palette_links argument can be any of the following:
#' A character vector of RGB or named colors. Examples: palette(), 
#' c("#000000", "#0000FF", "#FFFFFF"), topo.colors(10)
#' 
#' The name of an RColorBrewer palette, e.g. "BuPu" or "Greens".
#' 
#' The full name of a viridis palette: "viridis", "magma", "inferno", 
#' or "plasma".
#' 
#' A function that receives a single value between 0 and 1 and returns a color.
#'  Examples: colorRamp(c("#000000", "#FFFFFF"), interpolate = "spline").
#'  
#' @author Bernd Gruber -- Post to \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' require("dartR.data")
#' gl.map.interactive(bandicoot.gl)
#' cols <- c("red","blue","yellow")[as.numeric(pop(platypus.gl))]
#' gl.map.interactive(platypus.gl, ind.circle.cols=cols, ind.circle.cex=10, 
#' ind.circle.transparency=0.5)

gl.map.interactive <- function(x,
                               matrix = NULL,
                               standard = TRUE,
                               symmetric = TRUE,
                               pop.labels = TRUE,
                               pop.labels.cex = 12,
                               ind.circles = TRUE,
                               ind.circle.cols = NULL,
                               ind.circle.cex = 10,
                               ind.circle.transparency = 0.8,        
                               palette_links = NULL,
                               leg_title = NULL,
                               provider = "Esri.NatGeoWorldMap",
                               verbose = NULL) {
    
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
    pkg <- "leaflet"
    if (!(requireNamespace(pkg, quietly = TRUE))) {
      cat(error(
        "Package",
        pkg,
        " needed for this function to work. Please install it.\n"
      ))
      return(-1)
    }
    
    pkg <- "leaflet.minicharts"
    if (!(requireNamespace(pkg, quietly = TRUE))) {
      cat(error(
        "Package",
        pkg,
        " needed for this function to work. Please install it.\n"
      ))
      return(-1)
    } else {
       
        if (is.null(x@other$latlon)) {
            stop(error(
                "No valid coordinates are supplied at gl@other$latlon"
            ))
        }
        
        if (sum(colnames(x@other$latlon) %in% c("lat", "lon")) != 2) {
            stop(error(
 "Coordinates under gl@other$latlon are not named 'lat' and 'lon'."
            ))
        }
        
        if (!is.null(matrix)) {
            if (nrow(matrix) != nInd(x) & nrow(matrix) != nPop(x)) {
                stop(
                    error(
"The dimension of the provided matrix does neither match the number of 
individuals nor the number of populations."
                    )
                )
            }
        }
        
        if (is.null(ind.circle.cols))
            {
            cols <- rainbow(nPop(x))
            cols <- substr(cols, 1, 7)
            ic <- cols[as.numeric(pop(x))]
        } else ic <- ind.circle.cols
        
        
        df <- x@other$latlon
        centers <-
            apply(df, 2, function(xx)
                tapply(xx, pop(x), mean, na.rm = TRUE))
        # when there is just one population the output of centers is a vector 
        #the following lines fix this error
        if (nPop(x) == 1) {
            centers <- data.frame(lon = centers[1], lat = centers[2])
            row.names(centers) <- popNames(x)
        }
        # Add default OpenStreetMap map tiles
        m <- leaflet::leaflet() %>%
            leaflet::addTiles()
        
        if (ind.circles) {
            m <- m %>%
                leaflet::addCircles(
                    lng = df$lon,
                    lat = df$lat,
                    popup = indNames(x),
                    color = ic,
                    opacity = ind.circle.transparency,
                    weight = ind.circle.cex
                    
                )
        }
        
        if (pop.labels) {
            m <- m %>%
                leaflet::addLabelOnlyMarkers(
                    lng = centers[, "lon"],
                    lat = centers[, "lat"],
                    label = popNames(x),
                    labelOptions = leaflet::labelOptions(
                        noHide = T,
                        direction = "top",
                        textOnly = T,
                        textsize = paste0(pop.labels.cex, "px")
                    )
                )
        }
        
        if (!is.null(matrix)) {
          matrix <- matrix[order(indNames(x)),]
            # standardize
            if (standard) {
                matrix[, ] <-
                    ((matrix[, ] - min(matrix, na.rm = T)) / 
                       (max(matrix, na.rm = T) - 
                          min(matrix, na.rm = T))) * 9 + 1
            }
            
            if (nrow(matrix) == nPop(x)) {
                xys <- centers
            } else {
                xys <- df
            }
          
          if(is.null(palette_links)){
            palette_links <- 
          diverging_palette(length(unique(unlist(unname(as.vector(matrix))))))
          }
            
          qpal <- leaflet::colorNumeric(
            palette = palette_links,
            domain = unique(unlist(unname(as.vector(matrix)))))
          
            if (symmetric) {
                for (ii in 1:nrow(matrix)) {
                    for (i in ii:nrow(matrix)) {
                        if (!is.null(matrix[i, ii]) & matrix[i, ii] > 0 ){
                            m <- m %>%
                                leaflet::addPolylines(
                                    lng = c(xys[i, "lon"], xys[ii, "lon"]),
                                    lat = c(xys[i, "lat"], xys[ii, "lat"]),
                                    color = qpal(matrix[i,ii]),
                                    opacity = 1
                                )
                        }else{
                          next()
                        }
                    }
                }
              m <- m %>% leaflet::addLegend(
                pal = qpal, 
                values = unique(unlist(unname(as.vector(matrix)))), 
                group = "addPolylines", 
                position = "bottomleft",
                title = leg_title) 
       
            }
            
            if (!symmetric) {
                for (i in 1:nrow(matrix)) {
                    for (ii in 1:nrow(matrix)) {
                        if (abs((i - ii)) != 0) {
                            from <- xys[i, ]
                            to <- xys[ii, ]
                            if (!is.null(matrix[i, ii]) &
                                !is.null(matrix[ii, i])) {
                                if (matrix[i, ii] > matrix[ii, i]){
                                    lcols <-"#FFAA00"
                                }else{
                                    lcols <-"#00AAFF"
                                }
                                if (matrix[i, ii] == matrix[ii, i]){
                                    lcols <-"#00AA00"
                                }
                            } else{
                                lcols <-"#333333"
                            }
                            m <- m %>%
                                leaflet.minicharts::addFlows(
                                    lng0 = as.numeric(from["lon"]),
                                    lng1 = as.numeric(to["lon"]),
                                    lat0 = as.numeric(from["lat"]),
                                    lat1 = as.numeric(to["lat"]),
                                    flow = matrix[i, ii],
                                    color = lcols,
                                    maxThickness = 10,
                                    minThickness = 0,
                                    maxFlow = max(matrix,
                                                  na.rm = T),
                                    opacity = 0.8
                                )
                        }
                    }
                }
            }
        }
        # FLAG SCRIPT END
        
        if (verbose >= 1) {
            cat(report("Completed:", funname, "\n"))
        }
        
        # RETURN
        m %>%
            leaflet::addProviderTiles(provider)
        
    }
}
