#' @name gl.map.structure
#' @title Maps a STRUCTURE plot using a genlight object
#' @description
#' This function takes the output of plotstructure (the q matrix) and maps the
#' q-matrix across using the population centers from the genlight object that
#' was used to run the structure analysis via \code{\link{gl.run.structure}})
#' and plots the typical structure bar plots on a spatial map, providing a
#' barplot for each subpopulation. Therefore it requires coordinates from a
#'  genlight object. This kind of plots should support the interpretation of the
#'   spatial structure of a population, but in principle is not different from
#'   \code{\link{gl.plot.structure}}
#'
#' @param qmat Q-matrix from a structure run followed by a clumpp run object
#' [from \code{\link{gl.run.structure}} and \code{\link{gl.plot.structure}}]
#'  [required].
#' @param x Name of the genlight object containing the coordinates in the
#'  \code{\@other$latlon} slot to calculate the population centers [required].
#' @param provider Provider	passed to leaflet. Check \link[leaflet]{providers}
#' for a list of possible backgrounds [default "Esri.NatGeoWorldMap"].
#' @param scalex Scaling factor to determine the size of the bars in x direction 
#' [default 1].
#' @param scaley Scaling factor to determine the size of the bars in y direction
#'  [default 1].
#' @param movepops A two-dimensional data frame that allows to move the center of
#' the barplots manually in case they overlap. Often if populations are
#' horizontally close to each other. This needs to be a data.frame of the
#' dimensions [rows=number of populations, columns = 2 (lon/lat)]. For each
#' population you have to specify the x and y (lon and lat) units you want to
#' move the center of the plot, (see example for details) [default NULL].
#' @param pop.labels Switch for population labels below the parplots 
#' [default TRUE].
#' @param pop.labels.cex Size of population labels [default 12].
#' @return An interactive map that shows the structure plots broken down by 
#' population.
#' @author Bernd Gruber (Post to \url{https://groups.google.com/d/forum/dartr})
#' @details
#' Creates a mapped version of structure plots. For possible background maps
#' check as specified via the provider:
#' \url{http://leaflet-extras.github.io/leaflet-providers/preview/index.html}.
#' You may need to adjust scalex and scaley values [default 1], as the size
#' depends on the scale of the map and the position of the populations.
#' @return returns the map and a list of the qmat split into sorted matrices per
#'  population. This can be used to create your own map.
#' @examples
#' \dontrun{
#' #CLUMPP needs to be installed to be able to run the example
#' #bc <- bandicoot.gl[,1:100]
#' #sr <- gl.run.structure(bc, k.range = 2:5, num.k.rep = 3, exec = './structure.exe')
#' #ev <- gl.evanno(sr)
#' #ev
#' #qmat <- gl.plot.structure(sr, k=3, CLUMPP='d:/structure/')
#' #head(qmat)
#' #gl.map.structure(qmat, bc, scalex=1, scaley=0.5)
#' #move population 4 (out of 5) 0.5 degrees to the right and populations 1
#' #0.3 degree to the top of the map.
#' #mp <- data.frame(lon=c(0,0,0,0.5,0), lat=c(-0.3,0,0,0,0))
#' #gl.map.structure(qmat, bc, scalex=1, scaley=0.5, movepops=mp)
#' }
#' @export
#' @seealso \code{\link{gl.run.structure}},  \link[strataG]{clumpp},
#' \code{\link{gl.plot.structure}}
#' @references
#' \itemize{
#' \item Pritchard, J.K., Stephens, M., Donnelly, P. (2000) Inference of
#' population structure using multilocus genotype data. Genetics 155, 945-959.
#' \item Archer, F. I., Adams, P. E. and Schneiders, B. B. (2016) strataG: An R
#'  package for manipulating, summarizing and analysing population genetic data.
#'   Mol Ecol Resour. doi:10.1111/1755-0998.12559
#' \item Evanno, G., Regnaut, S., and J. Goudet. 2005. Detecting the number of
#' clusters of individuals using the software STRUCTURE: a simulation study.
#' Molecular Ecology 14:2611-2620.
#' \item Mattias Jakobsson and Noah A. Rosenberg. 2007. CLUMPP: a cluster
#' matching and permutation program for dealing with label switching and
#' multimodality in analysis of population structure. Bioinformatics
#' 23(14):1801-1806. Available at
#' \href{http://web.stanford.edu/group/rosenberglab/clumppDownload.html}{clumpp}
#' }

gl.map.structure <- function(qmat,
                             x,
                             provider = "Esri.NatGeoWorldMap",
                             scalex = 1,
                             scaley = 1,
                             movepops = NULL,
                             pop.labels = TRUE,
                             pop.labels.cex = 12) {
    ff <- qmat[, 4:(ncol(qmat))]
    
    df <- x@other$latlon
    centers <-
        apply(df, 2, function(xx)
            tapply(xx, pop(x), mean, na.rm = TRUE))
    
    if (!is.null(movepops)) {
        if (nrow(movepops) != nrow(centers))
            stop(
                error(
                    "The provided movepops data.frame has not the corret number of rows, please check. It needs to have the same numbers of rows as the number populations in your genlight object."
                )
            )
        centers[, 1] <- centers[, 1] + movepops[, 1]
        centers[, 2] <- centers[, 2] + movepops[, 2]
        
    }
    
    
    sc <-
        match(rownames(centers), levels(factor(qmat$orig.pop)))
    if (any(is.na(sc)))
        cat(
            error(
                "Population names (coordinates) in the genlight object do not match population in your q-matrix. Please check both."
            )
        )
    centers <- centers[sc, ]
    cx <- centers[, "lon"]
    cy <- centers[, "lat"]
    sx <- abs(diff(range(centers[, "lon"]))) / (100) * scalex
    sy <- 20 * sx * scaley
    #
    qmat$orig.pop <- factor(qmat$orig.pop)
    npops <- length(levels(qmat$orig.pop))
    ll <- data.frame(cbind(as.numeric(qmat$orig.pop), ff))
    zz <- do.call(order, unname(as.list(ll)))
    bb <- qmat[zz, ]
    bb$orig.pop <- factor(bb$orig.pop)
    ff <- bb[, 4:(ncol(bb))]
    
    out <- list()
    m1 <- leaflet::leaflet() %>%
        leaflet::addProviderTiles(provider = provider)
    for (p in 1:npops) {
        qmi <- ff[bb$orig.pop == levels(bb$orig.pop)[p], ]
        out[[p]] <- bb[bb$orig.pop == levels(bb$orig.pop)[p], ]
        names(out)[p] <- levels(bb$orig.pop)[p]
        qmi1 <- cbind(rep(0, nrow(qmi)), qmi)
        for (xx in 1:nrow(qmi1))
            qmi1[xx, ] <- cumsum(as.numeric(qmi1[xx, ]))
        
        
        for (ii in 1:nrow(qmi1)) {
            for (i in 1:(ncol(qmi1) - 1)) {
                oo <- (ii - nrow(qmi) / 2) * sx
                
                m1 <- m1 %>%
                    leaflet::addRectangles(
                        cx[p] + oo,
                        cy[p] + qmi1[ii, i] * sy,
                        cx[p] + oo + sx,
                        cy[p] + qmi1[ii, i + 1] * sy,
                        opacity = 0,
                        color = rainbow(ncol(ff))[i],
                        fillOpacity = 0.8
                    )
                
            }
        }
        
    }
    if (pop.labels)
        m1 <- m1 %>%
        leaflet::addLabelOnlyMarkers(
            lng = centers[, "lon"],
            lat = centers[, "lat"] - sy * 0.1,
            label = rownames(centers),
            labelOptions = leaflet::labelOptions(
                noHide = T,
                direction = "center",
                textOnly = T,
                textsize = paste0(pop.labels.cex, "px")
            )
        )
    
    
    print(m1)
    # mapshot(m1, file='./Rplot.png', remove_controls = TRUE)
    return(out)
    # %>% addLegend(labels=paste('Group',1:ncol(ff)), colors=rainbow(ncol(ff)),position ='topright' )
    
    # if (save) mapshot(m1, file='./Rplot.png', remove_controls = TRUE)
    
    
}
